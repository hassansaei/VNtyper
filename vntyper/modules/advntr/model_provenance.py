"""Identify the adVNTR model a run will use, and refuse the combinations that lie.

adVNTR derives its read-fetch window from the model's own content. Before #268 the
shipped hg38 model described 840 bp of an array GRCh38 carries at 3,525 bp, so adVNTR
fetched a quarter of the locus and every read outside it was invisible -- silently.
A v2 model records the array's genomic end explicitly, which is what corrects it.

That correction only holds if the model and the binary agree. adVNTR 2.0.3 selects the
legacy columns by name, so it would ignore a recorded end and reproduce the truncated
window with no error at all. Wrong answers are worse than a refusal, so this module
refuses, and it does so before the run rather than after.

It also produces the provenance record `pipeline_summary.json` carries, so a result can
be traced to the model that produced it.
"""

from __future__ import annotations

import hashlib
import logging
import re
import shlex
import sqlite3
import subprocess
import threading
import time
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, cast

logger = logging.getLogger(__name__)

#: adVNTR releases before this ignore a model's recorded genomic end.
SPAN_AWARE_ADVNTR = (2, 0, 4)

_BARE_VERSION_LINE = re.compile(r"^\s*(\d+)\.(\d+)\.(\d+)\s*$")
_LEGACY_VERSION_LINE = re.compile(r"^\s*adVNTR\s+(\d+)\.(\d+)\.(\d+):[^\r\n]*$")

# Measured under concurrent ``mamba run`` launches: mamba can exit zero while
# surrounding the tool's real output with these lock-contention diagnostics. That is
# a launcher failure, not malformed adVNTR output. Keep the classification narrower
# than either word on its own so an ordinary status-zero parse failure is never retried.
_MAMBA_RESOURCE_LOCK_MARKER = "Could not set lock (Resource temporarily unavailable)"
_MAMBA_PROCESS_LOCK = re.compile(r"Cannot lock '[^'\r\n]*/\.cache/mamba/proc'")
_MAMBA_PROCESS_WAIT_MARKER = "Waiting for other mamba process to finish"
_VERSION_PROBE_ATTEMPTS = 3
_VERSION_PROBE_RETRY_SECONDS = 0.1

#: Every statement this module runs, keyed by the two table names it accepts. No SQL is
#: assembled from a variable, so there is no identifier for anything to be injected into.
#: SQLite cannot parameterise an identifier, so a fixed literal is the only safe form.
_PRAGMA = {
    "vntrs": "PRAGMA table_info(vntrs)",
    "vntrs_v2": "PRAGMA table_info(vntrs_v2)",
}
_SELECT = {
    "vntrs": {
        False: "SELECT id, chromosome, ref_start, repeats FROM vntrs",
        True: "SELECT id, chromosome, ref_start, repeats, ref_end FROM vntrs",
    },
    "vntrs_v2": {
        False: "SELECT id, chromosome, ref_start, repeats FROM vntrs_v2",
        True: "SELECT id, chromosome, ref_start, repeats, ref_end FROM vntrs_v2",
    },
}

#: The MUC1 VNTR. A database carrying more rows is selected from, not rejected.
MUC1_VID = 25561

_REINSTALL = "Reinstall the reference bundle with `vntyper install-references -d reference --references hg38`."


class AdvntrModelError(RuntimeError):
    """The resolved adVNTR model cannot be used with the available adVNTR."""


class AdvntrProbeStatus(Enum):
    """Terminal classification of one bounded adVNTR version probe."""

    VERIFIED = "verified"
    VERSIONED_LAUNCH_FAILURE = "versioned_launch_failure"
    LAUNCH_FAILURE = "launch_failure"
    UNPARSEABLE_SUCCESS = "unparseable_success"
    TRANSIENT_EXHAUSTED = "transient_exhausted"


@dataclass(frozen=True)
class AdvntrVersionOutcome:
    """Typed terminal result carried from the version subprocess to the pipeline guard."""

    status: AdvntrProbeStatus
    version: tuple[int, int, int] | None = None
    message: str = ""

    def __post_init__(self) -> None:
        """Reject every state that contradicts the outcome discriminant.

        Raises:
            ValueError: If the status, version, and message do not form one valid outcome.
        """
        error: str | None = None
        if not isinstance(self.status, AdvntrProbeStatus):
            error = "Invalid adVNTR version outcome: status must be an AdvntrProbeStatus."
        elif self.status in {AdvntrProbeStatus.VERIFIED, AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE}:
            version = self.version
            if not (
                type(version) is tuple
                and len(version) == 3
                and all(type(part) is int and part >= 0 for part in version)
            ):
                error = f"Invalid {self.status.value} adVNTR version outcome: version must be exactly three non-negative integers."
            elif self.status is AdvntrProbeStatus.VERIFIED and self.message != "":
                error = "Invalid verified adVNTR version outcome: failure message must be empty."
            elif self.status is AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE and (
                not isinstance(self.message, str) or not self.message.startswith("adVNTR version launch failed")
            ):
                error = (
                    "Invalid versioned_launch_failure adVNTR version outcome: message must start with "
                    "'adVNTR version launch failed'."
                )
        elif self.version is not None:
            error = f"Invalid {self.status.value} adVNTR version outcome: version must be None."
        else:
            message_prefix = {
                AdvntrProbeStatus.LAUNCH_FAILURE: "adVNTR version launch failed",
                AdvntrProbeStatus.UNPARSEABLE_SUCCESS: "adVNTR version command succeeded",
                AdvntrProbeStatus.TRANSIENT_EXHAUSTED: "adVNTR version detection exhausted",
            }[self.status]
            if not isinstance(self.message, str) or not self.message.startswith(message_prefix):
                error = (
                    f"Invalid {self.status.value} adVNTR version outcome: message must start with '{message_prefix}'."
                )

        if error is not None:
            logger.error(error)
            raise ValueError(error)


class AdvntrVersionProbe:
    """Concurrency-safe successful-version cache scoped to one pipeline run.

    A pipeline creates one instance and discards it when the run ends. Versions parsed
    from a successful subprocess are cached by the configured command, including
    known-but-incompatible versions. Launch failures (even one carrying a legacy
    banner) and malformed output are never cached.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._versions: dict[str, tuple[int, int, int]] = {}

    def detect(self, command: str) -> AdvntrVersionOutcome:
        """Classify the configured command after bounded lock-contention retries.

        Args:
            command: Configured adVNTR command, possibly a multi-token ``mamba run`` command.

        Returns:
            A typed verified, versioned-launch-failure, unversioned-launch-failure,
            unparseable-success, or exhausted outcome.
        """
        with self._lock:
            cached = self._versions.get(command)
            if cached is not None:
                return AdvntrVersionOutcome(AdvntrProbeStatus.VERIFIED, version=cached)

            try:
                argv = shlex.split(command) + ["--version"]
            except ValueError as exc:
                return AdvntrVersionOutcome(
                    AdvntrProbeStatus.LAUNCH_FAILURE,
                    message=f"adVNTR version launch failed: invalid tools.advntr command: {exc}.",
                )
            attempt = 1
            while True:
                try:
                    result = subprocess.run(argv, capture_output=True, text=True, check=False)
                except FileNotFoundError:
                    return AdvntrVersionOutcome(
                        AdvntrProbeStatus.LAUNCH_FAILURE,
                        message=f"adVNTR version launch failed: command not found: {command}.",
                    )
                except PermissionError:
                    return AdvntrVersionOutcome(
                        AdvntrProbeStatus.LAUNCH_FAILURE,
                        message=f"adVNTR version launch failed: permission denied: {command}.",
                    )
                except OSError as exc:
                    return AdvntrVersionOutcome(
                        AdvntrProbeStatus.LAUNCH_FAILURE,
                        message=f"adVNTR version launch failed for {command}: {exc}.",
                    )

                combined_output = f"{result.stdout}\n{result.stderr}"
                version = (
                    _parse_probe_output(result.stdout, result.stderr)
                    if result.returncode == 0
                    else _parse_failed_probe_output(result.stdout, result.stderr)
                )
                if result.returncode == 0 and version is not None:
                    self._versions[command] = version
                    return AdvntrVersionOutcome(AdvntrProbeStatus.VERIFIED, version=version)

                if _is_transient_mamba_lock_failure(argv, combined_output):
                    if attempt == _VERSION_PROBE_ATTEMPTS:
                        msg = (
                            f"adVNTR version detection exhausted {_VERSION_PROBE_ATTEMPTS} attempts "
                            "after transient mamba launch failures."
                        )
                        return AdvntrVersionOutcome(AdvntrProbeStatus.TRANSIENT_EXHAUSTED, message=msg)
                    logger.warning(
                        f"Transient mamba launch failure while detecting the adVNTR version; "
                        f"retrying ({attempt}/{_VERSION_PROBE_ATTEMPTS})."
                    )
                    attempt += 1
                    time.sleep(_VERSION_PROBE_RETRY_SECONDS)
                    continue

                if result.returncode != 0:
                    if version is not None:
                        return AdvntrVersionOutcome(
                            AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE,
                            version=version,
                            message=f"adVNTR version launch failed: command exited with status {result.returncode}.",
                        )
                    return AdvntrVersionOutcome(
                        AdvntrProbeStatus.LAUNCH_FAILURE,
                        message=f"adVNTR version launch failed: command exited with status {result.returncode}.",
                    )
                return AdvntrVersionOutcome(
                    AdvntrProbeStatus.UNPARSEABLE_SUCCESS,
                    message="adVNTR version command succeeded but its response was unparseable or ambiguous.",
                )


def _is_transient_mamba_lock_failure(argv: list[str], output: str) -> bool:
    """Return whether output is the measured mamba process-lock contention failure."""
    measured_process_lock = bool(_MAMBA_PROCESS_LOCK.search(output)) and _MAMBA_PROCESS_WAIT_MARKER in output
    return (
        len(argv) >= 2
        and Path(argv[0]).name == "mamba"
        and argv[1] == "run"
        and "libmamba" in output
        and (_MAMBA_RESOURCE_LOCK_MARKER in output or measured_process_lock)
    )


def _parse_probe_output(stdout: str, stderr: str) -> tuple[int, int, int] | None:
    """Return one unambiguous version from strict adVNTR answer lines on either stream."""
    candidates: set[tuple[int, int, int]] = set()
    for stream in (stdout, stderr):
        for line in stream.splitlines():
            match = _BARE_VERSION_LINE.fullmatch(line) or _LEGACY_VERSION_LINE.fullmatch(line)
            if match is not None:
                candidates.add((int(match.group(1)), int(match.group(2)), int(match.group(3))))
    if len(candidates) != 1:
        return None
    return next(iter(candidates))


def _parse_failed_probe_output(stdout: str, stderr: str) -> tuple[int, int, int] | None:
    """Return one tagged legacy banner version from a failed subprocess."""
    candidates: set[tuple[int, int, int]] = set()
    for stream in (stdout, stderr):
        for line in stream.splitlines():
            match = _LEGACY_VERSION_LINE.fullmatch(line)
            if match is not None:
                candidates.add((int(match.group(1)), int(match.group(2)), int(match.group(3))))
    if len(candidates) != 1:
        return None
    return next(iter(candidates))


def parse_advntr_version(text: str | None) -> tuple[int, int, int] | None:
    """Extract one unambiguous version from strict adVNTR answer lines.

    Returns None when no version can be read. A strict tagged legacy banner can identify
    adVNTR 2.0.3 even though that release has no `--version` flag; unrelated diagnostic
    version tokens are never guessed to be the answer.
    """
    if not text:
        return None
    return _parse_probe_output(text, "")


def detect_advntr_version(config: dict[str, Any], *, probe: AdvntrVersionProbe | None = None) -> AdvntrVersionOutcome:
    """Return the typed terminal result of the configured adVNTR version probe.

    Args:
        config: Complete pipeline configuration containing ``tools.advntr``.
        probe: Run-scoped probe to share a successful answer. A standalone call gets
            an isolated probe and therefore no process-global cache.

    Returns:
        A typed verified, versioned-launch-failure, unversioned-launch-failure,
        unparseable-success, or exhausted outcome.
    """
    command = config.get("tools", {}).get("advntr")
    if not command:
        return AdvntrVersionOutcome(
            AdvntrProbeStatus.LAUNCH_FAILURE,
            message="adVNTR version launch failed: no command is configured.",
        )
    active_probe = probe if probe is not None else AdvntrVersionProbe()
    return active_probe.detect(command)


def require_verified_advntr_version(outcome: AdvntrVersionOutcome) -> tuple[int, int, int]:
    """Return a verified version or raise with the probe's class-specific terminal message.

    Args:
        outcome: Typed terminal result from :func:`detect_advntr_version`.

    Returns:
        The verified three-part adVNTR version.

    Raises:
        ValueError: If the outcome itself violates its discriminant invariant.
        RuntimeError: If the command did not produce a verified version.
    """
    validated = AdvntrVersionOutcome(outcome.status, version=outcome.version, message=outcome.message)
    if validated.status is AdvntrProbeStatus.VERIFIED:
        return cast(tuple[int, int, int], validated.version)
    logger.error(validated.message)
    raise RuntimeError(validated.message)


def describe_model(db_path: str | Path) -> dict[str, Any]:
    """Return the provenance record for a model database.

    Deliberately reports `window_bp` -- the interval adVNTR will actually fetch -- and
    never a bare sum of segment lengths presented as a genomic span. Conflating those
    two quantities is the defect itself.
    """
    path = Path(db_path)
    if not path.is_file():
        raise AdvntrModelError(f"adVNTR model not found: {path}. {_REINSTALL}")

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    connection = sqlite3.connect(str(path))
    try:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        # Statements are selected from a fixed map rather than built by interpolating a
        # name into SQL. SQLite cannot parameterise an identifier, so the only safe form
        # is a literal the caller never influences.
        if "vntrs_v2" in tables:
            table, schema_version = "vntrs_v2", "v2"
        elif "vntrs" in tables:
            table, schema_version = "vntrs", "v1"
        else:
            raise AdvntrModelError(
                f"{path} contains no adVNTR model table (expected `vntrs_v2` or `vntrs`). {_REINSTALL}"
            )

        columns = {row[1] for row in connection.execute(_PRAGMA[table])}
        rows = list(connection.execute(_SELECT[table]["ref_end" in columns]))
    finally:
        connection.close()

    if not rows:
        raise AdvntrModelError(f"{path} holds no VNTR rows. {_REINSTALL}")
    if len(rows) == 1:
        row = rows[0]
    else:
        # A general adVNTR database carries thousands of VNTRs. Select MUC1 rather than
        # assuming the file is the single-row bundle.
        muc1 = [candidate for candidate in rows if candidate[0] == MUC1_VID]
        if len(muc1) != 1:
            raise AdvntrModelError(
                f"{path} holds {len(rows)} VNTR rows and {len(muc1)} with VID "
                f"{MUC1_VID}; exactly one MUC1 model is required. {_REINSTALL}"
            )
        row = muc1[0]
    vid, chromosome, ref_start, repeats = row[0], row[1], int(row[2]), row[3]
    ref_end = int(row[4]) if len(row) > 4 and row[4] is not None else None
    segments = repeats.split(",")

    if ref_end is None:
        # No recorded end: adVNTR falls back to start + sum(len(units)), which is the
        # window it will really use. Report that, so the truncation is visible.
        ref_end = ref_start + sum(len(segment) for segment in segments)
    if ref_end <= ref_start:
        raise AdvntrModelError(f"{path} records ref_end={ref_end} at or before ref_start={ref_start}. {_REINSTALL}")

    return {
        "path": str(path),
        "sha256": digest,
        "schema_version": schema_version,
        "vid": vid,
        "genomic_interval": f"{chromosome}:{ref_start}-{ref_end}",
        "window_bp": ref_end - ref_start,
        "n_segments": len(segments),
        "n_distinct_segments": len(set(segments)),
        "max_segment_len": max(len(segment) for segment in segments),
    }


def require_compatible_advntr(model: dict[str, Any], advntr_version: tuple[int, int, int] | None) -> None:
    """Raise unless this model can be trusted with this adVNTR.

    Refuses two combinations, both of which would otherwise produce wrong results in
    silence: a v2 model read by an adVNTR that ignores its recorded end, and a v1
    model, which carries the truncated hg38 window regardless of the binary.
    """
    if model["schema_version"] == "v1":
        raise AdvntrModelError(
            f"{model['path']} is a v1 adVNTR model. On GRCh38 it describes "
            f"{model['window_bp']} bp of a 3,525 bp repeat array, so adVNTR would see "
            f"part of the locus and miss variants outside it "
            f"(hassansaei/VNtyper#268). {_REINSTALL}"
        )

    if advntr_version is None or advntr_version < SPAN_AWARE_ADVNTR:
        found = ".".join(str(part) for part in advntr_version) if advntr_version else "unknown"
        wanted = ".".join(str(part) for part in SPAN_AWARE_ADVNTR)
        raise AdvntrModelError(
            f"adVNTR {found} cannot read the recorded genomic end in "
            f"{model['path']}; it would fetch a truncated window and report "
            f"confident but incomplete results. Install adVNTR >= {wanted}."
        )


def require_compatible_advntr_outcome(model: dict[str, Any], outcome: AdvntrVersionOutcome) -> tuple[int, int, int]:
    """Return a verified compatible version, preserving a failed launch's diagnosis.

    A pre-2.0.4 adVNTR has no ``--version`` option: argparse exits nonzero but includes
    one tagged version banner. That remains a launch failure. Its version may only be
    used to produce the more actionable incompatibility refusal; a compatible banner
    on a failed subprocess never authorizes a run.

    Args:
        model: Parsed adVNTR model provenance.
        outcome: Typed version-probe outcome.

    Returns:
        The verified compatible version from a successful subprocess.

    Raises:
        AdvntrModelError: If the detected version and model are incompatible.
        RuntimeError: If the subprocess did not complete successfully.
        ValueError: If the outcome violates its typed invariant.
    """
    validated = AdvntrVersionOutcome(outcome.status, version=outcome.version, message=outcome.message)
    if validated.status is AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE:
        version = cast(tuple[int, int, int], validated.version)
        require_compatible_advntr(model, version)
        logger.error(validated.message)
        raise RuntimeError(validated.message)

    version = require_verified_advntr_version(validated)
    require_compatible_advntr(model, version)
    return version
