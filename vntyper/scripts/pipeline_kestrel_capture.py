"""Emit replayable Kestrel calibration evidence from a supported pipeline run.

Cutoff calibration replays the complete pre-filter candidate population rather than
re-running Kestrel once per threshold, so it needs that population preserved. A final
result table cannot stand in for it: by the time it is written, every candidate the
gates rejected is gone, and those are exactly the rows a lower cutoff would admit.

This module owns the run-side half of that contract. It builds the capture through
:func:`vntyper.scripts.calibration_kestrel_capture.build_kestrel_capture` and installs
it atomically, so a run either leaves one complete capture or none at all. A run whose
VCF held no indel still produced evidence -- a complete empty population -- and
:meth:`CaptureWriter.write_empty` records it, because treating the absent file as
"no evidence" would silently drop every true negative that had no candidate.
"""

from __future__ import annotations

import logging
import os
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, NoReturn

import pandas as pd

from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_kestrel_capture import build_kestrel_capture, kestrel_capture_document
from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.version import __version__

if TYPE_CHECKING:  # pragma: no cover - typing only
    from vntyper.scripts.identity_candidates import IdentityTranslationComponent
    from vntyper.scripts.kestrel_selection import KestrelSelection

logger = logging.getLogger(__name__)

CAPTURE_POLICY_SCHEMA = "vntyper-kestrel-recruitment-v1"

_RAW_COLUMNS = ("Motifs", "POS", "REF", "ALT", "Sample", "Motif_sequence", "Variant")


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


@dataclass(frozen=True)
class CaptureAssets:
    """Everything a capture must commit to besides the candidate rows.

    Attributes:
        reference_file: Kestrel VNTR reference FASTA actually used by the run.
        motif_reference_file: Motif annotation FASTA actually used by the run.
        kestrel_jar: Kestrel JAR actually invoked.
        baseline_policy: Complete caller policy the run's configuration represents.
        decision_profile_sha256: Digest of the run's resolved decision profile.
        capture_policy: Recruitment and assembly parameters that force recapture
            when they change.
    """

    reference_file: Path
    motif_reference_file: Path
    kestrel_jar: Path
    baseline_policy: CallerPolicyValues
    decision_profile_sha256: str
    capture_policy: dict[str, object]


def build_capture_policy(*, kmer_sizes: Sequence[int], threads: int, input_scope: str) -> dict[str, object]:
    """Record the recruitment parameters a replay may not silently differ from.

    Changing how reads are recruited changes which candidates exist at all, so a
    capture taken under one recruitment policy cannot be replayed under another. The
    replay compares this document's digest and demands recapture when it differs.

    Args:
        kmer_sizes: Distinct positive k-mer sizes Kestrel was run with, in order.
        threads: Positive thread budget the run used.
        input_scope: Non-empty description of the input the run consumed.

    Returns:
        A JSON-compatible capture policy document.

    Raises:
        ValueError: If any parameter is empty, non-positive or repeated.
    """
    sizes = list(kmer_sizes)
    if not sizes or len(set(sizes)) != len(sizes):
        _fail("Kestrel capture policy requires distinct non-empty k-mer sizes")
    if any(isinstance(size, bool) or not isinstance(size, int) or size <= 0 for size in sizes):
        _fail("Kestrel capture policy k-mer sizes must be positive integers")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0:
        _fail("Kestrel capture policy threads must be a positive integer")
    if not isinstance(input_scope, str) or not input_scope or input_scope != input_scope.strip():
        _fail("Kestrel capture policy input scope must be non-empty trimmed text")
    return {
        "schema_version": CAPTURE_POLICY_SCHEMA,
        "input_scope": input_scope,
        "kmer_sizes": sizes,
        "threads": threads,
        "vntyper_version": __version__,
    }


def _read(path: Path, label: str) -> bytes:
    try:
        return Path(path).read_bytes()
    except OSError as error:
        _fail(f"cannot read Kestrel capture {label} {path}: {error}")


class CaptureWriter:
    """Write exactly one complete Kestrel capture for one pipeline run."""

    def __init__(self, destination: Path, assets: CaptureAssets) -> None:
        """Bind a writer to its destination and committed assets.

        Args:
            destination: Final capture path; its parent is created on write.
            assets: Files and policy the capture commits to.
        """
        self._destination = Path(destination)
        self._assets = assets
        self._written = False

    def __call__(
        self,
        raw_frame: pd.DataFrame,
        motif_frame: pd.DataFrame,
        kestrel_config: dict[str, Any],
        selection: KestrelSelection,
        identity_component: IdentityTranslationComponent,
    ) -> None:
        """Capture the complete pre-scoring candidate population.

        Args:
            raw_frame: Exact ordered post-VCF, pre-scoring candidate frame.
            motif_frame: Exact ordered parsed motif annotation table.
            kestrel_config: Complete frozen production Kestrel component.
            selection: Frozen production filtering, ordering and dominance decisions.
            identity_component: Frozen production identity translation tables.

        Raises:
            ValueError: If a capture was already written for this run, an asset cannot
                be read, or the frames do not satisfy the capture contract.
        """
        self._write(raw_frame, motif_frame, kestrel_config, selection, identity_component)

    def write_empty(self, motif_frame: pd.DataFrame, kestrel_config: dict[str, Any]) -> None:
        """Capture a completed run whose candidate population was empty.

        A VCF with no indel never reaches the scoring boundary, so no observation
        arrives. The run still evaluated the sample, and a calibration cohort needs
        that recorded rather than inferred from a missing file.

        Args:
            motif_frame: Exact ordered parsed motif annotation table.
            kestrel_config: Complete frozen production Kestrel component.

        Raises:
            ValueError: If a capture was already written for this run.
        """
        from vntyper.scripts.identity_candidates import translation_component_from_config
        from vntyper.scripts.kestrel_genotyping import _resolve_selection
        from vntyper.scripts.nomenclature import nomenclature_config

        empty = pd.DataFrame({column: pd.Series(dtype=object) for column in _RAW_COLUMNS})
        self._write(
            empty,
            motif_frame,
            kestrel_config,
            _resolve_selection(kestrel_config),
            translation_component_from_config(nomenclature_config),
        )

    def _write(
        self,
        raw_frame: pd.DataFrame,
        motif_frame: pd.DataFrame,
        kestrel_config: dict[str, Any],
        selection: KestrelSelection,
        identity_component: IdentityTranslationComponent,
    ) -> None:
        if self._written:
            _fail("Kestrel calibration capture is written once per run")
        assets = self._assets
        capture = build_kestrel_capture(
            raw_frame,
            motif_frame,
            kestrel_config=kestrel_config,
            baseline_policy=assets.baseline_policy,
            selection=selection,
            identity_component=identity_component,
            decision_profile_sha256=assets.decision_profile_sha256,
            reference_file_bytes=_read(assets.reference_file, "reference"),
            motif_reference_file_bytes=_read(assets.motif_reference_file, "motif reference"),
            kestrel_jar_bytes=_read(assets.kestrel_jar, "Kestrel JAR"),
            capture_policy=assets.capture_policy,
        )
        payload = canonical_json_bytes(kestrel_capture_document(capture)) + b"\n"
        self._install(payload)
        self._written = True
        logger.info(
            "Wrote Kestrel calibration capture with %d candidate rows to %s",
            len(capture.rows),
            self._destination,
        )

    def _install(self, payload: bytes) -> None:
        destination = self._destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary: Path | None = None
        installed = False
        try:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".tmp",
                delete=False,
            ) as candidate:
                temporary = Path(candidate.name)
                os.fchmod(candidate.fileno(), 0o600)
                candidate.write(payload)
                candidate.flush()
                os.fsync(candidate.fileno())
            os.replace(temporary, destination)
            installed = True
        except OSError as error:
            _fail(f"cannot write Kestrel calibration capture {destination}: {error}")
        finally:
            if temporary is not None and not installed:
                temporary.unlink(missing_ok=True)


def kestrel_capture_writer(destination: Path, *, assets: CaptureAssets) -> CaptureWriter:
    """Build the one-shot capture writer for a run.

    Args:
        destination: Final capture path.
        assets: Files and policy the capture commits to.

    Returns:
        A writer usable directly as the Kestrel raw-capture observer.

    Raises:
        ValueError: If the assets are not a validated :class:`CaptureAssets`.
    """
    if not isinstance(assets, CaptureAssets):
        _fail("Kestrel capture writer requires validated CaptureAssets")
    return CaptureWriter(destination, assets)
