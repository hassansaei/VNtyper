"""Safe persistence of curated alignment-preflight failures."""

from __future__ import annotations

import json
import logging
import os
import secrets
import stat
from collections.abc import Iterator
from contextlib import contextmanager, suppress
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.alignment_contract import (
    ReferenceAttempt,
    TimedOutProbeReason,
    preflight_error_payload,
    unresolvable_reference_message,
)

logger = logging.getLogger(__name__)

PREFLIGHT_ERROR_FILENAME = "preflight_error.json"
_PAYLOAD_KEYS = {"code", "message", "candidates"}


@dataclass(frozen=True)
class PreflightFailurePhase:
    """Stable public failure semantics for one preflight operation."""

    code: str
    message: str


UNCLASSIFIED_FAILURE = PreflightFailurePhase(
    "alignment_preflight_failed",
    "Alignment preflight failed before processing; inspect the server logs for the job.",
)
REFERENCE_POLICY_FAILURE = PreflightFailurePhase(
    "reference_policy_invalid",
    "CRAM reference candidate policy is invalid; configure a list ending in exactly one terminal "
    "htslib-resolved candidate.",
)
OUTPUT_SAFETY_FAILURE = PreflightFailurePhase(
    "preflight_output_unsafe",
    "Alignment preflight rejected an unsafe output or log destination; use a dedicated output directory "
    "and remove conflicting entries.",
)
VIEW_INDEX_FAILURE = PreflightFailurePhase(
    "alignment_index_unusable",
    "Alignment preflight could not prepare a safe local view and fresh index; remove conflicting run-output "
    "entries and verify samtools can index the alignment.",
)
SCAN_SELECTION_FAILURE = PreflightFailurePhase(
    "unmapped_scan_invalid",
    "Alignment unmapped-read scan selection failed; use auto or stream mode to avoid losing placed-unmapped reads.",
)
REFERENCE_PROBE_FAILURE = PreflightFailurePhase(
    "reference_unresolved",
    "CRAM reference preflight could not prove the requested target; verify the reference FASTA and target contigs.",
)
BAM_PROBE_FAILURE = PreflightFailurePhase(
    "alignment_probe_failed",
    "BAM preflight could not retrieve the requested target; verify the index and target coordinates.",
)
HEADER_PREPARATION_FAILURE = PreflightFailurePhase(
    "alignment_header_invalid",
    "Alignment preflight rejected the alignment header or declared assembly; verify the input and "
    "--reference-assembly.",
)
TARGET_PREPARATION_FAILURE = PreflightFailurePhase(
    "alignment_target_invalid",
    "Alignment preflight could not prepare the requested target; verify the BED file or configured regions.",
)


@dataclass
class PreflightErrorContext:
    """Mutable diagnostics retained while the preflight advances through stages."""

    output_dir: str | Path
    phase: PreflightFailurePhase = UNCLASSIFIED_FAILURE
    payload: dict | None = None


def _public_attempt(attempt: ReferenceAttempt) -> ReferenceAttempt:
    source, path, reason = attempt
    public_path = path.replace("\\", "/").rsplit("/", 1)[-1] if path is not None else None
    if reason in {"not supplied", "reference FASTA not found"}:
        public_reason = reason
    elif reason.startswith("reference FASTA unreadable"):
        public_reason = "reference FASTA unreadable"
    else:
        # Only the command runner can construct this subtype. The tool's own
        # output may imitate the timeout sentence, so parsing text is unsafe.
        public_reason = (
            f"probe timed out after {reason.timeout_seconds:g} seconds"
            if isinstance(reason, TimedOutProbeReason)
            else "probe exited non-zero"
        )
    return source, public_path, public_reason


def _public_identifier(value: str, fallback: str = "unknown") -> str:
    """Retain a short printable identifier while rejecting path-like content."""
    if not value or len(value) > 200 or "/" in value or "\\" in value or not value.isprintable():
        return fallback
    return value


def public_reference_error_payload(
    contig: str,
    m5: str | None,
    attempts: tuple[ReferenceAttempt, ...] | list[ReferenceAttempt],
) -> dict:
    """Build the public reference-failure artifact without worker paths.

    Args:
        contig: Target contig that could not be decoded.
        m5: CRAM header checksum for that contig, when present.
        attempts: Full internal candidate diagnostics.

    Returns:
        The exact three-field error contract with candidate basenames and
        controlled reasons instead of raw worker paths or samtools output.
    """
    public_attempts = tuple(_public_attempt(attempt) for attempt in attempts)
    public_contig = _public_identifier(contig)
    public_m5 = None if m5 is None else _public_identifier(m5)
    message = unresolvable_reference_message("CRAM input", public_contig, public_m5, public_attempts)
    return preflight_error_payload("reference_unresolved", message, public_attempts)


def public_preflight_error_payload(phase: PreflightFailurePhase) -> dict:
    """Build a public payload from stable phase metadata, never exception prose.

    Args:
        phase: Explicit operation active when preflight failed.

    Returns:
        The exact public contract with an empty candidate list.
    """
    return preflight_error_payload(phase.code, phase.message, ())


def _clear_preflight_error(output_dir: str | Path) -> None:
    """Remove a prior run's artifact through an unfollowed directory handle."""
    directory_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    directory_flags |= getattr(os, "O_DIRECTORY", 0)
    try:
        directory_descriptor = os.open(output_dir, directory_flags)
    except FileNotFoundError:
        return
    try:
        with suppress(FileNotFoundError):
            os.unlink(PREFLIGHT_ERROR_FILENAME, dir_fd=directory_descriptor)
    finally:
        os.close(directory_descriptor)


@contextmanager
def persist_preflight_failure(context: PreflightErrorContext) -> Iterator[None]:
    """Persist one curated artifact and re-raise the original exception unchanged.

    Args:
        context: Current stage and any structured reference-failure payload.

    Yields:
        Control to the complete preflight operation.
    """
    try:
        _clear_preflight_error(context.output_dir)
        yield
    except Exception:
        payload = context.payload or public_preflight_error_payload(context.phase)
        try:
            write_preflight_error(context.output_dir, payload)
        except Exception as artifact_error:
            logger.error(f"Unable to write preflight error artifact safely: {artifact_error}")
        raise


def write_preflight_error(output_dir: str | Path, payload: dict) -> Path:
    """Atomically install ``preflight_error.json`` without following its entry.

    Args:
        output_dir: Existing run output directory that owns the artifact.
        payload: Curated three-field preflight error payload.

    Returns:
        The installed artifact path.

    Raises:
        OSError: If the temporary write or atomic replacement fails.
        ValueError: If the payload does not have the exact public contract.
    """
    if set(payload) != _PAYLOAD_KEYS:
        message = "preflight error payload must contain exactly code, message, candidates"
        logger.error(message)
        raise ValueError(message)

    output = Path(output_dir)
    final_path = output / PREFLIGHT_ERROR_FILENAME
    directory_flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    directory_flags |= getattr(os, "O_DIRECTORY", 0)
    try:
        # Preserve a descriptor-bound spelling such as
        # ``/proc/<pid>/fd/<descriptor>/.``. Converting it to ``Path`` first
        # removes the trailing ``/.`` and turns the final component back into
        # a procfs symlink, which O_NOFOLLOW correctly refuses.
        directory_descriptor = os.open(output_dir, directory_flags)
    except OSError as error:
        message = "preflight error output must be an existing real directory"
        logger.error(message)
        raise ValueError(message) from error

    descriptor = -1
    temporary_name: str | None = None
    try:
        try:
            metadata = os.stat(PREFLIGHT_ERROR_FILENAME, dir_fd=directory_descriptor, follow_symlinks=False)
        except FileNotFoundError:
            metadata = None
        if (
            metadata is not None
            and not stat.S_ISLNK(metadata.st_mode)
            and (not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1)
        ):
            message = "preflight error destination must be a safe final entry"
            logger.error(message)
            raise ValueError(message)

        for _ in range(100):
            candidate = f".{PREFLIGHT_ERROR_FILENAME}.{secrets.token_hex(8)}.tmp"
            try:
                descriptor = os.open(
                    candidate,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0),
                    0o600,
                    dir_fd=directory_descriptor,
                )
            except FileExistsError:
                continue
            temporary_name = candidate
            break
        if temporary_name is None:
            raise OSError("unable to reserve a preflight error temporary file")

        with os.fdopen(descriptor, "w", encoding="utf-8") as temporary_file:
            descriptor = -1
            json.dump(payload, temporary_file, separators=(",", ":"))
            temporary_file.write("\n")
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(
            temporary_name,
            PREFLIGHT_ERROR_FILENAME,
            src_dir_fd=directory_descriptor,
            dst_dir_fd=directory_descriptor,
        )
        temporary_name = None
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary_name is not None:
            with suppress(OSError):
                os.unlink(temporary_name, dir_fd=directory_descriptor)
        os.close(directory_descriptor)
    return final_path
