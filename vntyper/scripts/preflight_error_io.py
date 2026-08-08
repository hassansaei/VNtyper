"""Safe persistence of curated alignment-preflight failures."""

from __future__ import annotations

import json
import logging
import os
import re
import secrets
import stat
from collections.abc import Iterator
from contextlib import contextmanager, suppress
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.alignment_contract import ReferenceAttempt, preflight_error_payload, unresolvable_reference_message

logger = logging.getLogger(__name__)

PREFLIGHT_ERROR_FILENAME = "preflight_error.json"
_PAYLOAD_KEYS = {"code", "message", "candidates"}
_PLACED_UNMAPPED = re.compile(r"idxstats reports (\d+) placed-unmapped reads")


@dataclass
class PreflightErrorContext:
    """Mutable diagnostics retained while the preflight advances through stages."""

    output_dir: str | Path
    payload: dict | None = None


def _public_attempt(attempt: ReferenceAttempt) -> ReferenceAttempt:
    source, path, reason = attempt
    public_path = path.replace("\\", "/").rsplit("/", 1)[-1] if path is not None else None
    if reason in {"not supplied", "reference FASTA not found"}:
        public_reason = reason
    elif reason.startswith("reference FASTA unreadable"):
        public_reason = "reference FASTA unreadable"
    else:
        public_reason = "probe exited non-zero"
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


def public_preflight_error_payload(error: Exception) -> dict:
    """Build a path-free public payload for a non-reference preflight failure.

    Args:
        error: Original exception retained for CLI diagnostics and server logs.

    Returns:
        The exact public contract with an empty candidate list.
    """
    detail = str(error)
    if "reference candidate" in detail:
        code = "reference_policy_invalid"
        if "exactly one terminal" in detail:
            message = "CRAM reference candidate policy must contain exactly one terminal htslib-resolved candidate."
        elif "must end with terminal" in detail:
            message = "CRAM reference candidate policy must end with the terminal htslib-resolved candidate."
        elif "duplicate explicit" in detail:
            message = "CRAM reference candidate policy contains a duplicate explicit source."
        elif "unknown reference candidate source" in detail:
            message = "CRAM reference candidate policy contains an unknown source."
        else:
            message = "CRAM reference candidate policy must be a list."
    elif "cram.local_ref_path" in detail:
        code = "reference_policy_invalid"
        message = (
            "CRAM reference policy rejected a remote REF_PATH; configure a local cache or explicitly allow "
            "ambient resolution."
        )
    elif "output_name" in detail:
        code = "preflight_output_unsafe"
        message = "Alignment preflight rejected an unsafe output name; use a single basename."
    elif detail.startswith("No usable BAM index"):
        code = "alignment_index_unusable"
        message = "No usable BAM index found. Create one with samtools index before retrying."
    elif detail.startswith("No usable CRAM index"):
        code = "alignment_index_unusable"
        message = "No usable CRAM index found. Create one with samtools index before retrying."
    elif "unknown alignment format" in detail:
        code = "alignment_format_invalid"
        message = "Alignment format is not supported by preflight."
    elif (placed := _PLACED_UNMAPPED.search(detail)) is not None:
        code = "unmapped_scan_lossy"
        message = f"Indexed CRAM scanning would omit {placed.group(1)} placed-unmapped reads; use auto or stream mode."
    elif "invalid unmapped scan mode" in detail:
        code = "unmapped_scan_invalid"
        message = "CRAM unmapped-read scan policy is invalid; use auto, indexed, or stream."
    elif detail.startswith("BAM preflight probe failed"):
        code = "alignment_probe_failed"
        message = "BAM preflight could not retrieve the requested target; verify the index and target coordinates."
    elif detail.startswith("Alignment view"):
        code = "preflight_output_unsafe"
        message = "Alignment preflight rejected an unsafe alignment view; choose a safe output directory."
    elif detail.startswith("Log path") or "command log" in detail:
        code = "preflight_output_unsafe"
        message = "Alignment preflight rejected an unsafe log destination; remove the conflicting entry and retry."
    elif detail.startswith(("Generated index", "Generated-index")) or "stale view index" in detail:
        code = "preflight_output_unsafe"
        message = (
            "Alignment preflight rejected an unsafe index destination or ownership record; "
            "remove the conflicting entry and retry."
        )
    elif detail.startswith("Derived"):
        code = "preflight_output_unsafe"
        message = (
            "Alignment preflight rejected an unsafe alignment view or index destination; "
            "remove the conflicting entry and retry."
        )
    elif "output directory" in detail:
        code = "preflight_output_unsafe"
        message = "Alignment preflight rejected the output directory; use a real directory outside the input tree."
    else:
        code = "alignment_preflight_failed"
        message = "Alignment preflight failed before processing; inspect the server logs for the job."
    return preflight_error_payload(code, message, ())


@contextmanager
def persist_preflight_failure(context: PreflightErrorContext) -> Iterator[None]:
    """Persist one curated artifact and re-raise the original exception unchanged.

    Args:
        context: Current stage and any structured reference-failure payload.

    Yields:
        Control to the complete preflight operation.
    """
    try:
        yield
    except Exception as error:
        payload = context.payload or public_preflight_error_payload(error)
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
        directory_descriptor = os.open(output, directory_flags)
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
