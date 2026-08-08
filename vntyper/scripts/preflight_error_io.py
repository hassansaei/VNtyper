"""Safe persistence of curated alignment-preflight failures."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from contextlib import suppress
from pathlib import Path

from vntyper.scripts.alignment_contract import ReferenceAttempt, preflight_error_payload, unresolvable_reference_message

logger = logging.getLogger(__name__)

PREFLIGHT_ERROR_FILENAME = "preflight_error.json"
_PAYLOAD_KEYS = {"code", "message", "candidates"}


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
    message = unresolvable_reference_message("CRAM input", contig, m5, public_attempts)
    return preflight_error_payload("reference_unresolved", message, public_attempts)


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

    final_path = Path(output_dir) / PREFLIGHT_ERROR_FILENAME
    descriptor, temporary_name = tempfile.mkstemp(
        dir=final_path.parent,
        prefix=f".{final_path.name}.",
        suffix=".tmp",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as temporary_file:
            descriptor = -1
            json.dump(payload, temporary_file, separators=(",", ":"))
            temporary_file.write("\n")
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, final_path)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        with suppress(OSError):
            os.unlink(temporary_path)
    return final_path
