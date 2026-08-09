"""Pure data and diagnostics contracts for alignment preflight.

The preflight stage decides which index names and reference candidates are
usable before any alignment conversion.  This module contains only those
deterministic decisions and the data passed between the stage boundaries.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from vntyper.scripts.alignment_binding import AlignmentBinding

logger = logging.getLogger(__name__)

FORMAT_BAM = "bam"
FORMAT_CRAM = "cram"
INDEX_SUFFIXES: dict[str, tuple[str, ...]] = {
    FORMAT_BAM: ("bai", "csi"),
    FORMAT_CRAM: ("crai", "csi"),
}

ReferenceAttempt = tuple[str, str | None, str]


def missing_reference_contig(diagnostic: str, known_contigs: Iterable[str]) -> str | None:
    """Extract a known contig from pinned samtools reference failures.

    Args:
        diagnostic: Combined samtools output from a failed CRAM decode.
        known_contigs: Contig names declared by the CRAM header.

    Returns:
        The referenced header contig, or ``None`` when the diagnostic shape is
        unrecognized or names a contig outside the header.
    """
    known = tuple(known_contigs)
    for pattern in (r"Unable to fetch reference ([^\s:]+)", r"ref '([^']+)' not present"):
        match = re.search(pattern, diagnostic)
        if match is not None and match.group(1) in known:
            return match.group(1)
    if "MD5 checksum reference mismatch" in diagnostic:
        for contig in known:
            if re.search(rf"(?:SN:|reference\s+){re.escape(contig)}(?:\s|\t|:|$)", diagnostic):
                return contig
    return None


def index_candidate_names(in_path: str, file_format: str, *, bai_only: bool = False) -> tuple[str, ...]:
    """Return index paths in their conventional and stem-based spellings.

    Args:
        in_path: Alignment path whose index is sought.
        file_format: Either ``"bam"`` or ``"cram"``.
        bai_only: Restrict BAM candidates to BAI. CRAM still offers CRAI then CSI.

    Returns:
        Candidate index paths, with appended-file spelling before stem spelling.

    Raises:
        ValueError: If ``file_format`` is not a supported alignment format.
    """
    try:
        suffixes = INDEX_SUFFIXES[file_format]
    except KeyError:
        message = f"unknown alignment format: {file_format}"
        logger.error(message)
        raise ValueError(message) from None

    if bai_only and file_format == FORMAT_BAM:
        suffixes = ("bai",)

    alignment = Path(in_path)
    return tuple(
        candidate
        for suffix in suffixes
        for candidate in (f"{alignment}.{suffix}", str(alignment.with_suffix(f".{suffix}")))
    )


@dataclass(frozen=True)
class AlignmentPlan:
    """Resolved paths and decisions needed by downstream alignment stages."""

    input_path: str
    view_path: str
    file_format: str
    index_path: str
    reference_path: str | None
    reference_source: str
    uncovered_contigs: tuple[str, ...]
    unmapped_scan: str
    binding: AlignmentBinding | None = field(default=None, repr=False, compare=False)

    def close(self) -> None:
        """Release the run-local alignment descriptor after its final consumer."""
        if self.binding is not None:
            self.binding.close()


def missing_index_message(in_path: str, file_format: str, tried: Iterable[str]) -> str:
    """Build an actionable message when a fresh run-local index cannot be built.

    Args:
        in_path: Alignment path whose index could not be resolved.
        file_format: Alignment format used while looking for an index.
        tried: Supplied index names inspected and protected by preflight.

    Returns:
        Actionable diagnostic naming the input, candidates, and indexing command.
    """
    candidates = ", ".join(str(path) for path in tried)
    return (
        f"Unable to build a fresh run-local {file_format.upper()} index for {in_path}. "
        f"Supplied index candidates inspected: {candidates}. Verify the alignment is readable and "
        f"samtools index {in_path} succeeds."
    )


def unresolvable_reference_message(
    in_path: str,
    contig: str,
    m5: str | None,
    attempts: Iterable[ReferenceAttempt],
) -> str:
    """Build diagnostics for a CRAM contig whose reference cannot be resolved.

    Args:
        in_path: CRAM path whose reference could not be resolved.
        contig: Contig requiring a reference sequence.
        m5: Reference sequence checksum from the alignment header, if available.
        attempts: Candidate reference sources, paths, and resolution reasons.

    Returns:
        Diagnostic containing the contig identity and every candidate attempt.
    """
    checksum = m5 if m5 is not None else "no M5"
    details = "; ".join(
        f"source={source}, path={path if path is not None else 'no path'}, reason={reason}"
        for source, path, reason in attempts
    )
    return f"Unable to resolve reference for {in_path}: contig={contig}, M5={checksum}. Candidates: {details}"


def preflight_error_payload(code: str, message: str, attempts: Iterable[ReferenceAttempt]) -> dict:
    """Return the serializable contents of ``preflight_error.json``.

    Args:
        code: Stable machine-readable preflight error code.
        message: Human-readable diagnostic for the failed preflight.
        attempts: Candidate reference sources, paths, and resolution reasons.

    Returns:
        Payload containing only the error code, message, and candidate attempts.
    """
    return {"code": code, "message": message, "candidates": list(attempts)}
