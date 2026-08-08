"""Pure parsing and scan selection for ``samtools idxstats`` output."""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

SCAN_INDEXED = "indexed"
SCAN_STREAM = "stream"


def parse_idxstats(text: str) -> tuple[int, int] | None:
    """Parse an idxstats table, failing closed when its structure is invalid.

    Args:
        text: Complete tab-separated output from ``samtools idxstats``.

    Returns:
        A ``(placed_unmapped, unplaced)`` pair, or ``None`` for malformed output.
    """
    lines = text.splitlines()
    if not lines:
        return None

    placed_unmapped = 0
    unplaced = 0
    star_rows = 0
    for line_number, line in enumerate(lines, start=1):
        fields = line.split("\t")
        if len(fields) != 4:
            return None

        try:
            counts = tuple(int(field) for field in fields[1:])
        except ValueError:
            return None
        if any(count < 0 for count in counts):
            return None

        if fields[0] == "*":
            star_rows += 1
            if line_number != len(lines):
                return None
            unplaced = counts[2]
        else:
            placed_unmapped += counts[2]

    if star_rows != 1:
        return None
    return placed_unmapped, unplaced


def choose_scan(configured: str, idxstats_text: str | None, exit_ok: bool) -> tuple[str, str]:
    """Choose a lossless unmapped-read scan from idxstats evidence.

    Args:
        configured: Requested scan mode (``"auto"``, ``"indexed"`` or ``"stream"``).
        idxstats_text: Captured idxstats output, or ``None`` when unavailable.
        exit_ok: Whether the idxstats command exited successfully.

    Returns:
        The selected scan mode and a human-readable reason.

    Raises:
        ValueError: If indexed scanning was forced despite placed unmapped reads.
    """
    if configured == SCAN_STREAM:
        return SCAN_STREAM, "stream scan configured"
    if not exit_ok:
        return SCAN_STREAM, "idxstats failed; using lossless stream scan"

    counts = parse_idxstats(idxstats_text if idxstats_text is not None else "")
    if counts is None:
        return SCAN_STREAM, "idxstats output is malformed; using lossless stream scan"

    placed_unmapped, _ = counts
    if placed_unmapped == 0:
        return SCAN_INDEXED, "idxstats reports no placed-unmapped reads"

    reason = f"idxstats reports {placed_unmapped} placed-unmapped reads; using stream scan"
    if configured == SCAN_INDEXED:
        logger.error(reason)
        raise ValueError(reason)
    return SCAN_STREAM, reason
