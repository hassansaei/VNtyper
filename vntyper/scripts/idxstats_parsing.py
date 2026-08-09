"""Pure parsing and scan selection for ``samtools idxstats`` output."""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

SCAN_INDEXED = "indexed"
SCAN_STREAM = "stream"
_MAX_DIAGNOSTIC_EXCERPT_CHARACTERS = 160
_TRUNCATION_MARKER = "...<truncated>"


def _diagnostic_excerpt(line: str) -> str:
    """Render one escaped idxstats line within the operator-message bound."""
    rendered = repr(line)
    if len(rendered) <= _MAX_DIAGNOSTIC_EXCERPT_CHARACTERS:
        return rendered
    prefix_length = _MAX_DIAGNOSTIC_EXCERPT_CHARACTERS - len(_TRUNCATION_MARKER) - 1
    return f"{rendered[:prefix_length]}{_TRUNCATION_MARKER}{rendered[-1]}"


def _parse_idxstats_with_error(text: str) -> tuple[tuple[int, int] | None, str]:
    """Parse an idxstats table and retain the exact fail-closed diagnostic."""
    lines = text.splitlines()
    if not lines:
        return None, "idxstats output is empty"

    placed_unmapped = 0
    unplaced = 0
    star_rows = 0
    for line_number, line in enumerate(lines, start=1):
        fields = line.split("\t")
        if len(fields) != 4 or not all(field.isascii() and field.isdecimal() for field in fields[1:]):
            return None, f"idxstats output is malformed at line {line_number}: {_diagnostic_excerpt(line)}"

        try:
            counts = tuple(int(field) for field in fields[1:])
        except ValueError:
            return None, f"idxstats output is malformed at line {line_number}: {_diagnostic_excerpt(line)}"

        if fields[0] == "*":
            star_rows += 1
            if line_number != len(lines):
                return None, f"idxstats output is malformed at line {line_number}: {_diagnostic_excerpt(line)}"
            unplaced = counts[2]
        else:
            placed_unmapped += counts[2]

    if star_rows == 0:
        return None, "idxstats output is missing its terminal '*' row"
    return (placed_unmapped, unplaced), ""


def parse_idxstats(text: str) -> tuple[int, int] | None:
    """Parse an idxstats table, failing closed when its structure is invalid.

    Args:
        text: Complete tab-separated output from ``samtools idxstats``.

    Returns:
        A ``(placed_unmapped, unplaced)`` pair, or ``None`` for malformed output.
    """
    counts, _error = _parse_idxstats_with_error(text)
    return counts


def choose_scan(configured: str, idxstats_text: str | None, exit_ok: bool) -> tuple[str, str]:
    """Choose a lossless unmapped-read scan from idxstats evidence.

    Args:
        configured: Requested scan mode (``"auto"``, ``"indexed"`` or ``"stream"``).
        idxstats_text: Captured idxstats output, or ``None`` when unavailable.
        exit_ok: Whether the idxstats command exited successfully.

    Returns:
        The selected scan mode and a human-readable reason.

    Raises:
        ValueError: If the configured mode is invalid or indexed scanning was forced despite
            placed unmapped reads.
    """
    if configured not in {"auto", SCAN_INDEXED, SCAN_STREAM}:
        message = f"invalid unmapped scan mode: {configured}"
        logger.error(message)
        raise ValueError(message)
    if configured == SCAN_STREAM:
        return SCAN_STREAM, "stream scan configured"
    if not exit_ok:
        return SCAN_STREAM, "idxstats failed; using lossless stream scan"

    counts, parse_error = _parse_idxstats_with_error(idxstats_text if idxstats_text is not None else "")
    if counts is None:
        return SCAN_STREAM, f"{parse_error}; using lossless stream scan"

    placed_unmapped, _ = counts
    if placed_unmapped == 0:
        return SCAN_INDEXED, "idxstats reports no placed-unmapped reads"

    reason = f"idxstats reports {placed_unmapped} placed-unmapped reads; using stream scan"
    if configured == SCAN_INDEXED:
        logger.error(reason)
        raise ValueError(reason)
    return SCAN_STREAM, reason
