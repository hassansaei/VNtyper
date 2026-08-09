"""Pure parsing and scan selection for ``samtools idxstats`` output."""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

SCAN_INDEXED = "indexed"
SCAN_STREAM = "stream"
_MAX_DIAGNOSTIC_EXCERPT_CHARACTERS = 160
_MAX_COUNT_DIGITS = 20
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


def _parse_indexed_count(text: str) -> int | None:
    """Parse the single non-negative count emitted by ``samtools view -c``."""
    lines = text.splitlines()
    if len(lines) != 1 or len(lines[0]) > _MAX_COUNT_DIGITS or not lines[0].isascii() or not lines[0].isdecimal():
        return None
    try:
        return int(lines[0])
    except ValueError:
        return None


def _unsafe_indexed_result(configured: str, reason: str) -> tuple[str, str]:
    """Fail closed, raising when an operator explicitly forced the unsafe mode."""
    if configured == SCAN_INDEXED:
        logger.error(reason)
        raise ValueError(reason)
    return SCAN_STREAM, reason


def choose_scan(
    configured: str,
    idxstats_text: str | None,
    exit_ok: bool,
    *,
    indexed_count_text: str | None = None,
    indexed_count_exit_ok: bool = False,
) -> tuple[str, str]:
    """Choose a lossless unmapped-read scan from idxstats evidence.

    Args:
        configured: Requested scan mode (``"auto"``, ``"indexed"`` or ``"stream"``).
        idxstats_text: Captured idxstats output, or ``None`` when unavailable.
        exit_ok: Whether the idxstats command exited successfully.
        indexed_count_text: Captured count from the exact flag-4 literal-``'*'``
            indexed consumer, when idxstats found no placed-unmapped records.
        indexed_count_exit_ok: Whether that indexed count command exited successfully.

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

    placed_unmapped, unplaced = counts
    if placed_unmapped != 0:
        reason = f"idxstats reports {placed_unmapped} placed-unmapped reads; using stream scan"
        return _unsafe_indexed_result(configured, reason)

    if not indexed_count_exit_ok:
        reason = "indexed '*' count failed; using lossless stream scan"
        return _unsafe_indexed_result(configured, reason)
    rendered_count = indexed_count_text if indexed_count_text is not None else ""
    indexed_count = _parse_indexed_count(rendered_count)
    if indexed_count is None:
        reason = (
            f"indexed '*' count output is malformed: {_diagnostic_excerpt(rendered_count)}; using lossless stream scan"
        )
        return _unsafe_indexed_result(configured, reason)
    if indexed_count != unplaced:
        reason = (
            f"indexed '*' count {indexed_count} differs from idxstats unplaced count {unplaced}; "
            "using lossless stream scan"
        )
        return _unsafe_indexed_result(configured, reason)
    return SCAN_INDEXED, f"idxstats and indexed '*' count agree on {unplaced} unplaced reads"
