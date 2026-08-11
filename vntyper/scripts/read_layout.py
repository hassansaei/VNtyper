"""Classify and route the FASTQs emitted by alignment conversion."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

logger = logging.getLogger(__name__)

_FASTQ_KEYS = ("r1", "r2", "other", "single")
_LAYOUTS = {"paired", "single", "mixed", "invalid", "empty"}


def classify_layout(r1: int, r2: int, other: int, single: int) -> str:
    """Classify four FASTQ record counts without tolerating stranded reads.

    Args:
        r1: Records emitted to the READ1 output.
        r2: Records emitted to the READ2 output.
        other: Records whose READ1 and READ2 flags are both unset.
        single: Singleton records emitted by ``samtools fastq -s``.

    Returns:
        ``"paired"``, ``"single"``, ``"mixed"``, ``"invalid"`` or ``"empty"``.

    Raises:
        ValueError: If any record count is negative.
    """
    counts = (r1, r2, other, single)
    if any(count < 0 for count in counts):
        raise ValueError("FASTQ record counts must be non-negative.")
    if not any(counts):
        return "empty"
    if r1 != r2:
        return "invalid"
    if r1 > 0 and r1 == r2 and other == 0 and single == 0:
        return "paired"
    if r1 == 0 and r2 == 0 and (other > 0) != (single > 0):
        return "single"
    return "mixed"


def route_fastqs(
    layout: str,
    paths: Mapping[str, str | Path],
    counts: Mapping[str, int],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return consumed paths and every non-empty path left stranded.

    Args:
        layout: Verdict returned by :func:`classify_layout`.
        paths: FASTQ paths keyed by ``r1``, ``r2``, ``other`` and ``single``.
        counts: Record counts under the same four keys.

    Returns:
        A pair ``(consumed, stranded_nonempty)``. Paths retain the fixed
        R1/R2/other/single order.

    Raises:
        ValueError: If ``layout`` is not a recognised verdict.
    """
    if layout not in _LAYOUTS:
        raise ValueError(f"Unknown FASTQ layout: {layout}")

    routed_keys: tuple[str, ...]
    if layout == "paired":
        routed_keys = ("r1", "r2")
    elif layout == "single":
        routed_keys = ("other",) if counts["other"] > 0 else ("single",)
    elif layout == "mixed":
        routed_keys = tuple(key for key in _FASTQ_KEYS if counts[key] > 0)
    else:
        routed_keys = ()

    consumed = tuple(str(paths[key]) for key in routed_keys if counts[key] > 0)
    resolved_consumed = tuple(Path(path).resolve(strict=False) for path in consumed)
    if len(set(resolved_consumed)) != len(resolved_consumed):
        raise ValueError("FASTQ routing selected duplicate filesystem paths.")
    routed = set(routed_keys)
    stranded = tuple(str(paths[key]) for key in _FASTQ_KEYS if counts[key] > 0 and key not in routed)
    return consumed, stranded
