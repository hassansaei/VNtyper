"""Failure-aware cleanup for pipeline-owned alignment plans."""

from __future__ import annotations

import logging

from vntyper.scripts.alignment_contract import AlignmentPlan

logger = logging.getLogger(__name__)


def close_alignment_plan(plan: AlignmentPlan | None, *, preserve_primary: bool) -> None:
    """Close an alignment plan without replacing an already-active primary outcome.

    Args:
        plan: Pipeline-owned alignment plan, when preflight reached that boundary.
        preserve_primary: Suppress and log cleanup failure because another outcome is active.

    Raises:
        RuntimeError: The plan's cleanup failure when no primary outcome is active.
    """
    if plan is None:
        return
    try:
        plan.close()
    except Exception:
        if not preserve_primary:
            raise
        logger.exception("Alignment plan cleanup failed while preserving the primary pipeline outcome.")
