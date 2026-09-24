"""Fold-local axis inventories: the thresholds each outer fold may select, from its training samples only.

A data-derived axis is a function of the samples it is derived from: its breakpoints, its
endpoint sentinel and, under ``--max-breakpoints``, which breakpoints the rank-space cap
retains. Deriving the axis once from every sample and then restricting an outer fold to it
leaks the held-out samples into that fold's candidate set in two ways (spec §14.1):

* **The sentinel.** The reject-everything endpoint is computed from the extreme observed
  value. When a held-out sample holds that extreme, it chooses its own fold's sentinel --
  a threshold placed exactly at its own statistic.
* **The cap.** Rank-space subsampling depends on every value's rank, so a held-out value
  can shift the retained ranks and displace a breakpoint the training samples produced.

:func:`fold_axis` removes both by calling the same derive function once per outer fold on
that fold's training samples alone, so each fold computes its own breakpoints, sentinel
and cap. Every value admissible in fold f is then a function of fold f's training
statistics only. The replayed inventory is the union of the full-data axis and every fold
axis (:func:`~vntyper.scripts.calibration_cutoff_axes.merge_axis_values`), and
:func:`fold_candidate_inventories` maps each fold's values to the candidate IDs of that
union, which is what ``evaluate_cutoff_arms(fold_inventories=...)`` admits per fold.

The module is caller-agnostic: the Kestrel axes pass per-sample observed-value sequences,
the adVNTR axes per-sample decisive statistics.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from typing import NoReturn, TypeVar

from vntyper.scripts.calibration_cutoff_axes import AxisBreakpoints, merge_axis_values
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate

logger = logging.getLogger(__name__)

T = TypeVar("T")


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _log_refusals(axis: str, refused: Mapping[tuple[float, str], None]) -> None:
    """Report each dropped breakpoint once, whichever inventories dropped it; it is routine, not an error."""
    for value, reason in refused:
        logger.info("cutoff axis %s breakpoint %r refused: %s", axis, value, reason)


def fold_axis(
    derive: Callable[[Mapping[str, T]], AxisBreakpoints],
    per_sample: Mapping[str, T],
    assignments: Mapping[str, int],
) -> tuple[AxisBreakpoints, dict[int, frozenset[int | float]]]:
    """Derive one axis from all samples and, per outer fold, from that fold's training samples.

    Args:
        derive: Builds the axis from a per-sample mapping; it carries the axis name, the
            baseline and the cap, so every fold is derived by exactly the full-data rule.
        per_sample: Per-sample statistics. A sample without a statistic may be absent.
        assignments: Specimen key to outer fold, from ``outer_fold_assignments``. Empty
            when there are too few groups to cross-validate.

    Returns:
        The full-data axis merged with every fold's values, with ``fold_capped`` naming
        the folds whose own inventory was capped, and each fold's value set. Without
        folds, the full-data axis unchanged and ``{}``.

    Raises:
        ValueError: If ``per_sample`` names a sample that has no outer fold, or the derive
            function refuses the data.
    """
    full = derive(per_sample)
    refused = dict.fromkeys(full.rejected)
    if not assignments:
        _log_refusals(full.axis, refused)
        return full, {}
    unassigned = sorted(key for key in per_sample if key not in assignments)
    if unassigned:
        _fail(f"cutoff fold axis {full.axis} has samples without an outer fold: {unassigned}")
    inventories: dict[int, frozenset[int | float]] = {}
    capped: list[int] = []
    for fold in sorted(set(assignments.values())):
        training = {key: value for key, value in per_sample.items() if assignments[key] != fold}
        axis = derive(training)
        refused.update(dict.fromkeys(axis.rejected))
        inventories[fold] = frozenset(axis.values)
        if axis.capped:
            capped.append(fold)
    _log_refusals(full.axis, refused)
    extra = set().union(*inventories.values())
    return replace(merge_axis_values(full, extra), fold_capped=tuple(capped)), inventories


def fold_candidate_inventories(
    derived: Sequence[tuple[AxisBreakpoints, Sequence[CutoffCandidate]]],
    fold_values: Mapping[str, Mapping[int, frozenset[int | float]]],
) -> dict[int, frozenset[str]]:
    """Map every fold's training-derived values to candidate IDs, across all axes.

    Args:
        derived: Each merged axis with its candidates, in value order (from
            ``axis_candidates``).
        fold_values: Axis name to the per-fold value sets :func:`fold_axis` returned.

    Returns:
        Outer fold to the IDs of every candidate its training samples derived, on any
        axis; ``{}`` when no axis has folds.

    Raises:
        ValueError: If the axes of ``fold_values`` differ from ``derived``, the axes
            disagree on the folds, or a fold value has no candidate on its merged axis.
    """
    names = [axis.axis for axis, _ in derived]
    if len(set(names)) != len(names) or set(names) != set(fold_values):
        _fail(f"cutoff fold inventories cover axes {sorted(fold_values)}, not the derived axes {sorted(names)}")
    folds = {frozenset(values) for values in fold_values.values()}
    if len(folds) > 1:
        _fail("cutoff fold inventories disagree on the outer folds across axes")
    inventories: dict[int, set[str]] = {fold: set() for fold in next(iter(folds), frozenset())}
    for axis, candidates in derived:
        if len(candidates) != len(axis.values):
            _fail(f"cutoff fold axis {axis.axis} has {len(candidates)} candidates for {len(axis.values)} values")
        ids = {value: candidate.candidate_id for value, candidate in zip(axis.values, candidates, strict=True)}
        for fold, values in fold_values[axis.axis].items():
            missing = sorted(value for value in values if value not in ids)
            if missing:
                _fail(f"cutoff fold {fold} value(s) {missing} of axis {axis.axis} have no candidate")
            inventories[fold].update(ids[value] for value in values)
    return {fold: frozenset(ids) for fold, ids in sorted(inventories.items())}
