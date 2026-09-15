"""Assemble the ``calibration-cutoff-report-v1`` document from replayed cutoff evidence.

This module owns what the published record *says*; ``calibration_cutoff_optimize`` owns
the run that produces the evidence, and ``calibration_cutoff_report`` owns the bytes that
reach disk. Nothing here reads a file or replays anything: every number it publishes comes
from the counters and interval estimators that already computed it, so the report cannot
disagree with the evaluation by recomputing it differently.

Three of the sections need a word of explanation.

*Plateau.* A threshold sweep is a step function. Between two adjacent observed values every
threshold gives the identical partition of the cohort, so a selected value is a member of
an interval, not a measurement. :func:`_plateau` walks outward from the selected candidate
while the per-sample outcome vector is unchanged and publishes that run together with the
neighbouring values at which the outcome does change.

*Failed selection.* An objective whose constraints no candidate satisfies is an outcome,
not an error, and it is published with the constraint that could not be met and the best
value that was actually reachable -- otherwise the reader cannot tell "the cohort cannot
support this floor" from "the search was misconfigured".

*Joint points.* Rows from two axes are never concatenated into a curve. When more than one
axis was swept the report carries a labelled table beside the curves, ordered by label and
by nothing else, because no single parameter separates those rows.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Final, NoReturn

from vntyper.scripts.calibration_caller_metrics import CallerObservation, calculate_caller_metrics
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues
from vntyper.scripts.calibration_cohort_manifest import CohortSample
from vntyper.scripts.calibration_cohort_metrics import caller_metrics_document
from vntyper.scripts.calibration_cutoff_advntr import AdvntrCutoffGridResult
from vntyper.scripts.calibration_cutoff_axes import AxisBreakpoints, axis_document
from vntyper.scripts.calibration_cutoff_curves import (
    axis_curve_document,
    build_axis_curve,
    build_joint_points,
    joint_points_document,
)
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate
from vntyper.scripts.calibration_cutoff_kestrel import KestrelGridReplay
from vntyper.scripts.calibration_cutoff_report import SCHEMA_VERSION
from vntyper.scripts.calibration_cutoff_selection import SearchSpec, cutoff_counts, cutoff_counts_document

logger = logging.getLogger(__name__)

BASELINE_ID: Final[str] = "baseline"
PROFILE_NAME: Final[str] = "research-decision-profile.json"
USAGE_HINT: Final[str] = "vntyper pipeline --research-decision-profile <output>/" + PROFILE_NAME + " --bam <reads>"
LIMITATIONS: Final[str] = (
    "Research use only. These cutoffs were derived from a previously examined local cohort; they carry no "
    "independent validation and no deployment approval."
)
_DUPLICATE_REASON: Final[str] = "not-the-first-seen-representative-of-its-declared-group"
_PLATEAU_NOTE: Final[str] = (
    "A threshold sweep is a step function: every value listed here reproduces the selected outcome exactly, "
    "so quoting the selected number alone is false precision."
)
_NO_SELECTION: Final[str] = "no tested cutoff satisfied the declared objective and its constraints"

#: Per-axis candidates and their breakpoints, in ascending axis-value order.
DerivedAxis = tuple[AxisBreakpoints, tuple[CutoffCandidate, ...]]


@dataclass(frozen=True)
class CutoffReportInputs:
    """Everything the published record states, already measured.

    Attributes:
        objective: The declared objective and its rate floors, published verbatim.
        caller: ``kestrel``, ``advntr`` or ``both``.
        folds: Requested outer folds.
        seed: Seed for the grouped allocation.
        workers: Replay worker count.
        max_breakpoints: Optional per-axis cap on tested breakpoints.
        samples: Every declared cohort row, before duplicate control.
        primary: The scored group representatives.
        derived: Per-axis breakpoints with their complete replayed candidates.
        anchors: Axis name to the candidate reproducing the shipped policy.
        arms: Replayed observations keyed by policy id, including the baseline arm.
        replay: The grid replay, supplying the evidence digests.
        parity: The proven baseline-parity record.
        evaluation: Grouped cross-validation and the separate full-data selection.
        advntr_result: The native adVNTR grid result, when adVNTR was requested.
        comparisons: Axis name to the exact production comparator it sweeps.
        cohort_manifest_sha256: Digest of the declared cohort manifest.
        capture_manifest_sha256: Digest of the capture association manifest.
        generator_version: Identity of the generator that produced this record.
    """

    objective: SearchSpec
    caller: str
    folds: int
    seed: int
    workers: int
    max_breakpoints: int | None
    samples: Sequence[CohortSample]
    primary: Sequence[CohortSample]
    derived: Sequence[DerivedAxis]
    anchors: Mapping[str, CutoffCandidate]
    arms: Mapping[str, tuple[CallerObservation, ...]]
    replay: KestrelGridReplay
    parity: Mapping[str, Any]
    evaluation: Mapping[str, Any]
    advntr_result: AdvntrCutoffGridResult | None
    comparisons: Mapping[str, str]
    cohort_manifest_sha256: str
    capture_manifest_sha256: str
    generator_version: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _truth_set(samples: Sequence[CohortSample], primary: Sequence[CohortSample]) -> dict[str, Any]:
    """Composition of the scored roster, and every duplicate dropped to reach it."""
    representatives = {sample.group_id: sample.sample_id for sample in primary}
    genotypes: dict[str, int] = {}
    assemblies: dict[str, int] = {}
    for sample in primary:
        label = "unknown" if sample.genotype is None else ("positive" if sample.genotype else "negative")
        genotypes[label] = genotypes.get(label, 0) + 1
        assemblies[sample.assembly] = assemblies.get(sample.assembly, 0) + 1
    scored = {sample.sample_id for sample in primary}
    dropped = [
        {
            "sample_id": sample.sample_id,
            "group_id": sample.group_id,
            "representative": representatives[sample.group_id],
            "reason": _DUPLICATE_REASON,
        }
        for sample in samples
        if sample.sample_id not in scored
    ]
    return {
        "sample_count": len(samples),
        "primary_count": len(primary),
        "dropped_count": len(dropped),
        "dropped_duplicates": dropped,
        "by_genotype": dict(sorted(genotypes.items())),
        "by_assembly": dict(sorted(assemblies.items())),
        "unknown_truth_count": genotypes.get("unknown", 0),
        "truth_variant_count": sum(sample.truth_variant is not None for sample in primary),
    }


def _outcomes(rows: Sequence[CallerObservation]) -> tuple[bool | None, ...]:
    """The per-sample decision vector of one arm, in its fixed group order."""
    return tuple(row.called_positive for row in rows)


def _plateau(
    axis: AxisBreakpoints,
    candidates: Sequence[CutoffCandidate],
    arms: Mapping[str, Sequence[CallerObservation]],
    selected: str,
    comparison: str,
) -> dict[str, Any]:
    """The contiguous run of tested values that reproduce the selected outcome exactly.

    Args:
        axis: The axis the selected candidate belongs to.
        candidates: That axis's candidates, ascending in axis value.
        arms: Replayed observations keyed by candidate id.
        selected: The selected candidate id on this axis.
        comparison: The exact production comparator this axis sweeps.

    Returns:
        The plateau record, naming the tested values that are indistinguishable from the
        selected one and the neighbouring values at which the outcome changes.

    Raises:
        ValueError: If the selected candidate does not belong to the supplied axis.
    """
    positions = {candidate.candidate_id: index for index, candidate in enumerate(candidates)}
    if selected not in positions:
        _fail(f"cutoff report selected candidate {selected} is not a breakpoint of axis {axis.axis}")
    index = positions[selected]
    target = _outcomes(arms[selected])
    low = index
    while low > 0 and _outcomes(arms[candidates[low - 1].candidate_id]) == target:
        low -= 1
    high = index
    while high + 1 < len(candidates) and _outcomes(arms[candidates[high + 1].candidate_id]) == target:
        high += 1
    values = [float(value) for value in axis.values[low : high + 1]]
    return {
        "axis": axis.axis,
        "comparison": comparison,
        "selected_value": float(axis.values[index]),
        "equivalent_values": values,
        "interval_low": values[0],
        "interval_high": values[-1],
        "open_below": float(axis.values[low - 1]) if low > 0 else None,
        "open_above": float(axis.values[high + 1]) if high + 1 < len(candidates) else None,
        "width": values[-1] - values[0],
        "single_value": len(values) == 1,
        "note": _PLATEAU_NOTE,
    }


def _rate(events: int, total: int) -> Fraction | None:
    return Fraction(events, total) if total else None


def _infeasible(arms: Mapping[str, Sequence[CallerObservation]], spec: SearchSpec, reason: str) -> dict[str, Any]:
    """Name the constraint no candidate met, and the best value that was reachable."""
    best: dict[str, float | None] = {}
    unsatisfiable: list[str] = []
    achievable: dict[str, Fraction | None] = {"min_sensitivity": None, "min_specificity": None}
    for rows in arms.values():
        counts = cutoff_counts(rows)
        for name, value in (
            ("min_sensitivity", _rate(counts.true_positives, counts.positive_count)),
            ("min_specificity", _rate(counts.true_negatives, counts.negative_count)),
        ):
            current = achievable[name]
            if value is not None and (current is None or value > current):
                achievable[name] = value
    for name, floor in (("min_sensitivity", spec.min_sensitivity), ("min_specificity", spec.min_specificity)):
        if floor is None:
            continue
        reached = achievable[name]
        best[name] = None if reached is None else float(reached)
        if reached is None or reached < Fraction(str(floor)):
            unsatisfiable.append(name)
    if reason == "training-truth-class-missing":
        note = "the scored roster does not contain both truth classes, so no objective can be optimized"
    elif unsatisfiable:
        note = "no tested cutoff satisfies " + ", ".join(unsatisfiable)
    else:
        note = "each declared constraint is individually reachable, but no single tested cutoff satisfies them all"
    return {"unsatisfiable_constraints": unsatisfiable, "best_achievable": best, "note": note, "reason": reason}


def _old_versus_derived(baseline: CallerPolicyValues, selected: CallerPolicyValues) -> list[dict[str, Any]]:
    """Every decision pointer, its shipped value and its derived value, changed or not."""
    return [
        {
            "pointer": pointer,
            "baseline_value": baseline.values[pointer],
            "derived_value": selected.values[pointer],
            "changed": baseline.values[pointer] != selected.values[pointer],
        }
        for pointer in sorted(baseline.values)
    ]


def _cutoff_rows(
    arms: Mapping[str, Sequence[CallerObservation]], derived: Sequence[DerivedAxis]
) -> list[dict[str, Any]]:
    """One published row per tested cutoff, counts and exact intervals included."""
    located: dict[str, tuple[str, float, CutoffCandidate]] = {}
    for axis, candidates in derived:
        for index, candidate in enumerate(candidates):
            located[candidate.candidate_id] = (axis.axis, float(axis.values[index]), candidate)
    ordered = [BASELINE_ID, *[candidate.candidate_id for _, candidates in derived for candidate in candidates]]
    rows: list[dict[str, Any]] = []
    for policy_id in ordered:
        if policy_id not in arms:
            continue
        observations = tuple(arms[policy_id])
        found = located.get(policy_id)
        rows.append(
            {
                "policy_id": policy_id,
                "axis": None if found is None else found[0],
                "value": None if found is None else found[1],
                "policy_sha256": None if found is None else found[2].policy.sha256,
                "parameters": {} if found is None else dict(found[2].parameters),
                "counts": cutoff_counts_document(cutoff_counts(observations)),
                "metrics": caller_metrics_document(calculate_caller_metrics(observations)),
            }
        )
    return rows


def _boundary_support(curve: Mapping[str, Any]) -> dict[str, Any]:
    """The selected axis's boundary support, with an explicit warning when it is empty."""
    support = dict(curve["boundary_support"])
    support["warnings"] = [
        f"no {label}-truth sample lies inside the tested band of axis {curve['axis']}, so the selected cutoff "
        f"is not constrained by any observed {label}"
        for label, key in (("positive", "positives_within_band"), ("negative", "negatives_within_band"))
        if not support.get(key)
    ]
    return support


def _joint_points(inputs: CutoffReportInputs) -> dict[str, Any] | None:
    """The labelled multi-axis table, which exists only when more than one axis was swept."""
    if len(inputs.derived) < 2:
        return None
    labelled = {
        f"{axis.axis}={float(axis.values[index])}": (candidate, inputs.arms[candidate.candidate_id])
        for axis, candidates in inputs.derived
        for index, candidate in enumerate(candidates)
        if candidate.parameters
    }
    return joint_points_document(build_joint_points(labelled, inputs.anchors[inputs.derived[0][0].axis]))


def build_cutoff_report_document(inputs: CutoffReportInputs) -> dict[str, Any]:
    """Assemble the complete ``calibration-cutoff-report-v1`` document.

    Args:
        inputs: Everything already measured by the run, described field by field on
            :class:`CutoffReportInputs`.

    Returns:
        The complete report document. The profile section is deliberately left
        unavailable: only the run that actually writes and re-reads the profile may
        claim it round-tripped, so the caller replaces this section on success.

    Raises:
        ValueError: If a curve cannot be built because a candidate varies off its axis,
            if the arms disagree about the roster, or if the selected candidate is not a
            breakpoint of the axis it claims.
    """
    curves = [
        axis_curve_document(
            build_axis_curve(
                axis,
                candidates,
                inputs.arms,
                comparison=inputs.comparisons[axis.axis],
                phase="policy-selection",
            )
        )
        for axis, candidates in inputs.derived
    ]
    selected = inputs.evaluation["final_selection"]["policy_id"]
    by_candidate = {
        candidate.candidate_id: (axis, candidates) for axis, candidates in inputs.derived for candidate in candidates
    }
    first_axis, first_candidates = inputs.derived[0]
    anchor = inputs.anchors[first_axis.axis]
    resolved = selected if selected in by_candidate else anchor.candidate_id
    axis, candidates = by_candidate.get(resolved, (first_axis, first_candidates))
    policy = {candidate.candidate_id: candidate.policy for candidate in candidates}.get(resolved, anchor.policy)
    curve = next(entry for entry in curves if entry["axis"] == axis.axis)
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "selected" if selected is not None else "infeasible",
        "successful": selected is not None,
        "usage_hint": USAGE_HINT,
        "caller": inputs.caller,
        "objective": {
            "objective": inputs.objective.objective,
            "min_sensitivity": inputs.objective.min_sensitivity,
            "min_specificity": inputs.objective.min_specificity,
        },
        "folds_requested": inputs.folds,
        "seed": inputs.seed,
        "workers": inputs.workers,
        "max_breakpoints": inputs.max_breakpoints,
        "truth_set": _truth_set(inputs.samples, inputs.primary),
        "provenance": {
            "cohort_manifest_sha256": inputs.cohort_manifest_sha256,
            "capture_manifest_sha256": inputs.capture_manifest_sha256,
            "baseline_policy_sha256": inputs.replay.baseline_policy_sha256,
            "grid_replay_sha256": inputs.replay.sha256,
            "grid_input_sha256": inputs.replay.input_sha256,
            "capture_file_sha256": dict(inputs.replay.capture_file_sha256),
            "native_file_sha256": dict(inputs.replay.native_file_sha256),
            "generator_version": inputs.generator_version,
            "advntr": None if inputs.advntr_result is None else {"sha256": inputs.advntr_result.sha256},
        },
        "baseline_parity": dict(inputs.parity),
        "axes": [axis_document(axis) for axis, _ in inputs.derived],
        "cutoffs": _cutoff_rows(inputs.arms, inputs.derived),
        "curves": curves,
        "joint_points": _joint_points(inputs),
        "rejected_candidates": [
            {"axis": entry.axis, "value": value, "reason": reason}
            for entry, _ in inputs.derived
            for value, reason in entry.rejected
        ],
        "evaluation": dict(inputs.evaluation),
        "selection": {
            "policy_id": selected,
            "reason": inputs.evaluation["final_selection"]["reason"],
            "axis": axis.axis if selected is not None else None,
            "value": float(axis.values[[entry.candidate_id for entry in candidates].index(resolved)])
            if selected is not None
            else None,
            "policy_sha256": policy.sha256 if selected is not None else None,
            "plateau": _plateau(axis, candidates, inputs.arms, resolved, inputs.comparisons[axis.axis])
            if selected is not None
            else None,
            "infeasible": None
            if selected is not None
            else _infeasible(inputs.arms, inputs.objective, str(inputs.evaluation["final_selection"]["reason"])),
        },
        "old_versus_derived": _old_versus_derived(anchor.policy, policy),
        "boundary_support": _boundary_support(curve),
        "profile": {"status": "unavailable", "reason": _NO_SELECTION},
        "limitations": LIMITATIONS,
    }


__all__ = [
    "BASELINE_ID",
    "LIMITATIONS",
    "PROFILE_NAME",
    "USAGE_HINT",
    "CutoffReportInputs",
    "DerivedAxis",
    "build_cutoff_report_document",
]
