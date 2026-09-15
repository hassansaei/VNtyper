"""Complete-roster outer-fold evaluation of actual caller cutoff observations."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from vntyper.scripts.calibration_caller_metrics import (
    CallerObservation,
    calculate_caller_metrics,
    validate_caller_observations,
)
from vntyper.scripts.calibration_cohort_metrics import (
    caller_metrics_document,
    group_folds,
    paired_caller_differences,
)
from vntyper.scripts.calibration_cutoff_selection import (
    CutoffSelection,
    SearchSpec,
    cutoff_counts,
    cutoff_counts_document,
    select_cutoff_policy,
)


def _selection_document(selection: CutoffSelection) -> dict[str, Any]:
    return {
        "policy_id": selection.policy_id,
        "reason": selection.reason,
        "training_count": selection.training_count,
        "eligible_candidates": selection.eligible_candidates,
        "selected_metrics": (
            cutoff_counts_document(selection.selected_metrics) if selection.selected_metrics is not None else None
        ),
        "baseline_metrics": cutoff_counts_document(selection.baseline_metrics),
    }


def _metrics(rows: tuple[CallerObservation, ...]) -> dict[str, Any]:
    return {
        "counts": cutoff_counts_document(cutoff_counts(rows)),
        "exact": caller_metrics_document(calculate_caller_metrics(rows)),
    }


def _bound_arms(
    arms: Mapping[str, Sequence[CallerObservation]], baseline_id: str
) -> dict[str, tuple[CallerObservation, ...]]:
    if not isinstance(arms, Mapping) or baseline_id not in arms:
        raise ValueError("cutoff evaluation requires an explicit baseline arm")
    if any(not isinstance(name, str) or not name or name.strip() != name for name in arms):
        raise ValueError("cutoff evaluation policy IDs must be nonempty trimmed strings")
    validated = {name: validate_caller_observations(tuple(rows)) for name, rows in arms.items()}
    reference = tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in validated[baseline_id])
    for rows in validated.values():
        if tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in rows) != reference:
            raise ValueError("cutoff evaluation requires complete identical specimen, group and truth rosters")
    return validated


def evaluate_cutoff_arms(
    arms: Mapping[str, Sequence[CallerObservation]],
    *,
    baseline_id: str = "baseline",
    spec: SearchSpec,
    folds: int = 5,
    seed: int = 20260915,
) -> dict[str, Any]:
    """Compare a fixed native baseline with training-selected held-out predictions.

    Args:
        arms: Actual replay/native outcomes over the same complete primary roster.
        baseline_id: Unchanged policy retained in the candidate inventory.
        spec: Explicit training objective and optional class-rate constraints.
        folds: Requested number of outcome-independent outer folds, at least two.
        seed: Seed for group allocation and descriptive paired uncertainty.

    Returns:
        JSON-compatible counts, exact metric intervals, fold decisions and individual
        held-out predictions. Full-data selection and operating points are separate.
        An absent training class uses the unchanged baseline as an explicit fallback.
        No final selection or insufficient groups produces unavailable status.

    Raises:
        ValueError: For incomplete/inconsistent rosters, invalid folds or search spec.
    """
    validated = _bound_arms(arms, baseline_id)
    baseline = validated[baseline_id]
    assignments = group_folds({r.key: r.group_key for r in baseline}, folds=folds, seed=seed)
    indices = {name: {row.key: row for row in rows} for name, rows in validated.items()}
    heldout: dict[str, CallerObservation] = {}
    used_policy: dict[str, str] = {}
    fold_records = []
    for fold in sorted(set(assignments.values())):
        training = sorted(key for key, value in assignments.items() if value != fold)
        held = sorted(key for key, value in assignments.items() if value == fold)
        selection = select_cutoff_policy(validated, training, baseline_id, spec)
        policy_id = selection.policy_id if selection.policy_id is not None else baseline_id
        for key in held:
            heldout[key] = indices[policy_id][key]
            used_policy[key] = policy_id
        fold_records.append(
            {
                "fold": fold,
                "training_keys": training,
                "held_out_keys": held,
                "selection": _selection_document(selection),
                "used_policy": policy_id,
                "fallback_reason": selection.reason if selection.policy_id is None else None,
            }
        )
    final = select_cutoff_policy(validated, [row.key for row in baseline], baseline_id, spec)
    heldout_rows = tuple(heldout[row.key] for row in baseline) if assignments else ()
    return {
        "schema_version": "calibration-cutoff-evaluation-v1",
        "status": "available" if assignments and final.policy_id is not None else "unavailable",
        "cross_validation_available": bool(assignments),
        "status_reason": ("insufficient-cross-validation-groups" if not assignments else final.reason),
        "folds_requested": folds,
        "seed": seed,
        "objective": spec.objective,
        "min_sensitivity": spec.min_sensitivity,
        "min_specificity": spec.min_specificity,
        "final_selection": _selection_document(final),
        "folds": fold_records,
        "rows": [
            {
                "key": row.key,
                "group_key": row.group_key,
                "truth_positive": row.truth_positive,
                "baseline": row.called_positive,
                "held_out": heldout[row.key].called_positive if row.key in heldout else None,
                "selected_policy": used_policy.get(row.key),
                "fold": assignments.get(row.key),
            }
            for row in baseline
        ],
        "baseline": _metrics(baseline),
        "held_out": _metrics(heldout_rows) if heldout_rows else None,
        "full_data_operating_points": {
            name: cutoff_counts_document(cutoff_counts(rows)) for name, rows in validated.items()
        },
        "full_data_operating_points_scope": "descriptive searched-cohort points; separate from held-out evaluation",
        "paired_differences": paired_caller_differences(baseline, heldout_rows, seed=seed) if heldout_rows else None,
        "limitations": "Previously examined data or revised searches remain exploratory; no independent approval.",
    }
