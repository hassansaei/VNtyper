"""Complete-roster outer-fold evaluation of actual caller cutoff observations."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from typing import Any, NoReturn

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

logger = logging.getLogger(__name__)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


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
        _fail("cutoff evaluation requires an explicit baseline arm")
    if any(not isinstance(name, str) or not name or name.strip() != name for name in arms):
        _fail("cutoff evaluation policy IDs must be nonempty trimmed strings")
    validated = {name: validate_caller_observations(tuple(rows)) for name, rows in arms.items()}
    reference = tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in validated[baseline_id])
    for rows in validated.values():
        if tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in rows) != reference:
            _fail("cutoff evaluation requires complete identical specimen, group and truth rosters")
    return validated


def _bound_contributors(
    contributors: object, arms: Mapping[str, tuple[CallerObservation, ...]], baseline_id: str
) -> dict[str, frozenset[str]] | None:
    """Validate the candidate-to-contributing-samples map against the bound roster."""
    if contributors is None:
        return None
    if not isinstance(contributors, Mapping):
        _fail("cutoff evaluation contributors must map candidate IDs to sample keys")
    roster = {row.key for row in arms[baseline_id]}
    bound: dict[str, frozenset[str]] = {}
    for name, keys in contributors.items():
        if name not in arms or name == baseline_id:
            _fail(f"cutoff evaluation contributors name {name!r}, which is not a non-baseline candidate")
        if not isinstance(keys, frozenset) or not keys or not keys <= roster:
            _fail(f"cutoff evaluation contributors for {name!r} must be a nonempty frozenset of roster keys")
        bound[name] = keys
    return bound


def _bound_inventories(
    inventories: object,
    arms: Mapping[str, tuple[CallerObservation, ...]],
    baseline_id: str,
    folds: set[int],
) -> dict[int, frozenset[str]]:
    """Validate the fold-to-admissible-candidates map against the computed folds and arms."""
    if not isinstance(inventories, Mapping):
        _fail("cutoff evaluation fold inventories must map outer folds to candidate IDs")
    if set(inventories) != folds:
        _fail(f"cutoff evaluation fold inventories name folds {sorted(inventories)}, not the computed {sorted(folds)}")
    bound: dict[int, frozenset[str]] = {}
    for fold, ids in inventories.items():
        if not isinstance(ids, frozenset) or any(name not in arms or name == baseline_id for name in ids):
            _fail(f"cutoff evaluation fold inventory {fold} must be a frozenset of non-baseline candidate IDs")
        bound[fold] = ids
    return bound


def outer_fold_assignments(rows: Sequence[tuple[str, str, bool | None]], *, folds: int, seed: int) -> dict[str, int]:
    """The outer fold of every specimen, grouped and stratified by truth label.

    Rows are one representative per independent group, so each group has exactly one
    truth label. Stratifying on it keeps a class with at least ``folds`` groups present in
    every fold; the labels steer allocation only, never a fold's selection.

    Args:
        rows: ``(key, group_key, truth_positive)`` per specimen.
        folds: Requested number of outer folds, at least two.
        seed: Seed for group allocation.

    Returns:
        Specimen key to outer fold; empty when there are too few groups to cross-validate.

    Raises:
        ValueError: For invalid fold counts or seeds.
    """
    return group_folds(
        {key: group for key, group, _ in rows},
        folds=folds,
        seed=seed,
        strata={group: truth for _, group, truth in rows},
    )


def evaluate_cutoff_arms(
    arms: Mapping[str, Sequence[CallerObservation]],
    *,
    baseline_id: str = "baseline",
    spec: SearchSpec,
    folds: int = 5,
    seed: int = 20260915,
    contributors: Mapping[str, frozenset[str]] | None = None,
    fold_inventories: Mapping[int, frozenset[str]] | None = None,
    tie_keys: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Compare a fixed native baseline with training-selected held-out predictions.

    Args:
        arms: Actual replay/native outcomes over the same complete primary roster.
        baseline_id: Unchanged policy retained in the candidate inventory.
        spec: Explicit training objective and optional class-rate constraints.
        folds: Requested number of outer folds, at least two. Groups are allocated
            stratified by truth label, so a class with at least ``folds`` groups is
            present in every fold.
        seed: Seed for group allocation and descriptive paired uncertainty.
        contributors: Optional map from a data-derived candidate ID to the sample keys
            whose observed feature values produced its threshold. Inside an outer fold
            such a candidate is admissible only when at least one contributor is a
            training sample, so a held-out sample's feature values cannot create the
            breakpoint its own fold selects. Candidates absent from the map (the
            baseline, anchors reproducing it, endpoint sentinels) are always admissible.
            ``None`` keeps every candidate admissible in every fold. The full-data
            selection always searches every candidate.
        fold_inventories: Optional map from each outer fold to the candidate IDs its
            training samples alone derived, exclusive with ``contributors``. In fold f only
            the baseline and ``fold_inventories[f]`` are admissible. This is leak-free:
            every admissible value -- breakpoints, endpoint sentinel and cap alike -- is a
            function of fold f's training statistics only, so no held-out value can
            create, remove or displace a threshold its own fold may select. The keys must
            equal the computed folds; with too few groups to cross-validate that is ``{}``.
            The full-data selection always searches every candidate.
        tie_keys: Optional map from every non-baseline candidate ID to a content key,
            passed to every fold and full-data selection (``select_cutoff_policy``). A
            caller whose IDs encode a rank in an inventory that held-out samples helped
            build must pass content keys, or a held-out value could renumber the
            candidates and flip a tie its own fold resolves.

    Returns:
        JSON-compatible counts, exact metric intervals, fold decisions and individual
        held-out predictions. Full-data selection and operating points are separate.
        An absent training class uses the unchanged baseline as an explicit fallback.
        No final selection or insufficient groups produces unavailable status. Each
        fold records how many candidates were admissible to its selection.

    Raises:
        ValueError: For incomplete/inconsistent rosters, invalid folds or search spec,
            or contributors naming an unknown candidate, the baseline, or non-roster keys,
            fold inventories whose folds differ from the computed ones or that name an
            unknown candidate or the baseline, or both contributors and fold inventories.
    """
    if contributors is not None and fold_inventories is not None:
        _fail("cutoff evaluation takes contributors or fold inventories, not both")
    validated = _bound_arms(arms, baseline_id)
    derived = _bound_contributors(contributors, validated, baseline_id)
    baseline = validated[baseline_id]
    assignments = outer_fold_assignments(
        [(r.key, r.group_key, r.truth_positive) for r in baseline], folds=folds, seed=seed
    )
    inventories = (
        None
        if fold_inventories is None
        else _bound_inventories(fold_inventories, validated, baseline_id, set(assignments.values()))
    )
    indices = {name: {row.key: row for row in rows} for name, rows in validated.items()}
    heldout: dict[str, CallerObservation] = {}
    used_policy: dict[str, str] = {}
    fold_records = []
    for fold in sorted(set(assignments.values())):
        training = sorted(key for key, value in assignments.items() if value != fold)
        held = sorted(key for key, value in assignments.items() if value == fold)
        seen = set(training)
        if inventories is not None:
            admissible = {
                name: rows for name, rows in validated.items() if name == baseline_id or name in inventories[fold]
            }
        else:
            admissible = {
                name: rows
                for name, rows in validated.items()
                if derived is None or name not in derived or not derived[name].isdisjoint(seen)
            }
        selection = select_cutoff_policy(admissible, training, baseline_id, spec, tie_keys=tie_keys)
        policy_id = selection.policy_id if selection.policy_id is not None else baseline_id
        for key in held:
            heldout[key] = indices[policy_id][key]
            used_policy[key] = policy_id
        fold_records.append(
            {
                "fold": fold,
                "training_keys": training,
                "held_out_keys": held,
                "admissible_candidates": len(admissible),
                "selection": _selection_document(selection),
                "used_policy": policy_id,
                "fallback_reason": selection.reason if selection.policy_id is None else None,
            }
        )
    final = select_cutoff_policy(validated, [row.key for row in baseline], baseline_id, spec, tie_keys=tie_keys)
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
        "fold_admissibility": (
            "training-derived-inventories"
            if inventories is not None
            else "all-candidates"
            if derived is None
            else "training-observed-breakpoints"
        ),
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
