"""Training-only objective selection using fixed-denominator integer counts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from fractions import Fraction

from vntyper.scripts.calibration_caller_metrics import CallerObservation, validate_caller_observations


@dataclass(frozen=True)
class SearchSpec:
    """An explicitly chosen development objective and optional rate constraints."""

    objective: str = "balanced-accuracy"
    min_sensitivity: float | None = None
    min_specificity: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.objective, str) or self.objective not in {
            "balanced-accuracy",
            "sensitivity",
            "specificity",
        }:
            raise ValueError("unsupported cutoff search objective")
        for value in (self.min_sensitivity, self.min_specificity):
            if value is not None and (
                type(value) not in {int, float} or not 0 <= value <= 1 or not math.isfinite(value)
            ):
                raise ValueError("cutoff search minimum rates must be finite fractions in [0,1]")


@dataclass(frozen=True)
class CutoffCounts:
    """Binary outcomes with no-calls retained in positive/negative denominators."""

    eligible_count: int
    positive_count: int
    negative_count: int
    unknown_truth_count: int
    true_positives: int
    false_negatives: int
    true_negatives: int
    false_positives: int
    positive_no_calls: int
    negative_no_calls: int
    no_calls: int


@dataclass(frozen=True)
class CutoffSelection:
    """A training decision, including explicit failure to select any candidate."""

    policy_id: str | None
    reason: str
    training_count: int
    eligible_candidates: int
    selected_metrics: CutoffCounts | None
    baseline_metrics: CutoffCounts


def cutoff_counts(observations: Sequence[CallerObservation]) -> CutoffCounts:
    """Count validated independent observations without computing confidence bounds.

    Args:
        observations: Complete primary representatives; empty input is allowed.

    Returns:
        Exact counts, including known-class and overall no-calls.

    Raises:
        ValueError: For malformed observations or duplicate specimen/group identities.
    """
    rows = validate_caller_observations(tuple(observations)) if observations else ()
    return CutoffCounts(
        len(rows),
        sum(r.truth_positive is True for r in rows),
        sum(r.truth_positive is False for r in rows),
        sum(r.truth_positive is None for r in rows),
        sum(r.truth_positive is True and r.called_positive is True for r in rows),
        sum(r.truth_positive is True and r.called_positive is False for r in rows),
        sum(r.truth_positive is False and r.called_positive is False for r in rows),
        sum(r.truth_positive is False and r.called_positive is True for r in rows),
        sum(r.truth_positive is True and r.called_positive is None for r in rows),
        sum(r.truth_positive is False and r.called_positive is None for r in rows),
        sum(r.called_positive is None for r in rows),
    )


def _fraction(events: int, total: int) -> Fraction | None:
    return Fraction(events, total) if total else None


def _rates(counts: CutoffCounts) -> dict[str, Fraction | None]:
    sensitivity = _fraction(counts.true_positives, counts.positive_count)
    specificity = _fraction(counts.true_negatives, counts.negative_count)
    return {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "false_positive_rate": _fraction(counts.false_positives, counts.negative_count),
        "precision": _fraction(counts.true_positives, counts.true_positives + counts.false_positives),
        # This is the harmonic mean of precision and fixed-denominator sensitivity:
        # positive no-calls are unsuccessful positive outcomes, not dropped cases.
        "f1": _fraction(
            2 * counts.true_positives, counts.positive_count + counts.true_positives + counts.false_positives
        ),
        "balanced_accuracy": (sensitivity + specificity) / 2
        if sensitivity is not None and specificity is not None
        else None,
        "no_call_rate": _fraction(counts.no_calls, counts.eligible_count),
    }


def cutoff_counts_document(counts: CutoffCounts) -> dict[str, int | float | None]:
    """Serialize counts and bounded rates, with undefined denominators as null.

    Args:
        counts: Counts computed from the supplied observation roster.

    Returns:
        Counts plus sensitivity, specificity, FPR, precision, F1, balanced accuracy
        and no-call rate. This fast representation contains no confidence bounds.
    """
    return {
        **asdict(counts),
        **{key: float(value) if value is not None else None for key, value in _rates(counts).items()},
    }


def _training_rows(rows: Sequence[CallerObservation], keys: set[str]) -> tuple[CallerObservation, ...]:
    selected = tuple(row for row in rows if row.key in keys)
    if len(selected) != len(keys) or {row.key for row in selected} != keys:
        raise ValueError("cutoff training roster is missing or duplicated in an arm")
    return validate_caller_observations(selected) if selected else ()


def select_cutoff_policy(
    arms: Mapping[str, Sequence[CallerObservation]],
    training_keys: Sequence[str],
    baseline_id: str,
    spec: SearchSpec,
) -> CutoffSelection:
    """Optimize a declared objective using only the explicitly supplied training keys.

    Args:
        arms: Actual replay/native observation inventories keyed by candidate ID.
        training_keys: Unique keys allowed to contribute labels or calls to selection.
        baseline_id: Explicit unchanged policy inventory member.
        spec: Objective and optional minimum sensitivity/specificity constraints.

    Returns:
        Selection with exact training counts. Both known truth classes are required.
        Ties prefer fewer FP, more TP, fewer no-calls, baseline, then stable ID.

    Raises:
        ValueError: For missing/duplicate training rows, inconsistent truth/group
            bindings, absent baseline or invalid constraints.
    """
    if not isinstance(spec, SearchSpec):
        raise ValueError("cutoff selection requires a typed SearchSpec")
    spec.__post_init__()
    if (
        not isinstance(arms, Mapping)
        or baseline_id not in arms
        or any(not isinstance(name, str) or not name for name in arms)
    ):
        raise ValueError("cutoff selection requires named candidates and an explicit baseline")
    if any(not isinstance(key, str) or not key for key in training_keys) or len(set(training_keys)) != len(
        training_keys
    ):
        raise ValueError("cutoff training keys must be unique nonempty strings")
    keys = set(training_keys)
    baseline_rows = _training_rows(arms[baseline_id], keys)
    baseline = cutoff_counts(baseline_rows)
    if not baseline.positive_count or not baseline.negative_count:
        return CutoffSelection(None, "training-truth-class-missing", len(keys), 0, None, baseline)
    identity = tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in baseline_rows)
    eligible: list[tuple[Fraction, int, int, int, bool, str, CutoffCounts]] = []
    for name, observations in arms.items():
        rows = _training_rows(observations, keys)
        if tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in rows) != identity:
            raise ValueError("cutoff training truth and group identities differ across arms")
        counts = cutoff_counts(rows)
        rates = _rates(counts)
        sensitivity, specificity = rates["sensitivity"], rates["specificity"]
        assert sensitivity is not None and specificity is not None
        if (spec.min_sensitivity is not None and sensitivity < Fraction(str(spec.min_sensitivity))) or (
            spec.min_specificity is not None and specificity < Fraction(str(spec.min_specificity))
        ):
            continue
        objective = rates[spec.objective.replace("-", "_")]
        assert objective is not None
        eligible.append(
            (
                -objective,
                counts.false_positives,
                -counts.true_positives,
                counts.no_calls,
                name != baseline_id,
                name,
                counts,
            )
        )
    if not eligible:
        return CutoffSelection(None, "no-candidate-satisfies-constraints", len(keys), 0, None, baseline)
    chosen = min(eligible)
    return CutoffSelection(chosen[5], "selected", len(keys), len(eligible), chosen[6], baseline)
