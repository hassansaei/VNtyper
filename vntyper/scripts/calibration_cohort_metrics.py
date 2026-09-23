"""Group-preserving folds and descriptive development-cohort error summaries."""

from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import fields
from fractions import Fraction

from vntyper.scripts.calibration_caller_metrics import (
    CallerMetrics,
    CallerObservation,
    CallerRate,
    validate_caller_observations,
)

_STRATUM_TYPES = (bool, int, str, type(None))


def _stratified_order(ordered: list[str], strata: object, generator: random.Random) -> list[str]:
    """Groups shuffled within each stratum, strata concatenated in a fixed label order."""
    if not isinstance(strata, Mapping) or any(group not in strata for group in ordered):
        raise ValueError("cohort fold strata must map every group to its label")
    if any(type(strata[group]) not in _STRATUM_TYPES for group in ordered):
        raise ValueError("cohort fold strata labels must be booleans, integers, strings or None")
    labels = sorted({strata[group] for group in ordered}, key=lambda label: (type(label).__name__, repr(label)))
    result: list[str] = []
    for label in labels:
        members = [group for group in ordered if strata[group] == label and type(strata[group]) is type(label)]
        generator.shuffle(members)
        result.extend(members)
    return result


def group_folds(
    groups: Mapping[str, str], *, folds: int, seed: int, strata: Mapping[str, object] | None = None
) -> dict[str, int]:
    """Assign every declared biological group wholly to one reproducible fold.

    Args:
        groups: Sample to biological group mapping.
        folds: Requested folds, at least two; reduced to the available groups.
        seed: Explicit deterministic random seed.
        strata: Optional group to label mapping (for example truth class). Groups are
            then shuffled within each label and dealt round-robin across folds with one
            running position over all labels, so a label with at least as many groups as
            folds reaches every fold and fold sizes still differ by at most one. ``None``
            keeps the unstratified assignment exactly.

    Returns:
        Sample to zero-based fold mapping, or empty when fewer than two groups exist.

    Raises:
        ValueError: For invalid folds, seed, identities or strata.
    """
    if type(folds) is not int or folds < 2 or type(seed) is not int:
        raise ValueError("cohort folds must be an integer >=2 and seed must be an integer")
    if any(
        not isinstance(key, str) or not key or not isinstance(value, str) or not value for key, value in groups.items()
    ):
        raise ValueError("cohort fold identities must be nonempty strings")
    ordered = sorted(set(groups.values()))
    if len(ordered) < 2:
        return {}
    generator = random.Random(seed)
    if strata is None:
        generator.shuffle(ordered)
    else:
        ordered = _stratified_order(ordered, strata, generator)
    assignments = {group: index % min(folds, len(ordered)) for index, group in enumerate(ordered)}
    return {key: assignments[groups[key]] for key in sorted(groups)}


def _interval(differences: Sequence[float], groups: Sequence[str], seed: int) -> list[float] | None:
    buckets: dict[str, list[float]] = {group: [] for group in sorted(set(groups))}
    for value, group in zip(differences, groups, strict=True):
        buckets[group].append(value)
    if len(buckets) < 2:
        return None
    if len(set(differences)) == 1:
        return [differences[0], differences[0]]
    generator = random.Random(seed)
    keys = tuple(buckets)
    sampled = []
    for _ in range(10_000):
        values = [value for group in generator.choices(keys, k=len(keys)) for value in buckets[group]]
        sampled.append(sum(values) / len(values))
    sampled.sort()
    return [sampled[249], sampled[9749]]


def length_metrics(
    truth: Sequence[float],
    predictions: Sequence[float | None],
    baseline: Sequence[float | None],
    groups: Sequence[str],
    *,
    seed: int,
) -> dict[str, object]:
    """Summarize outer-fold predictions, retaining unavailable eligible observations.

    Args:
        truth: Exact total-repeat targets.
        predictions: Held-out predictions or explicit unavailability.
        baseline: Corresponding training-fold target means.
        groups: Biological groups for descriptive paired cluster bootstrap.
        seed: Prespecified bootstrap seed.

    Returns:
        Counts, complete-case errors and paired baseline differences with explicit denominators.

    Raises:
        ValueError: For mismatched rows or nonfinite numerical values.
    """
    if not len(truth) == len(predictions) == len(baseline) == len(groups):
        raise ValueError("cohort metric arrays must have identical lengths")
    if any(not math.isfinite(value) for value in truth) or any(
        value is not None and not math.isfinite(value) for value in (*predictions, *baseline)
    ):
        raise ValueError("cohort metrics require finite numerical values")
    rows = [(target, value) for target, value in zip(truth, predictions, strict=True) if value is not None]
    errors = [value - target for target, value in rows]
    paired = [
        (abs(value - target), abs(old - target), group)
        for target, value, old, group in zip(truth, predictions, baseline, groups, strict=True)
        if value is not None and old is not None
    ]
    differences = [current - old for current, old, _group in paired]
    denominator = sum((target - sum(row[0] for row in rows) / len(rows)) ** 2 for target, _ in rows) if rows else 0
    return {
        "eligible": len(truth),
        "predicted": len(rows),
        "unavailable": len(truth) - len(rows),
        "independent_groups": len(set(groups)),
        "availability": len(rows) / len(truth) if truth else None,
        "mae": sum(map(abs, errors)) / len(errors) if errors else None,
        "rmse": math.sqrt(sum(value * value for value in errors) / len(errors)) if errors else None,
        "bias": sum(errors) / len(errors) if errors else None,
        "r2": 1 - sum(value * value for value in errors) / denominator if denominator else None,
        "baseline_mae_paired": sum(old for _, old, _ in paired) / len(paired) if paired else None,
        "paired_count": len(paired),
        "paired_mae_delta": sum(differences) / len(differences) if paired else None,
        "paired_mae_delta_interval": _interval(differences, [group for _, _, group in paired], seed),
        "uncertainty": "descriptive group bootstrap; does not include model-selection uncertainty",
    }


def caller_metrics_document(metrics: CallerMetrics) -> dict[str, object]:
    """Convert exact caller rates to JSON values without dropping their denominators.

    Args:
        metrics: Existing caller metric result from fixed complete observations.

    Returns:
        JSON-compatible metric fields and rate numerator/denominator/interval objects.
    """
    result: dict[str, object] = {}
    for field in fields(metrics):
        value = getattr(metrics, field.name)
        if isinstance(value, CallerRate):
            result[field.name] = {
                part.name: (float(item) if isinstance(item, Fraction) else item)
                for part in fields(value)
                for item in [getattr(value, part.name)]
            }
        else:
            result[field.name] = float(value) if isinstance(value, Fraction) else value
    return result


def paired_caller_differences(
    baseline: Sequence[CallerObservation],
    candidate: Sequence[CallerObservation],
    *,
    seed: int,
) -> dict[str, object]:
    """Measure held-out paired rate changes with descriptive group bootstrap intervals.

    Args:
        baseline: Complete native baseline observations, one per independent group.
        candidate: Same observations evaluated using the fold-selected policy.
        seed: Prespecified bootstrap seed.

    Returns:
        Candidate-minus-baseline rate differences and explicit eligible counts.

    Raises:
        ValueError: If identity/truth rosters differ.
    """
    old = validate_caller_observations(tuple(baseline))
    new = validate_caller_observations(tuple(candidate))
    if len(old) != len(new) or any(
        (a.key, a.group_key, a.truth_positive, a.truth_variants)
        != (b.key, b.group_key, b.truth_positive, b.truth_variants)
        for a, b in zip(old, new, strict=True)
    ):
        raise ValueError("paired caller outcomes require identical truth and identity rosters")
    result: dict[str, object] = {}
    for name, truth_class, called in (
        ("sensitivity", True, True),
        ("specificity", False, False),
        ("false_positive_rate", False, True),
        ("no_call_rate", None, None),
    ):
        pairs = [
            (a, b) for a, b in zip(old, new, strict=True) if truth_class is None or a.truth_positive is truth_class
        ]
        differences = [float(int(b.called_positive is called) - int(a.called_positive is called)) for a, b in pairs]
        result[name] = {
            "eligible": len(pairs),
            "delta": sum(differences) / len(differences) if pairs else None,
            "interval": _interval(differences, [a.group_key for a, _ in pairs], seed),
        }
    result["uncertainty"] = "descriptive group bootstrap of held-out predictions; excludes selection uncertainty"
    return result
