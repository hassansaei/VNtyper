"""Deterministic exact-rate statistics for grouped calibration evidence."""

from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from types import MappingProxyType


@dataclass(frozen=True)
class PairedObservation:
    """One candidate-versus-baseline value inside a connected leakage group."""

    group_id: str
    stratum: str
    baseline: Fraction
    candidate: Fraction

    def __post_init__(self) -> None:
        """Validate identifiers and exact bounded paired values."""
        if not isinstance(self.group_id, str) or not self.group_id:
            raise ValueError("bootstrap group identifier must be a non-empty string")
        if not isinstance(self.stratum, str) or not self.stratum:
            raise ValueError("bootstrap stratum must be a non-empty string")
        for value, label in ((self.baseline, "baseline"), (self.candidate, "candidate")):
            if not isinstance(value, Fraction) or not 0 <= value <= 1:
                raise ValueError(f"bootstrap {label} must be an exact unit Fraction")


@dataclass(frozen=True)
class BootstrapInterval:
    """One-sided and two-sided percentile bounds from complete-group resampling."""

    estimate: Fraction
    one_sided_lower: Fraction
    two_sided_lower: Fraction
    two_sided_upper: Fraction
    one_sided_noninferiority_p_value: Fraction
    iterations: int
    resampling_unit_count: int


@dataclass(frozen=True)
class BinomialInterval:
    """Two-sided exact Clopper-Pearson binomial interval."""

    estimate: Fraction
    lower: Fraction
    upper: Fraction
    events: int
    total: int


@dataclass(frozen=True)
class CurvePoint:
    """One deterministic ROC and precision-recall threshold row."""

    threshold: Fraction
    true_positives: int
    false_positives: int
    true_negatives: int
    false_negatives: int
    true_positive_rate: Fraction
    false_positive_rate: Fraction
    precision: Fraction
    recall: Fraction


def paired_group_bootstrap(
    observations: Sequence[PairedObservation],
    *,
    required_strata: Sequence[str],
    iterations: int,
    seed: int,
) -> BootstrapInterval:
    """Bootstrap complete connected groups, stratified by the declared cross-product."""
    if iterations != 10_000:
        raise ValueError("calibration bootstrap iterations must be exactly 10000")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("calibration bootstrap seed must be a non-negative integer")
    if not isinstance(observations, Sequence) or not observations:
        raise ValueError("calibration bootstrap requires observations")
    if any(not isinstance(row, PairedObservation) for row in observations):
        raise ValueError("calibration bootstrap rows must be PairedObservation values")
    if (
        not isinstance(required_strata, Sequence)
        or isinstance(required_strata, (str, bytes))
        or not required_strata
        or any(not isinstance(stratum, str) or not stratum for stratum in required_strata)
        or len(required_strata) != len(set(required_strata))
    ):
        raise ValueError("calibration bootstrap requires unique non-empty declared strata")

    group_strata: dict[str, set[str]] = {}
    strata: dict[str, dict[str, list[Fraction]]] = {stratum: {} for stratum in required_strata}
    for row in observations:
        if row.stratum not in strata:
            raise ValueError(f"calibration observation uses undeclared bootstrap stratum: {row.stratum}")
        group_strata.setdefault(row.group_id, set()).add(row.stratum)
        strata[row.stratum].setdefault(row.group_id, []).append(row.candidate - row.baseline)
    if any(len(names) != 1 for names in group_strata.values()):
        raise ValueError("connected bootstrap groups must not span strata")
    empty_strata = tuple(stratum for stratum, groups in strata.items() if not groups)
    if empty_strata:
        raise ValueError(f"calibration bootstrap has empty declared bootstrap strata: {list(empty_strata)}")

    estimate = _macro_mean(strata)
    rng = random.Random(seed)
    samples: list[Fraction] = []
    ordered_strata = tuple(required_strata)
    for _ in range(iterations):
        stratum_means: list[Fraction] = []
        for stratum in ordered_strata:
            groups = strata[stratum]
            group_ids = tuple(sorted(groups))
            drawn = tuple(rng.choice(group_ids) for _ in group_ids)
            values = [value for group_id in drawn for value in groups[group_id]]
            stratum_means.append(sum(values, start=Fraction(0)) / len(values))
        samples.append(sum(stratum_means, start=Fraction(0)) / len(stratum_means))
    samples.sort()
    return BootstrapInterval(
        estimate=estimate,
        one_sided_lower=_percentile(samples, Fraction(5, 100)),
        two_sided_lower=_percentile(samples, Fraction(25, 1000)),
        two_sided_upper=_percentile(samples, Fraction(975, 1000)),
        one_sided_noninferiority_p_value=Fraction(1 + sum(sample < 0 for sample in samples), iterations + 1),
        iterations=iterations,
        resampling_unit_count=len(group_strata),
    )


def clopper_pearson_interval(events: int, total: int, *, confidence: Fraction = Fraction(95, 100)) -> BinomialInterval:
    """Compute a two-sided exact binomial interval without optional dependencies."""
    if isinstance(events, bool) or not isinstance(events, int):
        raise ValueError("binomial events must be an integer")
    if isinstance(total, bool) or not isinstance(total, int) or total <= 0 or not 0 <= events <= total:
        raise ValueError("binomial total must be positive and contain the event count")
    if not isinstance(confidence, Fraction) or not 0 < confidence < 1:
        raise ValueError("binomial confidence must be an exact Fraction between zero and one")
    alpha_tail = float((1 - confidence) / 2)
    lower_float = 0.0 if events == 0 else _bisect_probability(events, total, alpha_tail, upper_tail=True)
    upper_float = 1.0 if events == total else _bisect_probability(events, total, alpha_tail, upper_tail=False)
    return BinomialInterval(
        Fraction(events, total),
        Fraction(lower_float).limit_denominator(10**12),
        Fraction(upper_float).limit_denominator(10**12),
        events,
        total,
    )


def deterministic_curves(scores: Sequence[tuple[Fraction, bool]]) -> tuple[CurvePoint, ...]:
    """Return deterministic ROC/PR counts at every observed score threshold."""
    if not isinstance(scores, Sequence) or not scores:
        raise ValueError("calibration curves require scored labels")
    for score, label in scores:
        if not isinstance(score, Fraction) or not isinstance(label, bool):
            raise ValueError("calibration curve rows require Fraction scores and Boolean labels")
    positives = sum(label for _, label in scores)
    negatives = len(scores) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("calibration curves require both positive and negative labels")
    points: list[CurvePoint] = []
    for threshold in sorted({score for score, _ in scores}, reverse=True):
        true_positives = sum(score >= threshold and label for score, label in scores)
        false_positives = sum(score >= threshold and not label for score, label in scores)
        false_negatives = positives - true_positives
        true_negatives = negatives - false_positives
        predicted_positive = true_positives + false_positives
        points.append(
            CurvePoint(
                threshold,
                true_positives,
                false_positives,
                true_negatives,
                false_negatives,
                Fraction(true_positives, positives),
                Fraction(false_positives, negatives),
                Fraction(true_positives, predicted_positive),
                Fraction(true_positives, positives),
            )
        )
    return tuple(points)


def joint_surface(grid: Mapping[str, Sequence[object]]) -> tuple[tuple[tuple[str, object], ...], ...]:
    """Enumerate the complete deterministic cross-product of declared grid axes."""
    if not isinstance(grid, Mapping) or not grid:
        raise ValueError("calibration joint surface requires a non-empty grid")
    names = tuple(sorted(grid))
    axes: list[tuple[object, ...]] = []
    for name in names:
        values = grid[name]
        if not isinstance(values, Sequence) or isinstance(values, (str, bytes)) or not values:
            raise ValueError(f"calibration joint surface axis {name} must be non-empty")
        axes.append(tuple(values))
    return tuple(tuple(zip(names, combination, strict=True)) for combination in product(*axes))


def holm_adjust(p_values: Mapping[str, Fraction]) -> Mapping[str, Fraction]:
    """Apply deterministic Holm family-wise multiplicity correction."""
    if not isinstance(p_values, Mapping) or not p_values:
        raise ValueError("Holm correction requires named p-values")
    for name, value in p_values.items():
        if not isinstance(name, str) or not name or not isinstance(value, Fraction) or not 0 <= value <= 1:
            raise ValueError("Holm correction requires non-empty names and exact unit Fractions")
    ordered = sorted(p_values.items(), key=lambda item: (item[1], item[0]))
    adjusted: dict[str, Fraction] = {}
    previous = Fraction(0)
    family_size = len(ordered)
    for index, (name, value) in enumerate(ordered):
        current = min(Fraction(1), value * (family_size - index))
        previous = max(previous, current)
        adjusted[name] = previous
    return MappingProxyType({name: adjusted[name] for name in p_values})


def _macro_mean(strata: Mapping[str, Mapping[str, Sequence[Fraction]]]) -> Fraction:
    means = []
    for groups in strata.values():
        values = [value for group_values in groups.values() for value in group_values]
        means.append(sum(values, start=Fraction(0)) / len(values))
    return sum(means, start=Fraction(0)) / len(means)


def _percentile(values: Sequence[Fraction], probability: Fraction) -> Fraction:
    index = math.floor(float(probability * (len(values) - 1)))
    return values[index]


def _bisect_probability(events: int, total: int, tail: float, *, upper_tail: bool) -> float:
    low = 0.0
    high = 1.0
    for _ in range(100):
        midpoint = (low + high) / 2
        probability = (
            _binomial_probability_at_least(events, total, midpoint)
            if upper_tail
            else _binomial_probability_at_most(events, total, midpoint)
        )
        if upper_tail:
            if probability < tail:
                low = midpoint
            else:
                high = midpoint
        elif probability > tail:
            low = midpoint
        else:
            high = midpoint
    return (low + high) / 2


def _binomial_probability_at_least(events: int, total: int, probability: float) -> float:
    return sum(
        math.comb(total, count) * probability**count * (1 - probability) ** (total - count)
        for count in range(events, total + 1)
    )


def _binomial_probability_at_most(events: int, total: int, probability: float) -> float:
    return sum(
        math.comb(total, count) * probability**count * (1 - probability) ** (total - count)
        for count in range(events + 1)
    )
