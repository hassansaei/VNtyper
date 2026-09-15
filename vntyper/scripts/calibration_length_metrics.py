"""Total-length errors and uncertainty with fixed independent-group denominators."""

from __future__ import annotations

import logging
import math
import random
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_statistics import clopper_pearson_interval

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LengthObservation:
    """One predeclared primary observation per specimen and independent group."""

    key: str
    group_key: str
    stratum: str
    truth: float
    prediction: float | None
    baseline_prediction: float


@dataclass(frozen=True)
class LengthMetrics:
    """Errors on paired assessable cases, availability on all eligible cases."""

    eligible_count: int
    assessable_count: int
    mae: float | None
    median_absolute_error: float | None
    rmse: float | None
    bias: float | None
    r_squared: float | None
    baseline_mae: float | None
    relative_mae_improvement: float | None
    within_tolerance: Fraction
    availability: Fraction
    tolerance_lower: Fraction
    availability_lower: Fraction


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _positive(value: object, field: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not 0 < value <= sys.float_info.max:
        _fail(f"length observation {field} must be finite and positive")


def _observations(rows: Sequence[LengthObservation]) -> tuple[LengthObservation, ...]:
    if not isinstance(rows, (tuple, list)) or not rows:
        _fail("length metrics require non-empty typed observations")
    keys: set[str] = set()
    groups: set[str] = set()
    for row in rows:
        if not isinstance(row, LengthObservation):
            _fail("length metrics require LengthObservation values")
        for field in ("key", "group_key", "stratum"):
            value = getattr(row, field)
            if not isinstance(value, str) or not value or value != value.strip():
                _fail(f"length observation {field} must be non-empty text")
        _positive(row.truth, "truth")
        _positive(row.baseline_prediction, "baseline_prediction")
        if row.prediction is not None:
            _positive(row.prediction, "prediction")
        if row.key in keys or row.group_key in groups:
            _fail("primary length observations must have unique specimen and independent-group representatives")
        keys.add(row.key)
        groups.add(row.group_key)
    return tuple(sorted(rows, key=lambda row: row.group_key))


def one_sided_binomial_lower(successes: int, total: int, *, confidence: Fraction = Fraction(95, 100)) -> Fraction:
    """Return an exact one-sided Clopper-Pearson lower proportion bound.

    A central interval with confidence 2*c-1 has lower-tail probability 1-c,
    exactly the required one-sided bound. The existing central-interval helper
    remains unchanged and is not relabeled as a one-sided 95% interval.

    Args:
        successes: Number of successful independent group representatives.
        total: Number of eligible independent group representatives.
        confidence: Exact confidence fraction strictly between one half and one.

    Returns:
        Lower confidence bound using the established exact-binomial inversion.

    Raises:
        ValueError: If counts or confidence are invalid.
    """
    if not isinstance(confidence, Fraction) or not Fraction(1, 2) < confidence < 1:
        _fail("one-sided confidence must be an exact fraction between one half and one")
    return clopper_pearson_interval(successes, total, confidence=2 * confidence - 1).lower


def calculate_length_metrics(rows: Sequence[LengthObservation]) -> LengthMetrics:
    """Measure fixed predictions without fitting or dropping unavailable cases.

    Args:
        rows: Predeclared independent primary observations with exact total truth.

    Returns:
        Paired error metrics and eligible-population proportion bounds. Undefined
        errors and R-squared are None, never zero or nonfinite JSON values.

    Raises:
        ValueError: If observations are empty, duplicated or numerically invalid.
    """
    observations = _observations(rows)
    count = len(observations)
    paired = tuple(row for row in observations if row.prediction is not None)
    available = len(paired)
    errors = [row.prediction - row.truth for row in paired if row.prediction is not None]
    absolute = [abs(error) for error in errors]
    within = sum(abs(error) <= max(10.0, 0.1 * row.truth) for row, error in zip(paired, errors, strict=True))
    mae = bias = baseline_mae = relative = rmse = med = r_squared = None
    if available:
        mae = math.fsum(error / available for error in absolute)
        bias = math.fsum(error / available for error in errors)
        med = _percentile(sorted(absolute), 0.5)
        rmse = math.hypot(*(error / math.sqrt(available) for error in errors))
        baseline_mae = math.fsum(abs(row.baseline_prediction - row.truth) / available for row in paired)
        if baseline_mae > 0:
            relative = 1 - mae / baseline_mae
            if not math.isfinite(relative):
                relative = None
        mean_truth = math.fsum(row.truth / available for row in paired)
        truth_norm = math.hypot(*(row.truth - mean_truth for row in paired))
        if truth_norm > 0:
            ratio = math.hypot(*errors) / truth_norm
            value = 1 - ratio * ratio
            if math.isfinite(value):
                r_squared = value
    return LengthMetrics(
        eligible_count=count,
        assessable_count=available,
        mae=mae,
        median_absolute_error=med,
        rmse=rmse,
        bias=bias,
        r_squared=r_squared,
        baseline_mae=baseline_mae,
        relative_mae_improvement=relative,
        within_tolerance=Fraction(within, count),
        availability=Fraction(available, count),
        tolerance_lower=one_sided_binomial_lower(within, count),
        availability_lower=one_sided_binomial_lower(available, count),
    )


def stratified_length_metrics(rows: Sequence[LengthObservation]) -> Mapping[str, LengthMetrics]:
    """Apply the same denominator contract independently to each observed stratum.

    Args:
        rows: Unique primary observations, including unavailable predictions.

    Returns:
        Immutable metrics keyed in increasing stratum order. Missing mandatory
        strata must additionally be rejected by the study's acceptance gate.

    Raises:
        ValueError: If any observation or representative is invalid.
    """
    grouped: dict[str, list[LengthObservation]] = {}
    for row in _observations(rows):
        grouped.setdefault(row.stratum, []).append(row)
    return MappingProxyType({name: calculate_length_metrics(grouped[name]) for name in sorted(grouped)})


def _percentile(values: list[float], fraction: float) -> float:
    position = (len(values) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def paired_length_error_interval(rows: Sequence[LengthObservation], *, seed: int) -> tuple[float, float] | None:
    """Bootstrap paired candidate-minus-baseline absolute errors by group.

    Args:
        rows: One frozen representative per independent group.
        seed: Frozen non-negative protocol seed.

    Returns:
        The 95% percentile interval from 10000 seeded paired group resamples, or
        None if fewer than two independent assessable observations are present.
        It is a validation interval only when the supplied evidence is validation.

    Raises:
        ValueError: If the seed or observations are invalid.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        _fail("length bootstrap seed must be a non-negative integer")
    observations = _observations(rows)
    differences = [
        abs(row.prediction - row.truth) - abs(row.baseline_prediction - row.truth)
        for row in observations
        if row.prediction is not None
    ]
    count = len(differences)
    if count < 2:
        return None
    generator = random.Random(seed)
    means = sorted(math.fsum(value / count for value in generator.choices(differences, k=count)) for _ in range(10_000))
    return _percentile(means, 0.025), _percentile(means, 0.975)
