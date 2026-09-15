"""Selection-only ROC/PR operating points from comparable production decisions."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import NoReturn

from vntyper.scripts.calibration_caller_metrics import CallerMetrics, CallerObservation, calculate_caller_metrics

logger = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class CallerOperatingPoint:
    """One replayed scalar cutoff; fixed_policy_sha256 binds all other policy values."""

    candidate_id: str
    threshold: Fraction
    fixed_policy_sha256: str
    observations: tuple[CallerObservation, ...]


@dataclass(frozen=True)
class CallerCurvePoint:
    """One measured cutoff with eligible-population ROC/PR coordinates."""

    candidate_id: str
    threshold: Fraction
    true_positives: int
    false_positives: int
    true_negatives: int
    false_negatives: int
    false_positive_rate: Fraction
    sensitivity: Fraction
    precision: Fraction | None


@dataclass(frozen=True)
class CallerCurves:
    """Observed points ordered strict-to-loose, with fixed no-calls kept visible.

    No endpoints or AUC are invented outside the measured cutoff grid. With
    no-calls present these are eligible-denominator operating curves; they do
    not silently condition performance on the assessable subset.
    """

    comparison: str
    fixed_policy_sha256: str
    eligible_count: int
    unknown_truth_count: int
    no_call_count: int
    positive_no_calls: int
    negative_no_calls: int
    points: tuple[CallerCurvePoint, ...]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object) -> None:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        _fail("caller curve identities must be lowercase SHA256 digests")


def _curve_point(point: CallerOperatingPoint, metrics: CallerMetrics) -> CallerCurvePoint:
    if metrics.sensitivity.estimate is None or metrics.false_positive_rate.estimate is None:
        _fail("caller curves require both positive and negative independent truth groups")
    return CallerCurvePoint(
        point.candidate_id,
        point.threshold,
        metrics.true_positives,
        metrics.false_positives,
        metrics.true_negatives,
        metrics.false_negatives,
        metrics.false_positive_rate.estimate,
        metrics.sensitivity.estimate,
        metrics.precision.estimate,
    )


def build_caller_curves(points: Sequence[CallerOperatingPoint], *, comparison: str, phase: str) -> CallerCurves:
    """Build a comparable scalar family's curve without reimplementing its caller.

    Args:
        points: Production replay/rerun decisions at predeclared scalar cutoffs.
            The producer must prove capture sufficiency and bind the complete
            eligible roster before supplying observations. Every point binds
            all non-threshold policy fields to the same fixed-policy digest.
        comparison: Exact production threshold comparator: <, <=, >, or >=.
        phase: Must be policy-selection; validation shows a frozen point only.

    Returns:
        Deterministically ordered ROC/PR points with explicit no-call counts.
        The recorded comparator documents the producer's semantics. This helper
        never computes calls from rounded p-values or applies a surrogate score.
        Multidimensional policies require operating-point tables instead.

    Raises:
        ValueError: If points change truth, eligible members, another policy
            dimension, or assessability, or violate scalar decision monotonicity.
    """
    if phase != "policy-selection":
        _fail("caller cutoff curves require policy-selection evidence")
    if not isinstance(comparison, str) or comparison not in {"<", "<=", ">", ">="}:
        _fail("caller curves require an exact scalar threshold comparator")
    if not isinstance(points, (tuple, list)) or not points:
        _fail("caller curves require non-empty typed operating points")
    candidate_ids: set[str] = set()
    thresholds: set[Fraction] = set()
    metrics_by_id: dict[str, CallerMetrics] = {}
    for point in points:
        if not isinstance(point, CallerOperatingPoint) or not isinstance(point.observations, tuple):
            _fail("caller curves require immutable CallerOperatingPoint values")
        _digest(point.candidate_id)
        _digest(point.fixed_policy_sha256)
        if not isinstance(point.threshold, Fraction):
            _fail("caller curve thresholds must be exact Fractions")
        if point.candidate_id in candidate_ids or point.threshold in thresholds:
            _fail("caller curve candidate identities and thresholds must be unique")
        candidate_ids.add(point.candidate_id)
        thresholds.add(point.threshold)
        metrics_by_id[point.candidate_id] = calculate_caller_metrics(point.observations)
    ordered = sorted(points, key=lambda point: point.threshold, reverse=comparison in {">", ">="})
    first = ordered[0]
    identities = {row.key: (row.group_key, row.truth_positive, row.truth_variants) for row in first.observations}
    unavailable = {row.key for row in first.observations if row.called_positive is None}
    previous_positive: set[str] = set()
    curve_points: list[CallerCurvePoint] = []
    for point in ordered:
        if point.fixed_policy_sha256 != first.fixed_policy_sha256:
            _fail("caller scalar curves cannot combine different fixed policies")
        observed_identities = {
            row.key: (row.group_key, row.truth_positive, row.truth_variants) for row in point.observations
        }
        if observed_identities != identities:
            _fail("caller curve points must have identical eligible rosters and truth")
        if {row.key for row in point.observations if row.called_positive is None} != unavailable:
            _fail("caller scalar curves cannot hide changes in assessability")
        positive = {row.key for row in point.observations if row.called_positive is True}
        if not previous_positive <= positive:
            _fail("caller decisions violate the declared scalar threshold direction")
        previous_positive = positive
        curve_points.append(_curve_point(point, metrics_by_id[point.candidate_id]))
    metrics = metrics_by_id[first.candidate_id]
    return CallerCurves(
        comparison,
        first.fixed_policy_sha256,
        metrics.eligible_count,
        metrics.unknown_truth_count,
        metrics.no_calls,
        metrics.positive_no_calls,
        metrics.negative_no_calls,
        tuple(curve_points),
    )
