"""One-axis ROC/PR curves, and labelled tables for everything that is not one axis.

A ROC or precision-recall curve asserts one claim: *holding every other decision fixed,
this is what moving one threshold does*. That holds only when exactly one parameter axis
varies. Joining points from different axes -- or from policies differing in more than the
axis value -- draws a line no single parameter can trace, and an operating point read off
it recovers a policy that was never replayed. :func:`build_axis_curve` therefore refuses
any candidate differing from the baseline at a pointer the axis does not name;
:func:`build_joint_points` reports multi-parameter results as a label-ordered table.

Curve construction is not reimplemented: :func:`build_axis_curve` assembles
``CallerOperatingPoint`` values for ``build_caller_curves``, which independently re-checks
that every point binds the same fixed-policy digest, carries an identical eligible roster
and truth, hides no change in assessability, and respects the declared threshold direction.

Every rate carries the exact two-sided Clopper-Pearson interval ``calculate_caller_metrics``
computes through ``clopper_pearson_interval``, plus the one-sided upper false-positive-rate
bound it takes from ``one_sided_fpr_upper``: a few dozen negative-truth groups with zero
observed false positives still admit a materially non-zero rate, which a bare ``0.0`` hides.
All stay exact :class:`fractions.Fraction` values until the document functions, and unknown
truth and no-calls are never re-derived.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields
from fractions import Fraction
from types import MappingProxyType
from typing import Final, NoReturn

from vntyper.scripts.calibration_caller_curves import CallerCurves, CallerOperatingPoint, build_caller_curves
from vntyper.scripts.calibration_caller_metrics import (
    CallerMetrics,
    CallerObservation,
    CallerRate,
    calculate_caller_metrics,
    validate_caller_observations,
)
from vntyper.scripts.calibration_caller_policy import CallerPolicyScalar
from vntyper.scripts.calibration_cohort_metrics import caller_metrics_document
from vntyper.scripts.calibration_cutoff_axes import AxisBreakpoints, axis_document
from vntyper.scripts.calibration_cutoff_grid import CutoffCandidate
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_CURVE_SCHEMA: Final[str] = "calibration-cutoff-curve-v1"
_JOINT_SCHEMA: Final[str] = "calibration-cutoff-joint-v1"
_ID_DIGITS: Final[int] = 4


@dataclass(frozen=True)
class AxisCurve:
    """ROC/PR over one axis, with every other decision value held fixed and hashed.

    Attributes:
        axis: Stable axis name, as published by ``calibration_cutoff_axes``.
        statistic: The measured quantity the axis thresholds, from the axis document.
        pointers: The decision-profile JSON pointers this axis moves together.
        curves: Operating points ordered strict-to-loose by ``build_caller_curves``.
        baseline_threshold: The shipped value, from the candidate whose changed
            ``parameters`` are empty.
        interval_by_threshold: ``repr(float(threshold))`` -> exact sensitivity,
            specificity and precision intervals plus the one-sided upper FPR bound.
        boundary_support: Truth-labelled samples whose decisive statistic lies inside the
            band the tested thresholds span, with its exact ends. Zero means no sample of
            that class constrains the cutoff, and the report layer must say so.
        rejected: The axis document's rejected ``(value, reason)`` breakpoints.
    """

    axis: str
    statistic: str
    pointers: tuple[str, ...]
    curves: CallerCurves
    baseline_threshold: float
    interval_by_threshold: Mapping[str, Mapping[str, object]]
    boundary_support: Mapping[str, object]
    rejected: tuple[tuple[float, str], ...]


@dataclass(frozen=True)
class JointOperatingPoint:
    """One labelled multi-parameter policy: a table row, never a point on a curve.

    Attributes:
        label: Caller-supplied name; the table is ordered by it and by nothing else.
        values: Every pointer whose value differs from the baseline policy, sorted.
        metrics: Exact counts and intervals over this point's own replayed cohort.
    """

    label: str
    values: Mapping[str, object]
    metrics: CallerMetrics


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _threshold_key(threshold: Fraction) -> str:
    """One threshold's document key: the round-trip repr of its float value."""
    return repr(float(threshold))


def _intervals(metrics: CallerMetrics) -> Mapping[str, object]:
    """Exact two-sided Clopper-Pearson rates, plus the one-sided upper FPR bound."""
    rates: dict[str, object] = {
        name: MappingProxyType({f.name: getattr(getattr(metrics, name), f.name) for f in fields(CallerRate)})
        for name in ("sensitivity", "specificity", "precision")
    }
    rates["false_positive_rate_one_sided_upper"] = metrics.fpr_one_sided_upper
    return MappingProxyType(rates)


def _json(value: object) -> object:
    """Convert exact Fractions to float, once, at the document boundary."""
    if isinstance(value, Mapping):
        return {str(key): _json(item) for key, item in value.items()}
    return float(value) if isinstance(value, Fraction) else value


def _index(axis: AxisBreakpoints, candidate: object) -> int:
    """The breakpoint a candidate replays, from the ``<axis>-NNNN`` id convention."""
    if not isinstance(candidate, CutoffCandidate):
        _fail("cutoff axis curves require CutoffCandidate values")
    suffix = candidate.candidate_id.removeprefix(f"{axis.axis}-")
    if suffix == candidate.candidate_id or len(suffix) != _ID_DIGITS or not suffix.isdigit():
        _fail(f"cutoff candidate {candidate.candidate_id} is not a breakpoint of axis {axis.axis}")
    position = int(suffix)
    if position >= len(axis.values):
        _fail(f"cutoff candidate {candidate.candidate_id} is outside the breakpoints of axis {axis.axis}")
    return position


def _fixed(axis: AxisBreakpoints, candidate: CutoffCandidate) -> Mapping[str, CallerPolicyScalar]:
    """Everything this candidate holds fixed: the policy minus the axis's own pointers."""
    moved = set(axis.pointers)
    return {pointer: value for pointer, value in candidate.policy.values.items() if pointer not in moved}


def _roster(rows: Sequence[CallerObservation]) -> tuple[tuple[str, str, bool | None], ...]:
    """The (key, group, truth) roster of one already-validated, group-ordered arm."""
    return tuple((row.key, row.group_key, row.truth_positive) for row in rows)


def _boundary_support(points: Sequence[CallerOperatingPoint]) -> Mapping[str, object]:
    """Count truth-labelled samples whose decisive statistic lies inside the tested band.

    One axis partitions each sample by one statistic, so a call changes across the tested
    range exactly when that statistic lies in the band. No-call rows never change --
    ``build_caller_curves`` refuses any change in assessability.
    """
    thresholds = [point.threshold for point in points]
    calls: dict[str, set[bool | None]] = {}
    for point in points:
        for row in point.observations:
            calls.setdefault(row.key, set()).add(row.called_positive)
    truth = {row.key: row.truth_positive for row in points[0].observations}
    moved = [truth[key] for key, seen in calls.items() if len(seen) > 1]
    return MappingProxyType(
        {
            "positives_within_band": moved.count(True),
            "negatives_within_band": moved.count(False),
            "band_low": min(thresholds),
            "band_high": max(thresholds),
        }
    )


def build_axis_curve(
    axis: AxisBreakpoints,
    candidates: Sequence[CutoffCandidate],
    arms: Mapping[str, Sequence[CallerObservation]],
    *,
    comparison: str,
    phase: str,
) -> AxisCurve:
    """Build ROC/PR for one varying axis with every other decision value held fixed.

    Args:
        axis: Breakpoints from ``calibration_cutoff_axes``: the statistic, the pointers the
            axis may move, and the values it could not test.
        candidates: Complete replayed policies with this axis's ``<axis>-NNNN`` ids.
            Exactly one must reproduce the shipped baseline, so the curve is readable
            against the operating point production actually runs.
        arms: Replayed observations keyed by candidate id, one per independent group.
        comparison: The exact production comparator: ``<``, ``<=``, ``>`` or ``>=``.
        phase: Must be ``policy-selection``; ``build_caller_curves`` enforces this.

    Returns:
        The axis curve with exact per-threshold intervals and boundary support. Counts,
        no-calls and unknown truth are exactly what ``calculate_caller_metrics`` reports.

    Raises:
        ValueError: If the axis is forged; if a candidate is not a breakpoint of this
            axis, is repeated, has no replayed arm, declares a different pointer
            inventory, or varies a pointer the axis does not name; if the arms disagree on
            the (key, group, truth) roster; if the baseline candidate is absent or
            ambiguous; or if ``build_caller_curves`` refuses the points.
    """
    document = axis_document(axis)
    if isinstance(candidates, str) or not isinstance(candidates, Sequence) or not candidates:
        _fail(f"cutoff axis {axis.axis} curves require a nonempty candidate sequence")
    if not isinstance(arms, Mapping):
        _fail(f"cutoff axis {axis.axis} curves require observations keyed by candidate identifier")
    indexed = sorted(((_index(axis, candidate), candidate) for candidate in candidates), key=lambda item: item[0])
    if len({position for position, _ in indexed}) != len(indexed):
        _fail(f"cutoff axis {axis.axis} curve repeats a breakpoint candidate")
    anchors = [(position, candidate) for position, candidate in indexed if not candidate.parameters]
    if len(anchors) != 1:
        _fail(
            f"cutoff axis {axis.axis} curves need exactly one baseline candidate, found "
            f"{sorted(candidate.candidate_id for _, candidate in anchors)}"
        )
    anchor_position, anchor = anchors[0]
    held = _fixed(axis, anchor)
    # Bind every non-axis decision value into the digest ``build_caller_curves`` re-checks.
    digest = canonical_sha256(
        {
            "schema_version": "calibration-cutoff-fixed-policy-v1",
            "required_callers": list(anchor.policy.required_callers),
            "axis_pointers": list(axis.pointers),
            "values": {pointer: held[pointer] for pointer in sorted(held)},
        }
    )
    rows: dict[str, tuple[CallerObservation, ...]] = {}
    for _, candidate in indexed:
        if candidate.candidate_id not in arms:
            _fail(f"cutoff candidate {candidate.candidate_id} has no replayed observations")
        rows[candidate.candidate_id] = validate_caller_observations(tuple(arms[candidate.candidate_id]))
    roster = _roster(rows[anchor.candidate_id])
    points: list[CallerOperatingPoint] = []
    intervals: dict[str, Mapping[str, object]] = {}
    for position, candidate in indexed:
        observations = rows[candidate.candidate_id]
        values = _fixed(axis, candidate)
        if set(values) != set(held):
            _fail(f"cutoff candidate {candidate.candidate_id} declares different policy pointers than the baseline")
        for pointer in sorted(values):
            if values[pointer] != held[pointer]:
                _fail(
                    f"cutoff candidate {candidate.candidate_id} varies off-axis pointer {pointer}; a curve needs "
                    f"exactly one moving axis, so this belongs in a joint operating-point table"
                )
        if _roster(observations) != roster:
            _fail(f"cutoff candidate {candidate.candidate_id} replays a different roster than the baseline")
        threshold = Fraction(axis.values[position])
        points.append(CallerOperatingPoint(candidate.policy.sha256, threshold, digest, observations))
        intervals[_threshold_key(threshold)] = _intervals(calculate_caller_metrics(observations))
    curves = build_caller_curves(tuple(points), comparison=comparison, phase=phase)
    return AxisCurve(
        axis=axis.axis,
        statistic=str(document["statistic"]),
        pointers=axis.pointers,
        curves=curves,
        baseline_threshold=float(axis.values[anchor_position]),
        interval_by_threshold=MappingProxyType(intervals),
        boundary_support=_boundary_support(points),
        rejected=axis.rejected,
    )


def build_joint_points(
    labelled: Mapping[str, tuple[CutoffCandidate, Sequence[CallerObservation]]],
    baseline: CutoffCandidate,
) -> tuple[JointOperatingPoint, ...]:
    """Build labelled joint operating points. Never a curve.

    Multi-parameter policies have no single axis to plot against, so each row names every
    pointer it moves away from the baseline and carries its own explicit denominators.

    Args:
        labelled: Label -> (complete replayed policy, its observations).
        baseline: The shipped policy every row's ``values`` are stated against.

    Returns:
        One point per label, ordered by label.

    Raises:
        ValueError: If the baseline, a label, a pair, a candidate or its observations
            are malformed, or a candidate declares a different pointer inventory.
    """
    if not isinstance(baseline, CutoffCandidate):
        _fail("joint operating points require a CutoffCandidate baseline")
    if not isinstance(labelled, Mapping) or not labelled:
        _fail("joint operating points require a nonempty labelled mapping")
    for label in labelled:
        if not isinstance(label, str) or not label or label.strip() != label:
            _fail("joint operating point labels must be nonempty trimmed strings")
    result: list[JointOperatingPoint] = []
    for label in sorted(labelled):
        entry = labelled[label]
        if not isinstance(entry, tuple) or len(entry) != 2:
            _fail(f"joint operating point {label} requires a candidate and its observations")
        candidate, rows = entry
        if not isinstance(candidate, CutoffCandidate):
            _fail(f"joint operating point {label} requires a CutoffCandidate")
        if isinstance(rows, str) or not isinstance(rows, Sequence):
            _fail(f"joint operating point {label} requires a sequence of observations")
        if set(candidate.policy.values) != set(baseline.policy.values):
            _fail(f"joint operating point {label} declares different policy pointers than the baseline")
        base = baseline.policy.values
        moved = {pointer: value for pointer, value in candidate.policy.values.items() if base[pointer] != value}
        result.append(
            JointOperatingPoint(
                label,
                MappingProxyType({pointer: moved[pointer] for pointer in sorted(moved)}),
                calculate_caller_metrics(tuple(rows)),
            )
        )
    return tuple(result)


def axis_curve_document(curve: AxisCurve) -> dict[str, object]:
    """Project one axis curve as fresh canonical ``calibration-cutoff-curve-v1`` JSON.

    Args:
        curve: A curve from :func:`build_axis_curve`.

    Returns:
        One row per tested threshold, ordered strict to loose. Exact Fractions become
        floats here and only here. The no-call and unknown-truth counts repeat on every
        row because ``build_caller_curves`` refuses any point that changes them.

    Raises:
        ValueError: If the argument is not an :class:`AxisCurve`.
    """
    if not isinstance(curve, AxisCurve):
        _fail("cutoff curve documents require an AxisCurve")
    curves = curve.curves
    rows = [
        {
            "candidate_id": point.candidate_id,
            "threshold": float(point.threshold),
            "true_positives": point.true_positives,
            "false_positives": point.false_positives,
            "true_negatives": point.true_negatives,
            "false_negatives": point.false_negatives,
            "positive_no_calls": curves.positive_no_calls,
            "negative_no_calls": curves.negative_no_calls,
            "unknown_truth_count": curves.unknown_truth_count,
            "false_positive_rate": float(point.false_positive_rate),
            "sensitivity": float(point.sensitivity),
            "precision": None if point.precision is None else float(point.precision),
            "intervals": _json(curve.interval_by_threshold[_threshold_key(point.threshold)]),
        }
        for point in curves.points
    ]
    return {
        "schema_version": _CURVE_SCHEMA,
        "axis": curve.axis,
        "statistic": curve.statistic,
        "pointers": list(curve.pointers),
        "comparison": curves.comparison,
        "fixed_policy_sha256": curves.fixed_policy_sha256,
        "baseline_threshold": curve.baseline_threshold,
        "thresholds": [float(point.threshold) for point in curves.points],
        "eligible_count": curves.eligible_count,
        "unknown_truth_count": curves.unknown_truth_count,
        "no_call_count": curves.no_call_count,
        "positive_no_calls": curves.positive_no_calls,
        "negative_no_calls": curves.negative_no_calls,
        "boundary_support": _json(curve.boundary_support),
        "rejected": [{"value": value, "reason": reason} for value, reason in curve.rejected],
        "points": rows,
    }


def joint_points_document(points: Sequence[JointOperatingPoint]) -> dict[str, object]:
    """Project labelled joint points as canonical ``calibration-cutoff-joint-v1`` JSON.

    Args:
        points: Points from :func:`build_joint_points`, in their label order.

    Returns:
        A JSON-compatible table: each row names the pointers it moves and carries the full
        metric document for its own cohort. There is no threshold column and no curve
        geometry, because no single parameter separates these rows.

    Raises:
        ValueError: If the argument is not a sequence of :class:`JointOperatingPoint`.
    """
    if isinstance(points, str) or not isinstance(points, Sequence):
        _fail("joint operating point documents require a sequence of joint points")
    rows: list[dict[str, object]] = []
    for point in points:
        if not isinstance(point, JointOperatingPoint):
            _fail("joint operating point documents require JointOperatingPoint values")
        rows.append(
            {
                "label": point.label,
                "values": {pointer: point.values[pointer] for pointer in sorted(point.values)},
                "metrics": caller_metrics_document(point.metrics),
            }
        )
    return {"schema_version": _JOINT_SCHEMA, "point_count": len(rows), "points": rows}
