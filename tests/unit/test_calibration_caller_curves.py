"""Selection-only curves derived from explicit production-replay decisions."""

from dataclasses import replace
from fractions import Fraction

import pytest

from vntyper.scripts.calibration_caller_curves import CallerOperatingPoint, build_caller_curves
from vntyper.scripts.calibration_caller_metrics import CallerObservation

pytestmark = pytest.mark.unit


def _point(threshold: int, calls: tuple[bool | None, ...]) -> CallerOperatingPoint:
    truths = (True, True, False, False)
    rows = tuple(
        CallerObservation(str(i), str(i), truth, ("v",) if truth else (), call, ("v",) if call else (), ())
        for i, (truth, call) in enumerate(zip(truths, calls, strict=True))
    )
    return CallerOperatingPoint(str(threshold) * 64, Fraction(threshold, 10), "a" * 64, rows)


def test_curve_points_use_replayed_calls_and_do_not_invent_threshold_endpoints() -> None:
    strict = _point(9, (True, False, False, False))
    loose = _point(2, (True, True, True, False))
    result = build_caller_curves((loose, strict), comparison=">=", phase="policy-selection")
    assert result.comparison == ">="
    assert result.fixed_policy_sha256 == "a" * 64
    assert result.eligible_count == 4
    assert result.no_call_count == 0
    assert [point.threshold for point in result.points] == [Fraction(9, 10), Fraction(1, 5)]
    assert [(point.false_positive_rate, point.sensitivity, point.precision) for point in result.points] == [
        (Fraction(0), Fraction(1, 2), Fraction(1)),
        (Fraction(1, 2), Fraction(1), Fraction(2, 3)),
    ]
    assert [
        (point.true_positives, point.false_positives, point.true_negatives, point.false_negatives)
        for point in result.points
    ] == [(1, 0, 2, 1), (2, 1, 1, 0)]


@pytest.mark.parametrize("comparison", ["<", "<="])
def test_lower_score_direction_orders_cutoffs_by_increasing_laxity(comparison: str) -> None:
    points = (_point(2, (False, False, False, False)), _point(9, (True, True, False, False)))
    result = build_caller_curves(tuple(reversed(points)), comparison=comparison, phase="policy-selection")
    assert [p.threshold for p in result.points] == [Fraction(1, 5), Fraction(9, 10)]
    assert result.points[0].precision is None
    assert result.points[1].sensitivity == 1


def test_fixed_no_calls_remain_visible_in_curve_denominators() -> None:
    points = (_point(9, (False, None, False, None)), _point(2, (True, None, True, None)))
    result = build_caller_curves(points, comparison=">", phase="policy-selection")
    assert result.no_call_count == 2
    assert result.positive_no_calls == result.negative_no_calls == 1
    assert result.points[-1].sensitivity == result.points[-1].false_positive_rate == Fraction(1, 2)
    assert result.points[-1].false_negatives == result.points[-1].true_negatives == 0


@pytest.mark.parametrize("phase", ["training", "validation", "locked-heldout", "development-assessment", None])
def test_no_search_curve_from_validation_or_other_outcomes(phase: object) -> None:
    with pytest.raises(ValueError, match="policy-selection"):
        build_caller_curves((_point(9, (True, False, False, False)),), comparison=">=", phase=phase)  # type: ignore[arg-type]


@pytest.mark.parametrize("mode", ["reversed-decision", "availability", "truth", "roster", "policy", "threshold", "id"])
def test_incomparable_or_inconsistent_operating_points_fail(mode: str) -> None:
    strict = _point(9, (True, False, False, False))
    loose = _point(2, (True, True, True, False))
    if mode == "reversed-decision":
        loose = _point(2, (False, True, True, False))
    elif mode == "availability":
        loose = _point(2, (True, None, True, False))
    elif mode == "truth":
        loose = replace(
            loose, observations=(replace(loose.observations[0], truth_variants=("other",)), *loose.observations[1:])
        )
    elif mode == "roster":
        loose = replace(loose, observations=loose.observations[1:])
    elif mode == "policy":
        loose = replace(loose, fixed_policy_sha256="b" * 64)
    elif mode == "threshold":
        loose = replace(loose, threshold=strict.threshold)
    else:
        loose = replace(loose, candidate_id=strict.candidate_id)
    with pytest.raises(ValueError):
        build_caller_curves((strict, loose), comparison=">=", phase="policy-selection")


@pytest.mark.parametrize(
    "changes", [{"threshold": 0.5}, {"candidate_id": "bad"}, {"fixed_policy_sha256": "bad"}, {"observations": []}]
)
def test_malformed_point_contracts_fail(changes: dict[str, object]) -> None:
    point = replace(_point(9, (True, False, False, False)), **changes)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        build_caller_curves((point,), comparison=">=", phase="policy-selection")


@pytest.mark.parametrize(
    "points,comparison",
    [((), ">"), ("bad", ">"), ([object()], ">"), ((_point(9, (True, False, False, False)),), "approx")],
)
def test_invalid_collections_and_comparators_fail(points: object, comparison: str) -> None:
    with pytest.raises(ValueError):
        build_caller_curves(points, comparison=comparison, phase="policy-selection")  # type: ignore[arg-type]


def test_single_truth_class_is_explicitly_unassessable_for_roc() -> None:
    point = _point(9, (True, False, False, False))
    point = replace(point, observations=point.observations[:2])
    with pytest.raises(ValueError, match="both positive and negative"):
        build_caller_curves((point,), comparison=">=", phase="policy-selection")
