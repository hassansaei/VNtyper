"""Independent arithmetic examples for caller operating-point statistics."""

from dataclasses import replace
from fractions import Fraction

import pytest

from vntyper.scripts.calibration_caller_metrics import (
    CallerObservation,
    calculate_caller_metrics,
    one_sided_fpr_upper,
)

pytestmark = pytest.mark.unit


def _row(key: str, truth: bool | None, call: bool | None, /, **changes: object) -> CallerObservation:
    return replace(
        CallerObservation(
            key, key, truth, None if truth is None else (("v1",) if truth else ()), call, ("v1",) if call else (), ()
        ),
        **changes,
    )


def test_no_calls_stay_in_truth_denominators_and_are_not_negative_calls() -> None:
    rows = (
        _row("tp", True, True),
        _row("fn", True, False),
        _row("pn", True, None),
        _row("tn", False, False),
        _row("fp", False, True),
        _row("nn", False, None),
        _row("unknown", None, None),
    )
    result = calculate_caller_metrics(rows)
    assert (result.eligible_count, result.known_truth_count, result.unknown_truth_count) == (7, 6, 1)
    assert (result.true_positives, result.false_negatives, result.true_negatives, result.false_positives) == (
        1,
        1,
        1,
        1,
    )
    assert (result.positive_no_calls, result.negative_no_calls, result.no_calls) == (1, 1, 3)
    assert (
        result.sensitivity.estimate
        == result.specificity.estimate
        == result.false_positive_rate.estimate
        == Fraction(1, 3)
    )
    assert result.precision.estimate == Fraction(1, 2)
    assert result.no_call_rate.estimate == Fraction(3, 7)
    assert result.assessability.estimate == Fraction(4, 7)
    assert result.study_prevalence == Fraction(1, 2)
    assert result.exact_variant_recovery.estimate == Fraction(1, 3)


def test_variant_identity_is_scored_only_when_independently_known() -> None:
    result = calculate_caller_metrics(
        (
            _row("right", True, True),
            _row("wrong", True, True, called_variants=("v2",), tier_a_variants=("v2",)),
            _row("missing-truth-identity", True, True, truth_variants=None, called_variants=("v2",)),
            _row("negative", False, True, called_variants=("v2",), tier_a_variants=("v2",)),
        )
    )
    assert result.exact_variant_recovery.events == 1
    assert result.exact_variant_recovery.total == 2
    assert result.positive_truth_missing_identity == 1
    assert result.wrong_identity_groups == 2
    assert result.wrong_tier_a_identity_groups == 2
    assert result.identity_assessable_count == 3


def test_undefined_denominators_are_null_and_all_no_calls_remain_counted() -> None:
    result = calculate_caller_metrics((_row("one", None, None),))
    for rate in (
        result.sensitivity,
        result.specificity,
        result.false_positive_rate,
        result.precision,
        result.exact_variant_recovery,
    ):
        assert (rate.events, rate.total, rate.estimate, rate.lower, rate.upper) == (0, 0, None, None, None)
    assert result.fpr_one_sided_upper is None
    assert result.study_prevalence is None
    assert result.no_call_rate.estimate == 1


def test_fpr_bound_is_one_sided_and_straddles_the_one_percent_gate() -> None:
    assert float(one_sided_fpr_upper(0, 299)) == pytest.approx(1 - 0.05 ** (1 / 299), abs=1e-12)
    assert one_sided_fpr_upper(0, 299) < Fraction(1, 100) < one_sided_fpr_upper(0, 298)
    assert one_sided_fpr_upper(299, 299) == 1


@pytest.mark.parametrize(
    "events,total,alpha",
    [
        (True, 3, Fraction(1, 20)),
        (0, 0, Fraction(1, 20)),
        (4, 3, Fraction(1, 20)),
        (0, 3, 0.05),
        (0, 3, Fraction(0)),
        (0, 3, Fraction(1, 2)),
    ],
)
def test_invalid_fpr_bound_inputs_fail(events: object, total: object, alpha: object) -> None:
    with pytest.raises(ValueError):
        one_sided_fpr_upper(events, total, alpha=alpha)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "changes",
    [
        {"key": " "},
        {"group_key": ""},
        {"truth_positive": 1},
        {"called_positive": 0},
        {"truth_variants": []},
        {"truth_variants": ("v1", "v1")},
        {"truth_variants": ()},
        {"truth_variants": ("v2", "v1")},
        {"called_variants": ("",)},
        {"tier_a_variants": ("absent",)},
        {"called_positive": None},
        {"truth_positive": False},
        {"truth_positive": None},
    ],
)
def test_inconsistent_observations_fail(changes: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        calculate_caller_metrics((_row("one", True, True, **changes),))


@pytest.mark.parametrize(
    "rows",
    [
        (),
        "rows",
        [object()],
        (_row("same", True, True), _row("same", True, True)),
        (_row("one", True, True), _row("two", False, False, group_key="one")),
    ],
)
def test_invalid_or_repeated_primary_rosters_fail(rows: object) -> None:
    with pytest.raises(ValueError):
        calculate_caller_metrics(rows)  # type: ignore[arg-type]


def test_metrics_are_order_invariant_and_rate_uncertainty_is_central_95_percent() -> None:
    rows = (_row("a", True, True), _row("b", False, False))
    result = calculate_caller_metrics(rows)
    assert result == calculate_caller_metrics(tuple(reversed(rows)))
    assert float(result.sensitivity.lower) == pytest.approx(0.025)
    assert result.sensitivity.upper == 1
    assert float(result.fpr_one_sided_upper) == pytest.approx(0.95)
