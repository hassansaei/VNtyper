"""Independent arithmetic and denominator checks for total-length evaluation."""

from fractions import Fraction
from importlib import import_module

import pytest

pytestmark = pytest.mark.unit


def test_length_error_metrics_have_hand_computable_values():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    observations = (
        m.LengthObservation("sample-1", "group-1", "nominal", 100.0, 90.0, 110.0),
        m.LengthObservation("sample-2", "group-2", "nominal", 120.0, 130.0, 110.0),
    )
    result = m.calculate_length_metrics(observations)
    assert result.eligible_count == 2
    assert result.assessable_count == 2
    assert result.mae == 10
    assert result.median_absolute_error == 10
    assert result.rmse == 10
    assert result.bias == 0
    assert result.r_squared == 0
    assert result.baseline_mae == 10
    assert result.relative_mae_improvement == 0
    assert result.within_tolerance == Fraction(1)
    assert result.availability == Fraction(1)


def test_missing_prediction_counts_against_coverage_and_tolerance_not_as_zero():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    observations = (
        m.LengthObservation("sample-1", "group-1", "nominal", 100.0, 100.0, 110.0),
        m.LengthObservation("sample-2", "group-2", "nominal", 200.0, None, 110.0),
    )
    result = m.calculate_length_metrics(observations)
    assert result.eligible_count == 2
    assert result.assessable_count == 1
    assert result.mae == 0
    assert result.baseline_mae == 10  # Same assessable specimen, not a 50-repeat oracle advantage.
    assert result.relative_mae_improvement == 1
    assert result.availability == Fraction(1, 2)
    assert result.within_tolerance == Fraction(1, 2)
    assert result.r_squared is None


def test_all_missing_predictions_preserve_eligible_count_and_undefined_errors():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    result = m.calculate_length_metrics((m.LengthObservation("sample-1", "group-1", "nominal", 100, None, 90),))
    assert result.eligible_count == 1
    assert result.assessable_count == 0
    assert result.mae is result.median_absolute_error is result.rmse is result.bias is None
    assert result.r_squared is result.baseline_mae is result.relative_mae_improvement is None
    assert result.availability == result.within_tolerance == 0
    assert result.availability_lower == result.tolerance_lower == 0


def test_zero_baseline_error_never_divides_or_creates_improvement():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    result = m.calculate_length_metrics((m.LengthObservation("sample-1", "group-1", "nominal", 100, 100, 100),))
    assert result.baseline_mae == 0
    assert result.relative_mae_improvement is None
    assert result.r_squared is None


def test_tolerance_is_inclusive_and_uses_the_larger_absolute_or_relative_limit():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    values = [(50, 60), (200, 220), (200, 220.0001)]
    rows = tuple(
        m.LengthObservation(f"sample-{i}", f"group-{i}", "nominal", truth, prediction, 100)
        for i, (truth, prediction) in enumerate(values)
    )
    result = m.calculate_length_metrics(rows)
    assert result.within_tolerance == Fraction(2, 3)


def test_binomial_lower_bound_is_one_sided_not_two_sided():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    assert float(m.one_sided_binomial_lower(60, 60)) == pytest.approx(0.05 ** (1 / 60), abs=1e-10)
    assert m.one_sided_binomial_lower(0, 60) == 0
    # Independent beta-CDF closed form: for1 success of2, P(p)=1-(1-p)^2.
    assert float(m.one_sided_binomial_lower(1, 2)) == pytest.approx(1 - 0.95**0.5, abs=1e-10)


@pytest.mark.parametrize("confidence", [0.95, True, Fraction(1, 2), Fraction(1), Fraction(0)])
def test_binomial_confidence_requires_an_exact_valid_fraction(confidence):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    with pytest.raises(ValueError, match="confidence"):
        m.one_sided_binomial_lower(1, 2, confidence=confidence)


@pytest.mark.parametrize("successes,total", [(True, 2), (-1, 2), (3, 2), (0, 0), (1, True)])
def test_binomial_invalid_counts_are_rejected(successes, total):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    with pytest.raises(ValueError):
        m.one_sided_binomial_lower(successes, total)


@pytest.mark.parametrize(
    "field,value",
    [
        ("truth", None),
        ("truth", True),
        ("truth", 0),
        ("truth", float("nan")),
        ("prediction", float("inf")),
        ("prediction", 0),
        ("prediction", False),
        ("baseline_prediction", None),
        ("baseline_prediction", -1),
        ("key", ""),
        ("group_key", ""),
        ("stratum", ""),
    ],
)
def test_observations_refuse_invalid_values_before_metrics(field, value):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    row = {
        "key": "sample-1",
        "group_key": "group-1",
        "stratum": "nominal",
        "truth": 100,
        "prediction": 110,
        "baseline_prediction": 90,
    }
    row[field] = value
    with pytest.raises(ValueError, match=field):
        m.calculate_length_metrics((m.LengthObservation(**row),))


@pytest.mark.parametrize(
    "keys,groups",
    [(["sample-1", "sample-1"], ["group-1", "group-2"]), (["sample-1", "sample-2"], ["group-1", "group-1"])],
)
def test_primary_metrics_refuse_duplicate_specimens_or_group_representatives(keys, groups):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = tuple(
        m.LengthObservation(key, group, "nominal", 100, 110, 90) for key, group in zip(keys, groups, strict=True)
    )
    with pytest.raises(ValueError, match="unique"):
        m.calculate_length_metrics(rows)


def test_strata_preserve_all_eligible_observations_including_unavailable_ones():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (
        m.LengthObservation("sample-1", "group-1", "short", 100, 105, 90),
        m.LengthObservation("sample-2", "group-2", "long", 200, None, 90),
    )
    by_stratum = m.stratified_length_metrics(rows)
    assert tuple(by_stratum) == ("long", "short")
    assert by_stratum["long"].availability == 0
    assert by_stratum["long"].eligible_count == 1
    assert by_stratum["short"].mae == 5
    with pytest.raises(TypeError):
        by_stratum["other"] = by_stratum["short"]


def test_bootstrap_uses_paired_error_differences_and_is_order_independent():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = tuple(m.LengthObservation(f"sample-{i}", f"group-{i}", "nominal", 100, 101, 110) for i in range(3))
    interval = m.paired_length_error_interval(rows, seed=7)
    assert interval == (-9.0, -9.0)
    assert m.paired_length_error_interval(tuple(reversed(rows)), seed=7) == interval


def test_bootstrap_excludes_unavailable_pairs_without_calling_them_accurate():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (
        m.LengthObservation("sample-1", "group-1", "nominal", 100, 101, 110),
        m.LengthObservation("sample-2", "group-2", "nominal", 100, None, 110),
    )
    assert m.paired_length_error_interval(rows, seed=7) is None  # Need two assessable independent groups.
    assert m.calculate_length_metrics(rows).availability == Fraction(1, 2)


@pytest.mark.parametrize("seed", [True, -1, 1.5])
def test_bootstrap_rejects_invalid_seed(seed):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (m.LengthObservation("sample-1", "group-1", "nominal", 100, 101, 110),)
    with pytest.raises(ValueError, match="seed"):
        m.paired_length_error_interval(rows, seed=seed)


@pytest.mark.parametrize("rows", [(), [], None, ["not-observation"]])
def test_empty_or_untyped_observations_cannot_produce_successful_metrics(rows):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    with pytest.raises(ValueError):
        m.calculate_length_metrics(rows)


def test_bootstrap_retains_both_directions_of_paired_change():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (
        m.LengthObservation("sample-1", "group-1", "nominal", 100, 101, 110),
        m.LengthObservation("sample-2", "group-2", "nominal", 100, 110, 101),
    )
    # Two equally likely independent groups give means -9, 0, +9 with mass 1/4, 1/2, 1/4.
    assert m.paired_length_error_interval(rows, seed=7) == (-9, 9)
    assert m.paired_length_error_interval(tuple(reversed(rows)), seed=7) == (-9, 9)


def test_extreme_finite_errors_have_finite_median_without_intermediate_overflow():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = tuple(m.LengthObservation(f"sample-{i}", f"group-{i}", "nominal", 1, 1e308, 2) for i in range(2))
    result = m.calculate_length_metrics(rows)
    assert result.median_absolute_error == pytest.approx(1e308)
    assert result.mae == pytest.approx(1e308)
    assert result.rmse == pytest.approx(1e308)


def test_unrepresentable_integer_input_has_the_same_validation_error_as_infinity():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (m.LengthObservation("sample-1", "group-1", "nominal", 10**1000, 100, 100),)
    with pytest.raises(ValueError, match="truth"):
        m.calculate_length_metrics(rows)


def test_twenty_percent_improvement_is_not_rounded_below_its_gate():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (m.LengthObservation("sample-1", "group-1", "nominal", 100, 110, 112.5),)
    assert m.calculate_length_metrics(rows).relative_mae_improvement == 0.2


def test_tolerance_can_follow_a_separately_frozen_study_contract():
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (m.LengthObservation("sample-1", "group-1", "nominal", 100, 108, 112),)
    result = m.calculate_length_metrics(rows, tolerance_absolute=5, tolerance_relative=0.05)
    assert result.within_tolerance == 0
    assert result.availability == 1


@pytest.mark.parametrize("absolute,relative", [(True, 0.1), (0, 0.1), (10, float("nan")), (10, 2)])
def test_invalid_tolerance_limits_are_rejected(absolute, relative):
    m = import_module("vntyper.scripts.calibration_length_metrics")
    rows = (m.LengthObservation("sample-1", "group-1", "nominal", 100, 108, 112),)
    with pytest.raises(ValueError, match="tolerance"):
        m.calculate_length_metrics(rows, tolerance_absolute=absolute, tolerance_relative=relative)
