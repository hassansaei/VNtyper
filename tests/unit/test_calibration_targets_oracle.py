"""Invented-fixture comparisons against the independent calibration math oracle."""

from fractions import Fraction
from importlib import import_module
from pathlib import Path

import pytest

from tests.golden.oracle_import_guard import assert_independent_import_closure
from tests.unit.test_calibration_length import fit, training_row
from tests.unit.test_length_estimation import features
from vntyper.scripts.calibration_caller_curves import CallerOperatingPoint, build_caller_curves
from vntyper.scripts.calibration_caller_metrics import CallerObservation, calculate_caller_metrics, one_sided_fpr_upper
from vntyper.scripts.calibration_length_metrics import (
    LengthObservation,
    calculate_length_metrics,
    one_sided_binomial_lower,
)
from vntyper.scripts.calibration_statistics import clopper_pearson_interval
from vntyper.scripts.length_estimation import estimate_total_repeats

pytestmark = pytest.mark.unit


def oracle():
    return import_module("tests.golden.calibration_targets_oracle")


def test_oracle_recursive_import_closure_excludes_all_production_and_third_party_code():
    module = oracle()
    source = Path(module.__file__)
    assert assert_independent_import_closure(source, Path(__file__).parents[2]) == (source.resolve(),)


@pytest.mark.parametrize("feature", ["A", "F"])
@pytest.mark.parametrize("truth", [(61, 108, 161, 205), (80, 90, 120, 150)])
def test_real_length_fit_and_inference_match_exact_rational_normal_equations(feature, truth):
    x = (1, 2, 3, 4)
    expected = oracle().affine_fit(tuple(zip(x, truth, strict=True)))
    rows = tuple(
        training_row(str(i), value, total, feature=feature)
        for i, (value, total) in enumerate(zip(x, truth, strict=True))
    )
    result = fit(rows, f"affine-{feature}")
    model = result.outcomes[0].model
    assert model is not None
    assert model.intercept == pytest.approx(float(expected[0]))
    assert model.coefficients == pytest.approx((float(expected[1]),))
    assert result.baseline.mean_total_repeat_count == pytest.approx(float(expected[2]))
    # In-range held-out features with exact integer read depths, no refitting.
    for value in (1.5, 2.5, 3.5):
        depth = int(100 * value)
        depths = (100, depth, depth, 100, 100, 100) if feature == "A" else (100, depth, depth, depth, depth, 100)
        heldout = features(depths=depths, manifest_key="heldout")
        actual = estimate_total_repeats(heldout, model, evidence_domain="synthetic")
        assert actual.status == "estimated"
        assert actual.estimated_total_repeat_count == pytest.approx(float(expected[0] + expected[1] * Fraction(value)))


def test_oracle_affine_known_solution_and_singular_design_refusal():
    assert oracle().affine_fit(((1, 7), (2, 11), (3, 15))) == (3, 4, 11)
    for rows in ((), ((1, 7),), ((1, 7), (1, 8))):
        with pytest.raises(ValueError):
            oracle().affine_fit(rows)


@pytest.mark.parametrize(
    "values",
    [
        ((100, 80, 130), (120, 130, 130), (160, 164, 130), (300, None, 130)),
        ((100, None, 100),),
        ((100, 100, 100),),
        ((50, 60, 100), (200, 220, 100), (200, 220.0001, 100)),
    ],
)
def test_length_errors_and_baseline_use_same_assessable_rows_and_fixed_availability(values):
    expected = oracle().length_summary(values)
    actual = calculate_length_metrics(tuple(LengthObservation(str(i), str(i), *row) for i, row in enumerate(values)))
    for name, value in expected.items():
        observed = getattr(actual, name)
        assert observed is None if value is None else observed == pytest.approx(float(value)), name


def test_oracle_length_hand_arithmetic():
    result = oracle().length_summary(((100, 80, 130), (120, 130, 130), (160, 164, 130), (300, None, 130)))
    assert result["mae"] == Fraction(34, 3)
    assert result["bias"] == -2
    assert result["baseline_mae"] == Fraction(70, 3)
    assert result["availability"] == Fraction(3, 4)
    assert result["within_tolerance"] == Fraction(1, 2)
    with pytest.raises(ValueError):
        oracle().length_summary(())


def observations():
    # Three known-identity positives, one positive lacking identity, three negatives,
    # and unknown truth. No-calls remain in their eligible denominators.
    return (
        (True, True, ("a",), ("a",), ()),
        (True, True, ("a",), ("a", "b"), ("b",)),
        (True, None, ("a",), (), ()),
        (True, True, None, ("b",), ("b",)),
        (False, False, (), (), ()),
        (False, True, (), ("b",), ("b",)),
        (False, None, (), (), ()),
        (None, None, None, (), ()),
    )


def test_caller_fixed_denominators_exact_sets_and_missing_identity_match_oracle():
    raw = observations()
    expected = oracle().caller_summary(raw)
    actual = calculate_caller_metrics(
        tuple(
            CallerObservation(str(i), str(i), truth, ids, call, called, tier)
            for i, (truth, call, ids, called, tier) in enumerate(raw)
        )
    )
    for name, value in expected.items():
        observed = getattr(actual, name)
        if hasattr(observed, "estimate"):
            observed = observed.estimate
        assert observed == value, name
    assert expected["sensitivity"] == Fraction(3, 4)
    assert expected["specificity"] == expected["false_positive_rate"] == Fraction(1, 3)
    assert expected["exact_variant_recovery"] == Fraction(1, 3)
    assert expected["wrong_identity_groups"] == expected["wrong_tier_a_identity_groups"] == 2


def test_unknown_only_caller_population_has_undefined_truth_metrics():
    raw = ((None, None, None, (), ()),)
    expected = oracle().caller_summary(raw)
    assert expected["sensitivity"] is expected["specificity"] is expected["exact_variant_recovery"] is None
    assert expected["no_call_rate"] == 1
    with pytest.raises(ValueError):
        oracle().caller_summary(())


@pytest.mark.parametrize("comparison", ["<", "<=", ">", ">="])
def test_actual_curve_points_match_independent_cutoff_decisions_including_equality_and_nocalls(comparison):
    scores = (Fraction(1, 4), Fraction(1, 2), None, Fraction(1, 2), Fraction(3, 4), None)
    truths = (True, True, True, False, False, False)
    thresholds = (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
    expected = oracle().cutoff_points(truths, scores, thresholds, comparison)
    points = []
    # This adapter supplies invented frozen decisions to the actual curve engine.
    # It is a metrics/curve proof, not a claim of production read-capture replay.
    for i, (threshold, calls, summary) in enumerate(expected):
        rows = tuple(
            CallerObservation(str(j), str(j), truth, ("v",) if truth else (), call, ("v",) if call else (), ())
            for j, (truth, call) in enumerate(zip(truths, calls, strict=True))
        )
        points.append(CallerOperatingPoint(f"{i + 1:064x}", threshold, "a" * 64, rows))
        assert summary["no_calls"] == 2
    actual = build_caller_curves(tuple(reversed(points)), comparison=comparison, phase="policy-selection")
    assert len(actual.points) == len(thresholds)
    for point, (threshold, _, summary) in zip(actual.points, expected, strict=True):
        assert point.threshold == threshold
        assert point.sensitivity == summary["sensitivity"]
        assert point.false_positive_rate == summary["false_positive_rate"]
        assert point.precision == summary["precision"]


@pytest.mark.parametrize("events,total", [(0, 1), (1, 1), (0, 10), (1, 10), (5, 10), (9, 10), (10, 10), (3, 60)])
def test_native_binomial_bounds_match_direct_exact_polynomial_inversion(events, total):
    expected = oracle().binomial_interval(events, total, Fraction(1, 40))
    actual = clopper_pearson_interval(events, total)
    assert float(actual.lower) == pytest.approx(expected[0], abs=1e-10)
    assert float(actual.upper) == pytest.approx(expected[1], abs=1e-10)
    one_sided = oracle().binomial_interval(events, total, Fraction(1, 20))
    assert float(one_sided_fpr_upper(events, total)) == pytest.approx(one_sided[1], abs=1e-10)
    assert float(one_sided_binomial_lower(events, total)) == pytest.approx(one_sided[0], abs=1e-10)


def test_exact_polynomial_oracle_has_known_binomial_boundary_values():
    assert oracle().binomial_interval(0, 1, Fraction(1, 20)) == pytest.approx((0, 0.95))
    assert oracle().binomial_interval(1, 1, Fraction(1, 20)) == pytest.approx((0.05, 1))
    for events, total, alpha in (
        (-1, 2, Fraction(1, 20)),
        (3, 2, Fraction(1, 20)),
        (0, 0, Fraction(1, 20)),
        (1, 2, Fraction(0)),
    ):
        with pytest.raises(ValueError):
            oracle().binomial_interval(events, total, alpha)


@pytest.mark.parametrize(
    "truths,scores,comparison", [((True,), (), ">"), ((), (), ">"), ((True,), (Fraction(1),), "==")]
)
def test_cutoff_oracle_refuses_unmatched_or_empty_population_and_unsupported_comparator(truths, scores, comparison):
    with pytest.raises(ValueError):
        oracle().cutoff_points(truths, scores, (Fraction(1),), comparison)


@pytest.mark.parametrize("comparison,at_boundary", [("<", False), ("<=", True), (">", False), (">=", True)])
def test_cutoff_oracle_itself_pins_equality_instead_of_merely_agreeing_with_curve_adapter(comparison, at_boundary):
    points = oracle().cutoff_points(
        (True, False, True), (Fraction(1, 2), Fraction(1, 2), None), (Fraction(1, 2),), comparison
    )
    assert len(points) == 1
    threshold, calls, summary = points[0]
    assert threshold == Fraction(1, 2)
    assert calls == (at_boundary, at_boundary, None)
    assert summary["sensitivity"] == (Fraction(1, 2) if at_boundary else 0)
    assert summary["false_positive_rate"] == int(at_boundary)


@pytest.mark.parametrize(
    "events,total,tail",
    [(True, 2, Fraction(1, 20)), (1, True, Fraction(1, 20)), (1, 101, Fraction(1, 20)), (1, 2, 0.05)],
)
def test_binomial_oracle_refuses_unsupported_size_and_nonexact_inputs(events, total, tail):
    with pytest.raises(ValueError):
        oracle().binomial_interval(events, total, tail)
