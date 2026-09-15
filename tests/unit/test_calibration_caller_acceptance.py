"""Fixed caller safety gates and finite selection on invented independent groups."""

from dataclasses import replace
from fractions import Fraction

import pytest

from vntyper.scripts.calibration_caller_acceptance import (
    CallerGateRules,
    CallerSelectionEntry,
    evaluate_caller_acceptance,
    select_caller_candidate,
)
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_roster import decode_caller_eligible_roster

pytestmark = pytest.mark.unit


def _evaluate(candidate, baseline, roster, rules, *, phase="policy-selection"):
    return evaluate_caller_acceptance(candidate, baseline, roster, rules, phase=phase)


def _population():
    baseline = []
    candidate = []
    members = []
    for index in range(360):
        key = f"g{index:03}"
        truth = index < 60
        row = CallerObservation(key, key, truth, ("v",) if truth else (), truth, ("v",) if truth else (), ())
        candidate.append(row)
        # Baseline has two false positives, enough to meet the frozen benefit.
        baseline.append(replace(row, called_positive=True, called_variants=("v",)) if index in {60, 61} else row)
        members.append({"key": key, "group_key": key, "strata": ["nominal"]})
    return tuple(candidate), tuple(baseline), decode_caller_eligible_roster(members)


def _rules(**changes):
    return replace(CallerGateRules(("nominal",), 60, 299, 17), **changes)


def test_safe_candidate_passes_fixed_gates_and_demonstrates_selection_benefit():
    candidate, baseline, roster = _population()
    result = _evaluate(candidate, baseline, roster, _rules())
    assert result.status == "passed"
    assert result.selection_benefit is True
    assert result.pooled.sensitivity_interval == (Fraction(0), Fraction(0))
    assert result.pooled.no_call_increase == 0
    assert result.pooled.candidate.false_positive_rate.estimate == 0
    assert result.pooled.baseline.false_positive_rate.estimate == Fraction(1, 150)
    assert result.strata["nominal"].status == "passed"
    assert result.eligible_roster_sha256 == roster.sha256


def test_one_control_false_positive_fails_the_upper_bound_despite_high_sensitivity():
    candidate, baseline, roster = _population()
    candidate = (*candidate[:60], replace(candidate[60], called_positive=True, called_variants=("v",)), *candidate[61:])
    result = _evaluate(candidate, baseline, roster, _rules())
    assert result.status == "failed"
    assert "false_positive_upper_bound" in result.pooled.reasons


def test_suppress_all_cannot_pass_by_eliminating_false_positives():
    candidate, baseline, roster = _population()
    candidate = tuple(replace(row, called_positive=None, called_variants=()) for row in candidate)
    result = _evaluate(candidate, baseline, roster, _rules())
    assert result.status == "failed"
    assert {"sensitivity_noninferiority", "no_call_increase"} <= set(result.pooled.reasons)


def test_wrong_tier_a_identity_is_a_hard_failure_even_if_binary_detection_is_correct():
    candidate, baseline, roster = _population()
    candidate = (replace(candidate[0], called_variants=("wrong",), tier_a_variants=("wrong",)), *candidate[1:])
    result = _evaluate(candidate, baseline, roster, _rules())
    assert "wrong_tier_a_identity" in result.pooled.reasons
    assert result.status == "failed"


@pytest.mark.parametrize("truth_positive,index", [(True, 0), (None, 62)])
def test_candidate_tier_a_call_without_assessable_identity_is_insufficient_in_every_population(truth_positive, index):
    candidate, baseline, roster = _population()
    missing_truth = replace(
        candidate[index],
        truth_positive=truth_positive,
        truth_variants=None,
        called_positive=True,
        called_variants=("v",),
        tier_a_variants=("v",),
    )
    candidate = (*candidate[:index], missing_truth, *candidate[index + 1 :])
    baseline = (
        *baseline[:index],
        replace(
            baseline[index],
            truth_positive=truth_positive,
            truth_variants=None,
            called_positive=truth_positive is True,
            called_variants=("v",) if truth_positive is True else (),
            tier_a_variants=(),
        ),
        *baseline[index + 1 :],
    )

    result = _evaluate(candidate, baseline, roster, _rules())

    assert result.status == "insufficient-evidence"
    assert result.pooled.reasons == ("unassessable_tier_a_identity",)
    assert result.strata["nominal"].status == "insufficient-evidence"
    assert result.strata["nominal"].reasons == ("unassessable_tier_a_identity",)
    assert result.pooled.candidate.sensitivity.events == 60
    assert result.pooled.candidate.sensitivity.total == 60
    assert result.pooled.candidate.wrong_tier_a_identity_groups == 0


@pytest.mark.parametrize("truth_positive,index", [(True, 0), (None, 62)])
def test_missing_truth_identity_without_candidate_tier_a_call_remains_permitted(truth_positive, index):
    candidate, baseline, roster = _population()
    candidate = (
        *candidate[:index],
        replace(
            candidate[index],
            truth_positive=truth_positive,
            truth_variants=None,
            called_positive=True,
            called_variants=("v",),
        ),
        *candidate[index + 1 :],
    )
    baseline = (
        *baseline[:index],
        replace(
            baseline[index],
            truth_positive=truth_positive,
            truth_variants=None,
            called_positive=truth_positive is True,
            called_variants=("v",) if truth_positive is True else (),
        ),
        *baseline[index + 1 :],
    )

    result = _evaluate(candidate, baseline, roster, _rules())

    assert result.status == "passed"
    assert "unassessable_tier_a_identity" not in result.pooled.reasons


def test_baseline_only_unassessable_tier_a_call_does_not_penalize_candidate():
    candidate, baseline, roster = _population()
    candidate = (replace(candidate[0], truth_variants=None), *candidate[1:])
    baseline = (
        replace(baseline[0], truth_variants=None, tier_a_variants=("v",)),
        *baseline[1:],
    )

    result = _evaluate(candidate, baseline, roster, _rules())

    assert result.status == "passed"
    assert "unassessable_tier_a_identity" not in result.pooled.reasons


def test_unassessable_candidate_tier_a_identity_cannot_enter_selection():
    candidate, baseline, roster = _population()
    candidate = (
        replace(candidate[0], truth_variants=None, tier_a_variants=("v",)),
        *candidate[1:],
    )
    baseline = (replace(baseline[0], truth_variants=None), *baseline[1:])
    result = _evaluate(candidate, baseline, roster, _rules())

    assert result.selection_benefit is True
    assert result.status == "insufficient-evidence"
    assert select_caller_candidate((CallerSelectionEntry("a" * 64, 1, result),)) is None


def test_missing_or_sparse_strata_cannot_be_hidden_by_pooled_success():
    candidate, baseline, roster = _population()
    result = _evaluate(candidate, baseline, roster, _rules(required_strata=("absent", "nominal")))
    assert result.status == "insufficient-evidence"
    assert result.strata["absent"].reasons == ("missing_required_stratum",)
    result = _evaluate(candidate, baseline, roster, _rules(minimum_positive_groups=61))
    assert "insufficient_positive_groups" in result.pooled.reasons


def test_noop_candidate_is_not_selected_as_an_improvement():
    candidate, _, roster = _population()
    result = _evaluate(candidate, candidate, roster, _rules())
    assert result.status == "passed"
    assert result.selection_benefit is False
    assert select_caller_candidate((CallerSelectionEntry("a" * 64, 1, result),)) is None


def test_finite_selection_uses_safety_then_simplicity_and_stable_identity():
    candidate, baseline, roster = _population()
    result = _evaluate(candidate, baseline, roster, _rules())
    entries = (
        CallerSelectionEntry("c" * 64, 2, result),
        CallerSelectionEntry("b" * 64, 1, result),
        CallerSelectionEntry("a" * 64, 1, result),
    )
    assert select_caller_candidate(entries) == "a" * 64
    assert select_caller_candidate(tuple(reversed(entries))) == "a" * 64


def test_truth_or_population_changes_fail_before_statistical_evaluation():
    candidate, baseline, roster = _population()
    with pytest.raises(ValueError, match="roster"):
        _evaluate(candidate[:-1], baseline, roster, _rules())
    changed = (replace(baseline[0], truth_variants=("different",)), *baseline[1:])
    with pytest.raises(ValueError, match="truth"):
        _evaluate(candidate, changed, roster, _rules())


@pytest.mark.parametrize(
    "changes",
    [
        {"required_strata": ()},
        {"required_strata": ["nominal"]},
        {"minimum_positive_groups": True},
        {"minimum_negative_groups": 0},
        {"seed": -1},
        {"maximum_fpr_upper": 0.01},
        {"minimum_exact_benefit": Fraction(0)},
    ],
)
def test_acceptance_rules_are_strict_frozen_values(changes):
    with pytest.raises(ValueError):
        _rules(**changes)


def test_validation_outcomes_cannot_enter_policy_selection():
    candidate, baseline, roster = _population()
    result = _evaluate(candidate, baseline, roster, _rules(), phase="validation")
    with pytest.raises(ValueError, match="policy-selection"):
        select_caller_candidate((CallerSelectionEntry("a" * 64, 1, result),))


def test_invalid_phase_is_rejected_before_scoring():
    candidate, baseline, roster = _population()
    with pytest.raises(ValueError, match="phase"):
        _evaluate(candidate, baseline, roster, _rules(), phase="training")


def test_paired_interval_detects_sensitivity_loss_and_uses_inclusive_frozen_boundary():
    candidate, baseline, roster = _population()
    candidate = (replace(candidate[0], called_positive=False, called_variants=()), *candidate[1:])
    result = _evaluate(candidate, baseline, roster, _rules())
    assert result.pooled.sensitivity_interval == (Fraction(-3, 60), Fraction(0))
    assert "sensitivity_noninferiority" in result.pooled.reasons
    equal = _evaluate(candidate, baseline, roster, _rules(minimum_sensitivity_delta_lower=Fraction(-3, 60)))
    assert "sensitivity_noninferiority" not in equal.pooled.reasons


def test_no_call_increase_gate_uses_the_fixed_population_and_inclusive_limit():
    candidate, baseline, roster = _population()
    candidate = tuple(
        replace(row, called_positive=None, called_variants=()) if 62 <= index < 69 else row
        for index, row in enumerate(candidate)
    )
    result = _evaluate(candidate, baseline, roster, _rules(maximum_no_call_increase=Fraction(7, 360)))
    assert result.pooled.no_call_increase == Fraction(7, 360)
    assert "no_call_increase" not in result.pooled.reasons
    stricter = _evaluate(candidate, baseline, roster, _rules(maximum_no_call_increase=Fraction(6, 360)))
    assert "no_call_increase" in stricter.pooled.reasons


def test_missing_positive_and_negative_truth_rates_remain_insufficient():
    candidate, baseline, roster = _population()
    unknown = tuple(replace(row, truth_positive=None, truth_variants=None) for row in candidate)
    result = _evaluate(unknown, unknown, roster, _rules())
    assert result.status == "insufficient-evidence"
    assert "undefined_sensitivity_interval" in result.pooled.reasons
    assert "insufficient_negative_groups" in result.pooled.reasons
    assert result.selection_benefit is False


@pytest.mark.parametrize("mode", ["empty", "untyped", "id", "duplicate", "complexity", "bindings", "rules"])
def test_invalid_selection_contracts_fail(mode):
    candidate, baseline, roster = _population()
    result = _evaluate(candidate, baseline, roster, _rules())
    entry = CallerSelectionEntry("a" * 64, 1, result)
    entries = (entry,)
    if mode == "empty":
        entries = ()
    elif mode == "untyped":
        entries = (object(),)
    elif mode == "id":
        entries = (replace(entry, candidate_id="bad"),)
    elif mode == "duplicate":
        entries = (entry, entry)
    elif mode == "complexity":
        entries = (replace(entry, free_parameters=True),)
    elif mode == "bindings":
        entries = (
            entry,
            replace(entry, candidate_id="b" * 64, acceptance=replace(result, baseline_predictions_sha256="0" * 64)),
        )
    else:
        with pytest.raises(ValueError, match="CallerGateRules"):
            _evaluate(candidate, baseline, roster, object())
        return
    with pytest.raises(ValueError):
        select_caller_candidate(entries)
