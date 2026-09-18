"""Complete-roster cutoff evaluation separates selection from held-out outcomes."""

from dataclasses import replace

import pytest

from tests.unit.test_calibration_cutoff_selection import rows
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_cutoff_selection import SearchSpec

pytestmark = pytest.mark.unit


def arms():
    truth = (True, True, True, True, False, False, False, False)
    return {
        "baseline": rows((True, False, True, False, True, False, True, False), truth),
        "candidate": rows((True, True, True, True, False, False, False, False), truth),
    }


def test_complete_evaluation_retains_actual_training_and_held_out_predictions():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    result = evaluate_cutoff_arms(arms(), spec=SearchSpec("balanced-accuracy"), folds=2, seed=7)
    assert result["status"] == "available" and result["cross_validation_available"]
    assert result["final_selection"]["policy_id"] == "candidate"
    assert result["baseline"]["counts"]["balanced_accuracy"] == 0.5
    assert result["held_out"]["counts"]["balanced_accuracy"] == 1
    assert result["held_out"]["exact"]["sensitivity"]["estimate"] == 1
    assert len(result["rows"]) == 8
    for fold in result["folds"]:
        assert not set(fold["training_keys"]) & set(fold["held_out_keys"])
        assert len(fold["training_keys"]) + len(fold["held_out_keys"]) == 8
        assert fold["selection"]["training_count"] == len(fold["training_keys"])
        assert fold["selection"]["policy_id"] == "candidate"
        assert fold["fallback_reason"] is None
    assert result["full_data_operating_points"]["candidate"]["false_positives"] == 0
    assert result["paired_differences"]["sensitivity"]["delta"] == 0.5


def test_mutating_held_out_truth_cannot_change_its_fold_selection():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    original = arms()
    first = evaluate_cutoff_arms(original, spec=SearchSpec("balanced-accuracy"), folds=2, seed=7)
    fold = first["folds"][0]
    held = set(fold["held_out_keys"])
    changed = {
        name: tuple(
            replace(row, truth_positive=None, truth_variants=None) if row.key in held else row for row in values
        )
        for name, values in original.items()
    }
    second = evaluate_cutoff_arms(changed, spec=SearchSpec("balanced-accuracy"), folds=2, seed=7)
    assert second["folds"][0] == fold
    assert second["baseline"]["counts"]["unknown_truth_count"] == len(held)


def test_unknown_truth_and_no_calls_remain_in_roster_and_denominators():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    records = rows((True, None, False, None)) + (CallerObservation("x", "x", None, None, None, (), ()),)
    result = evaluate_cutoff_arms({"baseline": records}, spec=SearchSpec("balanced-accuracy"), folds=2, seed=7)
    assert len(result["rows"]) == 5
    counts = result["held_out"]["counts"]
    assert counts["eligible_count"] == 5 and counts["unknown_truth_count"] == 1
    assert counts["sensitivity"] == counts["specificity"] == 0.5
    assert counts["no_calls"] == 3
    assert result["held_out"]["exact"]["negative_no_calls"] == 1


def test_singleton_has_no_cv_and_no_exportable_selection():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    result = evaluate_cutoff_arms({"baseline": rows((True,), (True,))}, spec=SearchSpec("balanced-accuracy"))
    assert result["status"] == "unavailable"
    assert not result["cross_validation_available"]
    assert result["final_selection"]["policy_id"] is None
    assert result["held_out"] is None and result["paired_differences"] is None
    assert result["rows"][0]["held_out"] is None and result["rows"][0]["selected_policy"] is None


def test_training_class_absence_explicitly_falls_back_to_baseline():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    result = evaluate_cutoff_arms(
        {"baseline": rows((True, False), (True, False))}, spec=SearchSpec("balanced-accuracy"), folds=2
    )
    assert result["final_selection"]["policy_id"] == "baseline"
    assert all(fold["fallback_reason"] == "training-truth-class-missing" for fold in result["folds"])
    assert all(row["selected_policy"] == "baseline" for row in result["rows"])


def test_unsatisfied_full_data_constraints_do_not_yield_a_selected_export():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    result = evaluate_cutoff_arms(
        {"baseline": rows((True, False, True, False))}, spec=SearchSpec("balanced-accuracy", min_specificity=1), folds=2
    )
    assert result["status"] == "unavailable"
    assert result["final_selection"]["policy_id"] is None
    assert result["final_selection"]["reason"] == "no-candidate-satisfies-constraints"
    assert len(result["rows"]) == 4


@pytest.mark.parametrize("mutation", ["missing", "truth", "group", "duplicate", "identity"])
def test_full_roster_and_truth_binding_precede_any_fold_selection(monkeypatch, mutation):
    from vntyper.scripts import calibration_cutoff_evaluation as module

    values = arms()
    candidate = list(values["candidate"])
    if mutation == "missing":
        candidate.pop()
    elif mutation == "truth":
        candidate[0] = replace(candidate[0], truth_positive=False, truth_variants=())
    elif mutation == "group":
        candidate[0] = replace(candidate[0], group_key="other")
    elif mutation == "duplicate":
        candidate.append(candidate[0])
    else:
        candidate[0] = replace(candidate[0], truth_variants=("invented-variant",))
    values["candidate"] = tuple(candidate)

    def forbidden(*args, **kwargs):
        raise AssertionError("selection ran before roster validation")

    monkeypatch.setattr(module, "select_cutoff_policy", forbidden)
    with pytest.raises(ValueError):
        module.evaluate_cutoff_arms(values, spec=SearchSpec("balanced-accuracy"))


@pytest.mark.parametrize("values", [{}, {"other": rows((True, False), (True, False))}, {"baseline": ()}])
def test_missing_baseline_or_empty_population_is_refused(values):
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    with pytest.raises(ValueError):
        evaluate_cutoff_arms(values, spec=SearchSpec("balanced-accuracy"))
