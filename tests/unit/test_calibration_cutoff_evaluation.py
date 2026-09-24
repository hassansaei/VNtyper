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


def test_mutating_held_out_truth_cannot_change_its_fold_selection(monkeypatch):
    """Given a fold allocation, selection reads training truth only.

    Allocation is stratified by truth label, so relabelling a sample may legitimately move
    it to another fold. The allocation is therefore frozen at the original one here, which
    isolates the property under test: a held-out label never reaches its fold's selection.
    """
    from vntyper.scripts import calibration_cutoff_evaluation as module
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    original = arms()
    frozen = module.group_folds(
        {row.key: row.group_key for row in original["baseline"]},
        folds=2,
        seed=7,
        strata={row.group_key: row.truth_positive for row in original["baseline"]},
    )
    monkeypatch.setattr(module, "group_folds", lambda *args, **kwargs: dict(frozen))
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


def test_outer_folds_are_stratified_by_truth_so_no_fold_lacks_a_class():
    """27 negatives and 55 positives: every one of ten folds holds both classes."""
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    truth = (False,) * 27 + (True,) * 55
    records = {"baseline": rows(truth, truth)}
    for seed in (3, 7, 20260915):
        result = evaluate_cutoff_arms(records, spec=SearchSpec("balanced-accuracy"), folds=10, seed=seed)
        for fold in result["folds"]:
            held = set(fold["held_out_keys"])
            classes = {row.truth_positive for row in records["baseline"] if row.key in held}
            assert classes == {False, True}, (seed, fold["fold"])


def _leaky_arms():
    """``leaky`` is the best candidate everywhere, but only sample ``0`` produced its value."""
    values = arms()
    return {"baseline": values["baseline"], "leaky": values["candidate"]}


def test_a_candidate_contributed_only_by_held_out_samples_is_not_admissible_in_that_fold():
    """Held-out feature values must not create the breakpoints a fold selects from."""
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    result = evaluate_cutoff_arms(
        _leaky_arms(),
        spec=SearchSpec("balanced-accuracy"),
        folds=2,
        seed=7,
        contributors={"leaky": frozenset({"0"})},
    )
    by_holdout = {"0" in fold["held_out_keys"]: fold for fold in result["folds"]}

    assert by_holdout[True]["used_policy"] == "baseline"
    assert by_holdout[True]["admissible_candidates"] == 1
    assert by_holdout[True]["selection"]["eligible_candidates"] == 1
    assert by_holdout[False]["used_policy"] == "leaky"
    assert by_holdout[False]["admissible_candidates"] == 2
    # The full-data selection is not a fold, so it still searches the complete grid.
    assert result["final_selection"]["policy_id"] == "leaky"
    assert result["fold_admissibility"] == "training-observed-breakpoints"


def test_without_contributors_every_candidate_stays_admissible_in_every_fold():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    result = evaluate_cutoff_arms(_leaky_arms(), spec=SearchSpec("balanced-accuracy"), folds=2, seed=7)

    assert all(fold["used_policy"] == "leaky" for fold in result["folds"])
    assert all(fold["admissible_candidates"] == 2 for fold in result["folds"])
    assert result["fold_admissibility"] == "all-candidates"


def test_a_candidate_contributed_by_any_training_sample_stays_admissible():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    every = frozenset(str(index) for index in range(8))
    result = evaluate_cutoff_arms(
        _leaky_arms(), spec=SearchSpec("balanced-accuracy"), folds=2, seed=7, contributors={"leaky": every}
    )

    assert all(fold["used_policy"] == "leaky" for fold in result["folds"])


@pytest.mark.parametrize(
    "contributors",
    [
        {"unknown-policy": frozenset({"0"})},
        {"baseline": frozenset({"0"})},
        {"leaky": frozenset({"not-a-sample"})},
        {"leaky": frozenset()},
        {"leaky": ["0"]},
        [("leaky", frozenset({"0"}))],
    ],
)
def test_malformed_contributors_are_refused(contributors):
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    with pytest.raises(ValueError, match="contributors"):
        evaluate_cutoff_arms(
            _leaky_arms(), spec=SearchSpec("balanced-accuracy"), folds=2, seed=7, contributors=contributors
        )


def _fold_rows(values):
    return [(row.key, row.group_key, row.truth_positive) for row in values["baseline"]]


def test_outer_fold_assignments_are_the_folds_the_evaluation_publishes():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms, outer_fold_assignments

    values = _leaky_arms()
    for seed in (3, 7, 20260915):
        assignments = outer_fold_assignments(_fold_rows(values), folds=2, seed=seed)
        result = evaluate_cutoff_arms(values, spec=SearchSpec("balanced-accuracy"), folds=2, seed=seed)
        assert assignments == {row["key"]: row["fold"] for row in result["rows"]}


def test_outer_fold_assignments_are_empty_for_a_singleton_cohort():
    from vntyper.scripts.calibration_cutoff_evaluation import outer_fold_assignments

    assert outer_fold_assignments([("only", "only", True)], folds=2, seed=7) == {}


def _inventories_without_leaky_where(values, heldout_key, seed=7):
    """Admit ``leaky`` in every fold except the one holding ``heldout_key`` out."""
    from vntyper.scripts.calibration_cutoff_evaluation import outer_fold_assignments

    assignments = outer_fold_assignments(_fold_rows(values), folds=2, seed=seed)
    return {
        fold: frozenset() if fold == assignments[heldout_key] else frozenset({"leaky"})
        for fold in set(assignments.values())
    }


def test_a_candidate_absent_from_a_fold_inventory_is_never_used_in_that_fold():
    """``leaky`` scores best on every training set, yet the fold whose inventory lacks it cannot use it."""
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    values = _leaky_arms()
    result = evaluate_cutoff_arms(
        values,
        spec=SearchSpec("balanced-accuracy"),
        folds=2,
        seed=7,
        fold_inventories=_inventories_without_leaky_where(values, "0"),
    )
    by_holdout = {"0" in fold["held_out_keys"]: fold for fold in result["folds"]}

    assert by_holdout[True]["used_policy"] == "baseline"
    assert by_holdout[True]["admissible_candidates"] == 1
    assert by_holdout[False]["used_policy"] == "leaky"
    assert by_holdout[False]["admissible_candidates"] == 2
    held = set(by_holdout[True]["held_out_keys"])
    assert all(row["selected_policy"] == "baseline" for row in result["rows"] if row["key"] in held)
    # The full-data selection is descriptive and still searches every replayed candidate.
    assert result["final_selection"]["policy_id"] == "leaky"
    assert result["fold_admissibility"] == "training-derived-inventories"


def test_contributors_and_fold_inventories_are_mutually_exclusive():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    values = _leaky_arms()
    with pytest.raises(ValueError, match="cutoff evaluation takes contributors or fold inventories, not both"):
        evaluate_cutoff_arms(
            values,
            spec=SearchSpec("balanced-accuracy"),
            folds=2,
            seed=7,
            contributors={"leaky": frozenset({"0"})},
            fold_inventories=_inventories_without_leaky_where(values, "0"),
        )


@pytest.mark.parametrize(
    "mutate",
    [
        lambda inventories: {**inventories, 99: frozenset()},
        lambda inventories: dict(list(inventories.items())[:1]),
        lambda inventories: {},
        lambda inventories: {fold: frozenset({"baseline"}) for fold in inventories},
        lambda inventories: {fold: frozenset({"unknown-policy"}) for fold in inventories},
        lambda inventories: {fold: {"leaky"} for fold in inventories},
        lambda inventories: [(fold, ids) for fold, ids in inventories.items()],
    ],
)
def test_malformed_fold_inventories_are_refused(mutate):
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    values = _leaky_arms()
    inventories = mutate(_inventories_without_leaky_where(values, "0"))
    with pytest.raises(ValueError, match="fold inventor"):
        evaluate_cutoff_arms(
            values, spec=SearchSpec("balanced-accuracy"), folds=2, seed=7, fold_inventories=inventories
        )


def test_empty_fold_inventories_keep_the_unavailable_result_of_a_singleton_cohort():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms

    records = {"baseline": rows((True,), (True,))}
    plain = evaluate_cutoff_arms(records, spec=SearchSpec("balanced-accuracy"))
    inventoried = evaluate_cutoff_arms(records, spec=SearchSpec("balanced-accuracy"), fold_inventories={})

    assert inventoried["status"] == "unavailable" and not inventoried["cross_validation_available"]
    assert inventoried["fold_admissibility"] == "training-derived-inventories"
    assert {**inventoried, "fold_admissibility": None} == {**plain, "fold_admissibility": None}
    with pytest.raises(ValueError, match="fold inventor"):
        evaluate_cutoff_arms(records, spec=SearchSpec("balanced-accuracy"), fold_inventories={0: frozenset()})
