"""Cutoff objectives use training labels and fixed no-call denominators."""

from dataclasses import replace

import pytest

from vntyper.scripts.calibration_caller_metrics import CallerObservation

pytestmark = pytest.mark.unit


def rows(calls, truth=(True, True, False, False)):
    return tuple(
        CallerObservation(str(i), str(i), target, () if target is False else None, call, (), ())
        for i, (target, call) in enumerate(zip(truth, calls, strict=True))
    )


def test_counts_keep_no_calls_in_truth_denominators_and_zero_events_are_explicit():
    from vntyper.scripts.calibration_cutoff_selection import cutoff_counts, cutoff_counts_document

    result = cutoff_counts_document(cutoff_counts(rows((True, None, False, None))))
    assert result["true_positives"] == result["true_negatives"] == 1
    assert result["positive_no_calls"] == result["negative_no_calls"] == 1
    assert result["sensitivity"] == result["specificity"] == result["balanced_accuracy"] == 0.5
    assert result["false_positive_rate"] == 0
    assert result["precision"] == 1 and result["f1"] == pytest.approx(2 / 3)
    empty = cutoff_counts_document(cutoff_counts(()))
    assert empty["sensitivity"] is None and empty["precision"] is None


@pytest.mark.parametrize(
    "objective,expected", [("balanced-accuracy", "balanced"), ("sensitivity", "sensitive"), ("specificity", "specific")]
)
def test_objectives_choose_real_tradeoffs(objective, expected):
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    arms = {
        "baseline": rows((False, False, True, True)),
        "balanced": rows((True, True, True, False)),
        "sensitive": rows((True, True, True, True)),
        "specific": rows((False, False, False, False)),
    }
    # Sensitivity ties prefer fewer FP, so use an extra positive to make its benefit strict.
    arms = {
        name: row + (CallerObservation("4", "4", True, None, name == "sensitive", (), ()),)
        for name, row in arms.items()
    }
    choice = select_cutoff_policy(arms, ["0", "1", "2", "3", "4"], "baseline", SearchSpec(objective))
    assert choice.policy_id == expected


def test_constraints_exclude_ineligible_candidates_and_report_impossibility():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    arms = {"baseline": rows((True, False, False, False)), "sensitive": rows((True, True, True, False))}
    keys = ["0", "1", "2", "3"]
    assert (
        select_cutoff_policy(arms, keys, "baseline", SearchSpec("balanced-accuracy", min_specificity=1)).policy_id
        == "baseline"
    )
    assert (
        select_cutoff_policy(arms, keys, "baseline", SearchSpec("balanced-accuracy", min_sensitivity=1)).policy_id
        == "sensitive"
    )
    result = select_cutoff_policy(
        arms, keys, "baseline", SearchSpec("balanced-accuracy", min_sensitivity=1, min_specificity=1)
    )
    assert result.policy_id is None and result.reason == "no-candidate-satisfies-constraints"


def test_ties_use_counts_then_baseline_then_stable_id():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    keys = ["0", "1", "2", "3"]
    a = rows((True, False, False, False))
    assert (
        select_cutoff_policy(
            {"z": a, "baseline": a, "a": a}, keys, "baseline", SearchSpec("balanced-accuracy")
        ).policy_id
        == "baseline"
    )
    worse = rows((False, False, True, True))
    assert (
        select_cutoff_policy(
            {"z": a, "baseline": worse, "a": a}, keys, "baseline", SearchSpec("balanced-accuracy")
        ).policy_id
        == "a"
    )
    # Equal balanced accuracy, but fewer false positives precede more true positives.
    assert (
        select_cutoff_policy(
            {"baseline": rows((True, True, True, False)), "specific": a},
            keys,
            "baseline",
            SearchSpec("balanced-accuracy"),
        ).policy_id
        == "specific"
    )


def test_selection_cannot_see_held_out_truth_or_calls():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    arms = {"baseline": rows((True, False, False, True)), "candidate": rows((True, True, False, False))}
    first = select_cutoff_policy(arms, ["0", "2"], "baseline", SearchSpec("balanced-accuracy"))
    changed = {
        name: tuple(
            replace(r, truth_positive=None, truth_variants=None, called_positive=None) if r.key in {"1", "3"} else r
            for r in rr
        )
        for name, rr in arms.items()
    }
    assert select_cutoff_policy(changed, ["0", "2"], "baseline", SearchSpec("balanced-accuracy")) == first
    assert first.training_count == 2 and first.policy_id == "baseline"


@pytest.mark.parametrize("keys", [[], ["0", "1"], ["2", "3"]])
def test_missing_training_truth_class_returns_no_selection(keys):
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    result = select_cutoff_policy(
        {"baseline": rows((True, True, False, False))}, keys, "baseline", SearchSpec("balanced-accuracy")
    )
    assert result.policy_id is None and result.reason == "training-truth-class-missing"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"objective": "auc"},
        {"objective": "f1"},
        {"objective": "balanced-accuracy", "min_sensitivity": True},
        {"objective": "balanced-accuracy", "min_sensitivity": -0.1},
        {"objective": "balanced-accuracy", "min_specificity": 1.1},
        {"objective": "balanced-accuracy", "min_specificity": float("nan")},
        {"objective": "balanced-accuracy", "min_specificity": "0.5"},
    ],
)
def test_invalid_search_spec_is_refused(kwargs):
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec

    with pytest.raises(ValueError):
        SearchSpec(**kwargs)


def test_training_roster_truth_and_baseline_must_bind():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    base = rows((True, True, False, False))
    keys = [r.key for r in base]
    invalid = [
        ({"other": base}, keys),
        ({"baseline": base}, ["missing"]),
        ({"baseline": base}, ["0", "0"]),
        ({"baseline": base, "candidate": base[:-1]}, keys),
        (
            {"baseline": base, "candidate": (replace(base[0], truth_positive=False, truth_variants=()),) + base[1:]},
            keys,
        ),
    ]
    for arms, training in invalid:
        with pytest.raises(ValueError):
            select_cutoff_policy(arms, training, "baseline", SearchSpec("balanced-accuracy"))


def test_malformed_spec_types_and_huge_fraction_are_value_errors():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    for objective in ([], {}):
        with pytest.raises(ValueError):
            SearchSpec(objective=objective)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        SearchSpec("balanced-accuracy", min_sensitivity=10**400)
    with pytest.raises(ValueError, match="typed"):
        select_cutoff_policy({"baseline": rows((True, True, False, False))}, ["0", "2"], "baseline", None)  # type: ignore[arg-type]


def test_specificity_ties_prefer_more_true_positives_then_fewer_no_calls():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    keys = ["0", "1", "2", "3"]
    arms = {"baseline": rows((False, False, False, False)), "better": rows((True, False, False, False))}
    assert select_cutoff_policy(arms, keys, "baseline", SearchSpec("specificity")).policy_id == "better"
    arms = {"baseline": rows((True, None, False, False)), "complete": rows((True, False, False, False))}
    assert select_cutoff_policy(arms, keys, "baseline", SearchSpec("specificity")).policy_id == "complete"


def test_unknown_truth_contributes_availability_but_not_class_rates():
    from vntyper.scripts.calibration_cutoff_selection import cutoff_counts, cutoff_counts_document

    known = rows((True, True, False, False))
    unknown = CallerObservation("unknown", "unknown", None, None, None, (), ())
    metrics = cutoff_counts_document(cutoff_counts(known + (unknown,)))
    assert metrics["unknown_truth_count"] == 1 and metrics["eligible_count"] == 5
    assert metrics["sensitivity"] == metrics["specificity"] == 1
    assert metrics["no_call_rate"] == 0.2


def test_objective_has_no_default_and_must_be_declared_explicitly():
    from vntyper.scripts.calibration_cutoff_selection import OBJECTIVES, SearchSpec

    with pytest.raises(TypeError):
        SearchSpec()  # type: ignore[call-arg]
    assert OBJECTIVES == (
        "max-sensitivity-at-specificity",
        "youden-j",
        "max-f1",
        "balanced-accuracy",
        "sensitivity",
        "specificity",
    )
    for objective in OBJECTIVES:
        spec = SearchSpec(objective, min_specificity=0.5)
        assert spec.objective == objective


def test_max_sensitivity_at_specificity_requires_an_explicit_specificity_floor():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec

    with pytest.raises(ValueError, match="max-sensitivity-at-specificity.*min_specificity"):
        SearchSpec("max-sensitivity-at-specificity")
    with pytest.raises(ValueError, match="min_specificity"):
        SearchSpec("max-sensitivity-at-specificity", min_sensitivity=0.9)
    assert SearchSpec("max-sensitivity-at-specificity", min_specificity=0).min_specificity == 0


def test_max_sensitivity_at_specificity_maximizes_sensitivity_within_the_constraint():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    keys = ["0", "1", "2", "3"]
    arms = {
        "baseline": rows((False, False, False, False)),
        "greedy": rows((True, True, True, True)),
        "best": rows((True, True, True, False)),
        "mid": rows((True, False, False, False)),
    }
    spec = SearchSpec("max-sensitivity-at-specificity", min_specificity=0.5)
    assert select_cutoff_policy(arms, keys, "baseline", spec).policy_id == "best"
    strict = SearchSpec("max-sensitivity-at-specificity", min_specificity=1)
    assert select_cutoff_policy(arms, keys, "baseline", strict).policy_id == "mid"


def test_max_f1_and_balanced_accuracy_disagree_on_the_same_candidates():
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    keys = ["0", "1", "2", "3"]
    arms = {"baseline": rows((True, True, True, False)), "conservative": rows((True, False, False, False))}
    # Both arms score balanced accuracy 3/4, so the tie-break takes the fewer false positives.
    assert select_cutoff_policy(arms, keys, "baseline", SearchSpec("balanced-accuracy")).policy_id == "conservative"
    assert select_cutoff_policy(arms, keys, "baseline", SearchSpec("max-f1")).policy_id == "baseline"


def test_youden_j_is_reported_exactly_and_ranks_like_balanced_accuracy():
    from vntyper.scripts.calibration_cutoff_selection import (
        SearchSpec,
        cutoff_counts,
        cutoff_counts_document,
        select_cutoff_policy,
    )

    thirds = rows((True, False, None, True, True, False), truth=(True, True, True, False, False, False))
    document = cutoff_counts_document(cutoff_counts(thirds))
    # Fraction(1, 3) + Fraction(1, 3) - 1 is exactly -1/3; the float sum is -0.33333333333333337.
    assert document["sensitivity"] == document["specificity"] == 1 / 3
    assert document["youden_j"] == -1 / 3
    assert cutoff_counts_document(cutoff_counts(()))["youden_j"] is None

    keys = ["0", "1", "2", "3"]
    arms = {
        "baseline": rows((False, False, False, False)),
        "wide": rows((True, True, True, False)),
        "narrow": rows((True, False, False, False)),
    }
    for objective in ("youden-j", "balanced-accuracy"):
        assert select_cutoff_policy(arms, keys, "baseline", SearchSpec(objective)).policy_id == "narrow"
    arms["wide"] = rows((True, True, False, False))
    for objective in ("youden-j", "balanced-accuracy"):
        assert select_cutoff_policy(arms, keys, "baseline", SearchSpec(objective)).policy_id == "wide"


def test_tie_keys_break_ties_by_content_and_not_by_candidate_id():
    """A width-crossing ID (``-10001`` sorts before ``-9999``) must not decide a tie."""
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    arms = {
        "baseline": rows((False, False, False, False)),
        "axis-9999": rows((True, True, False, False)),
        "axis-10001": rows((True, True, False, False)),
    }
    keys, spec = ["0", "1", "2", "3"], SearchSpec("balanced-accuracy")

    assert select_cutoff_policy(arms, keys, "baseline", spec).policy_id == "axis-10001"
    tie_keys = {"axis-9999": "a" * 64, "axis-10001": "b" * 64}
    assert select_cutoff_policy(arms, keys, "baseline", spec, tie_keys=tie_keys).policy_id == "axis-9999"
    # The baseline keeps its preference over a tied candidate whatever the candidate's key.
    tied = {"baseline": arms["axis-9999"], "axis-9999": arms["axis-9999"]}
    assert select_cutoff_policy(tied, keys, "baseline", spec, tie_keys={"axis-9999": ""}).policy_id == "baseline"


@pytest.mark.parametrize("tie_keys", [{"axis-9999": "a"}, {"axis-9999": "a", "axis-10001": 3}, ["axis-9999"]])
def test_tie_keys_must_name_every_candidate_with_text(tie_keys):
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec, select_cutoff_policy

    arms = {"baseline": rows((False,) * 4), "axis-9999": rows((True,) * 4), "axis-10001": rows((True,) * 4)}
    with pytest.raises(ValueError, match="cutoff selection tie keys must map every candidate ID to text"):
        select_cutoff_policy(arms, ["0", "1", "2", "3"], "baseline", SearchSpec("sensitivity"), tie_keys=tie_keys)
