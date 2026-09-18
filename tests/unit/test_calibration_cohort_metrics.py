"""Independent arithmetic and leakage checks for exploratory cohort comparisons."""

import pytest

pytestmark = pytest.mark.unit


def test_group_folds_are_deterministic_and_never_split_relatives():
    from vntyper.scripts.calibration_cohort_metrics import group_folds

    groups = {"a": "family", "b": "family", "c": "other", "d": "third"}
    result = group_folds(groups, folds=3, seed=7)
    assert result["a"] == result["b"]
    assert len(set(result.values())) == 3
    assert result == group_folds(dict(reversed(tuple(groups.items()))), folds=3, seed=7)
    with pytest.raises(ValueError):
        group_folds(groups, folds=1, seed=7)


def test_length_metrics_preserve_missing_and_use_paired_training_baselines():
    from vntyper.scripts.calibration_cohort_metrics import length_metrics

    result = length_metrics([10.0, 20.0, 30.0], [12.0, 18.0, None], [15.0, 15.0, 15.0], ["a", "b", "c"], seed=7)
    assert result["eligible"] == 3
    assert result["predicted"] == 2
    assert result["availability"] == pytest.approx(2 / 3)
    assert result["mae"] == 2
    assert result["rmse"] == 2
    assert result["bias"] == 0
    assert result["baseline_mae_paired"] == 5
    assert result["paired_mae_delta"] == -3
    assert result["paired_mae_delta_interval"] == [-3, -3]
    assert result["r2"] == pytest.approx(0.84)


def test_empty_and_constant_length_targets_have_no_fabricated_r2():
    from vntyper.scripts.calibration_cohort_metrics import length_metrics

    assert length_metrics([], [], [], [], seed=1)["mae"] is None
    assert length_metrics([5.0, 5.0], [5.0, 5.0], [5.0, 5.0], ["a", "b"], seed=1)["r2"] is None
    with pytest.raises(ValueError):
        length_metrics([5.0], [float("nan")], [5.0], ["a"], seed=1)


def test_paired_caller_deltas_use_same_truth_denominators_and_keep_no_calls():
    from dataclasses import replace

    from vntyper.scripts.calibration_caller_metrics import CallerObservation
    from vntyper.scripts.calibration_cohort_metrics import paired_caller_differences

    baseline = tuple(CallerObservation(str(i), str(i), i < 2, None if i < 2 else (), False, (), ()) for i in range(4))
    candidate = tuple(replace(row, called_positive=True) if row.truth_positive else row for row in baseline)
    result = paired_caller_differences(baseline, candidate, seed=7)
    assert result["sensitivity"]["delta"] == 1
    assert result["sensitivity"]["interval"] == [1, 1]
    assert result["specificity"]["delta"] == 0
    assert result["no_call_rate"]["delta"] == 0
    with pytest.raises(ValueError):
        paired_caller_differences(baseline, candidate[:-1], seed=7)
