"""Native caller policy arms preserve no-calls and fold-local selection."""

from dataclasses import replace

import pytest

from vntyper.scripts.calibration_cohort_manifest import CohortSample

pytestmark = pytest.mark.unit


def sample(tmp_path, key="a", truth=True):
    return CohortSample(key, tmp_path / (key + ".bam"), "GRCh38", truth, None, "sample:" + key)


@pytest.mark.parametrize("advntr_samples", [(), ("a",), ("a", "b")])
def test_run_root_discovers_advntr_and_keeps_partial_roster_unavailable(tmp_path, advntr_samples):
    from vntyper.scripts.calibration_cohort_callers import read_policy_arms

    samples = (sample(tmp_path, "a"), sample(tmp_path, "b"))
    for row in samples:
        root = tmp_path / "runs" / row.sample_id
        (root / "kestrel").mkdir(parents=True)
        (root / "kestrel/kestrel_result.tsv").write_text("Confidence\nHigh_Precision\n")
        if row.sample_id in advntr_samples:
            (root / "advntr").mkdir()
            (root / "advntr/output_adVNTR_result.tsv").write_text("VID\tState\n25561\tI1_2_C\n")
    _baseline, arms, metadata = read_policy_arms(samples, None, tmp_path / "runs")
    expected = ["kestrel", "advntr"] if advntr_samples else ["kestrel"]
    assert metadata["baseline"]["required_callers"] == expected
    assert [row.called_positive for row in arms["baseline"]] == [
        True if not advntr_samples or row.sample_id in advntr_samples else None for row in samples
    ]


def test_native_final_positive_negative_and_missing_are_distinct(tmp_path):
    from vntyper.scripts.calibration_cohort_callers import read_native_observation

    positive = tmp_path / "positive.tsv"
    positive.write_text("Confidence\nHigh_Precision\n")
    negative = tmp_path / "negative.tsv"
    negative.write_text(
        "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
    )
    row = sample(tmp_path)
    assert read_native_observation(row, positive, None, ("kestrel",)).called_positive is True
    assert read_native_observation(row, negative, None, ("kestrel",)).called_positive is False
    assert read_native_observation(row, tmp_path / "missing", None, ("kestrel",)).called_positive is None
    assert read_native_observation(row, positive, None, ("kestrel", "advntr")).called_positive is None


def test_training_selection_cannot_see_held_out_truth(tmp_path):
    from vntyper.scripts.calibration_caller_metrics import CallerObservation
    from vntyper.scripts.calibration_cohort_callers import choose_policy

    a = CallerObservation("a", "a", True, None, False, (), ())
    b = CallerObservation("b", "b", False, (), False, (), ())
    arms = {"baseline": (a, b), "better": (replace(a, called_positive=True), b)}
    assert choose_policy(arms, ("a", "b"), "baseline") == "better"
    # With one truth class only, no safety-benefit claim/selection is possible.
    assert choose_policy(arms, ("a",), "baseline") == "baseline"


def test_missing_baseline_inventory_is_refused(tmp_path):
    from vntyper.scripts.calibration_cohort_callers import read_policy_arms

    path = tmp_path / "policies.json"
    path.write_text(
        '{"schema_version":"cohort-caller-policies-v1","baseline":"absent","required_callers":["kestrel"],"training_scope":"fixed-before-cohort","policies":[]}'
    )
    with pytest.raises(ValueError):
        read_policy_arms((sample(tmp_path),), path, None)


def test_actual_policy_inventory_drives_group_held_out_comparison(tmp_path):
    import json

    from vntyper.scripts.calibration_cohort_callers import compare_caller_arms

    positive = tmp_path / "positive.tsv"
    positive.write_text("Confidence\nHigh_Precision\n")
    negative = tmp_path / "negative.tsv"
    negative.write_text(
        "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
    )
    samples = tuple(sample(tmp_path, f"s{i}", i % 2 == 0) for i in range(8))
    document = {
        "schema_version": "cohort-caller-policies-v1",
        "baseline": "baseline",
        "required_callers": ["kestrel"],
        "training_scope": "fixed-before-cohort",
        "policies": [],
    }
    for name, cutoff in (("baseline", 0.01), ("candidate", 0.005)):
        document["policies"].append(
            {
                "policy_id": name,
                "cutoff": cutoff,
                "comparison": ">=",
                "samples": {
                    row.sample_id: {
                        "kestrel_result": "positive.tsv" if name == "candidate" and row.genotype else "negative.tsv",
                        "advntr_result": None,
                    }
                    for row in samples
                },
            }
        )
    path = tmp_path / "policies.json"
    path.write_text(json.dumps(document))
    result = compare_caller_arms(samples, path, None, folds=4, seed=7)
    assert result["baseline"]["sensitivity"]["estimate"] == 0
    assert result["selected_out_of_fold"]["sensitivity"]["estimate"] == 1
    assert result["selected_out_of_fold"]["specificity"]["estimate"] == 1
    assert result["selected_out_of_fold"]["fpr_one_sided_upper"] > 0
    assert len(result["operating_points"]) == 2
    assert all(name == "candidate" for name in result["selected_by_fold"].values())
    document["policies"][1]["samples"].pop("s0")
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="exact sample roster"):
        compare_caller_arms(samples, path, None, folds=4, seed=7)


def test_explicit_caller_target_with_no_truth_is_unavailable(tmp_path):
    from vntyper.scripts.calibration_cohort_callers import compare_caller_arms

    positive = tmp_path / "positive.tsv"
    positive.write_text("Confidence\nHigh_Precision\n")
    row = replace(sample(tmp_path), genotype=None, kestrel_result=positive)
    result = compare_caller_arms((row,), None, None, folds=5, seed=7)
    assert result["status"] == "unavailable"


def test_dual_caller_report_exposes_individual_calls_and_explicit_union(tmp_path):
    from vntyper.scripts.calibration_cohort_callers import compare_caller_arms

    kestrel = tmp_path / "kestrel.tsv"
    kestrel.write_text(
        "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
    )
    advntr = tmp_path / "advntr.tsv"
    advntr.write_text("VID\tState\n25561\tI1_3_A\n")
    rows = tuple(
        replace(sample(tmp_path, str(i), True), kestrel_result=kestrel, advntr_result=advntr) for i in range(3)
    )
    result = compare_caller_arms(rows, None, None, folds=3, seed=7)
    assert result["composition"] == "native-positive union (OR); requires all declared callers assessable"
    assert result["per_caller"]["kestrel"]["baseline"]["true_positives"] == 0
    assert result["per_caller"]["advntr"]["baseline"]["true_positives"] == 3
    assert result["baseline"]["true_positives"] == 3
    assert len(result["predictions"]) == 3
    assert result["paired_differences"]["sensitivity"]["delta"] == 0
