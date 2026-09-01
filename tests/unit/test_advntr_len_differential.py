import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import scripts.advntr_len_differential as differential

pytestmark = pytest.mark.unit


def _clean_sweep_result() -> dict[str, Any]:
    return {
        "probes": 1,
        "identical": 1,
        "differing": 0,
        "unchanged_classes": {
            "names": list(differential.UNCHANGED_CLASSES),
            "probes": 1,
            "identical": 1,
            "violations": [],
        },
        "by_class": {"probes": {"no_len_token": 1}, "differing": {}},
        "trailing_ampersand": {"probes": 0, "differing": 0},
        "reporting_delta": {"gained": [], "lost": []},
        "sign_fix_delta": {
            "moved": 0,
            "lost": [],
            "gained": [],
            "residue_of_lost": {},
            "pure_states_moved": [],
        },
        "oracle": {
            "predicted": 0,
            "differed_but_not_predicted": [],
            "predicted_but_identical": [],
        },
        "model_problems": [],
        "differences": [],
    }


@pytest.mark.parametrize(
    ("state", "category", "differs"),
    [
        ("D50_2", "no_len_token", False),
        ("I9_2_A_LEN3", "terminal_len_single_part", False),
        ("D50_2&I9_2_A_LEN3", "terminal_len_compound", False),
        ("I9_2_A_LEN3&D50_2", "material_after_first_len", True),
        ("I9_2_A_LENX&I50_2_A_LEN3", "stray_len_literal", True),
        ("I9_2_A_LEN0&D50_2", "material_after_first_len", False),
    ],
)
def test_classification_and_oracle(state: str, category: str, differs: bool) -> None:
    assert differential.classify(state) == category
    assert differential.oracle_predicts_difference(state) is differs


def test_historic_insertion_length_only_parses_a_terminal_len_token() -> None:
    states = pd.Series(["I9_2_A_LEN3", "I9_2_A_LEN3&D50_2", "D50_2"])
    assert differential.historic_insertion_len(states).tolist() == [3, 0, 0]


def test_frame_series_comes_from_the_resolved_packaged_decision_profile() -> None:
    ins_frame, del_frame = differential.accepted_frames()

    assert len(ins_frame) == len(del_frame) == 100
    assert {"1", "4", "298"} <= ins_frame
    assert {"2", "5", "299"} <= del_frame


def test_signed_survival_keeps_insertion_deletion_and_zero_nets_disjoint() -> None:
    insertion_length = pd.Series([4, 0, 1])
    deletion_length = pd.Series([0, 2, 1])
    kept_ins, kept_del = differential.survival(insertion_length, deletion_length)
    assert kept_ins.tolist() == [True, False, False]
    assert kept_del.tolist() == [False, True, False]
    assert not (kept_ins & kept_del).any()


def test_absolute_survival_reproduces_the_retired_mixed_state_result() -> None:
    insertion_length = pd.Series([3])
    deletion_length = pd.Series([1])
    kept_ins, kept_del = differential.absolute_frame_survival(insertion_length, deletion_length)
    assert kept_ins.tolist() == [False]
    assert kept_del.tolist() == [True]


def test_real_generators_cover_every_changed_and_unchanged_state_family() -> None:
    states = differential.generate_states()
    unchanged_states = differential.unchanged_class_states()

    assert len(states) == len(set(states)) == 52_511
    assert len(unchanged_states) == len(set(unchanged_states)) == 13_277
    assert Counter(map(differential.classify, unchanged_states)) == {
        "no_len_token": 3_905,
        "terminal_len_single_part": 12,
        "terminal_len_compound": 9_360,
    }
    assert {
        "I9_2_A_LEN9&",
        "I9_2_A_LEN9&D50_2",
        "I9_2_A_LEN9&I50_2",
        "I9_2_A_LEN9&I50_2_A_LEN3",
        "I9_2_A_LENX&I50_2_A_LEN3",
    }.issubset(states)

    variants = pd.Series(states)
    old_len = differential.historic_insertion_len(variants)
    new_len = variants.map(differential.advntr.sum_insertion_lengths).astype(int)
    classes = variants.map(differential.classify)
    assert Counter(classes[old_len != new_len]) == {
        "material_after_first_len": 38_941,
        "stray_len_literal": 2,
    }
    unchanged_mask = classes.isin(differential.UNCHANGED_CLASSES)
    assert (old_len[unchanged_mask] == new_len[unchanged_mask]).all()


def test_production_cross_check_reports_both_mismatch_directions(monkeypatch: pytest.MonkeyPatch) -> None:
    frame = differential.state_frame(["I9_2_A_LEN3", "D50_2"])
    modelled_ins = pd.Series([True, False])
    modelled_del = pd.Series([False, False])

    monkeypatch.setattr(differential.advntr, "advntr_processing_ins", lambda batch: batch.loc[[1]])
    monkeypatch.setattr(differential.advntr, "advntr_processing_del", lambda batch: batch.iloc[0:0])

    assert differential.cross_check_against_production(frame, modelled_ins, modelled_del) == [
        "advntr_processing_ins: model and production disagree (model-only [0], production-only [1])"
    ]


def test_bounded_sweep_preserves_exhaustive_counts_and_caps_each_class(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    states = [
        "D50_2",
        "I9_2_A_LEN3",
        "D50_2&I9_2_A_LEN3",
        "I9_2_A_LEN3&D50_2",
        "I9_2_A_LENX&I50_2_A_LEN3",
        "I9_2_A_LEN0&D50_2",
    ]
    monkeypatch.setattr(differential, "generate_states", lambda: states)

    result = differential.sweep(max_examples=1)

    assert result["probes"] == 6
    assert result["differing"] == 2
    assert result["oracle"]["predicted"] == 2
    assert result["by_class"]["probes"] == {
        "no_len_token": 1,
        "terminal_len_single_part": 1,
        "terminal_len_compound": 1,
        "material_after_first_len": 2,
        "stray_len_literal": 1,
    }
    recorded_by_class = pd.Series([record["class"] for record in result["differences"]]).value_counts()
    assert recorded_by_class.le(1).all()


def test_bounded_sweep_records_exactly_the_cap_for_same_class_differences(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    states = [
        "I22_2_G_LEN1&D50_2",
        "I80_2_A_LEN2&D50_2",
        "I14_2_G_LEN14&D50_2",
    ]
    monkeypatch.setattr(differential, "generate_states", lambda: states)

    result = differential.sweep(max_examples=2)

    assert result["differing"] == 3
    assert [record["state"] for record in result["differences"]] == states[:2]


def test_main_writes_sorted_parseable_json_for_a_clean_sweep(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    result = _clean_sweep_result()
    output = tmp_path / "result.json"
    monkeypatch.setattr(differential, "sweep", lambda max_examples: result)

    assert differential.main(["--out", str(output)]) == 0
    parsed = json.loads(output.read_text())
    assert list(parsed) == sorted(result)
    assert parsed == result


def test_main_reports_every_regression_class(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    result = _clean_sweep_result()
    result["model_problems"] = ["production mismatch"]
    result["unchanged_classes"]["violations"] = ["unchanged"]
    result["oracle"]["differed_but_not_predicted"] = ["unexpected"]
    result["oracle"]["predicted_but_identical"] = ["missing"]
    result["sign_fix_delta"]["pure_states_moved"] = ["pure"]
    result["sign_fix_delta"]["residue_of_lost"] = {1: 1}
    monkeypatch.setattr(differential, "sweep", lambda max_examples: result)

    assert differential.main([]) == 1
    output = capsys.readouterr().out
    assert "REGRESSION" in output
    assert "production mismatch" in output
    assert "unchanged class differed: 'unchanged'" in output
    assert "differed but not predicted: 'unexpected'" in output
    assert "predicted but identical: 'missing'" in output
    assert "sign fix moved a pure state: 'pure'" in output
    assert "sign fix dropped a pathogenic-frame state (residue 1)" in output
