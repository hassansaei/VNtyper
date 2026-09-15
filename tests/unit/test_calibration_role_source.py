"""Role source metadata is checked before any truth or result file is opened."""

from copy import deepcopy
from importlib import import_module
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_calibration_target_runs import run_document
from vntyper.scripts.calibration_manifest import connected_leakage_groups
from vntyper.scripts.calibration_target_contract import decode_target_study
from vntyper.scripts.calibration_target_runs import decode_target_runs

pytestmark = pytest.mark.unit


def source_fixture():
    study = decode_target_study(study_document())
    runs_raw = run_document()
    runs_raw["runs"][0]["manifest_key"] = "artifact-0"
    runs = decode_target_runs(runs_raw)
    raw = {
        "schema_version": "calibration-role-source-v2",
        "study_sha256": study.sha256,
        "role": "training",
        "roster": [
            {
                "key": "artifact-0",
                "group_key": connected_leakage_groups(study.partitions)["artifact-0"],
                "strata": ["all"],
            }
        ],
        "excluded": [],
        "evidence_domains": {"artifact-0": "synthetic"},
        "identities_by_key": {
            "artifact-0": [
                {"namespace": "physical-readset", "sha256": "b" * 64},
                {"namespace": "specimen", "sha256": "c" * 64},
            ]
        },
        "truth_asset": {"path": "/sealed/training-truth.json", "sha256": "d" * 64, "size_bytes": 17},
        "run_manifest_sha256": runs.sha256,
    }
    return study, runs, raw


def test_metadata_only_role_binding_and_roundtrip():
    module = import_module("vntyper.scripts.calibration_role_source")
    study, runs, raw = source_fixture()
    with patch("builtins.open", side_effect=AssertionError("outcome opened")):
        source = module.decode_role_source(raw, study=study, runs=runs, expected_role="training")
    assert module.role_source_document(source, study=study, runs=runs) == raw
    assert source.keys == ("artifact-0",)
    assert source.identities == (("physical-readset", "b" * 64), ("specimen", "c" * 64))
    assert source.truth_asset.sha256 == "d" * 64


@pytest.mark.parametrize(
    "change", ["role", "study", "runs", "group", "key", "physical", "domain", "missing-exclusion", "extra"]
)
def test_source_cannot_move_members_groups_or_commitments(change):
    module = import_module("vntyper.scripts.calibration_role_source")
    study, runs, raw = source_fixture()
    if change == "role":
        raw["role"] = "validation"
    elif change == "study":
        raw["study_sha256"] = "0" * 64
    elif change == "runs":
        raw["run_manifest_sha256"] = "0" * 64
    elif change == "group":
        raw["roster"][0]["group_key"] = "forged-group"
    elif change == "key":
        raw["roster"][0]["key"] = "artifact-2"
    elif change == "physical":
        raw["identities_by_key"]["artifact-0"][0]["sha256"] = "0" * 64
    elif change == "domain":
        raw["evidence_domains"]["artifact-0"] = "validated"
    elif change == "missing-exclusion":
        raw["excluded"] = [{"key": "artifact-0", "reason": "truth-unavailable"}]
    else:
        raw["extra"] = 1
    with pytest.raises(ValueError):
        module.decode_role_source(raw, study=study, runs=runs, expected_role="training")


def test_failed_measurement_stays_in_roster_and_is_not_a_truth_exclusion():
    module = import_module("vntyper.scripts.calibration_role_source")
    study, runs, raw = source_fixture()
    changed = deepcopy(raw)
    changed["excluded"] = [{"key": "another", "reason": "low-coverage"}]
    with pytest.raises(ValueError):
        module.decode_role_source(changed, study=study, runs=runs, expected_role="training")
