"""Closed run commitments preserve technical arms without opening their outcomes."""

from __future__ import annotations

from collections.abc import MutableMapping
from copy import deepcopy
from dataclasses import replace
from typing import cast

import pytest

pytestmark = pytest.mark.unit


def run_document(*, target: str = "length") -> dict:
    role = "length_features" if target == "length" else "kestrel_capture"
    return {
        "schema_version": "calibration-runs-v2",
        "target": target,
        "runs": [
            {
                "manifest_key": "artifact-001",
                "policy_sha256": "a" * 64,
                "input_sha256": "b" * 64,
                "producer_sha256": "c" * 64,
                "baseline_assets_sha256": "f" * 64,
                "capture_policy_sha256": "a" * 64,
                "execution_kind": "measurement" if target == "length" else "baseline-rerun",
                "exit_code": 0,
                **({"vntr_ids": []} if target == "callers" else {}),
                "assets": {role: {"path": "/invented inputs/feature.json", "sha256": "d" * 64, "size_bytes": 17}},
            }
        ],
    }


def test_manifest_is_immutable_and_keeps_same_artifact_policy_arms() -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs, target_runs_document

    raw = run_document()
    second = deepcopy(raw["runs"][0])
    second["policy_sha256"] = "e" * 64
    raw["runs"].append(second)
    decoded = decode_target_runs(raw)
    assert len(decoded.runs) == 2
    assert decoded.runs[0].input_sha256 == decoded.runs[1].input_sha256
    assert target_runs_document(decoded) == raw
    raw["runs"][0]["assets"]["length_features"]["sha256"] = "0" * 64
    assert decoded.runs[0].assets["length_features"].sha256 == "d" * 64
    with pytest.raises(TypeError):
        cast(MutableMapping, decoded.runs[0].assets)["other"] = decoded.runs[0].assets["length_features"]


@pytest.mark.parametrize(
    "field,value", [("target", "dominance"), ("schema_version", "calibration-runs-v1"), ("extra", 1)]
)
def test_wrong_schema_target_and_extra_root_fields_fail(field: str, value: object) -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs

    raw = run_document()
    raw[field] = value
    with pytest.raises(ValueError):
        decode_target_runs(raw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("manifest_key", ""),
        ("policy_sha256", "x" * 64),
        ("producer_sha256", None),
        ("exit_code", True),
        ("exit_code", -1),
    ],
)
def test_run_fields_cannot_hide_missing_identity_or_completion(field: str, value: object) -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs

    raw = run_document()
    raw["runs"][0][field] = value
    with pytest.raises(ValueError):
        decode_target_runs(raw)


def test_duplicates_reordering_and_changed_artifact_bytes_fail() -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs

    raw = run_document()
    raw["runs"].append(deepcopy(raw["runs"][0]))
    with pytest.raises(ValueError):
        decode_target_runs(raw)
    raw["runs"][1]["policy_sha256"] = "e" * 64
    raw["runs"][1]["input_sha256"] = "f" * 64
    with pytest.raises(ValueError):
        decode_target_runs(raw)
    raw["runs"][1]["input_sha256"] = "b" * 64
    raw["runs"].reverse()
    with pytest.raises(ValueError):
        decode_target_runs(raw)


@pytest.mark.parametrize(
    "field,value",
    [("path", "relative.json"), ("path", "/input/../file"), ("size_bytes", True), ("size_bytes", -1), ("sha256", "")],
)
def test_asset_commitments_are_closed_and_exact(field: str, value: object) -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs

    raw = run_document()
    raw["runs"][0]["assets"]["length_features"][field] = value
    with pytest.raises(ValueError):
        decode_target_runs(raw)


def test_target_and_required_caller_assets_are_not_inferred_from_available_files() -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs, require_run_assets

    raw = run_document(target="callers")
    run = decode_target_runs(raw).runs[0]
    with pytest.raises(ValueError):
        require_run_assets("callers", run, ("advntr", "kestrel"))
    raw["target"] = "length"
    with pytest.raises(ValueError):
        decode_target_runs(raw)
    raw = run_document()
    raw["runs"][0]["exit_code"] = 1
    run = decode_target_runs(raw).runs[0]
    with pytest.raises(ValueError, match="complete"):
        require_run_assets("length", run, ())


def test_role_selection_requires_every_requested_key_and_policy_without_io() -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs, select_target_runs

    manifest = decode_target_runs(run_document())
    assert select_target_runs(manifest, ("artifact-001",), ("a" * 64,)) == manifest.runs
    for keys, policies in ((("missing",), ("a" * 64,)), (("artifact-001",), ("e" * 64,)), ((), ("a" * 64,))):
        with pytest.raises(ValueError):
            select_target_runs(manifest, keys, policies)


def test_completed_assets_match_the_explicit_caller_and_execution_contract() -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs, require_run_assets

    raw = run_document(target="callers")
    row = raw["runs"][0]
    asset = row["assets"]["kestrel_capture"]
    row["assets"]["kestrel_result"] = deepcopy(asset)
    require_run_assets("callers", decode_target_runs(raw).runs[0], ("kestrel",))
    row["vntr_ids"] = [25561]
    for name in ("advntr_capture", "advntr_result", "advntr_model"):
        row["assets"][name] = deepcopy(asset)
    require_run_assets("callers", decode_target_runs(raw).runs[0], ("advntr", "kestrel"))
    with pytest.raises(ValueError, match="Kestrel-only"):
        require_run_assets("callers", decode_target_runs(raw).runs[0], ("kestrel",))
    row["execution_kind"] = "scalar-replay"
    with pytest.raises(ValueError, match="missing required"):
        require_run_assets("callers", decode_target_runs(raw).runs[0], ("advntr", "kestrel"))
    for name in ("advntr_replay_manifest", "advntr_replay_policy", "advntr_replay_result"):
        row["assets"][name] = deepcopy(asset)
    require_run_assets("callers", decode_target_runs(raw).runs[0], ("advntr", "kestrel"))
    require_run_assets("length", decode_target_runs(run_document()).runs[0], ())
    with pytest.raises(ValueError, match="no caller set"):
        require_run_assets("length", decode_target_runs(run_document()).runs[0], ("kestrel",))


def test_typed_digest_and_unauthorized_roster_changes_are_rejected() -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs, select_target_runs, target_runs_document

    manifest = decode_target_runs(run_document())
    with pytest.raises(ValueError, match="digest differs"):
        target_runs_document(replace(manifest, sha256="0" * 64))
    for keys, policies in (
        (("artifact-001", "artifact-001"), ("a" * 64,)),
        (("artifact-001",), ("a" * 64, "a" * 64)),
        (("artifact-001",), ("bad",)),
        ((" ",), ("a" * 64,)),
    ):
        with pytest.raises(ValueError):
            select_target_runs(manifest, keys, policies)


@pytest.mark.parametrize("ids", [[True], [0], [2, 1], [1, 1], "25561"])
def test_caller_vntr_roster_rejects_coercion_and_ambiguity(ids: object) -> None:
    from vntyper.scripts.calibration_target_runs import decode_target_runs

    raw = run_document(target="callers")
    raw["runs"][0]["vntr_ids"] = ids
    with pytest.raises(ValueError, match="VNTR roster"):
        decode_target_runs(raw)
