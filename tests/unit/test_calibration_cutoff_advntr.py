"""Tests for native adVNTR cutoff-grid replay orchestration."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import pytest

from tests.unit.test_advntr_calibration_policy import capabilities, capture_policy
from vntyper.modules.advntr.advntr_calibration_policy import advntr_canonical_sha256, decode_advntr_capabilities
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, decode_caller_policy_values
from vntyper.scripts.calibration_cutoff_advntr import evaluate_advntr_cutoff_grid

pytestmark = pytest.mark.unit


def _policy(*, cutoff: float = 0.001, support: int = 3, kestrel_floor: float = 0.5) -> CallerPolicyValues:
    values: dict[str, object] = {
        "/components/kestrel/alt_filtering/gg_depth_score_threshold": 0.5,
        "/components/kestrel/confidence_assignment/reporting_floor": kestrel_floor,
        "/components/kestrel/confidence_assignment/var_active_region_threshold": 1,
        "/components/kestrel/confidence_assignment/depth_score_thresholds/low": 0.2,
        "/components/kestrel/confidence_assignment/depth_score_thresholds/high": 0.8,
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 1,
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 2,
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": 3,
        "/components/advntr/calibrated_calling/mode": "legacy",
        "/components/advntr/calibrated_calling/cutoff": cutoff,
        "/components/advntr/calibrated_calling/minimum_read_support": support,
        "/components/advntr/calibrated_calling/rare_unit_fraction": None,
        "/components/advntr/calibrated_calling/adapter_filter": False,
        "/components/advntr/calibrated_calling/minimum_read_match_ratio": 0.6,
        "/components/advntr/calibrated_calling/prune_reverse": False,
    }
    return decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": ["advntr", "kestrel"],
            "values": values,
        }
    )


def _json_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii") + b"\n"


def _capture(
    path: Path, *, vntr_id: int, baseline: CallerPolicyValues, producer: dict[str, object] | None = None
) -> None:
    values = baseline.values
    prefix = "/components/advntr/calibrated_calling/"
    model_locus: dict[str, object] = {"id": vntr_id}
    document: dict[str, object] = {
        "schema_version": "advntr-frameshift-capture-v2",
        "completion": "completed-vntr",
        "producer": producer
        or {
            "package_version": "2.4.0",
            "build_id": "a" * 64,
            "source_revision": "b" * 40,
        },
        "assets": {
            "model_sha256": "c" * 64,
            "model_locus_sha256": advntr_canonical_sha256(model_locus),
            "background_sha256": None,
            "loaded_background_sha256": None,
        },
        "model_locus": model_locus,
        "loaded_background": None,
        "capture_policy": capture_policy(),
        "caller_policy": {
            "schema_version": "advntr-frameshift-policy-v1",
            "mode": values[f"{prefix}mode"],
            "cutoff": values[f"{prefix}cutoff"],
            "minimum_read_support": values[f"{prefix}minimum_read_support"],
        },
        "locus": {"vntr_id": vntr_id, "read_length": 150, "is_haploid": False, "selected_read_count": 10},
        "unit_geometry": {},
        "reference_order": [],
        "flank_boundaries": {},
        "warnings": [],
        "occurrences": [],
        "spans": [],
        "evidence_rows": [],
        "candidate_traversal": [],
        "decision_visits": [],
    }
    path.write_bytes(_json_bytes(document))


class _ReplayTool:
    def __init__(
        self,
        *,
        audit_key: str | None = None,
        baseline_parity: bool = True,
        replay_returncode: int = 0,
    ) -> None:
        self.replay_calls: list[tuple[str, ...]] = []
        self.audit_key = audit_key
        self.baseline_parity = baseline_parity
        self.replay_returncode = replay_returncode

    def __call__(self, argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if argv[-2:] == ("capabilities", "--json"):
            return subprocess.CompletedProcess(argv, 0, json.dumps(capabilities()), "")
        self.replay_calls.append(argv)
        if self.replay_returncode:
            return subprocess.CompletedProcess(argv, self.replay_returncode, "", "synthetic failure")
        manifest_path = Path(argv[argv.index("--manifest") + 1])
        policy_path = Path(argv[argv.index("--policy") + 1])
        capture_root = Path(argv[argv.index("--capture-root") + 1])
        output = Path(argv[argv.index("--output") + 1])
        manifest_raw = manifest_path.read_bytes()
        policy_raw = policy_path.read_bytes()
        manifest = json.loads(manifest_raw)
        policy = json.loads(policy_raw)
        results = []
        for row in manifest["captures"]:
            capture_raw = (capture_root / row["filename"]).read_bytes()
            records = []
            for line in capture_raw.splitlines():
                record = json.loads(line)
                vntr_id = record["locus"]["vntr_id"]
                audit = row["key"] == self.audit_key
                calls = [] if row["key"] != "sample-positive" else [{"state": "I22_2_G_LEN1"}]
                records.append(
                    {
                        "vntr_id": vntr_id,
                        "result": {
                            "schema_version": "advntr-frameshift-replay-result-v1",
                            "vntr_id": vntr_id,
                            "capture_record_sha256": advntr_canonical_sha256(record),
                            "policy_sha256": advntr_canonical_sha256(policy),
                            "capture_producer": {
                                "package_version": "2.4.0",
                                "build_id": "a" * 64,
                                "source_revision": "b" * 40,
                            },
                            "capture_assets": record["assets"],
                            "loaded_background_sha256": None,
                            "baseline_parity": self.baseline_parity,
                            "decision_visits": [],
                            "calls": calls,
                            "warnings": [],
                            "capture_audit": {
                                "attribution_outside_trials": ["synthetic-state"] if audit else [],
                                "calibrated_policy_domain_errors": [],
                            },
                        },
                    }
                )
            results.append({"key": row["key"], "capture_sha256": row["sha256"], "vntrs": records})
        document = {
            "schema_version": "advntr-frameshift-replay-output-v1",
            "manifest_file_sha256": hashlib.sha256(manifest_raw).hexdigest(),
            "manifest_sha256": advntr_canonical_sha256(manifest),
            "policy_file_sha256": hashlib.sha256(policy_raw).hexdigest(),
            "policy_sha256": advntr_canonical_sha256(policy),
            "background_file_sha256": None,
            "replay_producer": {
                "package_version": "2.4.0",
                "build_id": "a" * 64,
                "source_revision": "b" * 40,
            },
            "results": results,
        }
        output.mkdir()
        (output / "replay.json").write_bytes(_json_bytes(document))
        return subprocess.CompletedProcess(argv, 0, "", "")


def test_replays_each_unique_advntr_policy_once_and_keeps_full_results(tmp_path: Path) -> None:
    baseline = _policy()
    capture_root = tmp_path / "source captures"
    capture_root.mkdir()
    paths = {"sample-negative": capture_root / "negative.jsonl", "sample-positive": capture_root / "positive.jsonl"}
    for index, path in enumerate(paths.values(), start=17):
        _capture(path, vntr_id=index, baseline=baseline)
    tool = _ReplayTool()
    output = tmp_path / "private output"

    result = evaluate_advntr_cutoff_grid(
        paths,
        {
            "baseline": baseline,
            "kestrel-only-difference": _policy(kestrel_floor=0.4),
            "relaxed-cutoff": _policy(cutoff=0.005),
        },
        baseline_policy_id="baseline",
        executable_path=tmp_path / "tool with spaces" / "advntr",
        output=output,
        runner=tool,
    )

    assert len(tool.replay_calls) == 2
    by_policy = {row.policy_id: row for row in result.policies}
    assert by_policy["baseline"].execution_id == by_policy["kestrel-only-difference"].execution_id
    assert by_policy["relaxed-cutoff"].execution_id != by_policy["baseline"].execution_id
    assert {row.key: row.called_positive for row in by_policy["baseline"].samples} == {
        "sample-negative": False,
        "sample-positive": True,
    }
    assert json.loads(by_policy["baseline"].samples[1].loci[0].raw_json)["calls"] == [{"state": "I22_2_G_LEN1"}]
    assert (output / "grid-result.json").is_file()
    assert all(path.stat().st_mode & 0o077 == 0 for path in output.rglob("*"))


def test_audit_failure_is_unassessable_and_never_negative(tmp_path: Path) -> None:
    baseline = _policy()
    capture = tmp_path / "capture.jsonl"
    _capture(capture, vntr_id=17, baseline=baseline)

    result = evaluate_advntr_cutoff_grid(
        {"sample-audit": capture},
        {"baseline": baseline},
        baseline_policy_id="baseline",
        executable_path=tmp_path / "advntr",
        output=tmp_path / "result",
        runner=_ReplayTool(audit_key="sample-audit"),
    )

    sample = result.policies[0].samples[0]
    assert sample.assessable is False
    assert sample.called_positive is None


def test_distinct_samples_may_have_identical_complete_capture_bytes(tmp_path: Path) -> None:
    baseline = _policy()
    first = tmp_path / "first.jsonl"
    second = tmp_path / "second.jsonl"
    _capture(first, vntr_id=17, baseline=baseline)
    second.write_bytes(first.read_bytes())

    result = evaluate_advntr_cutoff_grid(
        {"sample-a": first, "sample-b": second},
        {"baseline": baseline},
        baseline_policy_id="baseline",
        executable_path=tmp_path / "advntr",
        output=tmp_path / "result",
        runner=_ReplayTool(),
    )

    assert [(row.key, row.called_positive) for row in result.policies[0].samples] == [
        ("sample-a", False),
        ("sample-b", False),
    ]


def test_refuses_capture_baseline_or_producer_mismatch_before_tool_execution(tmp_path: Path) -> None:
    baseline = _policy()
    capture = tmp_path / "capture.jsonl"
    _capture(capture, vntr_id=17, baseline=_policy(cutoff=0.005))
    tool = _ReplayTool()

    with pytest.raises(ValueError, match="baseline caller policy"):
        evaluate_advntr_cutoff_grid(
            {"sample": capture},
            {"baseline": baseline},
            baseline_policy_id="baseline",
            executable_path=tmp_path / "advntr",
            output=tmp_path / "result",
            runner=tool,
        )
    assert tool.replay_calls == []
    assert not (tmp_path / "result").exists()

    _capture(capture, vntr_id=17, baseline=baseline)
    second = tmp_path / "second.jsonl"
    _capture(
        second,
        vntr_id=17,
        baseline=baseline,
        producer={"package_version": "2.4.0", "build_id": "d" * 64, "source_revision": "b" * 40},
    )
    with pytest.raises(ValueError, match="producer"):
        evaluate_advntr_cutoff_grid(
            {"one": capture, "two": second},
            {"baseline": baseline},
            baseline_policy_id="baseline",
            executable_path=tmp_path / "advntr",
            output=tmp_path / "other-result",
            runner=tool,
        )


def test_refuses_nonlegacy_or_capture_changing_candidates(tmp_path: Path) -> None:
    baseline = _policy()
    capture = tmp_path / "capture.jsonl"
    _capture(capture, vntr_id=17, baseline=baseline)
    changed = _policy()
    changed_values = dict(changed.values)
    changed_values["/components/advntr/calibrated_calling/prune_reverse"] = True
    changed_document: dict[str, object] = {
        "schema_version": "calibration-caller-policy-values-v1",
        "required_callers": list(changed.required_callers),
        "values": changed_values,
    }
    changed_policy = decode_caller_policy_values(changed_document)

    with pytest.raises(ValueError, match="only legacy cutoff and support"):
        evaluate_advntr_cutoff_grid(
            {"sample": capture},
            {"baseline": baseline, "changed-capture": changed_policy},
            baseline_policy_id="baseline",
            executable_path=tmp_path / "advntr",
            output=tmp_path / "result",
            runner=_ReplayTool(),
        )


@pytest.mark.parametrize(
    ("tool", "error_type", "message"),
    [
        (_ReplayTool(baseline_parity=False), ValueError, "baseline parity"),
        (_ReplayTool(replay_returncode=9), RuntimeError, "process failed"),
    ],
)
def test_failed_native_replay_never_publishes_partial_results(
    tmp_path: Path,
    tool: _ReplayTool,
    error_type: type[Exception],
    message: str,
) -> None:
    baseline = _policy()
    capture = tmp_path / "capture.jsonl"
    _capture(capture, vntr_id=17, baseline=baseline)
    output = tmp_path / "result"

    with pytest.raises(error_type, match=message):
        evaluate_advntr_cutoff_grid(
            {"sample": capture},
            {"baseline": baseline},
            baseline_policy_id="baseline",
            executable_path=tmp_path / "advntr",
            output=output,
            runner=tool,
        )

    assert len(tool.replay_calls) == 1
    assert not output.exists()


def test_the_advntr_signature_ignores_kestrel_values_and_names_every_advntr_value() -> None:
    """One public signature decides which policies share a native execution (grid and guard alike)."""
    from vntyper.scripts.calibration_cutoff_advntr import advntr_signature

    baseline = _policy()

    assert advntr_signature(baseline) == advntr_signature(_policy(kestrel_floor=0.25))
    assert advntr_signature(baseline) != advntr_signature(_policy(cutoff=0.002))
    assert advntr_signature(baseline) != advntr_signature(_policy(support=4))
    assert len(advntr_signature(baseline)) == 7


def _grid(executions: dict[str, str]) -> Any:
    from vntyper.scripts.calibration_cutoff_advntr import AdvntrCutoffGridResult, AdvntrCutoffPolicyResult

    rows = tuple(AdvntrCutoffPolicyResult(pid, "c" * 64, execution, ()) for pid, execution in executions.items())
    return AdvntrCutoffGridResult(
        Path("/nonexistent"), "a", "b" * 64, decode_advntr_capabilities(capabilities()), rows, "d" * 64
    )


def test_the_execution_guard_binds_each_signature_to_exactly_one_execution() -> None:
    from vntyper.scripts.calibration_cutoff_advntr import require_one_execution_per_signature

    policies = {"a": _policy(), "b": _policy(kestrel_floor=0.25), "c": _policy(cutoff=0.002)}

    assert require_one_execution_per_signature(_grid({"a": "x", "b": "x", "c": "y"}), policies) == 2
    with pytest.raises(ValueError, match=r"did not replay the policies \['c'\]"):
        require_one_execution_per_signature(_grid({"a": "x", "b": "x"}), policies)
    with pytest.raises(ValueError, match="replayed by 2 executions"):
        require_one_execution_per_signature(_grid({"a": "x", "b": "z", "c": "y"}), policies)
    with pytest.raises(ValueError, match="execution x served 2 distinct adVNTR policies"):
        require_one_execution_per_signature(_grid({"a": "x", "b": "x", "c": "x"}), policies)
