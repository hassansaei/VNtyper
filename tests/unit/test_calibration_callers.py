"""Finite caller grid evaluation over complete role evidence."""

from copy import deepcopy
from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_caller_policy import policy_document
from tests.unit.test_calibration_caller_protocol import changed_policy, protocol_document
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
from vntyper.scripts.calibration_caller_protocol import decode_caller_protocol
from vntyper.scripts.calibration_caller_roster import decode_caller_eligible_roster
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def _protocol():
    baseline = decode_caller_policy_values(policy_document())
    raw = protocol_document(baseline)
    acceptance = raw["acceptance"]
    assert isinstance(acceptance, dict)
    acceptance.update(
        minimum_positive_groups=1,
        minimum_negative_groups=1,
        maximum_fpr_upper=1.0,
        minimum_sensitivity_delta_lower=-1.0,
    )
    return decode_caller_protocol(raw, baseline_policy=baseline)


def _roster():
    return decode_caller_eligible_roster(
        [
            {"key": "specimen-a", "group_key": "group-a", "strata": ["fallback", "nominal"]},
            {"key": "specimen-b", "group_key": "group-b", "strata": ["fallback", "nominal"]},
            {"key": "specimen-c", "group_key": "group-c", "strata": ["fallback", "nominal"]},
            {"key": "specimen-d", "group_key": "group-d", "strata": ["fallback", "nominal"]},
        ]
    )


def _evidence(protocol=None, *, phase: str = "policy-selection"):
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    if protocol is None:
        protocol = _protocol()
    rows = []
    for key, group, truth in (
        ("specimen-a", "group-a", True),
        ("specimen-b", "group-b", True),
        ("specimen-c", "group-c", False),
        ("specimen-d", "group-d", False),
    ):
        rows.append(
            {
                "key": key,
                "group_key": group,
                "truth_positive": truth,
                "truth_variants": ["variant-a"] if truth else [],
                "called_positive": truth if key != "specimen-b" else False,
                "called_variants": ["variant-a"] if truth and key != "specimen-b" else [],
                "tier_a_variants": [],
                "disposition": "called" if truth and key != "specimen-b" else "zero-candidate",
                "source_evidence_sha256": (key[-1] * 64),
            }
        )
    baseline_policy = {
        "candidate_id": protocol.baseline_policy_sha256,
        "policy_sha256": protocol.baseline_policy_sha256,
        "execution_kind": "baseline-rerun",
        "capture_policy_sha256": "a" * 64,
        "rows": rows,
    }
    policies = [baseline_policy]
    for index, candidate in enumerate(protocol.candidates):
        candidate_rows = [dict(row) for row in rows]
        if index == 0:
            candidate_rows[1].update(
                called_positive=True,
                called_variants=["variant-a"],
                disposition="called",
            )
        policies.append(
            {
                "candidate_id": candidate.candidate_id,
                "policy_sha256": candidate.policy.sha256,
                "execution_kind": "scalar-replay",
                "capture_policy_sha256": "a" * 64,
                "rows": candidate_rows,
            }
        )
    policies.sort(key=lambda row: str(row["candidate_id"]))
    raw = {
        "schema_version": "calibration-caller-role-evidence-v1",
        "phase": phase,
        "protocol_sha256": protocol.sha256,
        "eligible_roster_sha256": _roster().sha256,
        "run_manifest_sha256": "b" * 64,
        "baseline_assets_sha256": "c" * 64,
        "policies": policies,
        "replay_equivalence": [],
    }
    evidence = artifacts.decode_caller_role_evidence(raw)
    baseline = next(row for row in evidence.policies if row.execution_kind == "baseline-rerun")
    replay_rows = [
        {key: value for key, value in row.items() if key not in {"disposition", "source_evidence_sha256"}}
        for row in rows
    ]
    raw["replay_equivalence"] = [
        {
            "capture_policy_sha256": "a" * 64,
            "baseline_policy_sha256": protocol.baseline_policy_sha256,
            "baseline_rerun_sha256": baseline.sha256,
            "baseline_replay_sha256": canonical_sha256(replay_rows),
            "baseline_replay_rows": replay_rows,
        }
    ]
    return artifacts.decode_caller_role_evidence(raw)


def test_policy_selection_evaluates_complete_grid_and_selects_benefit() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    protocol = _protocol()
    evidence = _evidence(protocol)

    result = callers.evaluate_caller_grid(protocol, _roster(), evidence)

    assert len(result.candidates) == len(protocol.candidates)
    assert result.selection.status == "selected"
    assert result.selection.selected_candidate_id == protocol.candidates[0].candidate_id
    assert result.candidates[0].acceptance.pooled.candidate is not None
    assert result.candidates[0].acceptance.pooled.candidate.no_calls == 0


def test_zero_candidate_and_no_call_rows_are_retained_in_denominators() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol = _protocol()
    raw = artifacts.caller_role_evidence_document(_evidence(protocol))
    policies = raw["policies"]
    assert isinstance(policies, list)
    for policy in policies:
        assert isinstance(policy, dict)
        rows = policy["rows"]
        assert isinstance(rows, list)
        if policy["execution_kind"] != "baseline-rerun":
            rows[3].update(called_positive=None, disposition="no-call")
    evidence = artifacts.decode_caller_role_evidence({key: value for key, value in raw.items() if key != "sha256"})

    result = callers.evaluate_caller_grid(protocol, _roster(), evidence)

    assert all(candidate.acceptance.pooled.candidate.no_calls == 1 for candidate in result.candidates)
    assert len(evidence.policies[0].rows) == 4


def test_selection_refuses_missing_policy_or_replay_attestation() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol = _protocol()
    evidence = _evidence(protocol)

    raw = artifacts.caller_role_evidence_document(evidence)
    raw.pop("sha256")
    raw["policies"] = raw["policies"][:-1]
    with pytest.raises(ValueError, match="policy roster"):
        callers.evaluate_caller_grid(protocol, _roster(), artifacts.decode_caller_role_evidence(raw))


def test_scalar_replay_requires_actual_baseline_observation_equivalence() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol = _protocol()
    raw = artifacts.caller_role_evidence_document(_evidence(protocol))
    raw.pop("sha256")
    equivalence = raw["replay_equivalence"]
    assert isinstance(equivalence, list) and isinstance(equivalence[0], dict)
    rows = equivalence[0]["baseline_replay_rows"]
    assert isinstance(rows, list) and isinstance(rows[0], dict)
    rows[0].update(called_positive=None, called_variants=[], disposition=None)
    rows[0].pop("disposition")
    equivalence[0]["baseline_replay_sha256"] = canonical_sha256(rows)

    evidence = artifacts.decode_caller_role_evidence(raw)
    with pytest.raises(ValueError, match="not equivalent"):
        callers.evaluate_caller_grid(protocol, _roster(), evidence)
    raw = artifacts.caller_role_evidence_document(evidence)
    raw.pop("sha256")
    raw["replay_equivalence"] = []
    with pytest.raises(ValueError, match="equivalence"):
        callers.evaluate_caller_grid(protocol, _roster(), artifacts.decode_caller_role_evidence(raw))


def test_nonselection_requires_only_frozen_candidate_and_never_selects() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol = _protocol()
    selected = protocol.candidates[0].candidate_id
    raw = artifacts.caller_role_evidence_document(_evidence(protocol))
    raw.pop("sha256")
    raw["phase"] = "validation"
    raw["policies"] = [
        row for row in raw["policies"] if row["candidate_id"] in {selected, protocol.baseline_policy_sha256}
    ]
    evidence = artifacts.decode_caller_role_evidence(raw)

    result = callers.evaluate_caller_grid(protocol, _roster(), evidence, fixed_candidate_id=selected)

    assert len(result.candidates) == 1
    assert result.selection.status == "not-applicable"
    assert result.selection.selected_candidate_id == selected
    with pytest.raises(ValueError, match="fixed_candidate_id"):
        callers.evaluate_caller_grid(protocol, _roster(), evidence)


def test_no_beneficial_candidate_is_a_completed_failed_selection() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol = _protocol()
    raw = artifacts.caller_role_evidence_document(_evidence(protocol))
    raw.pop("sha256")
    policies = raw["policies"]
    assert isinstance(policies, list)
    baseline = next(row for row in policies if row["execution_kind"] == "baseline-rerun")
    for policy in policies:
        if policy["execution_kind"] != "baseline-rerun":
            policy["rows"] = deepcopy(baseline["rows"])
    evidence = artifacts.decode_caller_role_evidence(raw)

    result = callers.evaluate_caller_grid(protocol, _roster(), evidence)

    assert result.selection.status == "no-feasible-candidate"
    assert result.selection.selected_candidate_id is None


@pytest.mark.parametrize(
    "pointer,value",
    [
        ("/components/advntr/calibrated_calling/mode", "legacy"),
        ("/components/advntr/calibrated_calling/cutoff", 0.025),
        ("/components/advntr/calibrated_calling/minimum_read_support", 5),
    ],
)
def test_advntr_native_decision_policy_changes_use_verified_scalar_replay(pointer: str, value: object) -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    baseline = decode_caller_policy_values(policy_document())
    candidate = changed_policy(**{pointer: value})
    raw_protocol = protocol_document(baseline, [(candidate, 1)])
    acceptance = raw_protocol["acceptance"]
    assert isinstance(acceptance, dict)
    acceptance.update(
        minimum_positive_groups=1,
        minimum_negative_groups=1,
        maximum_fpr_upper=1.0,
        minimum_sensitivity_delta_lower=-1.0,
    )
    protocol = decode_caller_protocol(raw_protocol, baseline_policy=baseline)

    result = callers.evaluate_caller_grid(protocol, _roster(), _evidence(protocol))

    assert result.candidates[0].execution_kind == "scalar-replay"


@pytest.mark.parametrize(
    "pointer,value",
    [
        ("/components/advntr/calibrated_calling/rare_unit_fraction", 0.2),
        ("/components/advntr/calibrated_calling/adapter_filter", True),
        ("/components/advntr/calibrated_calling/minimum_read_match_ratio", 0.8),
        ("/components/advntr/calibrated_calling/prune_reverse", True),
    ],
)
def test_advntr_capture_policy_change_is_recapture_and_needs_no_scalar_attestation(
    pointer: str, value: object
) -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    baseline = decode_caller_policy_values(policy_document())
    candidate = changed_policy(**{pointer: value})
    raw_protocol = protocol_document(baseline, [(candidate, 1)])
    acceptance = raw_protocol["acceptance"]
    assert isinstance(acceptance, dict)
    acceptance.update(
        minimum_positive_groups=1,
        minimum_negative_groups=1,
        maximum_fpr_upper=1.0,
        minimum_sensitivity_delta_lower=-1.0,
    )
    protocol = decode_caller_protocol(raw_protocol, baseline_policy=baseline)
    raw = artifacts.caller_role_evidence_document(_evidence(protocol))
    raw.pop("sha256")
    policies = raw["policies"]
    assert isinstance(policies, list)
    candidate_row = next(row for row in policies if row["candidate_id"] == candidate.sha256)
    candidate_row.update(execution_kind="recapture", capture_policy_sha256="e" * 64)
    raw["replay_equivalence"] = []

    result = callers.evaluate_caller_grid(protocol, _roster(), artifacts.decode_caller_role_evidence(raw))

    assert result.candidates[0].execution_kind == "recapture"


def test_public_boundary_revalidates_directly_forged_typed_evidence() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    protocol = _protocol()
    evidence = _evidence(protocol)

    with pytest.raises(ValueError, match="canonical content or digest"):
        callers.evaluate_caller_grid(protocol, _roster(), replace(evidence, sha256="f" * 64))
