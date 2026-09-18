"""Strict caller role-evidence, result, and report adapters."""

from copy import deepcopy
from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_callers import _evidence, _protocol, _roster

pytestmark = pytest.mark.unit


def test_role_evidence_roundtrip_binds_complete_rows_and_digest() -> None:
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    evidence = _evidence()

    document = artifacts.caller_role_evidence_document(evidence)

    assert artifacts.decode_caller_role_evidence(document) == evidence
    assert document["sha256"] == evidence.sha256


@pytest.mark.parametrize(
    ("field", "value"),
    [("called_positive", 1), ("truth_positive", 0), ("called_variants", ["z", "a"]), ("unexpected", True)],
)
def test_role_evidence_rejects_bool_integer_hazards_order_and_unknown_fields(field: str, value: object) -> None:
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    raw = artifacts.caller_role_evidence_document(_evidence())
    raw.pop("sha256")
    policies = raw["policies"]
    assert isinstance(policies, list) and isinstance(policies[0], dict)
    rows = policies[0]["rows"]
    assert isinstance(rows, list) and isinstance(rows[0], dict)
    rows[0][field] = value

    with pytest.raises(ValueError):
        artifacts.decode_caller_role_evidence(raw)


def test_disposition_must_match_call_and_policy_rows_are_unique() -> None:
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    raw = artifacts.caller_role_evidence_document(_evidence())
    raw.pop("sha256")
    policies = raw["policies"]
    assert isinstance(policies, list) and isinstance(policies[0], dict)
    rows = policies[0]["rows"]
    assert isinstance(rows, list) and isinstance(rows[0], dict)
    rows[0]["disposition"] = "no-call"
    with pytest.raises(ValueError, match="disposition"):
        artifacts.decode_caller_role_evidence(raw)

    duplicate = artifacts.caller_role_evidence_document(_evidence())
    duplicate.pop("sha256")
    duplicate["policies"] = [duplicate["policies"][0], duplicate["policies"][0]]
    with pytest.raises(ValueError, match="unique"):
        artifacts.decode_caller_role_evidence(duplicate)


def test_replay_digest_must_commit_the_actual_observations() -> None:
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    raw = artifacts.caller_role_evidence_document(_evidence())
    raw.pop("sha256")
    equivalence = raw["replay_equivalence"]
    assert isinstance(equivalence, list) and isinstance(equivalence[0], dict)
    equivalence[0]["baseline_replay_sha256"] = "f" * 64

    with pytest.raises(ValueError, match="actual observations"):
        artifacts.decode_caller_role_evidence(raw)


def test_result_roundtrip_recomputes_and_rejects_tampering() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol, roster = _protocol(), _roster()
    evidence = _evidence(protocol)
    result = callers.evaluate_caller_grid(protocol, roster, evidence)

    document = artifacts.caller_evaluation_document(result)

    assert artifacts.decode_caller_evaluation(document, protocol=protocol, roster=roster, evidence=evidence) == result
    forged = deepcopy(document)
    forged["selection"]["selected_candidate_id"] = protocol.candidates[1].candidate_id
    with pytest.raises(ValueError, match="recomputed"):
        artifacts.decode_caller_evaluation(forged, protocol=protocol, roster=roster, evidence=evidence)
    with pytest.raises(ValueError, match="canonical content or digest"):
        artifacts.caller_evaluation_document(replace(result, sha256="f" * 64))


def test_report_adapter_recomputes_and_does_not_render_specimen_keys() -> None:
    callers = import_module("vntyper.scripts.calibration_callers")
    artifacts = import_module("vntyper.scripts.calibration_caller_artifacts")
    protocol, roster = _protocol(), _roster()
    evidence = _evidence(protocol)
    result = callers.evaluate_caller_grid(protocol, roster, evidence)

    html = artifacts.render_caller_evaluation(result, protocol=protocol, roster=roster, evidence=evidence)

    assert "Eligible independent groups" in html
    assert result.selection.selected_candidate_id in html
    assert "specimen-a" not in html
