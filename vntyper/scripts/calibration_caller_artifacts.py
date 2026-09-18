"""Strict canonical artifacts and reports for finite caller evaluation."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import replace
from typing import NoReturn, cast

from vntyper.scripts.calibration_caller_curves import CallerCurves
from vntyper.scripts.calibration_caller_metrics import CallerObservation, validate_caller_observations
from vntyper.scripts.calibration_caller_protocol import CallerProtocol
from vntyper.scripts.calibration_caller_report import CallerReportCandidate, render_caller_comparison
from vntyper.scripts.calibration_caller_roster import CallerEligibleRoster
from vntyper.scripts.calibration_callers import (
    CallerEvaluationResult,
    CallerEvidenceRow,
    CallerPolicyEvidence,
    CallerReplayEquivalence,
    CallerRoleEvidence,
    EvaluationPhase,
    EvidenceDisposition,
    ExecutionKind,
    caller_evaluation_payload,
    evaluate_caller_grid,
)
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_PHASES = {"policy-selection", "validation", "locked-heldout", "development-assessment"}
_KINDS = {"baseline-rerun", "scalar-replay", "recapture"}
_DISPOSITIONS = {"called", "negative", "no-call", "zero-candidate", "unsupported"}
_ROOT_FIELDS = {
    "schema_version",
    "phase",
    "protocol_sha256",
    "eligible_roster_sha256",
    "run_manifest_sha256",
    "baseline_assets_sha256",
    "policies",
    "replay_equivalence",
}
_POLICY_FIELDS = {
    "candidate_id",
    "policy_sha256",
    "execution_kind",
    "capture_policy_sha256",
    "rows",
}
_ROW_FIELDS = {
    "key",
    "group_key",
    "truth_positive",
    "truth_variants",
    "called_positive",
    "called_variants",
    "tier_a_variants",
    "disposition",
    "source_evidence_sha256",
}
_EQUIVALENCE_FIELDS = {
    "capture_policy_sha256",
    "baseline_policy_sha256",
    "baseline_rerun_sha256",
    "baseline_replay_sha256",
    "baseline_replay_rows",
}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"caller {label} fields differ from the closed contract")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"caller {label} must be a lowercase SHA256 digest")
    return value


def _strings(value: object, label: str, *, nullable: bool = False) -> tuple[str, ...] | None:
    if nullable and value is None:
        return None
    if not isinstance(value, list) or any(
        not isinstance(item, str) or not item or item.strip() != item for item in value
    ):
        _fail(f"caller {label} must be a list of non-empty trimmed strings")
    if value != sorted(set(value)):
        _fail(f"caller {label} must be sorted and unique")
    return tuple(value)


def _decision(value: object, label: str) -> bool | None:
    if value is not None and not isinstance(value, bool):
        _fail(f"caller {label} must be Boolean or null")
    return cast(bool | None, value)


def _decode_observation(value: object, *, evidence: bool) -> tuple[CallerObservation, str | None, str | None]:
    fields = _ROW_FIELDS if evidence else _ROW_FIELDS - {"disposition", "source_evidence_sha256"}
    raw = _object(value, fields, "evidence row" if evidence else "baseline replay row")
    key = raw["key"]
    group = raw["group_key"]
    if any(not isinstance(item, str) or not item or item.strip() != item for item in (key, group)):
        _fail("caller evidence row identities must be non-empty trimmed strings")
    truth_variants = _strings(raw["truth_variants"], "truth variants", nullable=True)
    called_variants = _strings(raw["called_variants"], "called variants")
    tier_a = _strings(raw["tier_a_variants"], "tier-A variants")
    observation = CallerObservation(
        cast(str, key),
        cast(str, group),
        _decision(raw["truth_positive"], "truth"),
        truth_variants,
        _decision(raw["called_positive"], "call"),
        cast(tuple[str, ...], called_variants),
        cast(tuple[str, ...], tier_a),
    )
    validate_caller_observations((observation,))
    if not evidence:
        return observation, None, None
    disposition = raw["disposition"]
    if not isinstance(disposition, str) or disposition not in _DISPOSITIONS:
        _fail("caller evidence disposition is unsupported")
    expected_calls = {
        "called": True,
        "negative": False,
        "zero-candidate": False,
        "no-call": None,
        "unsupported": None,
    }
    if observation.called_positive is not expected_calls[disposition]:
        _fail("caller evidence disposition differs from the recorded call")
    return observation, disposition, _digest(raw["source_evidence_sha256"], "source evidence")


def _observation_document(row: CallerObservation) -> dict[str, object]:
    validate_caller_observations((row,))
    return {
        "key": row.key,
        "group_key": row.group_key,
        "truth_positive": row.truth_positive,
        "truth_variants": None if row.truth_variants is None else list(row.truth_variants),
        "called_positive": row.called_positive,
        "called_variants": list(row.called_variants),
        "tier_a_variants": list(row.tier_a_variants),
    }


def _policy_payload(policy: CallerPolicyEvidence) -> dict[str, object]:
    return {
        "candidate_id": policy.candidate_id,
        "policy_sha256": policy.policy_sha256,
        "execution_kind": policy.execution_kind,
        "capture_policy_sha256": policy.capture_policy_sha256,
        "rows": [
            {
                **_observation_document(row.observation),
                "disposition": row.disposition,
                "source_evidence_sha256": row.source_evidence_sha256,
            }
            for row in policy.rows
        ],
    }


def _decode_policy(value: object) -> CallerPolicyEvidence:
    raw = _object(value, _POLICY_FIELDS, "policy evidence")
    kind = raw["execution_kind"]
    if not isinstance(kind, str) or kind not in _KINDS:
        _fail("caller policy evidence execution kind is unsupported")
    rows_raw = raw["rows"]
    if not isinstance(rows_raw, list) or not rows_raw:
        _fail("caller policy evidence rows must be a non-empty list")
    rows: list[CallerEvidenceRow] = []
    for value_row in rows_raw:
        observation, disposition, source = _decode_observation(value_row, evidence=True)
        rows.append(
            CallerEvidenceRow(
                observation,
                cast(EvidenceDisposition, disposition),
                cast(str, source),
            )
        )
    observations = validate_caller_observations(tuple(row.observation for row in rows))
    if observations != tuple(row.observation for row in rows):
        _fail("caller policy evidence rows must be sorted by independent group")
    candidate_id = _digest(raw["candidate_id"], "candidate identity")
    policy_sha256 = _digest(raw["policy_sha256"], "policy identity")
    if candidate_id != policy_sha256:
        _fail("caller candidate identity must equal its full policy digest")
    result = CallerPolicyEvidence(
        candidate_id,
        policy_sha256,
        cast(ExecutionKind, kind),
        _digest(raw["capture_policy_sha256"], "capture policy"),
        tuple(rows),
        "",
    )
    return replace(result, sha256=canonical_sha256(_policy_payload(result)))


def _equivalence_document(value: CallerReplayEquivalence) -> dict[str, object]:
    return {
        "capture_policy_sha256": value.capture_policy_sha256,
        "baseline_policy_sha256": value.baseline_policy_sha256,
        "baseline_rerun_sha256": value.baseline_rerun_sha256,
        "baseline_replay_sha256": value.baseline_replay_sha256,
        "baseline_replay_rows": [_observation_document(row) for row in value.baseline_replay_rows],
    }


def _decode_equivalence(value: object) -> CallerReplayEquivalence:
    raw = _object(value, _EQUIVALENCE_FIELDS, "replay equivalence")
    rows_raw = raw["baseline_replay_rows"]
    if not isinstance(rows_raw, list) or not rows_raw:
        _fail("caller baseline replay observations must be a non-empty list")
    rows = tuple(_decode_observation(item, evidence=False)[0] for item in rows_raw)
    if validate_caller_observations(rows) != rows:
        _fail("caller baseline replay observations must be sorted by independent group")
    replay_sha256 = _digest(raw["baseline_replay_sha256"], "equivalence baseline replay")
    if replay_sha256 != canonical_sha256([_observation_document(row) for row in rows]):
        _fail("caller baseline replay digest differs from its actual observations")
    return CallerReplayEquivalence(
        _digest(raw["capture_policy_sha256"], "equivalence capture policy"),
        _digest(raw["baseline_policy_sha256"], "equivalence baseline policy"),
        _digest(raw["baseline_rerun_sha256"], "equivalence baseline rerun"),
        replay_sha256,
        rows,
    )


def _role_payload(evidence: CallerRoleEvidence) -> dict[str, object]:
    return {
        "schema_version": "calibration-caller-role-evidence-v1",
        "phase": evidence.phase,
        "protocol_sha256": evidence.protocol_sha256,
        "eligible_roster_sha256": evidence.eligible_roster_sha256,
        "run_manifest_sha256": evidence.run_manifest_sha256,
        "baseline_assets_sha256": evidence.baseline_assets_sha256,
        "policies": [_policy_payload(policy) for policy in evidence.policies],
        "replay_equivalence": [_equivalence_document(item) for item in evidence.replay_equivalence],
    }


def decode_caller_role_evidence(value: object) -> CallerRoleEvidence:
    """Decode complete role outcomes and their replay-equivalence observations."""
    if not isinstance(value, Mapping) or set(value) not in (_ROOT_FIELDS, {*_ROOT_FIELDS, "sha256"}):
        _fail("caller role evidence fields differ from the closed contract")
    raw = value
    if raw["schema_version"] != "calibration-caller-role-evidence-v1":
        _fail("caller role evidence schema is unsupported")
    phase = raw["phase"]
    if not isinstance(phase, str) or phase not in _PHASES:
        _fail("caller role evidence phase is unsupported")
    policies_raw = raw["policies"]
    equivalence_raw = raw["replay_equivalence"]
    if not isinstance(policies_raw, list) or not policies_raw or not isinstance(equivalence_raw, list):
        _fail("caller role evidence policies and equivalence must be lists")
    policies = tuple(_decode_policy(item) for item in policies_raw)
    ids = tuple(policy.candidate_id for policy in policies)
    if ids != tuple(sorted(set(ids))):
        _fail("caller policy evidence identities must be sorted and unique")
    row_shapes = {tuple((row.observation.key, row.observation.group_key) for row in policy.rows) for policy in policies}
    if len(row_shapes) != 1:
        _fail("caller policies must retain the identical complete role-row roster")
    equivalence = tuple(_decode_equivalence(item) for item in equivalence_raw)
    keys = tuple((item.capture_policy_sha256, item.baseline_policy_sha256) for item in equivalence)
    if keys != tuple(sorted(set(keys))):
        _fail("caller replay equivalence identities must be sorted and unique")
    evidence = CallerRoleEvidence(
        cast(EvaluationPhase, phase),
        _digest(raw["protocol_sha256"], "protocol"),
        _digest(raw["eligible_roster_sha256"], "eligible roster"),
        _digest(raw["run_manifest_sha256"], "run manifest"),
        _digest(raw["baseline_assets_sha256"], "baseline assets"),
        policies,
        equivalence,
        "",
    )
    evidence = replace(evidence, sha256=canonical_sha256(_role_payload(evidence)))
    if "sha256" in raw and raw["sha256"] != evidence.sha256:
        _fail("caller role evidence digest differs from its canonical content")
    return evidence


def caller_role_evidence_document(evidence: CallerRoleEvidence) -> dict[str, object]:
    """Project role evidence after revalidating typed content and its digest."""
    if (
        not isinstance(evidence, CallerRoleEvidence)
        or not isinstance(evidence.policies, tuple)
        or not isinstance(evidence.replay_equivalence, tuple)
        or any(
            not isinstance(policy, CallerPolicyEvidence) or not isinstance(policy.rows, tuple)
            for policy in evidence.policies
        )
        or any(not isinstance(item, CallerReplayEquivalence) for item in evidence.replay_equivalence)
    ):
        _fail("caller role evidence projection requires typed evidence")
    payload = _role_payload(evidence)
    decoded = decode_caller_role_evidence(payload)
    if decoded != evidence:
        _fail("caller role evidence differs from its canonical content or digest")
    return {**payload, "sha256": evidence.sha256}


def caller_evaluation_document(result: CallerEvaluationResult) -> dict[str, object]:
    """Project a recomputable caller evaluation result with its canonical digest."""
    if not isinstance(result, CallerEvaluationResult):
        _fail("caller evaluation projection requires a typed result")
    payload = caller_evaluation_payload(result)
    if canonical_sha256(payload) != result.sha256:
        _fail("caller evaluation differs from its canonical content or digest")
    return {**payload, "sha256": result.sha256}


def decode_caller_evaluation(
    value: object,
    *,
    protocol: CallerProtocol,
    roster: CallerEligibleRoster,
    evidence: CallerRoleEvidence,
) -> CallerEvaluationResult:
    """Decode a result by recomputing it from opened protocol and complete evidence."""
    if not isinstance(value, Mapping) or set(value) != {*caller_evaluation_payload_placeholder(), "sha256"}:
        _fail("caller evaluation fields differ from the closed contract")
    selected_raw = value.get("selection")
    if not isinstance(selected_raw, Mapping):
        _fail("caller evaluation selection must be an object")
    fixed = selected_raw.get("selected_candidate_id") if evidence.phase != "policy-selection" else None
    recomputed = evaluate_caller_grid(protocol, roster, evidence, fixed_candidate_id=cast(str | None, fixed))
    if dict(value) != caller_evaluation_document(recomputed):
        _fail("caller evaluation differs from recomputed evidence")
    return recomputed


def caller_evaluation_payload_placeholder() -> set[str]:
    """Return exact caller evaluation payload keys for strict decoding."""
    return {
        "schema_version",
        "phase",
        "protocol_sha256",
        "eligible_roster_sha256",
        "role_evidence_sha256",
        "baseline_policy_sha256",
        "candidates",
        "selection",
    }


def render_caller_evaluation(
    result: CallerEvaluationResult,
    *,
    protocol: CallerProtocol,
    roster: CallerEligibleRoster,
    evidence: CallerRoleEvidence,
    curves: tuple[CallerCurves, ...] = (),
) -> str:
    """Render an evaluation after recomputing decisions from complete outcomes.

    Args:
        result: Hash-bound result produced by the finite-grid engine.
        protocol: Opened outcome-independent caller protocol.
        roster: Opened complete independent role roster.
        evidence: Opened complete production role evidence.
        curves: Optional verified selection-only comparable scalar curves.

    Returns:
        Escaped standalone HTML without specimen or group identifiers.

    Raises:
        ValueError: If any opened artifact, result, or curve is inconsistent.
    """
    checked = decode_caller_evaluation(
        caller_evaluation_document(result), protocol=protocol, roster=roster, evidence=evidence
    )
    by_id = {policy.candidate_id: policy for policy in evidence.policies}
    candidates = [
        CallerReportCandidate(
            protocol.baseline_policy_sha256,
            tuple(row.observation for row in by_id[protocol.baseline_policy_sha256].rows),
            (),
        )
    ]
    for evaluation in checked.candidates:
        reasons = tuple(
            dict.fromkeys(
                (
                    *evaluation.acceptance.pooled.reasons,
                    *(reason for gate in evaluation.acceptance.strata.values() for reason in gate.reasons),
                )
            )
        )
        candidates.append(
            CallerReportCandidate(
                evaluation.candidate_id,
                tuple(row.observation for row in by_id[evaluation.candidate_id].rows),
                reasons,
            )
        )
    return render_caller_comparison(
        tuple(candidates),
        roster,
        phase=checked.phase,
        baseline_id=checked.baseline_policy_sha256,
        selected_id=checked.selection.selected_candidate_id,
        protocol_sha256=checked.protocol_sha256,
        evidence_sha256=checked.role_evidence_sha256,
        required_strata=protocol.required_strata,
        curves=curves,
    )
