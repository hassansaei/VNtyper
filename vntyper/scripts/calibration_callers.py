"""Finite caller-policy evaluation over complete immutable role evidence."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, replace
from typing import Literal, NoReturn

from vntyper.scripts.calibration_caller_acceptance import (
    CallerAcceptance,
    CallerSelectionEntry,
    evaluate_caller_acceptance,
    select_caller_candidate,
)
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_protocol import CallerProtocol, caller_protocol_document
from vntyper.scripts.calibration_caller_roster import CallerEligibleRoster, caller_eligible_roster_document
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

EvaluationPhase = Literal["policy-selection", "validation", "locked-heldout", "development-assessment"]
ExecutionKind = Literal["baseline-rerun", "scalar-replay", "recapture"]
EvidenceDisposition = Literal["called", "negative", "no-call", "zero-candidate", "unsupported"]
SelectionStatus = Literal["selected", "no-feasible-candidate", "not-applicable"]
_NON_SELECTION_PHASES = {"validation", "locked-heldout", "development-assessment"}
_ADVNTR_PREFIX = "/components/advntr/calibrated_calling/"
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class CallerEvidenceRow:
    """One complete role member outcome and its source-artifact commitment."""

    observation: CallerObservation
    disposition: EvidenceDisposition
    source_evidence_sha256: str


@dataclass(frozen=True)
class CallerPolicyEvidence:
    """Complete production outcomes for one frozen caller policy."""

    candidate_id: str
    policy_sha256: str
    execution_kind: ExecutionKind
    capture_policy_sha256: str
    rows: tuple[CallerEvidenceRow, ...]
    sha256: str


@dataclass(frozen=True)
class CallerReplayEquivalence:
    """Actual baseline rerun and baseline scalar-replay observations."""

    capture_policy_sha256: str
    baseline_policy_sha256: str
    baseline_rerun_sha256: str
    baseline_replay_sha256: str
    baseline_replay_rows: tuple[CallerObservation, ...]


@dataclass(frozen=True)
class CallerRoleEvidence:
    """Hash-bound complete policy outcomes for one authorized role."""

    phase: EvaluationPhase
    protocol_sha256: str
    eligible_roster_sha256: str
    run_manifest_sha256: str
    baseline_assets_sha256: str
    policies: tuple[CallerPolicyEvidence, ...]
    replay_equivalence: tuple[CallerReplayEquivalence, ...]
    sha256: str


@dataclass(frozen=True)
class CallerCandidateEvaluation:
    """One protocol candidate evaluated against the fixed baseline."""

    candidate_id: str
    free_parameters: int
    execution_kind: ExecutionKind
    acceptance: CallerAcceptance


@dataclass(frozen=True)
class CallerSelection:
    """Selection decision; non-selection phases only identify their frozen point."""

    status: SelectionStatus
    selected_candidate_id: str | None
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class CallerEvaluationResult:
    """Hash-bound local caller evaluation; it grants no promotion authority."""

    phase: EvaluationPhase
    protocol_sha256: str
    eligible_roster_sha256: str
    role_evidence_sha256: str
    baseline_policy_sha256: str
    candidates: tuple[CallerCandidateEvaluation, ...]
    selection: CallerSelection
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _expected_execution(protocol: CallerProtocol, candidate_id: str) -> ExecutionKind:
    candidate = next(item for item in protocol.candidates if item.candidate_id == candidate_id)
    changed = {
        pointer
        for pointer, value in candidate.policy.values.items()
        if value != protocol.baseline_policy.values[pointer]
    }
    return "recapture" if any(pointer.startswith(_ADVNTR_PREFIX) for pointer in changed) else "scalar-replay"


def _observations(policy: CallerPolicyEvidence) -> tuple[CallerObservation, ...]:
    return tuple(row.observation for row in policy.rows)


def _require_replay_equivalence(
    candidate: CallerPolicyEvidence,
    baseline: CallerPolicyEvidence,
    attestations: tuple[CallerReplayEquivalence, ...],
) -> None:
    matching = tuple(
        item
        for item in attestations
        if item.capture_policy_sha256 == candidate.capture_policy_sha256
        and item.baseline_policy_sha256 == baseline.policy_sha256
    )
    if len(matching) != 1:
        _fail("caller scalar replay requires one actual baseline equivalence evidence pair")
    attestation = matching[0]
    if attestation.baseline_rerun_sha256 != baseline.sha256:
        _fail("caller replay equivalence baseline rerun digest differs from role evidence")
    if attestation.baseline_replay_rows != _observations(baseline):
        _fail("caller baseline production and replay observations are not equivalent")


def _result_payload(result: CallerEvaluationResult) -> dict[str, object]:
    return {
        "schema_version": "calibration-caller-evaluation-v1",
        "phase": result.phase,
        "protocol_sha256": result.protocol_sha256,
        "eligible_roster_sha256": result.eligible_roster_sha256,
        "role_evidence_sha256": result.role_evidence_sha256,
        "baseline_policy_sha256": result.baseline_policy_sha256,
        "candidates": [
            {
                "candidate_id": candidate.candidate_id,
                "free_parameters": candidate.free_parameters,
                "execution_kind": candidate.execution_kind,
                "acceptance": {
                    "status": candidate.acceptance.status,
                    "selection_benefit": candidate.acceptance.selection_benefit,
                    "pooled": {
                        "status": candidate.acceptance.pooled.status,
                        "reasons": list(candidate.acceptance.pooled.reasons),
                    },
                    "strata": {
                        name: {"status": gate.status, "reasons": list(gate.reasons)}
                        for name, gate in sorted(candidate.acceptance.strata.items())
                    },
                },
            }
            for candidate in result.candidates
        ],
        "selection": {
            "status": result.selection.status,
            "selected_candidate_id": result.selection.selected_candidate_id,
            "reasons": list(result.selection.reasons),
        },
    }


def evaluate_caller_grid(
    protocol: CallerProtocol,
    roster: CallerEligibleRoster,
    evidence: CallerRoleEvidence,
    *,
    fixed_candidate_id: str | None = None,
) -> CallerEvaluationResult:
    """Evaluate a complete finite caller grid or one already frozen policy.

    Args:
        protocol: Outcome-independent full-policy candidate declaration.
        roster: Exact independent representatives for this role.
        evidence: Complete baseline and candidate production outcomes.
        fixed_candidate_id: Required outside policy selection and forbidden during it.

    Returns:
        Immutable acceptance and deterministic selection evidence.

    Raises:
        ValueError: If bindings, policy coverage, replay parity, or phase use differ.
    """
    caller_protocol_document(protocol)
    caller_eligible_roster_document(roster)
    from vntyper.scripts.calibration_caller_artifacts import caller_role_evidence_document

    caller_role_evidence_document(evidence)
    if evidence.protocol_sha256 != protocol.sha256 or evidence.eligible_roster_sha256 != roster.sha256:
        _fail("caller role evidence differs from its protocol or eligible roster")
    by_id = {policy.candidate_id: policy for policy in evidence.policies}
    baseline = by_id.get(protocol.baseline_policy_sha256)
    if baseline is None or baseline.policy_sha256 != protocol.baseline_policy_sha256:
        _fail("caller role evidence lacks the frozen baseline policy")
    if baseline.execution_kind != "baseline-rerun":
        _fail("caller baseline evidence must be a production baseline-rerun")

    if evidence.phase == "policy-selection":
        if fixed_candidate_id is not None:
            _fail("caller policy-selection must not receive fixed_candidate_id")
        expected = {protocol.baseline_policy_sha256, *(item.candidate_id for item in protocol.candidates)}
        if set(by_id) != expected:
            _fail("caller policy-selection evidence must match the complete protocol policy roster")
        candidates = protocol.candidates
    elif evidence.phase in _NON_SELECTION_PHASES:
        candidate_ids = {item.candidate_id for item in protocol.candidates}
        if not isinstance(fixed_candidate_id, str) or fixed_candidate_id not in candidate_ids:
            _fail("caller nonselection evaluation requires a valid fixed_candidate_id")
        if set(by_id) != {protocol.baseline_policy_sha256, fixed_candidate_id}:
            _fail("caller nonselection evidence requires only baseline and fixed candidate policies")
        candidates = tuple(item for item in protocol.candidates if item.candidate_id == fixed_candidate_id)
    else:
        _fail("caller role evidence phase is unsupported")

    baseline_rows = _observations(baseline)
    expected_equivalence = {
        (by_id[candidate.candidate_id].capture_policy_sha256, protocol.baseline_policy_sha256)
        for candidate in candidates
        if _expected_execution(protocol, candidate.candidate_id) == "scalar-replay"
    }
    actual_equivalence = {
        (item.capture_policy_sha256, item.baseline_policy_sha256) for item in evidence.replay_equivalence
    }
    if actual_equivalence != expected_equivalence:
        _fail("caller role evidence must contain exactly the required replay equivalence pairs")
    evaluated: list[CallerCandidateEvaluation] = []
    for candidate in candidates:
        policy = by_id[candidate.candidate_id]
        if policy.policy_sha256 != candidate.policy.sha256:
            _fail("caller evidence policy digest differs from its full protocol policy")
        expected_kind = _expected_execution(protocol, candidate.candidate_id)
        if policy.execution_kind != expected_kind:
            _fail("caller evidence replay or recapture classification differs from policy changes")
        if expected_kind == "scalar-replay":
            if policy.capture_policy_sha256 != baseline.capture_policy_sha256:
                _fail("caller scalar replay must use the baseline capture policy")
            _require_replay_equivalence(policy, baseline, evidence.replay_equivalence)
        acceptance = evaluate_caller_acceptance(
            _observations(policy), baseline_rows, roster, protocol.gate_rules, phase=evidence.phase
        )
        evaluated.append(
            CallerCandidateEvaluation(candidate.candidate_id, candidate.free_parameters, expected_kind, acceptance)
        )

    evaluated_tuple = tuple(evaluated)
    if evidence.phase == "policy-selection":
        selected = select_caller_candidate(
            tuple(
                CallerSelectionEntry(row.candidate_id, row.free_parameters, row.acceptance) for row in evaluated_tuple
            )
        )
        selection = (
            CallerSelection("selected", selected, ())
            if selected is not None
            else CallerSelection("no-feasible-candidate", None, ("no_candidate_passed_acceptance_with_benefit",))
        )
    else:
        selection = CallerSelection("not-applicable", fixed_candidate_id, ("selection_not_allowed_for_phase",))
    result = CallerEvaluationResult(
        evidence.phase,
        protocol.sha256,
        roster.sha256,
        evidence.sha256,
        protocol.baseline_policy_sha256,
        evaluated_tuple,
        selection,
        "",
    )
    return replace(result, sha256=canonical_sha256(_result_payload(result)))


def caller_evaluation_payload(result: CallerEvaluationResult) -> dict[str, object]:
    """Return the canonical digest payload after validating typed result content."""
    if (
        not isinstance(result, CallerEvaluationResult)
        or result.phase not in {"policy-selection", "validation", "locked-heldout", "development-assessment"}
        or any(
            not isinstance(getattr(result, field), str) or _SHA256.fullmatch(getattr(result, field)) is None
            for field in (
                "protocol_sha256",
                "eligible_roster_sha256",
                "role_evidence_sha256",
                "baseline_policy_sha256",
            )
        )
        or not isinstance(result.candidates, tuple)
        or any(
            not isinstance(candidate, CallerCandidateEvaluation)
            or not isinstance(candidate.acceptance, CallerAcceptance)
            or isinstance(candidate.free_parameters, bool)
            or not isinstance(candidate.free_parameters, int)
            or candidate.free_parameters < 0
            for candidate in result.candidates
        )
        or not isinstance(result.selection, CallerSelection)
        or not isinstance(result.selection.reasons, tuple)
    ):
        _fail("caller evaluation result must use the immutable typed contract")
    return _result_payload(result)
