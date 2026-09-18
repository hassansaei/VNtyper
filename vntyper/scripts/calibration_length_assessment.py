"""Fixed-model development assessment for total-length research."""

from __future__ import annotations

import hashlib
import logging
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import NoReturn, cast

from vntyper.scripts.calibration_length_artifacts import (
    LengthTrainingProfile,
    validate_length_training_profile,
)
from vntyper.scripts.calibration_length_evaluation import (
    LengthEvaluationResult,
    evaluate_length_hypotheses,
    length_evaluation_document,
)
from vntyper.scripts.calibration_length_evidence import LengthRoleEvidence, length_role_evidence_document
from vntyper.scripts.calibration_length_metrics import LengthEligibleRoster
from vntyper.scripts.calibration_length_protocol import LengthProtocol
from vntyper.scripts.calibration_length_report import render_length_evaluation
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_estimation import EvidenceDomain, estimate_total_repeats
from vntyper.scripts.length_model import LengthModel

logger = logging.getLogger(__name__)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class LengthExploratoryPrediction:
    """Domain-relaxed numerical result excluded from every primary metric."""

    key: str
    group_key: str
    prediction: float | None
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class LengthDevelopmentAssessment:
    """Local assessment evidence that is permanently ineligible for promotion."""

    evidence_role: str
    promotion_eligible: bool
    candidate_id: str
    model_sha256: str
    study_sha256: str
    protocol_sha256: str
    eligible_roster_sha256: str
    baseline_sha256: str
    input_evidence_sha256: str
    partition_sha256: str
    run_manifest_sha256: str
    evaluation: LengthEvaluationResult
    selection_status: str
    exploratory: tuple[LengthExploratoryPrediction, ...]
    report_sha256: str
    sha256: str


@dataclass(frozen=True)
class LengthAssessmentOutput:
    """Assessment document content and its separately hashed offline report."""

    assessment: LengthDevelopmentAssessment
    report_html: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length assessment {label} must be a lowercase SHA256 digest")
    return value


def _exploratory(
    evidence: LengthRoleEvidence, model: LengthModel, result: LengthEvaluationResult
) -> tuple[LengthExploratoryPrediction, ...]:
    evaluated = next(candidate for candidate in result.candidates if candidate.model_sha256 == model.sha256)
    primary = {prediction.key: prediction for prediction in evaluated.predictions}
    rows: list[LengthExploratoryPrediction] = []
    for row in evidence.rows:
        prediction = primary[row.key]
        if prediction.availability_reasons != ("unsupported_evidence_domain",):
            rows.append(
                LengthExploratoryPrediction(row.key, row.group_key, None, ("domain_relaxation_not_applicable",))
            )
            continue
        relaxed = estimate_total_repeats(
            row.features,
            model,
            evidence_domain=cast(EvidenceDomain, model.applicability.domain),
        )
        rows.append(
            LengthExploratoryPrediction(
                row.key,
                row.group_key,
                relaxed.estimated_total_repeat_count,
                relaxed.reasons,
            )
        )
    return tuple(rows)


def _assessment_payload(assessment: LengthDevelopmentAssessment) -> dict[str, object]:
    return {
        "schema_version": "calibration-assessment-v1",
        "target": "length",
        "evidence_role": assessment.evidence_role,
        "promotion_eligible": assessment.promotion_eligible,
        "candidate_id": assessment.candidate_id,
        "model_sha256": assessment.model_sha256,
        "study_sha256": assessment.study_sha256,
        "protocol_sha256": assessment.protocol_sha256,
        "eligible_roster_sha256": assessment.eligible_roster_sha256,
        "baseline_sha256": assessment.baseline_sha256,
        "input_evidence_sha256": assessment.input_evidence_sha256,
        "partition_sha256": assessment.partition_sha256,
        "run_manifest_sha256": assessment.run_manifest_sha256,
        "evaluation": length_evaluation_document(assessment.evaluation),
        "selection_status": assessment.selection_status,
        "exploratory": [
            {
                "key": row.key,
                "group_key": row.group_key,
                "out_of_domain_prediction": row.prediction,
                "reasons": list(row.reasons),
            }
            for row in assessment.exploratory
        ],
        "report_sha256": assessment.report_sha256,
    }


def length_development_assessment_document(assessment: LengthDevelopmentAssessment) -> dict[str, object]:
    """Project one immutable, permanently non-promotable assessment artifact.

    Args:
        assessment: Result returned by :func:`assess_length_candidate`.

    Returns:
        Closed canonical JSON-compatible local assessment document.

    Raises:
        ValueError: If content, identities or the non-promotion state differ.
    """
    if (
        not isinstance(assessment, LengthDevelopmentAssessment)
        or assessment.evidence_role != "development-assessment"
        or assessment.promotion_eligible is not False
        or not isinstance(assessment.exploratory, tuple)
        or any(
            not isinstance(row, LengthExploratoryPrediction) or not isinstance(row.reasons, tuple)
            for row in assessment.exploratory
        )
    ):
        _fail("length development assessment requires immutable non-promotable content")
    for name in (
        "model_sha256",
        "study_sha256",
        "protocol_sha256",
        "eligible_roster_sha256",
        "baseline_sha256",
        "input_evidence_sha256",
        "partition_sha256",
        "run_manifest_sha256",
        "report_sha256",
        "sha256",
    ):
        _digest(getattr(assessment, name), name)
    if assessment.selection_status != "not-applicable" or assessment.evaluation.selection.status != "not-applicable":
        _fail("length development assessment cannot contain a selection decision")
    if (
        assessment.evaluation.phase != "development-assessment"
        or assessment.evaluation.protocol_sha256 != assessment.protocol_sha256
        or assessment.evaluation.eligible_roster_sha256 != assessment.eligible_roster_sha256
        or assessment.evaluation.baseline_sha256 != assessment.baseline_sha256
    ):
        _fail("length development assessment evaluation bindings differ")
    evaluated = tuple(
        candidate
        for candidate in assessment.evaluation.candidates
        if candidate.candidate_id == assessment.candidate_id and candidate.status == "evaluated"
    )
    if len(evaluated) != 1 or evaluated[0].model_sha256 != assessment.model_sha256:
        _fail("length development assessment candidate differs from its evaluation")
    expected_rows = tuple((row.key, row.group_key) for row in evaluated[0].predictions)
    if tuple((row.key, row.group_key) for row in assessment.exploratory) != expected_rows:
        _fail("length development assessment exploratory rows differ from evaluated rows")
    for row in assessment.exploratory:
        if row.prediction is not None and (
            isinstance(row.prediction, bool)
            or not isinstance(row.prediction, (int, float))
            or not math.isfinite(row.prediction)
            or row.prediction <= 0
            or row.reasons
        ):
            _fail("length development assessment exploratory prediction is invalid")
        if row.prediction is None and not row.reasons:
            _fail("length development assessment unavailable exploratory result requires reasons")
    payload = _assessment_payload(assessment)
    if canonical_sha256(payload) != assessment.sha256:
        _fail("length development assessment digest differs from its canonical content")
    return {**payload, "sha256": assessment.sha256}


def assess_length_candidate(
    profile: LengthTrainingProfile,
    evidence: LengthRoleEvidence,
    roster: LengthEligibleRoster,
    protocol: LengthProtocol,
    *,
    fixed_candidate_id: str,
) -> LengthAssessmentOutput:
    """Evaluate one frozen model on previously examined evidence without fitting.

    Args:
        profile: Context-verified research fit, baseline, study and training roster.
        evidence: Complete development-assessment feature and truth evidence.
        roster: Exact assessment eligible roster.
        protocol: Frozen role-specific assessment policy.
        fixed_candidate_id: Existing fitted candidate to evaluate.

    Returns:
        Immutable non-promotable assessment and its standalone HTML report.

    Raises:
        ValueError: If profile, evidence, role, candidate or bindings differ.
    """
    validate_length_training_profile(profile, protocol)
    training = profile.artifact
    length_role_evidence_document(evidence)
    if evidence.phase != "development-assessment":
        _fail("length assessment evidence must use the development-assessment role")
    if evidence.protocol_sha256 != protocol.sha256 or evidence.eligible_roster_sha256 != roster.sha256:
        _fail("length assessment evidence differs from its protocol or roster")
    if not isinstance(fixed_candidate_id, str) or not fixed_candidate_id:
        _fail("length development assessment requires a fixed candidate ID")
    matching = tuple(outcome for outcome in training.outcomes if outcome.candidate_id == fixed_candidate_id)
    if len(matching) != 1 or matching[0].status != "fitted" or matching[0].model is None:
        _fail("length development assessment fixed candidate must identify one fitted model")
    result = evaluate_length_hypotheses(
        training.outcomes,
        training.baseline,
        roster,
        protocol,
        evidence.rows,
        fixed_candidate_id=fixed_candidate_id,
    )
    report = render_length_evaluation(result, roster, protocol, training.baseline)
    model = matching[0].model
    assessment = LengthDevelopmentAssessment(
        "development-assessment",
        False,
        fixed_candidate_id,
        model.sha256,
        training.study_sha256,
        protocol.sha256,
        roster.sha256,
        training.baseline.sha256,
        evidence.sha256,
        evidence.partition_sha256,
        evidence.run_manifest_sha256,
        result,
        result.selection.status,
        _exploratory(evidence, model, result),
        hashlib.sha256(report.encode("utf-8")).hexdigest(),
        "",
    )
    assessment = replace(assessment, sha256=canonical_sha256(_assessment_payload(assessment)))
    length_development_assessment_document(assessment)
    return LengthAssessmentOutput(assessment, report)


def decode_length_development_assessment(
    value: object,
    *,
    profile: LengthTrainingProfile,
    evidence: LengthRoleEvidence,
    roster: LengthEligibleRoster,
    protocol: LengthProtocol,
    fixed_candidate_id: str,
) -> LengthAssessmentOutput:
    """Recompute and decode a development assessment from opened evidence.

    Args:
        value: Parsed local assessment document.
        profile: Context-verified research training profile.
        evidence: Exact opened development evidence.
        roster: Exact assessment eligible roster.
        protocol: Frozen assessment protocol.
        fixed_candidate_id: Existing fitted candidate to evaluate.

    Returns:
        Recomputed assessment and deterministic offline report.

    Raises:
        ValueError: If stored content differs from recomputed evidence.
    """
    if not isinstance(value, Mapping):
        _fail("length development assessment document must be an object")
    recomputed = assess_length_candidate(
        profile,
        evidence,
        roster,
        protocol,
        fixed_candidate_id=fixed_candidate_id,
    )
    if dict(value) != length_development_assessment_document(recomputed.assessment):
        _fail("length development assessment differs from recomputed evidence")
    return recomputed
