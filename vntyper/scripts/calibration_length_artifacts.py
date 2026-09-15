"""Bound training artifacts and development-only assessment for length models.

The row evidence in this module is a local projection over already opened typed
features. It does not authorize access to a role, confer custody, or make a model
portable for the production pipeline.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_length import (
    LengthBaseline,
    LengthFitOutcome,
    LengthTrainingMetadata,
    decode_length_baseline,
    decode_length_training_metadata,
    fit_length_hypotheses,
    length_baseline_document,
    length_training_metadata_document,
)
from vntyper.scripts.calibration_length_evidence import (
    LengthTrainingEvidence,
    length_training_evidence_document,
)
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    length_eligible_roster_document,
)
from vntyper.scripts.calibration_length_protocol import LengthProtocol, length_protocol_document
from vntyper.scripts.calibration_target_contract import (
    LengthBaselinePlan,
    TargetStudy,
    target_study_document,
)
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_model import LengthModel, decode_length_model, encode_length_model

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_TRAINING_FIELDS = {
    "schema_version",
    "study_sha256",
    "protocol_sha256",
    "training_roster_sha256",
    "metadata",
    "baseline",
    "outcomes",
    "sha256",
}


@dataclass(frozen=True)
class LengthTrainingArtifact:
    """Strict local training result bound to its study baseline plan."""

    study_sha256: str
    protocol_sha256: str
    training_roster_sha256: str
    metadata: LengthTrainingMetadata
    baseline: LengthBaseline
    outcomes: tuple[LengthFitOutcome, ...]
    sha256: str


@dataclass(frozen=True)
class LengthTrainingProfile:
    """Opened training artifact with the study and roster needed to verify it."""

    artifact: LengthTrainingArtifact
    study: TargetStudy
    training_roster: LengthEligibleRoster


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"length artifact {label} fields differ from the closed contract")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"length artifact {label} must be a lowercase SHA256 digest")
    return value


def _study_bindings(study: TargetStudy, metadata: LengthTrainingMetadata) -> tuple[LengthProtocol, LengthBaselinePlan]:
    target_study_document(study)
    length_training_metadata_document(metadata)
    if study.target != "length" or not isinstance(study.protocol, LengthProtocol):
        _fail("length training artifact requires a length target study")
    if not isinstance(study.baseline, LengthBaselinePlan):
        _fail("length training artifact requires a length baseline plan")
    baseline = study.baseline
    if metadata.study_sha256 != study.sha256:
        _fail("length training metadata study digest differs from the target study")
    if metadata.annotation_sha256 != baseline.annotation_sha256:
        _fail("length training metadata annotation differs from the study baseline plan")
    if metadata.counting_policy_sha256 != baseline.counting_policy_sha256:
        _fail("length training metadata counting policy differs from the study baseline plan")
    if metadata.maximum_condition_number != baseline.maximum_condition_number:
        _fail("length training metadata condition limit differs from the study baseline plan")
    if metadata.producer != baseline.producer:
        _fail("length training metadata producer differs from the study baseline plan")
    if metadata.applicability != study.applicability:
        _fail("length training metadata applicability differs from the target study")
    if dict(metadata.qc) != dict(study.protocol.qc):
        _fail("length training metadata QC differs from the target study protocol")
    return study.protocol, baseline


def _outcome_document(outcome: LengthFitOutcome) -> dict[str, object]:
    return {
        "candidate_id": outcome.candidate_id,
        "model_kind": outcome.model_kind,
        "status": outcome.status,
        "reasons": list(outcome.reasons),
        "model": None if outcome.model is None else encode_length_model(outcome.model),
    }


def _decode_outcomes(
    value: object, protocol: LengthProtocol, metadata: LengthTrainingMetadata
) -> tuple[LengthFitOutcome, ...]:
    if not isinstance(value, list) or len(value) != len(protocol.candidates):
        _fail("length training artifact outcomes must match the complete protocol roster")
    outcomes: list[LengthFitOutcome] = []
    for raw_value, hypothesis in zip(value, protocol.candidates, strict=True):
        raw = _object(raw_value, {"candidate_id", "model_kind", "status", "reasons", "model"}, "outcome")
        if raw["candidate_id"] != hypothesis.candidate_id or raw["model_kind"] != hypothesis.model_kind:
            _fail("length training artifact outcome differs from its protocol hypothesis")
        reasons_value = raw["reasons"]
        if not isinstance(reasons_value, list) or any(
            not isinstance(reason, str) or not reason or reason != reason.strip() for reason in reasons_value
        ):
            _fail("length training artifact outcome reasons must be non-empty strings")
        reasons = tuple(reasons_value)
        status = raw["status"]
        model_value = raw["model"]
        model: LengthModel | None = None
        if status == "fitted":
            if reasons or model_value is None:
                _fail("length training fitted outcome must have one model and no reasons")
            model = decode_length_model(model_value)
            if (
                model.model_kind != hypothesis.model_kind
                or model.study_sha256 != metadata.study_sha256
                or model.training_evidence_sha256 != metadata.training_evidence_sha256
                or model.annotation_sha256 != metadata.annotation_sha256
                or model.counting_policy_sha256 != metadata.counting_policy_sha256
                or model.applicability != metadata.applicability
                or model.producer != metadata.producer
                or any(model.qc[name] != metadata.qc[name] for name in metadata.qc)
            ):
                _fail("length training fitted model differs from its metadata bindings")
        elif status == "ineligible":
            if not reasons or model_value is not None:
                _fail("length training ineligible outcome must have reasons and no model")
        else:
            _fail("length training artifact outcome status is unsupported")
        outcomes.append(
            LengthFitOutcome(
                hypothesis.candidate_id,
                hypothesis.model_kind,
                cast(Literal["fitted", "ineligible"], status),
                reasons,
                model,
            )
        )
    return tuple(outcomes)


def _training_payload(artifact: LengthTrainingArtifact) -> dict[str, object]:
    return {
        "schema_version": "calibration-length-training-artifact-v1",
        "study_sha256": artifact.study_sha256,
        "protocol_sha256": artifact.protocol_sha256,
        "training_roster_sha256": artifact.training_roster_sha256,
        "metadata": length_training_metadata_document(artifact.metadata),
        "baseline": length_baseline_document(artifact.baseline),
        "outcomes": [_outcome_document(outcome) for outcome in artifact.outcomes],
    }


def _validate_training_artifact(
    artifact: LengthTrainingArtifact, study: TargetStudy, training_roster: LengthEligibleRoster
) -> None:
    protocol, _ = _study_bindings(study, artifact.metadata)
    length_eligible_roster_document(training_roster)
    if (
        artifact.study_sha256 != study.sha256
        or artifact.protocol_sha256 != protocol.sha256
        or artifact.training_roster_sha256 != training_roster.sha256
        or artifact.metadata.training_roster_sha256 != training_roster.sha256
    ):
        _fail("length training artifact study, protocol, or roster binding differs")
    baseline = decode_length_baseline(length_baseline_document(artifact.baseline))
    if (
        baseline.study_sha256 != study.sha256
        or baseline.training_evidence_sha256 != artifact.metadata.training_evidence_sha256
        or baseline.training_roster_sha256 != training_roster.sha256
        or baseline.independent_group_count != len(training_roster.members)
    ):
        _fail("length training artifact baseline differs from its metadata or roster")
    decoded_outcomes = _decode_outcomes(
        [_outcome_document(outcome) for outcome in artifact.outcomes], protocol, artifact.metadata
    )
    if decoded_outcomes != artifact.outcomes:
        _fail("length training artifact outcomes differ from their canonical content")
    if canonical_sha256(_training_payload(artifact)) != artifact.sha256:
        _fail("length training artifact digest differs from its canonical content")


def build_length_training_artifact(
    study: TargetStudy,
    training_roster: LengthEligibleRoster,
    metadata: LengthTrainingMetadata,
    evidence: LengthTrainingEvidence,
) -> LengthTrainingArtifact:
    """Fit all hypotheses after binding metadata to the study baseline plan.

    Args:
        study: Opened length target study declared before outcomes.
        training_roster: Exact separately frozen training representatives.
        metadata: Training evidence, applicability, QC and numerical bindings.
        evidence: Complete hash-bound training features and exact matched truth.

    Returns:
        Immutable local fit artifact containing every declared outcome.

    Raises:
        ValueError: If the study, metadata, roster, rows or fitted output differ.
    """
    protocol, _ = _study_bindings(study, metadata)
    length_training_evidence_document(evidence)
    if (
        evidence.study_sha256 != study.sha256
        or evidence.training_roster_sha256 != training_roster.sha256
        or metadata.training_evidence_sha256 != evidence.sha256
    ):
        _fail("length training evidence differs from its study, roster, or metadata")
    result = fit_length_hypotheses(evidence.rows, training_roster, protocol, metadata)
    artifact = LengthTrainingArtifact(
        study.sha256,
        result.protocol_sha256,
        result.training_roster_sha256,
        metadata,
        result.baseline,
        result.outcomes,
        "",
    )
    artifact = replace(artifact, sha256=canonical_sha256(_training_payload(artifact)))
    _validate_training_artifact(artifact, study, training_roster)
    return artifact


def length_training_artifact_document(
    artifact: LengthTrainingArtifact, *, study: TargetStudy, training_roster: LengthEligibleRoster
) -> dict[str, object]:
    """Project a typed training artifact after contextual validation.

    Args:
        artifact: Local training artifact to project.
        study: Exact opened target study that authorized fitting.
        training_roster: Exact opened training roster.

    Returns:
        Closed canonical JSON-compatible artifact document.

    Raises:
        ValueError: If content, hashes or study bindings differ.
    """
    if not isinstance(artifact, LengthTrainingArtifact):
        _fail("length training artifact projection requires a typed artifact")
    _validate_training_artifact(artifact, study, training_roster)
    payload = _training_payload(artifact)
    return {**payload, "sha256": artifact.sha256}


def decode_length_training_artifact(
    value: object, *, study: TargetStudy, training_roster: LengthEligibleRoster
) -> LengthTrainingArtifact:
    """Decode a local training artifact against its opened study and roster.

    Args:
        value: Parsed closed training artifact document.
        study: Exact opened target study that authorized fitting.
        training_roster: Exact opened training roster.

    Returns:
        Immutable verified local training artifact.

    Raises:
        ValueError: If schema, content, hashes or contextual bindings differ.
    """
    raw = _object(value, _TRAINING_FIELDS, "training root")
    if raw["schema_version"] != "calibration-length-training-artifact-v1":
        _fail("length training artifact schema is unsupported")
    metadata = decode_length_training_metadata(raw["metadata"], training_roster)
    protocol, _ = _study_bindings(study, metadata)
    artifact = LengthTrainingArtifact(
        _digest(raw["study_sha256"], "study digest"),
        _digest(raw["protocol_sha256"], "protocol digest"),
        _digest(raw["training_roster_sha256"], "training roster digest"),
        metadata,
        decode_length_baseline(raw["baseline"]),
        _decode_outcomes(raw["outcomes"], protocol, metadata),
        _digest(raw["sha256"], "training artifact digest"),
    )
    _validate_training_artifact(artifact, study, training_roster)
    return artifact


def bind_length_training_profile(
    artifact: LengthTrainingArtifact, study: TargetStudy, training_roster: LengthEligibleRoster
) -> LengthTrainingProfile:
    """Bind an opened artifact to the study and roster needed to verify it.

    Args:
        artifact: Decoded local training artifact.
        study: Exact opened target study that authorized fitting.
        training_roster: Exact opened training roster.

    Returns:
        Immutable contextual profile for later fixed-model evaluation.

    Raises:
        ValueError: If any artifact content or contextual binding differs.
    """
    if not isinstance(artifact, LengthTrainingArtifact):
        _fail("length training profile requires a typed training artifact")
    _validate_training_artifact(artifact, study, training_roster)
    return LengthTrainingProfile(artifact, study, training_roster)


def validate_length_training_profile(profile: LengthTrainingProfile, protocol: LengthProtocol) -> None:
    """Validate an opened profile for a role-specific compatible protocol.

    Args:
        profile: Context-bound training profile.
        protocol: Predeclared role-specific protocol used for evaluation.

    Raises:
        ValueError: If profile context, candidates or shared QC differ.
    """
    if not isinstance(profile, LengthTrainingProfile):
        _fail("length training profile requires a contextual typed profile")
    artifact = profile.artifact
    _validate_training_artifact(artifact, profile.study, profile.training_roster)
    length_protocol_document(protocol)
    if dict(artifact.metadata.qc) != dict(protocol.qc):
        _fail("length training profile bindings differ")
    outcomes = _decode_outcomes(
        [_outcome_document(outcome) for outcome in artifact.outcomes], protocol, artifact.metadata
    )
    if outcomes != artifact.outcomes or canonical_sha256(_training_payload(artifact)) != artifact.sha256:
        _fail("length training profile differs from its canonical content")
