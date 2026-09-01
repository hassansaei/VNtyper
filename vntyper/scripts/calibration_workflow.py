"""Calibration extraction, finite-grid fitting, and validation orchestration."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import cast

from vntyper.scripts.calibration_baseline import build_baseline
from vntyper.scripts.calibration_contract import EvidenceManifest
from vntyper.scripts.calibration_custody import (
    ConsumptionReceipt,
    open_locked_payload,
    retire_candidate,
    write_precommit,
)
from vntyper.scripts.calibration_features import (
    FeatureArtifact,
    LabelArtifact,
    decode_feature_artifact,
    validate_artifact_alignment,
)
from vntyper.scripts.calibration_manifest import StudyDeclaration, require_operation_roles
from vntyper.scripts.calibration_objective import CandidateEvaluation, count_free_parameters, select_candidate
from vntyper.scripts.calibration_profiles import build_generated_profile, validate_generated_allowlist
from vntyper.scripts.calibration_run_extraction import (
    decode_run_artifact_declaration,
    extract_completed_run,
)
from vntyper.scripts.calibration_statistics import holm_adjust
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.decision_profile import ResolvedDecisionProfile, load_packaged_decision_profile
from vntyper.version import __version__


@dataclass(frozen=True)
class ExtractedEvidence:
    """Immutable replay evidence with complete study and run hashes."""

    study: StudyDeclaration
    features: FeatureArtifact
    labels: LabelArtifact
    baseline: Mapping[str, object]
    run_hashes: Mapping[str, Mapping[str, str]]
    study_sha256: str
    dataset_sha256: str


@dataclass(frozen=True)
class FittedCandidate:
    """Generated profile selected from the protocol-frozen family."""

    profile: ResolvedDecisionProfile
    accessed_roles: tuple[str, ...]
    baseline_reproduced: bool
    evaluation: CandidateEvaluation


@dataclass(frozen=True)
class WorkflowAttestation:
    """Profile-bound result for one non-selection validation role."""

    role: str
    profile_sha256: str
    evidence_sha256: str
    accessed_roles: tuple[str, ...]


@dataclass(frozen=True)
class LockedEvaluation:
    """One profile-bound evaluation performed on externally held evidence."""

    role: str
    result: object
    profile_sha256: str
    evidence_sha256: str
    receipt: ConsumptionReceipt


def extract_evidence(
    study: StudyDeclaration,
    labels: LabelArtifact,
    runs: Mapping[str, object],
) -> ExtractedEvidence:
    """Derive features and shipped replay from complete immutable run artifacts."""
    if not isinstance(study, StudyDeclaration):
        raise ValueError("calibration extraction requires a StudyDeclaration")
    if not isinstance(runs, Mapping) or set(runs) != {member.key for member in study.partitions.members}:
        raise ValueError("calibration run roots must match the complete partition manifest")
    if not isinstance(labels, LabelArtifact):
        raise ValueError("calibration extraction requires a decoded label artifact")
    label_keys = tuple(row.manifest_key for row in labels.rows)
    manifest_keys = tuple(member.key for member in study.partitions.members)
    if label_keys != manifest_keys:
        raise ValueError("calibration label rows must align exactly with the partition manifest")
    extractions = []
    run_hashes: dict[str, Mapping[str, str]] = {}
    members = {member.key: member for member in study.partitions.members}
    for key in sorted(runs):
        declaration = decode_run_artifact_declaration(runs[key])
        extracted = extract_completed_run(key, members[key].assay_class, declaration)
        extractions.append(extracted)
        run_hashes[key] = extracted.artifact_sha256
    features = decode_feature_artifact(
        {
            "schema_version": "calibration-features-v1",
            "rows": [
                {
                    "feature_key": f"run-{extraction.manifest_key}",
                    "manifest_key": extraction.manifest_key,
                    "features": dict(extraction.features),
                }
                for extraction in extractions
            ],
        }
    )
    validate_artifact_alignment(features, labels, study.partitions)
    baseline = build_baseline(extractions, {row.manifest_key: row for row in labels.rows})
    baseline_copy = cast(Mapping[str, object], _freeze(baseline))
    dataset_sha256 = canonical_sha256(
        {
            "study_sha256": study.sha256,
            "features_sha256": features.sha256,
            "labels_sha256": labels.sha256,
            "baseline_sha256": canonical_sha256(baseline),
            "run_artifact_sha256": {key: dict(value) for key, value in run_hashes.items()},
        }
    )
    return ExtractedEvidence(
        study,
        features,
        labels,
        baseline_copy,
        MappingProxyType({key: MappingProxyType(dict(value)) for key, value in run_hashes.items()}),
        study.sha256,
        dataset_sha256,
    )


def fit_candidate(
    evidence: ExtractedEvidence,
    *,
    objective: str,
    evaluator: Callable[[ResolvedDecisionProfile], CandidateEvaluation],
) -> FittedCandidate:
    """Replay the baseline, evaluate the frozen grid, and emit the selected profile."""
    if not isinstance(evidence, ExtractedEvidence):
        raise ValueError("calibration fitting requires ExtractedEvidence")
    if objective != evidence.study.protocol.objective:
        raise ValueError("calibration fit objective differs from the snapshotted protocol")
    if not callable(evaluator):
        raise ValueError("calibration fit requires a candidate evaluator")
    accessed = require_operation_roles(evidence.study.partitions, "fit")
    _validate_baseline_replay(evidence.baseline)
    candidates: list[tuple[ResolvedDecisionProfile, CandidateEvaluation]] = []
    for candidate in evidence.study.protocol.candidates:
        component = candidate.as_component()
        if count_free_parameters(component) > evidence.study.protocol.maximum_free_parameters:
            continue
        profile = build_generated_profile(
            component,
            dataset_manifest_hash=evidence.dataset_sha256,
            partition_manifest_hash=evidence.study.partitions.sha256,
            seed=evidence.study.protocol.seed,
            objective=objective,
            generator_version=__version__,
        )
        evaluation = evaluator(profile)
        if (
            not isinstance(evaluation, CandidateEvaluation)
            or evaluation.metrics.candidate_profile_sha256 != profile.digest
        ):
            raise ValueError("calibration evaluator result does not bind the candidate profile hash")
        candidates.append((profile, evaluation))
    if not candidates:
        raise ValueError("calibration protocol contains no candidate within its free-parameter limit")
    adjusted = holm_adjust({profile.digest: evaluation.holm_adjusted_p_value for profile, evaluation in candidates})
    corrected = tuple(
        (profile, replace(evaluation, holm_adjusted_p_value=adjusted[profile.digest]))
        for profile, evaluation in candidates
    )
    selected = select_candidate(tuple(evaluation for _, evaluation in corrected), evidence.study.protocol)
    if selected is None:
        raise ValueError("calibration fitting found no admissible candidate")
    profile = next(profile for profile, evaluation in corrected if evaluation is selected)
    return FittedCandidate(profile, accessed, True, selected)


def validate_candidate(profile: ResolvedDecisionProfile, evidence: ExtractedEvidence) -> WorkflowAttestation:
    """Validate one fixed profile on only the declared validation role."""
    if not isinstance(evidence, ExtractedEvidence):
        raise ValueError("calibration validation requires ExtractedEvidence")
    packaged = load_packaged_decision_profile()
    validate_generated_allowlist(profile, packaged)
    accessed = require_operation_roles(evidence.study.partitions, "validate")
    return WorkflowAttestation("validation", profile.digest, evidence.dataset_sha256, accessed)


def evaluate_locked_candidate(
    profile: ResolvedDecisionProfile,
    payload_path: Path,
    evidence: EvidenceManifest,
    custody_dir: Path,
    *,
    evaluator: Callable[[bytes], object],
) -> LockedEvaluation:
    """Precommit and consume externally held-out evidence exactly once."""
    if not isinstance(evidence, EvidenceManifest):
        raise ValueError("locked calibration evaluation requires an EvidenceManifest")
    if evidence.role != "locked-heldout" or evidence.provenance != "external-custodian":
        raise ValueError("locked held-out evaluation requires external custodian evidence")
    if not callable(evaluator):
        raise ValueError("locked calibration evaluation requires an evaluator")
    packaged = load_packaged_decision_profile()
    validate_generated_allowlist(profile, packaged)
    precommit = write_precommit(
        custody_dir,
        profile.digest,
        evidence.protocol_sha256,
        evidence.features_sha256,
    )
    try:
        opened = open_locked_payload(payload_path, precommit, custody_dir)
        result = evaluator(opened.payload)
    except Exception as error:
        retire_candidate(custody_dir, profile.digest, evidence.features_sha256, str(error) or type(error).__name__)
        raise
    return LockedEvaluation(evidence.role, result, profile.digest, evidence.features_sha256, opened.receipt)


def _validate_baseline_replay(value: Mapping[str, object]) -> None:
    if set(value) != {"schema_version", "expected", "observed"}:
        raise ValueError("calibration baseline replay fields differ from the closed contract")
    if value["schema_version"] != "calibration-baseline-replay-v1":
        raise ValueError("calibration baseline replay schema version is unsupported")
    expected = value["expected"]
    observed = value["observed"]
    if not isinstance(expected, Mapping) or not isinstance(observed, Mapping):
        raise ValueError("calibration baseline replay projections must be objects")
    required = {"aggregate", "per_tier", "rows"}
    if set(expected) != required or set(observed) != required:
        raise ValueError("calibration baseline replay must contain aggregate, per-tier, and row projections")
    if expected != observed:
        raise ValueError("calibration fitting requires exact shipped baseline reproduction")
    rows = expected["rows"]
    if not isinstance(rows, tuple) or not rows:
        raise ValueError("calibration baseline replay rows must be a non-empty immutable sequence")
    row_fields = {
        "manifest_key",
        "order",
        "canonical_identity",
        "name",
        "confidence",
        "flag",
        "tier",
        "support",
        "tie",
        "abstention",
    }
    if any(not isinstance(row, Mapping) or set(row) != row_fields for row in rows):
        raise ValueError("calibration baseline replay rows lack required decision fields")


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(child) for child in value)
    return value
