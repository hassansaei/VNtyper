"""Calibration extraction, finite-grid fitting, and validation orchestration."""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import cast

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
    validate_artifact_alignment,
)
from vntyper.scripts.calibration_manifest import StudyDeclaration, require_operation_roles
from vntyper.scripts.calibration_objective import CandidateEvaluation, count_free_parameters, select_candidate
from vntyper.scripts.calibration_profiles import build_generated_profile, validate_generated_allowlist
from vntyper.scripts.calibration_statistics import holm_adjust
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import ResolvedDecisionProfile, load_packaged_decision_profile
from vntyper.scripts.identity_candidate_persistence import IDENTITY_CAPTURE_COLUMNS
from vntyper.scripts.nomenclature_bam_replay import decode_bam_replay_artifact
from vntyper.version import __version__

_REQUIRED_RUN_ARTIFACTS = (
    Path("pipeline_summary.json"),
    Path("kestrel/kestrel_pre_result.tsv"),
    Path("kestrel/bam_identity_replay.v1.json"),
)


@dataclass(frozen=True)
class ExtractedEvidence:
    """Immutable replay evidence with complete study and run hashes."""

    study: StudyDeclaration
    features: FeatureArtifact
    labels: LabelArtifact
    baseline: Mapping[str, object]
    run_hashes: Mapping[str, str]
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
    features: FeatureArtifact,
    labels: LabelArtifact,
    runs: Mapping[str, Path],
    *,
    baseline: Mapping[str, object],
) -> ExtractedEvidence:
    """Hash complete replay artifacts without opening BAMs or rerunning callers."""
    if not isinstance(study, StudyDeclaration):
        raise ValueError("calibration extraction requires a StudyDeclaration")
    validate_artifact_alignment(features, labels, study.partitions)
    if not isinstance(runs, Mapping) or set(runs) != {member.key for member in study.partitions.members}:
        raise ValueError("calibration run roots must match the complete partition manifest")
    if not isinstance(baseline, Mapping) or not baseline:
        raise ValueError("calibration extraction requires the shipped decision projection baseline")
    run_hashes: dict[str, str] = {}
    for key in sorted(runs):
        root = runs[key]
        if not isinstance(root, Path):
            raise ValueError("calibration run roots must be Path values")
        run_hashes[key] = _hash_replay_run(root)
    baseline_copy = cast(Mapping[str, object], _freeze(dict(baseline)))
    dataset_sha256 = canonical_sha256(
        {
            "features_sha256": features.sha256,
            "labels_sha256": labels.sha256,
            "run_hashes": run_hashes,
        }
    )
    return ExtractedEvidence(
        study,
        features,
        labels,
        baseline_copy,
        MappingProxyType(run_hashes),
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


def _hash_replay_run(root: Path) -> str:
    artifact_hashes: dict[str, str] = {}
    for relative in _REQUIRED_RUN_ARTIFACTS:
        path = root / relative
        try:
            raw = path.read_bytes()
        except OSError as error:
            raise ValueError(f"calibration replay artifact is missing or unreadable: {path}") from error
        artifact_hashes[relative.as_posix()] = hashlib.sha256(raw).hexdigest()
    summary = load_strict_json_object((root / _REQUIRED_RUN_ARTIFACTS[0]).read_bytes())
    if summary.get("schema_version") != 3:
        raise ValueError("calibration replay requires pipeline summary schema 3")
    pre_result_lines = (root / _REQUIRED_RUN_ARTIFACTS[1]).read_text(encoding="utf-8").splitlines()
    if not pre_result_lines:
        raise ValueError("calibration replay pre-result is empty")
    pre_result_header = pre_result_lines[0].split("\t")
    if IDENTITY_CAPTURE_COLUMNS[5] not in pre_result_header:
        raise ValueError("calibration replay pre-result lacks complete PR-A identity artifacts")
    replay = load_strict_json_object((root / _REQUIRED_RUN_ARTIFACTS[2]).read_bytes())
    decode_bam_replay_artifact(replay)
    return canonical_sha256(artifact_hashes)


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
    row_fields = {"order", "name", "confidence", "flag", "tier", "support", "tie", "abstention"}
    if any(not isinstance(row, Mapping) or set(row) != row_fields for row in rows):
        raise ValueError("calibration baseline replay rows lack required decision fields")


def _freeze(value: object) -> object:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(child) for child in value)
    return value
