"""Verified local research profile loading for length calibration."""

from __future__ import annotations

import hashlib
import json
import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_artifact_io import load_object, verify_checksums
from vntyper.scripts.calibration_candidate import CandidateEnvelope, decode_candidate, validate_candidate_payload
from vntyper.scripts.calibration_length_artifacts import (
    LengthTrainingProfile,
    bind_length_training_profile,
    decode_length_training_artifact,
)
from vntyper.scripts.calibration_length_metrics import decode_length_eligible_roster
from vntyper.scripts.calibration_payload import PayloadManifest, decode_payload_manifest, validate_payload_observations
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.calibration_target_contract import LengthBaselinePlan, TargetStudy, decode_target_study
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation, encode_length_annotation
from vntyper.scripts.length_model import LengthModel, decode_length_model, encode_length_model

logger = logging.getLogger(__name__)

_RESEARCH_FILES = {
    "candidate.json",
    "checksums.json",
    "evaluation.json",
    "payload",
    "payload-manifest.json",
    "report.html",
    "selection-evidence.json",
    "study.json",
    "training-artifact.json",
    "training-evidence.json",
    "training-roster.json",
}


@dataclass(frozen=True)
class LengthResearchProfile:
    """Verified local research candidate and its training context."""

    candidate: CandidateEnvelope
    study: TargetStudy
    training_profile: LengthTrainingProfile
    selected_protocol_candidate_id: str
    payload_manifest: PayloadManifest
    payload_files: Mapping[str, bytes]
    model: LengthModel
    annotation: LengthAnnotation
    projection_sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _checked_local_document(path: Path, label: str) -> dict[str, object]:
    raw = read_regular_path(path)
    value = load_strict_json_object(raw)
    if canonical_json_bytes(value) != raw:
        _fail(f"{label} must use canonical JSON bytes")
    return value


def _checked_local_json(path: Path, label: str) -> object:
    raw = read_regular_path(path)
    try:
        value = json.loads(raw)
        if canonical_json_bytes(value) != raw:
            _fail(f"{label} must use canonical JSON bytes")
        return value
    except (UnicodeDecodeError, ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{label} is invalid") from error


def build_length_payload(model: LengthModel, annotation: LengthAnnotation) -> tuple[PayloadManifest, bytes, bytes]:
    """Build the exact two-file research payload and its byte manifest.

    Args:
        model: Selected fitted research model.
        annotation: Exact measurement geometry bound by the model.

    Returns:
        Typed payload manifest, annotation bytes, and model bytes.
    """
    model_raw = canonical_json_bytes(encode_length_model(model))
    annotation_raw = canonical_json_bytes(encode_length_annotation(annotation))
    manifest = decode_payload_manifest(
        [
            {
                "path": "length-annotation.json",
                "size_bytes": len(annotation_raw),
                "sha256": hashlib.sha256(annotation_raw).hexdigest(),
            },
            {
                "path": "length-model.json",
                "size_bytes": len(model_raw),
                "sha256": hashlib.sha256(model_raw).hexdigest(),
            },
        ]
    )
    return manifest, annotation_raw, model_raw


def load_length_research_profile(profile_dir: Path) -> LengthResearchProfile:
    """Load one research-only fit profile and verify its exact local bindings.

    Args:
        profile_dir: Fit output containing the candidate, fit evidence, and payload.

    Returns:
        Immutable typed research profile with verified original payload bytes.

    Raises:
        ValueError: If inventory, checksums, schemas, hashes, or model bindings differ.
    """
    if not isinstance(profile_dir, Path) or not profile_dir.is_dir() or profile_dir.is_symlink():
        _fail("length research profile must be a nonsymlink directory")
    if {path.name for path in profile_dir.iterdir()} != _RESEARCH_FILES:
        _fail("length research profile inventory differs")
    verify_checksums(profile_dir)
    checksums = load_object(profile_dir / "checksums.json", "length research checksums")
    direct_files = _RESEARCH_FILES - {"checksums.json", "payload"}
    checksum_files = checksums.get("files")
    if not isinstance(checksum_files, Mapping) or set(checksum_files) != direct_files:
        _fail("length research profile checksum inventory differs")
    study = decode_target_study(_checked_local_document(profile_dir / "study.json", "length research study"))
    if study.target != "length" or not isinstance(study.baseline, LengthBaselinePlan):
        _fail("length research profile requires a length target study")
    roster = decode_length_eligible_roster(_checked_local_json(profile_dir / "training-roster.json", "length roster"))
    artifact = decode_length_training_artifact(
        _checked_local_document(profile_dir / "training-artifact.json", "length training artifact"),
        study=study,
        training_roster=roster,
    )
    training_profile = bind_length_training_profile(artifact, study, roster)
    training_evidence = _checked_local_document(profile_dir / "training-evidence.json", "length training evidence")
    if training_evidence.get("sha256") != artifact.metadata.training_evidence_sha256 or canonical_sha256(
        {key: value for key, value in training_evidence.items() if key != "sha256"}
    ) != training_evidence.get("sha256"):
        _fail("length research training evidence digest differs")

    payload_root = profile_dir / "payload"
    if not payload_root.is_dir() or payload_root.is_symlink():
        _fail("length research payload must be a nonsymlink directory")
    payload_names = {"length-annotation.json", "length-model.json"}
    if {path.name for path in payload_root.iterdir()} != payload_names:
        _fail("length research payload inventory differs")
    payload_files = {name: read_regular_path(payload_root / name) for name in sorted(payload_names)}
    manifest = decode_payload_manifest(
        _checked_local_json(profile_dir / "payload-manifest.json", "length research payload manifest")
    )
    observed = {name: (len(raw), hashlib.sha256(raw).hexdigest()) for name, raw in payload_files.items()}
    validate_payload_observations(manifest, observed)
    candidate = decode_candidate(_checked_local_document(profile_dir / "candidate.json", "length candidate"))
    validate_candidate_payload(candidate, manifest, expected_target="length")
    annotation = decode_length_annotation(load_strict_json_object(payload_files["length-annotation.json"]))
    model = decode_length_model(load_strict_json_object(payload_files["length-model.json"]))
    evaluation = _checked_local_document(profile_dir / "evaluation.json", "length selection evaluation")
    evaluation_sha = evaluation.get("sha256")
    if (
        evaluation_sha != candidate.selection_evidence_sha256
        or canonical_sha256({key: value for key, value in evaluation.items() if key != "sha256"}) != evaluation_sha
    ):
        _fail("length research selection evidence digest differs")
    selection = evaluation.get("selection")
    if not isinstance(selection, Mapping) or set(selection) != {
        "status",
        "selected_candidate_id",
        "selected_model_sha256",
        "reasons",
    }:
        _fail("length research selection fields differ")
    selected_id = selection["selected_candidate_id"]
    if (
        selection["status"] != "selected"
        or not isinstance(selected_id, str)
        or selection["selected_model_sha256"] != model.sha256
        or selection["reasons"] != []
    ):
        _fail("length research selection does not identify the payload model")
    fitted = tuple(outcome for outcome in artifact.outcomes if outcome.candidate_id == selected_id)
    if len(fitted) != 1 or fitted[0].model != model:
        _fail("length research payload model differs from its fitted outcome")
    if (
        candidate.study_sha256 != study.sha256
        or candidate.baseline_sha256 != artifact.baseline.sha256
        or candidate.partition_sha256 != study.partitions.sha256
        or candidate.training_evidence_sha256 != artifact.metadata.training_evidence_sha256
        or candidate.applicability != model.applicability
        or candidate.producer != model.producer
        or model.annotation_sha256 != annotation.sha256
    ):
        _fail("length research candidate, model, annotation, or training bindings differ")
    projection = {
        "schema_version": "length-research-profile-projection-v1",
        "candidate_sha256": candidate.sha256,
        "study_sha256": study.sha256,
        "training_artifact_sha256": artifact.sha256,
        "selected_protocol_candidate_id": selected_id,
        "payload_sha256": manifest.sha256,
        "selection_evidence_sha256": evaluation_sha,
    }
    return LengthResearchProfile(
        candidate,
        study,
        training_profile,
        selected_id,
        manifest,
        MappingProxyType(payload_files),
        model,
        annotation,
        canonical_sha256(projection),
    )
