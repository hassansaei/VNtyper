"""Strict no-follow loader for approved portable VNTR length model bundles."""

from __future__ import annotations

import hashlib
import json
import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from vntyper.scripts.calibration_candidate import CandidateEnvelope, decode_candidate, validate_candidate_payload
from vntyper.scripts.calibration_payload import (
    PayloadManifest,
    decode_payload_manifest,
    validate_payload_observations,
)
from vntyper.scripts.calibration_portable_projection import (
    PortableApproval,
    decode_portable_approval,
    validate_portable_approval_candidate,
)
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation
from vntyper.scripts.length_model import LengthModel, decode_length_model

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_BUNDLE_FILES = {
    "candidate.json",
    "payload-manifest.json",
    "portable-approval.json",
    "length-model.json",
    "length-annotation.json",
    "checksums.json",
}
_PAYLOAD_FILES = {"length-annotation.json", "length-model.json"}


@dataclass(frozen=True)
class LengthModelBundle:
    """Approved runtime model, annotation, and their transitive identities."""

    candidate: CandidateEnvelope
    payload: PayloadManifest
    approval: PortableApproval
    model: LengthModel
    annotation: LengthAnnotation
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _canonical_object(raw: bytes, label: str) -> dict[str, object]:
    try:
        value = load_strict_json_object(raw)
        if canonical_json_bytes(value) != raw:
            _fail(f"runtime length bundle {label} must use canonical JSON bytes")
        return value
    except (UnicodeDecodeError, ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"runtime length bundle {label} is invalid: {error}") from error


def _reject_nonfinite(token: str) -> NoReturn:
    _fail(f"runtime length bundle payload manifest contains non-finite JSON constant {token}")


def _canonical_manifest(raw: bytes) -> object:
    try:
        value = json.loads(raw, parse_constant=_reject_nonfinite)
        if canonical_json_bytes(value) != raw:
            _fail("runtime length bundle payload manifest must use canonical JSON bytes")
        return value
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"runtime length bundle payload manifest is invalid: {error}") from error


def _validate_checksums(raw_files: Mapping[str, bytes]) -> str:
    document = _canonical_object(raw_files["checksums.json"], "checksums")
    if set(document) != {"schema_version", "files"} or document["schema_version"] != "calibration-checksums-v1":
        _fail("runtime length bundle checksum fields or schema differ")
    checksums = document["files"]
    expected_names = _BUNDLE_FILES - {"checksums.json"}
    if not isinstance(checksums, Mapping) or set(checksums) != expected_names:
        _fail("runtime length bundle checksum inventory differs")
    for name in expected_names:
        expected = checksums[name]
        if not isinstance(expected, str) or _SHA256.fullmatch(expected) is None:
            _fail("runtime length bundle checksum must be lowercase SHA-256")
        if hashlib.sha256(raw_files[name]).hexdigest() != expected:
            _fail(f"runtime length bundle checksum differs for {name}")
    return canonical_sha256(document)


def _decode_bundle(raw: Mapping[str, bytes], bundle_sha256: str) -> LengthModelBundle:
    manifest = decode_payload_manifest(_canonical_manifest(raw["payload-manifest.json"]))
    if {item.path for item in manifest.files} != _PAYLOAD_FILES:
        _fail("runtime length payload must contain exactly its annotation and model")
    validate_payload_observations(
        manifest,
        {
            name: (len(raw[name]), hashlib.sha256(raw[name]).hexdigest())
            for name in ("length-annotation.json", "length-model.json")
        },
    )
    candidate = decode_candidate(_canonical_object(raw["candidate.json"], "candidate"))
    approval = decode_portable_approval(_canonical_object(raw["portable-approval.json"], "portable approval"))
    model = decode_length_model(_canonical_object(raw["length-model.json"], "model"))
    annotation = decode_length_annotation(_canonical_object(raw["length-annotation.json"], "annotation"))
    validate_candidate_payload(candidate, manifest, expected_target="length")
    validate_portable_approval_candidate(approval, candidate)
    _validate_model_bindings(candidate, model, annotation)
    return LengthModelBundle(candidate, manifest, approval, model, annotation, bundle_sha256)


def _validate_model_bindings(candidate: CandidateEnvelope, model: LengthModel, annotation: LengthAnnotation) -> None:
    if model.annotation_sha256 != annotation.sha256:
        _fail("runtime length model annotation binding differs from the opened annotation")
    if model.boundary_definition != annotation.boundary_definition:
        _fail("runtime length model and annotation target boundary definitions differ")
    if annotation.assembly not in model.applicability.assemblies:
        _fail("runtime length annotation assembly is outside model applicability")
    if (
        model.study_sha256 != candidate.study_sha256
        or model.training_evidence_sha256 != candidate.training_evidence_sha256
        or model.applicability != candidate.applicability
        or model.producer != candidate.producer
    ):
        _fail("runtime length model provenance differs from its candidate envelope")


def load_length_model_bundle(path: Path) -> LengthModelBundle:
    """Load one exact approved length model bundle through pinned descriptors.

    Args:
        path: Directory containing the fixed six-file runtime inventory.

    Returns:
        Immutable approved model, annotation, and transitive identities.

    Raises:
        ValueError: If the directory, bytes, schemas, hashes, or semantic bindings differ.
    """
    try:
        with SecureDirectoryReader.open(path, _BUNDLE_FILES) as reader:
            raw = reader.read_files(tuple(sorted(_BUNDLE_FILES)))
    except ValueError as error:
        raise ValueError(f"runtime length model bundle is invalid: {error}") from error
    bundle_sha256 = _validate_checksums(raw)
    return _decode_bundle(raw, bundle_sha256)
