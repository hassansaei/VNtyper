"""Immutable research candidates for caller and total-length calibration targets."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_payload import PayloadManifest
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

CalibrationTarget = Literal["callers", "length"]
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_REVISION = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_DIGEST_FIELDS = (
    "study_sha256",
    "baseline_sha256",
    "partition_sha256",
    "training_evidence_sha256",
    "selection_evidence_sha256",
    "payload_sha256",
)
_FIELDS = {
    "schema_version",
    "target",
    "candidate_id",
    *_DIGEST_FIELDS,
    "applicability",
    "producer",
    "status",
}


@dataclass(frozen=True)
class CandidateApplicability:
    """Predeclared measurement domain; caller membership is target-specific."""

    domain: str
    assemblies: tuple[str, ...]
    assay_classes: tuple[str, ...]
    input_scopes: tuple[str, ...]
    preprocessing_ids: tuple[str, ...]
    required_callers: tuple[str, ...]


@dataclass(frozen=True)
class CandidateProducer:
    """Build and feature provenance without input locations or specimen labels."""

    name: str
    version: str
    source_revision: str | None
    tool_versions: Mapping[str, str]
    feature_schema_sha256: str


@dataclass(frozen=True)
class CandidateEnvelope:
    """Hash-bound research artifact; construction does not confer promotion approval."""

    target: CalibrationTarget
    candidate_id: str
    study_sha256: str
    baseline_sha256: str
    partition_sha256: str
    training_evidence_sha256: str
    selection_evidence_sha256: str
    payload_sha256: str
    applicability: CandidateApplicability
    producer: CandidateProducer
    status: Literal["research-only"]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"{label} fields differ from the closed contract")
    return value


def _text(value: object, field: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        _fail(f"{field} must be non-empty text without surrounding whitespace")
    return value


def _digest(value: object, field: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"{field} must be a lowercase SHA-256 digest")
    return value


def _strings(value: object, field: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not value:
        _fail(f"{field} must be a non-empty sorted unique list")
    items = tuple(_text(item, field) for item in value)
    if items != tuple(sorted(set(items))):
        _fail(f"{field} must be a non-empty sorted unique list")
    return items


def _applicability(value: object, target: CalibrationTarget) -> CandidateApplicability:
    fields = {"domain", "assemblies", "assay_classes", "input_scopes", "preprocessing_ids"}
    if target == "callers":
        fields.add("required_callers")
    raw = _object(value, fields, "candidate applicability")
    domain = _text(raw["domain"], "domain")
    if domain not in {"synthetic", "external"}:
        _fail("domain must be synthetic or external")
    scopes = _strings(raw["input_scopes"], "input_scopes")
    if not set(scopes) <= {"full", "regional"}:
        _fail("input_scopes must contain only full or regional")
    callers = _strings(raw["required_callers"], "required_callers") if target == "callers" else ()
    if target == "callers" and ("kestrel" not in callers or not set(callers) <= {"kestrel", "advntr"}):
        _fail("required_callers must include kestrel and optionally advntr")
    return CandidateApplicability(
        domain,
        _strings(raw["assemblies"], "assemblies"),
        _strings(raw["assay_classes"], "assay_classes"),
        scopes,
        _strings(raw["preprocessing_ids"], "preprocessing_ids"),
        callers,
    )


def _producer(value: object) -> CandidateProducer:
    raw = _object(
        value,
        {"name", "version", "source_revision", "tool_versions", "feature_schema_sha256"},
        "candidate producer",
    )
    revision = raw["source_revision"]
    if revision is not None and (not isinstance(revision, str) or _REVISION.fullmatch(revision) is None):
        _fail("source_revision must be a lowercase full Git revision or null")
    tools = raw["tool_versions"]
    if not isinstance(tools, Mapping) or not tools:
        _fail("tool_versions must be a non-empty string mapping")
    versions = {_text(key, "tool_versions"): _text(value, "tool_versions") for key, value in tools.items()}
    return CandidateProducer(
        _text(raw["name"], "name"),
        _text(raw["version"], "version"),
        cast(str | None, revision),
        MappingProxyType(versions),
        _digest(raw["feature_schema_sha256"], "feature_schema_sha256"),
    )


def decode_candidate(value: object) -> CandidateEnvelope:
    """Decode a research candidate and verify its deterministic content identity.

    Args:
        value: Strictly decoded JSON candidate object.

    Returns:
        Immutable target and payload binding. It is not a validation attestation.

    Raises:
        ValueError: If the schema, domain, provenance or identity is invalid.
    """
    raw = _object(value, _FIELDS, "candidate")
    if raw["schema_version"] != "calibration-candidate-v2":
        _fail("candidate schema_version must be calibration-candidate-v2")
    if raw["target"] not in ("callers", "length"):
        _fail("candidate target must be callers or length")
    target = cast(CalibrationTarget, raw["target"])
    if raw["status"] != "research-only":
        _fail("candidate status must remain research-only; approval requires attestations")
    digests = {field: _digest(raw[field], field) for field in _DIGEST_FIELDS}
    applicability = _applicability(raw["applicability"], target)
    producer = _producer(raw["producer"])
    candidate_id = _digest(raw["candidate_id"], "candidate_id")
    if candidate_id != canonical_sha256({key: child for key, child in raw.items() if key != "candidate_id"}):
        _fail("candidate_id differs from the candidate's canonical bound content")
    return CandidateEnvelope(
        target=target,
        candidate_id=candidate_id,
        study_sha256=digests["study_sha256"],
        baseline_sha256=digests["baseline_sha256"],
        partition_sha256=digests["partition_sha256"],
        training_evidence_sha256=digests["training_evidence_sha256"],
        selection_evidence_sha256=digests["selection_evidence_sha256"],
        payload_sha256=digests["payload_sha256"],
        applicability=applicability,
        producer=producer,
        status="research-only",
        sha256=canonical_sha256(raw),
    )


def candidate_document(candidate: CandidateEnvelope) -> dict[str, object]:
    """Return a fresh canonicalizable document for a decoded candidate.

    Args:
        candidate: Immutable candidate binding.

    Returns:
        Independent JSON-compatible object with the original content identity.

    Raises:
        ValueError: If candidate is not a CandidateEnvelope.
    """
    if not isinstance(candidate, CandidateEnvelope):
        _fail("candidate must be a CandidateEnvelope")
    applicable = candidate.applicability
    applicability: dict[str, object] = {
        "domain": applicable.domain,
        "assemblies": list(applicable.assemblies),
        "assay_classes": list(applicable.assay_classes),
        "input_scopes": list(applicable.input_scopes),
        "preprocessing_ids": list(applicable.preprocessing_ids),
    }
    if candidate.target == "callers":
        applicability["required_callers"] = list(applicable.required_callers)
    return {
        "schema_version": "calibration-candidate-v2",
        "target": candidate.target,
        "candidate_id": candidate.candidate_id,
        **{field: getattr(candidate, field) for field in _DIGEST_FIELDS},
        "applicability": applicability,
        "producer": {
            "name": candidate.producer.name,
            "version": candidate.producer.version,
            "source_revision": candidate.producer.source_revision,
            "tool_versions": dict(candidate.producer.tool_versions),
            "feature_schema_sha256": candidate.producer.feature_schema_sha256,
        },
        "status": candidate.status,
    }


def validate_candidate_payload(
    candidate: CandidateEnvelope, payload: PayloadManifest, *, expected_target: CalibrationTarget
) -> None:
    """Check target and complete payload identity before reading model semantics.

    Args:
        candidate: Decoded immutable candidate.
        payload: Decoded full payload manifest.
        expected_target: Target selected by the consuming operation.

    Raises:
        ValueError: If types, target or payload binding differ.
    """
    if not isinstance(candidate, CandidateEnvelope):
        _fail("candidate must be a CandidateEnvelope")
    if not isinstance(payload, PayloadManifest):
        _fail("payload must be a PayloadManifest")
    if candidate.target != expected_target:
        _fail("candidate target differs from the consuming operation")
    if candidate.payload_sha256 != payload.sha256:
        _fail("candidate payload digest differs from the opened manifest")
