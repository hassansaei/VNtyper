"""Immutable research candidates for caller and total-length calibration targets."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_payload import (
    CALLER_BUNDLE_DESCRIPTOR_PATH,
    CallerBundleDescriptor,
    PayloadManifest,
    caller_bundle_descriptor_document,
    payload_manifest_document,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256

logger = logging.getLogger(__name__)

CalibrationTarget = Literal["callers", "length"]
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_REVISION = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_MAPPING_PROXY_TYPE: type[object] = type(MappingProxyType({}))
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
    aligner_name: str | None
    aligner_version: str | None
    aligner_arguments_sha256: str | None
    primary_secondary_marking: str | None
    counting_policy_sha256: str | None


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
    else:
        fields.update(
            {
                "aligner_name",
                "aligner_version",
                "aligner_arguments_sha256",
                "primary_secondary_marking",
                "counting_policy_sha256",
            }
        )
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
    aligner_name = _text(raw["aligner_name"], "aligner_name") if target == "length" else None
    aligner_version = _text(raw["aligner_version"], "aligner_version") if target == "length" else None
    aligner_arguments_sha256 = (
        _digest(raw["aligner_arguments_sha256"], "aligner_arguments_sha256") if target == "length" else None
    )
    marking = _text(raw["primary_secondary_marking"], "primary_secondary_marking") if target == "length" else None
    counting_policy_sha256 = (
        _digest(raw["counting_policy_sha256"], "counting_policy_sha256") if target == "length" else None
    )
    return CandidateApplicability(
        domain,
        _strings(raw["assemblies"], "assemblies"),
        _strings(raw["assay_classes"], "assay_classes"),
        scopes,
        _strings(raw["preprocessing_ids"], "preprocessing_ids"),
        callers,
        aligner_name,
        aligner_version,
        aligner_arguments_sha256,
        marking,
        counting_policy_sha256,
    )


def decode_candidate_applicability(value: object, *, target: CalibrationTarget) -> CandidateApplicability:
    """Decode the shared closed applicability contract for one calibration target.

    Args:
        value: Parsed JSON applicability object.
        target: Caller or length target selecting its exact fields.

    Returns:
        Immutable validated applicability.

    Raises:
        ValueError: If the target or applicability content is invalid.
    """
    if target not in ("callers", "length"):
        _fail("candidate applicability target must be callers or length")
    return _applicability(value, target)


def _applicability_document(applicable: CandidateApplicability, target: CalibrationTarget) -> dict[str, object]:
    raw: dict[str, object] = {
        "domain": applicable.domain,
        "assemblies": list(applicable.assemblies),
        "assay_classes": list(applicable.assay_classes),
        "input_scopes": list(applicable.input_scopes),
        "preprocessing_ids": list(applicable.preprocessing_ids),
    }
    if target == "callers":
        raw["required_callers"] = list(applicable.required_callers)
    else:
        raw.update(
            {
                "aligner_name": applicable.aligner_name,
                "aligner_version": applicable.aligner_version,
                "aligner_arguments_sha256": applicable.aligner_arguments_sha256,
                "primary_secondary_marking": applicable.primary_secondary_marking,
                "counting_policy_sha256": applicable.counting_policy_sha256,
            }
        )
    return raw


def candidate_applicability_document(
    applicable: CandidateApplicability, *, target: CalibrationTarget
) -> dict[str, object]:
    """Project and revalidate shared immutable applicability.

    Args:
        applicable: Typed applicability returned by its decoder.
        target: Caller or length target selecting its exact fields.

    Returns:
        Fresh JSON-compatible applicability object.

    Raises:
        ValueError: If typed content, target, or immutable collections are invalid.
    """
    if not isinstance(applicable, CandidateApplicability):
        _fail("candidate applicability must be a CandidateApplicability")
    if not all(
        isinstance(items, tuple)
        for items in (
            applicable.assemblies,
            applicable.assay_classes,
            applicable.input_scopes,
            applicable.preprocessing_ids,
            applicable.required_callers,
        )
    ):
        _fail("candidate applicability must use decoded immutable collections")
    raw = _applicability_document(applicable, target)
    if decode_candidate_applicability(raw, target=target) != applicable:
        _fail("candidate applicability differs from its decoded contract")
    return raw


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


def decode_candidate_producer(value: object) -> CandidateProducer:
    """Decode the shared closed producer provenance contract.

    Args:
        value: Parsed JSON producer object.

    Returns:
        Immutable validated producer provenance.

    Raises:
        ValueError: If producer fields or values are invalid.
    """
    return _producer(value)


def _producer_document(producer: CandidateProducer) -> dict[str, object]:
    return {
        "name": producer.name,
        "version": producer.version,
        "source_revision": producer.source_revision,
        "tool_versions": dict(producer.tool_versions),
        "feature_schema_sha256": producer.feature_schema_sha256,
    }


def candidate_producer_document(producer: CandidateProducer) -> dict[str, object]:
    """Project and revalidate shared immutable producer provenance.

    Args:
        producer: Typed producer returned by its decoder.

    Returns:
        Fresh JSON-compatible producer object.

    Raises:
        ValueError: If typed content or immutable mappings are invalid.
    """
    if not isinstance(producer, CandidateProducer):
        _fail("candidate producer must be a CandidateProducer")
    if not isinstance(producer.tool_versions, _MAPPING_PROXY_TYPE):
        _fail("candidate producer must use a decoded immutable tool_versions mapping")
    raw = _producer_document(producer)
    if decode_candidate_producer(raw) != producer:
        _fail("candidate producer differs from its decoded contract")
    return raw


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


def _candidate_document(candidate: CandidateEnvelope) -> dict[str, object]:
    applicability = _applicability_document(candidate.applicability, candidate.target)
    return {
        "schema_version": "calibration-candidate-v2",
        "target": candidate.target,
        "candidate_id": candidate.candidate_id,
        **{field: getattr(candidate, field) for field in _DIGEST_FIELDS},
        "applicability": applicability,
        "producer": _producer_document(candidate.producer),
        "status": candidate.status,
    }


def _require_candidate(candidate: object) -> CandidateEnvelope:
    if not isinstance(candidate, CandidateEnvelope):
        _fail("candidate must be a CandidateEnvelope")
    applicable = candidate.applicability
    if (
        not isinstance(applicable, CandidateApplicability)
        or not all(
            isinstance(items, tuple)
            for items in (
                applicable.assemblies,
                applicable.assay_classes,
                applicable.input_scopes,
                applicable.preprocessing_ids,
                applicable.required_callers,
            )
        )
        or not isinstance(candidate.producer, CandidateProducer)
        or not isinstance(candidate.producer.tool_versions, _MAPPING_PROXY_TYPE)
    ):
        _fail("candidate must use decoded immutable content")
    decoded = decode_candidate(_candidate_document(candidate))
    if decoded != candidate:
        _fail("candidate differs from its canonical content or digest")
    return candidate


def candidate_document(candidate: CandidateEnvelope) -> dict[str, object]:
    """Return a fresh canonicalizable document for a decoded candidate.

    Args:
        candidate: Immutable candidate binding.

    Returns:
        Independent JSON-compatible object with the original content identity.

    Raises:
        ValueError: If candidate is not a CandidateEnvelope.
    """
    return _candidate_document(_require_candidate(candidate))


def validate_candidate_payload(
    candidate: CandidateEnvelope,
    payload: PayloadManifest,
    *,
    expected_target: CalibrationTarget,
    caller_descriptor: CallerBundleDescriptor | None = None,
) -> None:
    """Check target and complete payload identity before reading model semantics.

    Args:
        candidate: Decoded immutable candidate.
        payload: Decoded full payload manifest.
        expected_target: Target selected by the consuming operation.
        caller_descriptor: Required closed composition descriptor for caller payloads.

    Raises:
        ValueError: If types, target or payload binding differ.
    """
    candidate = _require_candidate(candidate)
    payload_manifest_document(payload)
    if candidate.target != expected_target:
        _fail("candidate target differs from the consuming operation")
    if candidate.payload_sha256 != payload.sha256:
        _fail("candidate payload digest differs from the opened manifest")
    if candidate.target == "length":
        if caller_descriptor is not None:
            _fail("length candidate cannot accept a caller descriptor")
        return
    if caller_descriptor is None:
        _fail("caller candidate requires its caller descriptor")
    descriptor_document = caller_bundle_descriptor_document(caller_descriptor)
    if caller_descriptor.required_callers != candidate.applicability.required_callers:
        _fail("caller descriptor required_callers differ from candidate applicability")
    expected_hashes = {
        CALLER_BUNDLE_DESCRIPTOR_PATH: caller_descriptor.sha256,
        "decision-profile.json": caller_descriptor.decision_profile_sha256,
    }
    if caller_descriptor.advntr_policy_sha256 is not None:
        expected_hashes["advntr-policy.json"] = caller_descriptor.advntr_policy_sha256
    if caller_descriptor.background_sha256 is not None:
        expected_hashes["background.json"] = caller_descriptor.background_sha256
    observed_hashes = {item.path: item.sha256 for item in payload.files}
    if set(observed_hashes) != set(expected_hashes):
        _fail("caller payload file set differs from its descriptor")
    if observed_hashes != expected_hashes:
        _fail("caller payload component digest differs from its descriptor")
    descriptor_file = next(item for item in payload.files if item.path == CALLER_BUNDLE_DESCRIPTOR_PATH)
    if descriptor_file.size_bytes != len(canonical_json_bytes(descriptor_document)):
        _fail("caller descriptor size differs from its canonical bytes")
