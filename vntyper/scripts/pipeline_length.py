"""Strict preflight contracts for optional pipeline VNTR length measurement."""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal, NoReturn, cast

from vntyper.scripts.calibration_candidate import candidate_document, decode_candidate, validate_candidate_payload
from vntyper.scripts.calibration_payload import decode_payload_manifest, payload_manifest_document
from vntyper.scripts.calibration_portable_projection import (
    decode_portable_approval,
    portable_approval_document,
    validate_portable_approval_candidate,
)
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation, encode_length_annotation
from vntyper.scripts.length_feature_provenance import (
    LengthFeatureContext,
    decode_length_feature_context,
    encode_length_feature_context,
)
from vntyper.scripts.length_model import LengthModel, encode_length_model
from vntyper.scripts.length_model_bundle import LengthModelBundle, load_length_model_bundle

EvidenceDomain = Literal["synthetic", "external"]
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_CONFIG_FIELDS = {
    "schema_version",
    "measurement_enabled",
    "annotation_sha256",
    "model_sha256",
    "model_bundle_sha256",
    "portable_approval_sha256",
    "candidate_id",
    "context_sha256",
    "counting_policy_sha256",
    "evidence_domain",
}


@dataclass(frozen=True)
class LengthPipelineContext:
    """Operator-supplied provenance whose observed claims remain independently checked."""

    evidence_domain: EvidenceDomain
    measurement_context: LengthFeatureContext
    sha256: str


@dataclass(frozen=True)
class LengthPipelineConfiguration:
    """Resolved length inputs and path-free identity for one pipeline run."""

    measurement_enabled: bool
    annotation: LengthAnnotation | None
    model: LengthModel | None
    measurement_context: LengthFeatureContext | None
    evidence_domain: EvidenceDomain | None
    model_bundle_sha256: str | None
    portable_approval_sha256: str | None
    candidate_id: str | None
    context_sha256: str | None
    sha256: str


def _fail(message: str) -> NoReturn:
    raise ValueError(message)


def decode_length_pipeline_context(value: object) -> LengthPipelineContext:
    """Decode the exact canonicalizable pipeline context sidecar content.

    Args:
        value: Parsed JSON-compatible sidecar object.

    Returns:
        Immutable context with a canonical content digest.

    Raises:
        ValueError: If fields, evidence domain, or measurement context are invalid.
    """
    if not isinstance(value, dict) or set(value) != {"schema_version", "evidence_domain", "measurement_context"}:
        _fail("length pipeline context fields differ from the closed contract")
    if value["schema_version"] != "length-pipeline-context-v1":
        _fail("length pipeline context schema_version must be length-pipeline-context-v1")
    domain = value["evidence_domain"]
    if not isinstance(domain, str) or domain not in {"synthetic", "external"}:
        _fail("length pipeline evidence_domain must be synthetic or external")
    measurement_context = decode_length_feature_context(value["measurement_context"])
    document = {
        "schema_version": "length-pipeline-context-v1",
        "evidence_domain": domain,
        "measurement_context": encode_length_feature_context(measurement_context),
    }
    return LengthPipelineContext(cast(EvidenceDomain, domain), measurement_context, canonical_sha256(document))


def encode_length_pipeline_context(context: LengthPipelineContext) -> dict[str, object]:
    """Project and revalidate a pipeline context.

    Args:
        context: Context returned by :func:`decode_length_pipeline_context`.

    Returns:
        Fresh closed JSON-compatible content.

    Raises:
        ValueError: If typed content or its digest was replaced or forged.
    """
    if not isinstance(context, LengthPipelineContext):
        _fail("length pipeline context must be a LengthPipelineContext")
    document: dict[str, object] = {
        "schema_version": "length-pipeline-context-v1",
        "evidence_domain": context.evidence_domain,
        "measurement_context": encode_length_feature_context(context.measurement_context),
    }
    if decode_length_pipeline_context(document) != context:
        _fail("length pipeline context digest does not match canonical content")
    return document


def _annotation_positions(annotation: LengthAnnotation) -> set[int]:
    intervals = [*(annotation.core or ()), *(annotation.invariant or ())]
    intervals.extend(
        value for value in (annotation.array, annotation.left_flank, annotation.right_flank) if value is not None
    )
    return {position for interval in intervals for position in range(interval.start, interval.end)}


def _context_positions(context: LengthFeatureContext) -> set[int]:
    return {
        position
        for interval in context.counting_policy.queried_intervals
        for position in range(interval.start, interval.end)
    }


def _validate_annotation_context(annotation: LengthAnnotation, context: LengthFeatureContext) -> None:
    decoded_annotation = decode_length_annotation(encode_length_annotation(annotation))
    if decoded_annotation != annotation:
        _fail("length pipeline annotation digest does not match canonical content")
    encode_length_feature_context(context)
    if context.annotation_sha256 != annotation.sha256:
        _fail("length pipeline context annotation digest differs")
    if context.assembly != annotation.assembly:
        _fail("length pipeline context assembly differs")
    if context.reference_fasta_sha256 != annotation.reference_fasta_sha256:
        _fail("length pipeline context reference digest differs")
    if context.original_contig not in annotation.accepted_contigs:
        _fail("length pipeline context contig is outside annotation aliases")
    if _context_positions(context) != _annotation_positions(annotation):
        _fail("length pipeline context queried interval union differs from annotation")


def _validated_bundle(bundle: LengthModelBundle) -> LengthModelBundle:
    if not isinstance(bundle, LengthModelBundle):
        _fail("length pipeline model must come from an approved runtime bundle")
    candidate = decode_candidate(candidate_document(bundle.candidate))
    payload = decode_payload_manifest(payload_manifest_document(bundle.payload))
    approval = decode_portable_approval(portable_approval_document(bundle.approval))
    model_document = encode_length_model(bundle.model)
    model = bundle.model
    annotation = decode_length_annotation(encode_length_annotation(bundle.annotation))
    if canonical_sha256(model_document) != model.sha256 or annotation != bundle.annotation:
        _fail("length pipeline bundle model or annotation digest differs")
    validate_candidate_payload(candidate, payload, expected_target="length")
    validate_portable_approval_candidate(approval, candidate)
    if not isinstance(bundle.sha256, str) or _SHA256.fullmatch(bundle.sha256) is None:
        _fail("length pipeline model bundle digest must be lowercase SHA-256")
    if model.annotation_sha256 != annotation.sha256 or model.boundary_definition != annotation.boundary_definition:
        _fail("length pipeline bundle model and annotation bindings differ")
    return bundle


def build_length_pipeline_configuration(
    *,
    annotation: LengthAnnotation,
    pipeline_context: LengthPipelineContext,
    bundle: LengthModelBundle | None,
) -> LengthPipelineConfiguration:
    """Build enabled path-free configuration from already decoded inputs.

    Args:
        annotation: Explicit or bundle-supplied annotation.
        pipeline_context: Explicit pipeline provenance sidecar.
        bundle: Approved runtime model bundle, or ``None`` for measurement only.

    Returns:
        Immutable enabled configuration.

    Raises:
        ValueError: If typed inputs, geometry, provenance, or model applicability differ.
    """
    if not isinstance(annotation, LengthAnnotation) or not isinstance(pipeline_context, LengthPipelineContext):
        _fail("enabled length pipeline configuration requires typed annotation and context")
    encode_length_pipeline_context(pipeline_context)
    context = pipeline_context.measurement_context
    _validate_annotation_context(annotation, context)
    model = None
    bundle_sha256 = approval_sha256 = candidate_id = None
    if bundle is not None:
        bundle = _validated_bundle(bundle)
        if bundle.annotation != annotation:
            _fail("length pipeline annotation differs from runtime bundle annotation")
        model = bundle.model
        applicability = model.applicability
        observed = (
            pipeline_context.evidence_domain,
            context.assembly,
            context.assay_class,
            context.input_scope,
            context.preprocessing_id,
            context.aligner.name,
            context.aligner.version,
            context.aligner.arguments_sha256,
            context.aligner.primary_secondary_marking,
            context.counting_policy_sha256,
        )
        if (
            observed[0] != applicability.domain
            or observed[1] not in applicability.assemblies
            or observed[2] not in applicability.assay_classes
            or observed[3] not in applicability.input_scopes
            or observed[4] not in applicability.preprocessing_ids
            or observed[5:]
            != (
                applicability.aligner_name,
                applicability.aligner_version,
                applicability.aligner_arguments_sha256,
                applicability.primary_secondary_marking,
                applicability.counting_policy_sha256,
            )
        ):
            _fail("length pipeline context is outside model applicability")
        bundle_sha256 = bundle.sha256
        approval_sha256 = bundle.approval.sha256
        candidate_id = bundle.candidate.candidate_id
    configuration = LengthPipelineConfiguration(
        measurement_enabled=True,
        annotation=annotation,
        model=model,
        measurement_context=context,
        evidence_domain=pipeline_context.evidence_domain,
        model_bundle_sha256=bundle_sha256,
        portable_approval_sha256=approval_sha256,
        candidate_id=candidate_id,
        context_sha256=pipeline_context.sha256,
        sha256="",
    )
    return replace(configuration, sha256=canonical_sha256(_configuration_document(configuration)))


def _configuration_document(configuration: LengthPipelineConfiguration) -> dict[str, object]:
    return {
        "schema_version": "length-pipeline-configuration-v1",
        "measurement_enabled": configuration.measurement_enabled,
        "annotation_sha256": configuration.annotation.sha256 if configuration.annotation is not None else None,
        "model_sha256": configuration.model.sha256 if configuration.model is not None else None,
        "model_bundle_sha256": configuration.model_bundle_sha256,
        "portable_approval_sha256": configuration.portable_approval_sha256,
        "candidate_id": configuration.candidate_id,
        "context_sha256": configuration.context_sha256,
        "counting_policy_sha256": (
            configuration.measurement_context.counting_policy_sha256
            if configuration.measurement_context is not None
            else None
        ),
        "evidence_domain": configuration.evidence_domain,
    }


def encode_length_pipeline_configuration(configuration: LengthPipelineConfiguration) -> dict[str, object]:
    """Project and revalidate the path-free configuration identity.

    Args:
        configuration: Resolved length pipeline configuration.

    Returns:
        Fresh closed JSON-compatible identity document.

    Raises:
        ValueError: If the object is inconsistent or its digest is stale.
    """
    if not isinstance(configuration, LengthPipelineConfiguration):
        _fail("length pipeline configuration must be typed")
    document = _configuration_document(configuration)
    if set(document) != _CONFIG_FIELDS:
        _fail("length pipeline configuration fields differ")
    if configuration.measurement_enabled:
        if configuration.annotation is None or configuration.measurement_context is None:
            _fail("enabled length pipeline configuration requires annotation and context")
        _validate_annotation_context(configuration.annotation, configuration.measurement_context)
        if configuration.model is None:
            if any(
                value is not None
                for value in (
                    configuration.model_bundle_sha256,
                    configuration.portable_approval_sha256,
                    configuration.candidate_id,
                )
            ):
                _fail("measurement-only length configuration must not contain model approval fields")
        else:
            encode_length_model(configuration.model)
            if any(
                value is None
                for value in (
                    configuration.model_bundle_sha256,
                    configuration.portable_approval_sha256,
                    configuration.candidate_id,
                )
            ):
                _fail("model length configuration requires approval fields")
    elif any(
        value is not None for key, value in document.items() if key not in {"schema_version", "measurement_enabled"}
    ):
        _fail("disabled length pipeline configuration must have null identities")
    if canonical_sha256(document) != configuration.sha256:
        _fail("length pipeline configuration digest does not match canonical content")
    return document


def _read_canonical_object(path: Path, label: str) -> dict[str, object]:
    try:
        raw = read_regular_path(path)
        value = load_strict_json_object(raw)
        if canonical_json_bytes(value) != raw:
            _fail(f"length pipeline {label} must use canonical JSON bytes")
        return value
    except (UnicodeDecodeError, ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"length pipeline {label} is invalid: {error}") from error


def _disabled_configuration() -> LengthPipelineConfiguration:
    configuration = LengthPipelineConfiguration(False, None, None, None, None, None, None, None, None, "")
    return replace(configuration, sha256=canonical_sha256(_configuration_document(configuration)))


def resolve_length_pipeline_configuration(
    *,
    measurement_enabled: bool,
    model_path: Path | None,
    annotation_path: Path | None,
    context_path: Path | None,
) -> LengthPipelineConfiguration:
    """Resolve optional paths before alignment I/O and freeze their scientific identity.

    Args:
        measurement_enabled: Explicit measurement request; a model also enables measurement.
        model_path: Approved runtime bundle directory.
        annotation_path: Explicit annotation for measurement-only mode.
        context_path: Explicit complete measurement and evidence provenance sidecar.

    Returns:
        Disabled or enabled immutable configuration.

    Raises:
        ValueError: If flags, paths, canonical bytes, or scientific bindings differ.
    """
    if not isinstance(measurement_enabled, bool):
        _fail("length measurement_enabled must be boolean")
    for value, label in ((model_path, "model"), (annotation_path, "annotation"), (context_path, "context")):
        if value is not None and not isinstance(value, Path):
            _fail(f"length pipeline {label} path must be a Path")
    if model_path is None and not measurement_enabled:
        if annotation_path is not None or context_path is not None:
            _fail("length annotation or context requires length measurement")
        return _disabled_configuration()
    if context_path is None:
        _fail("enabled length measurement requires --length-context")
    if model_path is not None and annotation_path is not None:
        _fail("--length-model supplies its annotation; --length-annotation must be omitted")
    pipeline_context = decode_length_pipeline_context(_read_canonical_object(context_path, "context"))
    if model_path is not None:
        bundle = load_length_model_bundle(model_path)
        return build_length_pipeline_configuration(
            annotation=bundle.annotation,
            pipeline_context=pipeline_context,
            bundle=bundle,
        )
    if annotation_path is None:
        _fail("measurement-only length mode requires --length-annotation")
    annotation = decode_length_annotation(_read_canonical_object(annotation_path, "annotation"))
    return build_length_pipeline_configuration(annotation=annotation, pipeline_context=pipeline_context, bundle=None)
