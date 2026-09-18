"""Explicit portable caller configuration and frozen applicability checks."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from vntyper.modules.advntr.advntr_calibration_policy import (
    CapturePolicy,
    capture_policy_document,
    capture_policy_for_caller,
    decode_capture_policy,
)
from vntyper.scripts.calibration_caller_bundle import CallerModelBundle, load_caller_model_bundle
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)
_FIELDS = {
    "schema_version",
    "evidence_domain",
    "assembly",
    "assay_class",
    "input_scope",
    "preprocessing_id",
    "advntr_capture_policy",
}


@dataclass(frozen=True)
class CallerRuntimeContext:
    """Minimal operator applicability metadata and exact native execution policy."""

    evidence_domain: str
    assembly: str
    assay_class: str
    input_scope: str
    preprocessing_id: str
    advntr_capture_policy: CapturePolicy | None
    sha256: str


@dataclass(frozen=True)
class CallerPipelineConfiguration:
    """One approved portable bundle and its validated run applicability."""

    bundle: CallerModelBundle
    context: CallerRuntimeContext
    operator_paths: tuple[Path, Path]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def decode_caller_runtime_context(value: object) -> CallerRuntimeContext:
    """Decode closed context fields without inferring any build or reference identity.

    Args:
        value: Parsed runtime context document.

    Returns:
        Immutable context with a canonical digest.

    Raises:
        ValueError: If schema, fields, text, or native capture policy is invalid.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("caller runtime context fields differ from the closed contract")
    if value["schema_version"] != "caller-runtime-context-v1":
        _fail("caller runtime context schema differs")
    text = {}
    for field in ("evidence_domain", "assembly", "assay_class", "input_scope", "preprocessing_id"):
        item = value[field]
        if not isinstance(item, str) or not item or item != item.strip():
            _fail("caller runtime context values must be nonempty text without surrounding whitespace")
        text[field] = item
    if text["evidence_domain"] not in {"synthetic", "external"} or text["input_scope"] not in {"full", "regional"}:
        _fail("caller runtime context domain or input scope differs")
    policy = None if value["advntr_capture_policy"] is None else decode_capture_policy(value["advntr_capture_policy"])
    return CallerRuntimeContext(**text, advntr_capture_policy=policy, sha256=canonical_sha256(value))


def caller_runtime_context_document(context: CallerRuntimeContext) -> dict[str, object]:
    """Project a context and reject replaced typed content or digest.

    Args:
        context: Previously decoded context.

    Returns:
        Fresh closed JSON-compatible content.

    Raises:
        ValueError: If the typed content no longer matches its canonical digest.
    """
    if not isinstance(context, CallerRuntimeContext):
        _fail("caller runtime context must be decoded")
    document = {field: getattr(context, field) for field in _FIELDS - {"schema_version", "advntr_capture_policy"}}
    document.update(
        schema_version="caller-runtime-context-v1",
        advntr_capture_policy=None
        if context.advntr_capture_policy is None
        else capture_policy_document(context.advntr_capture_policy),
    )
    if decode_caller_runtime_context(document) != context:
        _fail("caller runtime context differs from canonical content")
    return document


def resolve_caller_pipeline_configuration(bundle_path: Path, context_path: Path) -> CallerPipelineConfiguration:
    """Load an approved bundle and check its existing applicability contract.

    Args:
        bundle_path: Closed portable runtime directory.
        context_path: Explicit regular-file applicability context.

    Returns:
        Resolved configuration with operator paths retained only for I/O protection.

    Raises:
        ValueError: If approval, applicability, or selected native capture binding differs.
    """
    bundle = load_caller_model_bundle(bundle_path)
    context = decode_caller_runtime_context(load_strict_json_object(read_regular_path(context_path)))
    applicability = bundle.candidate.applicability
    if (
        context.evidence_domain != applicability.domain
        or context.assembly not in applicability.assemblies
        or context.assay_class not in applicability.assay_classes
        or context.input_scope not in applicability.input_scopes
        or context.preprocessing_id not in applicability.preprocessing_ids
    ):
        _fail("caller runtime context falls outside frozen candidate applicability")
    native = bundle.advntr_policy
    capture = context.advntr_capture_policy
    if (native is None) != (capture is None):
        _fail("caller runtime native capture policy must be present exactly for dual-caller bundles")
    if native is not None and capture is not None:
        if (
            native.capture_policy_sha256 != capture.sha256
            or capture_policy_for_caller(capture, bundle.caller_policy) != capture
        ):
            _fail("caller runtime capture policy differs from frozen sidecar or selected caller values")
        versions = bundle.candidate.producer.tool_versions
        if "advntr" not in versions or "advntr_build_id" not in versions:
            _fail("caller runtime producer requires explicit advntr and advntr_build_id identities")
    return CallerPipelineConfiguration(bundle, context, (bundle_path, context_path))


def validate_caller_run_options(
    configuration: CallerPipelineConfiguration,
    *,
    assembly: str,
    extra_modules: Sequence[str],
    threads: int,
    additional_commands: str,
) -> None:
    """Check actual resolved command options before native assets or reads are opened.

    Args:
        configuration: Resolved approved caller configuration.
        assembly: Actual pipeline assembly.
        extra_modules: Explicitly selected optional modules.
        threads: Effective native thread count.
        additional_commands: Effective native extension flags.

    Raises:
        ValueError: If command assembly, enablement, threads, or extension flags differ.
    """
    caller_runtime_context_document(configuration.context)
    if assembly != configuration.context.assembly:
        _fail("caller runtime assembly differs from the actual pipeline assembly")
    capture = configuration.context.advntr_capture_policy
    if capture is not None:
        if "advntr" not in extra_modules:
            _fail("dual-caller calibration requires explicit --extra-modules advntr")
        if type(threads) is not int or threads != capture.threads:
            _fail("caller runtime effective threads differ from the frozen native capture policy")
        if additional_commands != "":
            _fail("caller runtime calibration requires empty additional_commands")
