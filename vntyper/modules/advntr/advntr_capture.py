"""Shell-free adapter for complete adVNTR v2 calibration captures."""

from __future__ import annotations

import hashlib
import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from vntyper.modules.advntr.advntr_calibration_policy import (
    AdvntrCapabilities,
    AdvntrToolPin,
    CapturePolicy,
    ProcessRunner,
    advntr_canonical_sha256,
    advntr_capabilities_document,
    capture_policy_document,
    capture_policy_for_caller,
    probe_advntr_capabilities,
)
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, caller_policy_values_document
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)

_CAPTURE_FIELDS = {
    "schema_version",
    "completion",
    "producer",
    "assets",
    "model_locus",
    "loaded_background",
    "capture_policy",
    "caller_policy",
    "locus",
    "unit_geometry",
    "reference_order",
    "flank_boundaries",
    "warnings",
    "occurrences",
    "spans",
    "evidence_rows",
    "candidate_traversal",
    "decision_visits",
}
_PRODUCER_FIELDS = {"package_version", "build_id", "source_revision"}
_ASSET_FIELDS = {"model_sha256", "model_locus_sha256", "background_sha256", "loaded_background_sha256"}
_LOCUS_FIELDS = {"vntr_id", "read_length", "is_haploid", "selected_read_count"}


@dataclass(frozen=True)
class CapturePlan:
    """One exact production invocation and its expected evidence bindings."""

    argv_prefix: tuple[str, ...]
    argv: tuple[str, ...]
    alignment_path: Path
    model_path: Path
    result_path: Path
    sink_path: Path
    working_directory: Path
    reference_path: Path | None
    vntr_ids: tuple[int, ...]
    capture_policy: CapturePolicy
    caller_policy: CallerPolicyValues
    model_sha256: str
    background_path: Path | None
    background_sha256: str | None


@dataclass(frozen=True)
class CaptureRecord:
    """One complete opaque upstream record with only transport bindings projected."""

    vntr_id: int
    raw_json: bytes
    record_sha256: str


@dataclass(frozen=True)
class CaptureResult:
    """Complete process-success capture roster and unchanged raw record bytes."""

    sink_sha256: str
    result_sha256: str
    capabilities: AdvntrCapabilities
    records: tuple[CaptureRecord, ...]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _absolute(path: object, label: str) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        _fail(f"{label} must be an absolute Path")
    return path


def _digest(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        _fail(f"{label} must be a lowercase SHA256 digest")
    return value


def _caller_document(caller: CallerPolicyValues) -> dict[str, object]:
    document = caller_policy_values_document(caller)
    values = document["values"]
    if not isinstance(values, dict) or "advntr" not in caller.required_callers:
        _fail("adVNTR capture requires a complete adVNTR caller policy")
    prefix = "/components/advntr/calibrated_calling/"
    return {
        "schema_version": "advntr-frameshift-policy-v1",
        "mode": values[f"{prefix}mode"],
        "cutoff": values[f"{prefix}cutoff"],
        "minimum_read_support": values[f"{prefix}minimum_read_support"],
    }


def caller_policy_document(caller: CallerPolicyValues) -> dict[str, object]:
    """Project the three upstream native decision values from a full caller policy."""
    return _caller_document(caller)


def _number(value: int | float) -> str:
    return str(value)


def calibrated_policy_argv(policy: CapturePolicy, caller: CallerPolicyValues, background: Path | None) -> list[str]:
    """Render the exact representable native policy for capture or runtime calls.

    Args:
        policy: Full resolved native capture policy.
        caller: Complete selected caller values.
        background: Exact-mode background path, otherwise None.

    Returns:
        Explicit genotype policy arguments including the thread count.

    Raises:
        ValueError: If policy values conflict or cannot be represented by the CLI.
    """
    expected = capture_policy_for_caller(policy, caller)
    if expected != policy:
        _fail("adVNTR capture policy differs from the full calibrated caller values")
    if (
        policy.platform != "illumina"
        or not policy.frameshift_mode
        or policy.is_haploid
        or policy.maximum_error_rate != 0.05
        or policy.legacy_error_rate != 0.01
        or policy.mapq_cutoff != 0
        or policy.base_quality_cutoff != 20
        or policy.maximum_low_quality_fraction != 0.1
        or not policy.enhanced_hmm
        or policy.trained_hmms
    ):
        _fail("adVNTR capture policy contains values the published genotype CLI cannot represent exactly")
    values = caller.values
    prefix = "/components/advntr/calibrated_calling/"
    argv = [
        "-t",
        str(policy.threads),
        "--frameshift-pvalue-cutoff",
        _number(values[f"{prefix}cutoff"]),  # type: ignore[arg-type]
        "--min-frameshift-read-support",
        str(values[f"{prefix}minimum_read_support"]),
    ]
    if policy.minimum_read_length is not None:
        argv.extend(("--min_read_length", str(policy.minimum_read_length)))
    if policy.prune_reverse:
        argv.append("--prune-reverse")
    if policy.filter_adapter_readthrough:
        argv.append("--filter-adapter-readthrough")
    if policy.minimum_read_match_ratio is not None:
        argv.extend(("--min-read-match-ratio", _number(policy.minimum_read_match_ratio)))
    if policy.minimum_relative_ru_coverage is not None:
        argv.extend(("--rare-unit-coverage-guard", _number(policy.minimum_relative_ru_coverage)))
    if not policy.use_reference_alignment:
        argv.append("--noref_aln")
    if policy.fully_covered_ru_only:
        argv.append("--fullru")
    if policy.caller_mode == "exact":
        if background is None:
            _fail("adVNTR exact capture requires a background")
        argv.extend(("--exact-frameshift-caller", "--frameshift-background", str(background)))
    elif background is not None:
        _fail("adVNTR legacy capture forbids a background")
    return argv


def build_capture_plan(
    argv_prefix: tuple[str, ...],
    *,
    alignment_path: Path,
    model_path: Path,
    result_path: Path,
    sink_path: Path,
    working_directory: Path,
    vntr_ids: tuple[int, ...],
    capture_policy: CapturePolicy,
    caller_policy: CallerPolicyValues,
    model_sha256: str,
    background_path: Path | None = None,
    background_sha256: str | None = None,
    reference_path: Path | None = None,
    additional_commands: str = "",
) -> CapturePlan:
    """Build one typed shell-free production capture invocation."""
    if (
        not isinstance(argv_prefix, tuple)
        or not argv_prefix
        or any(not isinstance(token, str) or not token or "\x00" in token for token in argv_prefix)
    ):
        _fail("adVNTR executable prefix must be a nonempty tuple of command tokens")
    paths = tuple(
        _absolute(path, label)
        for path, label in (
            (alignment_path, "adVNTR capture alignment"),
            (model_path, "adVNTR capture model"),
            (result_path, "adVNTR production result"),
            (sink_path, "adVNTR capture sink"),
            (working_directory, "adVNTR working directory"),
        )
    )
    alignment_path, model_path, result_path, sink_path, working_directory = paths
    if reference_path is not None:
        reference_path = _absolute(reference_path, "adVNTR capture reference")
    if background_path is not None:
        background_path = _absolute(background_path, "adVNTR capture background")
    if not isinstance(additional_commands, str) or additional_commands:
        _fail("calibration capture forbids advntr additional_commands")
    if (
        not isinstance(vntr_ids, tuple)
        or not vntr_ids
        or tuple(sorted(set(vntr_ids))) != vntr_ids
        or any(isinstance(value, bool) or not isinstance(value, int) or value <= 0 for value in vntr_ids)
    ):
        _fail("adVNTR capture VNTR identifiers must be sorted unique positive integers")
    _digest(model_sha256, "adVNTR model")
    if (background_path is None) != (background_sha256 is None):
        _fail("adVNTR background path and SHA256 must be present together")
    if background_sha256 is not None:
        _digest(background_sha256, "adVNTR background")
    policy_args = calibrated_policy_argv(capture_policy, caller_policy, background_path)
    argv = [
        *argv_prefix,
        "genotype",
        "-fs",
        "-vid",
        ",".join(str(value) for value in vntr_ids),
        "--alignment_file",
        str(alignment_path),
        "-o",
        str(result_path),
        "--outfmt",
        "vcf",
        "-m",
        str(model_path),
        "--working_directory",
        str(working_directory),
        "--frameshift-capture-version",
        "2",
        "--frameshift-calibration-out",
        str(sink_path),
    ]
    if reference_path is not None:
        argv.extend(("--reference_filename", str(reference_path)))
    argv.extend(policy_args)
    return CapturePlan(
        argv_prefix,
        tuple(argv),
        alignment_path,
        model_path,
        result_path,
        sink_path,
        working_directory,
        reference_path,
        vntr_ids,
        capture_policy,
        caller_policy,
        model_sha256,
        background_path,
        background_sha256,
    )


def _require_plan(plan: CapturePlan) -> CapturePlan:
    if not isinstance(plan, CapturePlan):
        _fail("adVNTR capture execution requires CapturePlan")
    rebuilt = build_capture_plan(
        plan.argv_prefix,
        alignment_path=plan.alignment_path,
        model_path=plan.model_path,
        result_path=plan.result_path,
        sink_path=plan.sink_path,
        working_directory=plan.working_directory,
        vntr_ids=plan.vntr_ids,
        capture_policy=plan.capture_policy,
        caller_policy=plan.caller_policy,
        model_sha256=plan.model_sha256,
        background_path=plan.background_path,
        background_sha256=plan.background_sha256,
        reference_path=plan.reference_path,
    )
    if rebuilt != plan:
        _fail("adVNTR capture plan differs from its typed argv")
    return plan


def _mapping(value: object, fields: set[str], label: str) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != fields:
        _fail(f"{label} fields differ from the closed upstream contract")
    return value


def _expected_producer(capabilities: AdvntrCapabilities) -> dict[str, object]:
    return {
        "package_version": capabilities.package_version,
        "build_id": capabilities.build_id,
        "source_revision": capabilities.source_revision,
    }


def capture_record_document(record: CaptureRecord) -> dict[str, object]:
    """Decode a fresh full upstream record after rechecking its immutable bytes."""
    if not isinstance(record, CaptureRecord) or hashlib.sha256(record.raw_json).hexdigest() != record.record_sha256:
        _fail("adVNTR capture record differs from its immutable bytes")
    document = load_strict_json_object(record.raw_json)
    if (
        advntr_canonical_sha256(document) != record.record_sha256
        or document.get("locus", {}).get("vntr_id") != record.vntr_id
    ):
        _fail("adVNTR capture record differs from its canonical content")
    return document


def decode_capture_sink(raw: bytes, plan: CapturePlan, capabilities: AdvntrCapabilities) -> CaptureResult:
    """Decode a complete newline-terminated sink while preserving full record content."""
    _require_plan(plan)
    advntr_capabilities_document(capabilities)
    if not isinstance(raw, bytes) or not raw or not raw.endswith(b"\n"):
        _fail("adVNTR capture sink must be nonempty complete newline-terminated JSONL")
    expected_capture = capture_policy_document(plan.capture_policy)
    expected_caller = _caller_document(plan.caller_policy)
    records: list[CaptureRecord] = []
    identifiers: list[int] = []
    for line in raw.splitlines():
        if not line:
            _fail("adVNTR capture sink contains an empty record")
        document = load_strict_json_object(line)
        _mapping(document, _CAPTURE_FIELDS, "adVNTR capture record")
        if document["schema_version"] != "advntr-frameshift-capture-v2" or document["completion"] != "completed-vntr":
            _fail("adVNTR capture record is not a completed v2 locus")
        digest = advntr_canonical_sha256(document)
        if hashlib.sha256(line).hexdigest() != digest:
            _fail("adVNTR capture record bytes must use the canonical upstream JSON encoding")
        if _mapping(document["producer"], _PRODUCER_FIELDS, "adVNTR capture producer") != _expected_producer(
            capabilities
        ):
            _fail("adVNTR capture producer differs from the capability preflight")
        assets = _mapping(document["assets"], _ASSET_FIELDS, "adVNTR capture assets")
        if assets["model_sha256"] != plan.model_sha256 or assets["background_sha256"] != plan.background_sha256:
            _fail("adVNTR capture model or background digest differs from the plan")
        if assets["model_locus_sha256"] != advntr_canonical_sha256(document["model_locus"]):
            _fail("adVNTR loaded model locus digest differs from its canonical content")
        loaded_background = document["loaded_background"]
        if plan.background_sha256 is None:
            if loaded_background is not None or assets["loaded_background_sha256"] is not None:
                _fail("adVNTR legacy capture must not bind a loaded background")
        elif not isinstance(loaded_background, dict) or assets["loaded_background_sha256"] != advntr_canonical_sha256(
            loaded_background
        ):
            _fail("adVNTR loaded background digest differs from its canonical content")
        if document["capture_policy"] != expected_capture or document["caller_policy"] != expected_caller:
            _fail("adVNTR capture policy differs from the plan")
        locus = _mapping(document["locus"], _LOCUS_FIELDS, "adVNTR capture locus")
        vntr_id = locus["vntr_id"]
        if isinstance(vntr_id, bool) or not isinstance(vntr_id, int) or vntr_id <= 0:
            _fail("adVNTR capture VNTR identifier must be a positive integer")
        identifiers.append(vntr_id)
        records.append(CaptureRecord(vntr_id, line, digest))
    if tuple(sorted(identifiers)) != plan.vntr_ids or len(identifiers) != len(set(identifiers)):
        _fail("adVNTR capture sink does not exactly match the planned VNTR roster")
    ordered = tuple(record for _, record in sorted(zip(identifiers, records, strict=True), key=lambda pair: pair[0]))
    return CaptureResult(hashlib.sha256(raw).hexdigest(), "", capabilities, ordered)


def _hash_file(path: Path, expected: str, label: str) -> bytes:
    content = read_regular_path(path)
    if hashlib.sha256(content).hexdigest() != expected:
        _fail(f"{label} bytes differ from the planned SHA256")
    return content


def run_capture(
    plan: CapturePlan,
    pin: AdvntrToolPin,
    *,
    runner: ProcessRunner = subprocess.run,
) -> CaptureResult:
    """Capability-preflight, run once, and decode only a successful complete capture."""
    plan = _require_plan(plan)
    capabilities = probe_advntr_capabilities(plan.argv_prefix, pin, runner=runner)
    if not plan.alignment_path.is_file() or not plan.working_directory.is_dir():
        _fail("adVNTR capture alignment or working directory is unavailable")
    if os.path.lexists(plan.sink_path) or os.path.lexists(plan.result_path):
        _fail("adVNTR capture outputs must be new paths")
    _hash_file(plan.model_path, plan.model_sha256, "adVNTR model")
    if plan.background_path is not None and plan.background_sha256 is not None:
        _hash_file(plan.background_path, plan.background_sha256, "adVNTR background")
    completed = runner(
        plan.argv,
        cwd=plan.working_directory,
        capture_output=True,
        text=True,
        check=False,
    )
    if not isinstance(completed, subprocess.CompletedProcess) or completed.returncode != 0:
        raise RuntimeError("adVNTR capture process failed")
    raw = read_regular_path(plan.sink_path)
    production_result = read_regular_path(plan.result_path)
    _hash_file(plan.model_path, plan.model_sha256, "adVNTR model")
    if plan.background_path is not None and plan.background_sha256 is not None:
        _hash_file(plan.background_path, plan.background_sha256, "adVNTR background")
    decoded = decode_capture_sink(raw, plan, capabilities)
    return CaptureResult(
        decoded.sink_sha256,
        hashlib.sha256(production_result).hexdigest(),
        capabilities,
        decoded.records,
    )


__all__ = [
    "CapturePlan",
    "CaptureRecord",
    "CaptureResult",
    "build_capture_plan",
    "caller_policy_document",
    "capture_policy_document",
    "capture_policy_for_caller",
    "capture_record_document",
    "decode_capture_sink",
    "run_capture",
]
