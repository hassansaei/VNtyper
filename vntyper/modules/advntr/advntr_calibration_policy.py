"""Strict VNtyper projection of the packaged adVNTR calibration ABI."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import re
import subprocess
from collections.abc import Callable, Mapping
from dataclasses import dataclass, fields, replace
from typing import NoReturn, cast

from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, caller_policy_values_document
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)

ADVNTR_CAPABILITIES = (
    "fit-background-installed-v1",
    "frameshift-calibration-capture-v1",
    "frameshift-run-local-policy-v1",
    "frameshift-calibration-capture-v2",
    "frameshift-replay-v1",
)
ADVNTR_CAPTURE_SCHEMA_VERSIONS = (1, 2)
ADVNTR_POLICY_SCHEMA_VERSIONS = (
    "advntr-frameshift-policy-v1",
    "advntr-frameshift-replay-policy-v1",
)
ADVNTR_BACKGROUND_RECIPE_IDS = ("recipe-v1",)
_CAPABILITY_FIELDS = {
    "schema_version",
    "package_version",
    "build_id",
    "source_revision",
    "capabilities",
    "capture_schema_versions",
    "policy_schema_versions",
    "background_recipe_ids",
}
_CAPTURE_FIELDS = {
    "platform",
    "frameshift_mode",
    "is_haploid",
    "caller_mode",
    "threads",
    "minimum_read_length",
    "prune_reverse",
    "filter_adapter_readthrough",
    "minimum_read_match_ratio",
    "minimum_relative_ru_coverage",
    "use_reference_alignment",
    "fully_covered_ru_only",
    "maximum_error_rate",
    "legacy_error_rate",
    "mapq_cutoff",
    "base_quality_cutoff",
    "maximum_low_quality_fraction",
    "enhanced_hmm",
    "trained_hmms",
}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_REVISION = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_ADVNTR_PREFIX = "/components/advntr/calibrated_calling/"


@dataclass(frozen=True)
class AdvntrCapabilities:
    """Installed adVNTR build identity and exact finite feature inventory."""

    package_version: str
    build_id: str
    source_revision: str | None
    capabilities: tuple[str, ...]
    capture_schema_versions: tuple[int, ...]
    policy_schema_versions: tuple[str, ...]
    background_recipe_ids: tuple[str, ...]
    sha256: str


@dataclass(frozen=True)
class AdvntrToolPin:
    """Exact installed adVNTR identity authorized by a caller study."""

    package_version: str
    build_id: str
    source_revision: str

    def __post_init__(self) -> None:
        _text(self.package_version, "adVNTR pinned package version")
        _digest(self.build_id, "adVNTR pinned build")
        if not isinstance(self.source_revision, str) or _REVISION.fullmatch(self.source_revision) is None:
            _fail("adVNTR pinned source revision must be a lowercase 40- or 64-character hash")


@dataclass(frozen=True)
class CapturePolicy:
    """Complete upstream raw capture policy; calibrated profiles narrow its domains."""

    platform: str
    frameshift_mode: bool
    is_haploid: bool
    caller_mode: str
    threads: int
    minimum_read_length: int | None
    prune_reverse: bool
    filter_adapter_readthrough: bool
    minimum_read_match_ratio: float | None
    minimum_relative_ru_coverage: float | None
    use_reference_alignment: bool
    fully_covered_ru_only: bool
    maximum_error_rate: float
    legacy_error_rate: float
    mapq_cutoff: int
    base_quality_cutoff: int
    maximum_low_quality_fraction: float
    enhanced_hmm: bool
    trained_hmms: bool
    sha256: str


ProcessRunner = Callable[..., subprocess.CompletedProcess[str]]


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value or "\x00" in value:
        _fail(f"{label} must be nonempty trimmed text")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail(f"{label} must be a lowercase SHA256 digest")
    return value


def _object(value: object, expected: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != expected:
        _fail(f"{label} fields differ from the closed contract")
    return value


def _exact_list(value: object, expected: tuple[object, ...], label: str) -> tuple[object, ...]:
    if not isinstance(value, list) or tuple(value) != expected:
        _fail(f"{label} differs from the published adVNTR ABI")
    return tuple(value)


def decode_advntr_capabilities(value: object) -> AdvntrCapabilities:
    """Decode the exact installed adVNTR 2.4 capability response."""
    raw = _object(value, _CAPABILITY_FIELDS, "adVNTR capabilities")
    if raw["schema_version"] != "advntr-capabilities-v1":
        _fail("adVNTR capability schema is unsupported")
    package_version = _text(raw["package_version"], "adVNTR package version")
    build_id = _digest(raw["build_id"], "adVNTR build ID")
    revision = raw["source_revision"]
    if revision is not None and (not isinstance(revision, str) or _REVISION.fullmatch(revision) is None):
        _fail("adVNTR source revision must be null or a lowercase 40- or 64-character hash")
    capabilities = _exact_list(raw["capabilities"], ADVNTR_CAPABILITIES, "adVNTR capabilities")
    capture_versions = _exact_list(
        raw["capture_schema_versions"], ADVNTR_CAPTURE_SCHEMA_VERSIONS, "adVNTR capture schema versions"
    )
    if any(isinstance(item, bool) or not isinstance(item, int) for item in capture_versions):
        _fail("adVNTR capture schema versions must be integers")
    policy_versions = _exact_list(
        raw["policy_schema_versions"], ADVNTR_POLICY_SCHEMA_VERSIONS, "adVNTR policy schema versions"
    )
    recipe_ids = _exact_list(raw["background_recipe_ids"], ADVNTR_BACKGROUND_RECIPE_IDS, "adVNTR background recipes")
    decoded = AdvntrCapabilities(
        package_version,
        build_id,
        revision,
        tuple(str(item) for item in capabilities),
        tuple(cast(int, item) for item in capture_versions),
        tuple(str(item) for item in policy_versions),
        tuple(str(item) for item in recipe_ids),
        "",
    )
    return replace(decoded, sha256=advntr_canonical_sha256(_capabilities_document(decoded)))


def _capabilities_document(capabilities: AdvntrCapabilities) -> dict[str, object]:
    return {
        "schema_version": "advntr-capabilities-v1",
        "package_version": capabilities.package_version,
        "build_id": capabilities.build_id,
        "source_revision": capabilities.source_revision,
        "capabilities": list(capabilities.capabilities),
        "capture_schema_versions": list(capabilities.capture_schema_versions),
        "policy_schema_versions": list(capabilities.policy_schema_versions),
        "background_recipe_ids": list(capabilities.background_recipe_ids),
    }


def advntr_capabilities_document(capabilities: AdvntrCapabilities) -> dict[str, object]:
    """Project and revalidate an immutable capability contract."""
    if not isinstance(capabilities, AdvntrCapabilities):
        _fail("adVNTR capability projection requires AdvntrCapabilities")
    raw = _capabilities_document(capabilities)
    if decode_advntr_capabilities(raw) != capabilities:
        _fail("adVNTR capabilities differ from their canonical content")
    return raw


def require_advntr_capabilities(capabilities: AdvntrCapabilities, pin: AdvntrToolPin) -> AdvntrCapabilities:
    """Require an identified installed build matching the study's exact pin."""
    advntr_capabilities_document(capabilities)
    if not isinstance(pin, AdvntrToolPin):
        _fail("adVNTR capability requirement needs AdvntrToolPin")
    AdvntrToolPin(pin.package_version, pin.build_id, pin.source_revision)
    if (capabilities.package_version, capabilities.build_id, capabilities.source_revision) != (
        pin.package_version,
        pin.build_id,
        pin.source_revision,
    ):
        _fail("installed adVNTR identity differs from the study pin")
    return capabilities


def _argv_prefix(value: object) -> tuple[str, ...]:
    if not isinstance(value, tuple) or not value:
        _fail("adVNTR executable prefix must be a nonempty tuple")
    for token in value:
        _text(token, "adVNTR executable token")
    return value


def probe_advntr_capabilities(
    argv_prefix: tuple[str, ...],
    pin: AdvntrToolPin,
    *,
    runner: ProcessRunner = subprocess.run,
) -> AdvntrCapabilities:
    """Probe the installed executable with a shell-free capability command."""
    prefix = _argv_prefix(argv_prefix)
    completed = runner((*prefix, "capabilities", "--json"), capture_output=True, text=True, check=False)
    if not isinstance(completed, subprocess.CompletedProcess) or completed.returncode != 0:
        raise RuntimeError("adVNTR capability probe failed")
    try:
        document = load_strict_json_object(completed.stdout)
    except (TypeError, ValueError) as error:
        raise ValueError("adVNTR capability probe returned invalid strict JSON") from error
    return require_advntr_capabilities(decode_advntr_capabilities(document), pin)


def advntr_canonical_sha256(value: object) -> str:
    """Hash the exact compact ASCII JSON representation used by adVNTR 2.4."""
    try:
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode(
            "ascii"
        )
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise ValueError("adVNTR document is not finite JSON-compatible content") from error
    return hashlib.sha256(encoded).hexdigest()


def _integer(value: object, label: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _fail(f"{label} must be an integer >= {minimum}")
    return value


def _fraction(value: object, label: str, maximum: float | None = 1.0) -> float:
    if (
        not isinstance(value, float)
        or not math.isfinite(value)
        or value < 0
        or (maximum is not None and value > maximum)
    ):
        _fail(f"{label} must be a finite nonnegative float in its declared domain")
    return value


def _boolean(value: object, label: str) -> bool:
    if not isinstance(value, bool):
        _fail(f"{label} must be boolean")
    return value


def _capture_document(policy: CapturePolicy) -> dict[str, object]:
    parameters = {field.name: getattr(policy, field.name) for field in fields(CapturePolicy) if field.name != "sha256"}
    return {"schema_version": "advntr-runtime-capture-policy-v1", "parameters": parameters}


def decode_capture_policy(value: object) -> CapturePolicy:
    """Decode all upstream raw capture fields without inventing defaults."""
    raw = _object(value, {"schema_version", "parameters"}, "adVNTR capture policy")
    if raw["schema_version"] != "advntr-runtime-capture-policy-v1":
        _fail("adVNTR capture policy schema is unsupported")
    if not isinstance(raw["parameters"], Mapping) or set(raw["parameters"]) != _CAPTURE_FIELDS:
        _fail("adVNTR capture policy parameters must contain every field exactly once")
    parameters = raw["parameters"]
    platform = parameters["platform"]
    caller_mode = parameters["caller_mode"]
    if platform not in {"illumina", "pacbio", "nanopore"} or not isinstance(platform, str):
        _fail("adVNTR capture platform is unsupported")
    if caller_mode not in {"legacy", "exact"} or not isinstance(caller_mode, str):
        _fail("adVNTR capture caller mode is unsupported")
    boolean_fields = (
        "frameshift_mode",
        "is_haploid",
        "prune_reverse",
        "filter_adapter_readthrough",
        "use_reference_alignment",
        "fully_covered_ru_only",
        "enhanced_hmm",
        "trained_hmms",
    )
    booleans = {name: _boolean(parameters[name], name) for name in boolean_fields}
    if caller_mode == "exact" and not booleans["frameshift_mode"]:
        _fail("adVNTR exact caller mode requires frameshift capture")
    if not booleans["enhanced_hmm"] or booleans["trained_hmms"]:
        _fail("adVNTR capture requires the published enhanced untrained HMM mode")
    minimum_length = parameters["minimum_read_length"]
    if minimum_length is not None:
        minimum_length = _integer(minimum_length, "minimum_read_length", 1)
    match_ratio = parameters["minimum_read_match_ratio"]
    if match_ratio is not None:
        match_ratio = _fraction(match_ratio, "minimum_read_match_ratio")
    rare_fraction = parameters["minimum_relative_ru_coverage"]
    if rare_fraction is not None:
        rare_fraction = _fraction(rare_fraction, "minimum_relative_ru_coverage", None)
    parsed = CapturePolicy(
        platform=platform,
        frameshift_mode=booleans["frameshift_mode"],
        is_haploid=booleans["is_haploid"],
        caller_mode=caller_mode,
        threads=_integer(parameters["threads"], "threads", 1),
        minimum_read_length=minimum_length,
        prune_reverse=booleans["prune_reverse"],
        filter_adapter_readthrough=booleans["filter_adapter_readthrough"],
        minimum_read_match_ratio=match_ratio,
        minimum_relative_ru_coverage=rare_fraction,
        use_reference_alignment=booleans["use_reference_alignment"],
        fully_covered_ru_only=booleans["fully_covered_ru_only"],
        maximum_error_rate=_fraction(parameters["maximum_error_rate"], "maximum_error_rate"),
        legacy_error_rate=_fraction(parameters["legacy_error_rate"], "legacy_error_rate"),
        mapq_cutoff=_integer(parameters["mapq_cutoff"], "mapq_cutoff", 0),
        base_quality_cutoff=_integer(parameters["base_quality_cutoff"], "base_quality_cutoff", 0),
        maximum_low_quality_fraction=_fraction(
            parameters["maximum_low_quality_fraction"], "maximum_low_quality_fraction"
        ),
        enhanced_hmm=booleans["enhanced_hmm"],
        trained_hmms=booleans["trained_hmms"],
        sha256="",
    )
    return replace(parsed, sha256=advntr_canonical_sha256(_capture_document(parsed)))


def capture_policy_document(policy: CapturePolicy) -> dict[str, object]:
    """Project and revalidate a complete immutable raw capture policy."""
    if not isinstance(policy, CapturePolicy):
        _fail("adVNTR capture policy projection requires CapturePolicy")
    raw = _capture_document(policy)
    if decode_capture_policy(raw) != policy:
        _fail("adVNTR capture policy differs from its canonical content or digest")
    return raw


def capture_policy_for_caller(baseline: CapturePolicy, caller: CallerPolicyValues) -> CapturePolicy:
    """Project the seven calibrated values onto a complete baseline raw policy."""
    capture_policy_document(baseline)
    caller_policy_values_document(caller)
    if "advntr" not in caller.required_callers:
        _fail("calibrated adVNTR capture requires adVNTR in required_callers")
    values = caller.values
    match_ratio = cast(float, values[f"{_ADVNTR_PREFIX}minimum_read_match_ratio"])
    rare_fraction = cast(float | None, values[f"{_ADVNTR_PREFIX}rare_unit_fraction"])
    projected = replace(
        baseline,
        caller_mode=str(values[f"{_ADVNTR_PREFIX}mode"]),
        prune_reverse=bool(values[f"{_ADVNTR_PREFIX}prune_reverse"]),
        filter_adapter_readthrough=bool(values[f"{_ADVNTR_PREFIX}adapter_filter"]),
        minimum_read_match_ratio=match_ratio,
        minimum_relative_ru_coverage=rare_fraction,
        sha256="",
    )
    return decode_capture_policy(_capture_document(projected))


def replay_policy_document(capture: CapturePolicy, caller: CallerPolicyValues) -> dict[str, object]:
    """Build the closed upstream replay policy without implementing its evaluator."""
    capture_policy_document(capture)
    caller_policy_values_document(caller)
    if "advntr" not in caller.required_callers:
        _fail("adVNTR replay policy requires adVNTR in required_callers")
    values = caller.values
    mode = values[f"{_ADVNTR_PREFIX}mode"]
    if capture.caller_mode != mode:
        _fail("adVNTR capture and replay caller modes differ")
    document: dict[str, object] = {
        "schema_version": "advntr-frameshift-replay-policy-v1",
        "capture_policy": capture_policy_document(capture),
        "caller_policy": {
            "schema_version": "advntr-frameshift-policy-v1",
            "mode": mode,
            "cutoff": values[f"{_ADVNTR_PREFIX}cutoff"],
            "minimum_read_support": values[f"{_ADVNTR_PREFIX}minimum_read_support"],
        },
    }
    advntr_canonical_sha256(document)
    return document
