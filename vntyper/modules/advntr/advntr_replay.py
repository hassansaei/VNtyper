"""Shell-free adapter for packaged adVNTR production frameshift replay."""

from __future__ import annotations

import hashlib
import json
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
    capture_policy_document,
    probe_advntr_capabilities,
    replay_policy_document,
)
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, caller_policy_values_document
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader, read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)

_OUTPUT_FIELDS = {
    "schema_version",
    "manifest_file_sha256",
    "manifest_sha256",
    "policy_file_sha256",
    "policy_sha256",
    "background_file_sha256",
    "replay_producer",
    "results",
}
_RESULT_FIELDS = {"key", "capture_sha256", "vntrs"}
_VNTR_FIELDS = {"vntr_id", "result"}
_LOCUS_RESULT_FIELDS = {
    "schema_version",
    "vntr_id",
    "capture_record_sha256",
    "policy_sha256",
    "capture_producer",
    "capture_assets",
    "loaded_background_sha256",
    "baseline_parity",
    "decision_visits",
    "calls",
    "warnings",
    "capture_audit",
}
_AUDIT_FIELDS = {"attribution_outside_trials", "calibrated_policy_domain_errors"}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _digest(value: object, label: str, *, nullable: bool = False) -> str | None:
    if nullable and value is None:
        return None
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        _fail(f"{label} must be a lowercase SHA256 digest")
    return value


def _positive_ids(value: object) -> tuple[int, ...]:
    if (
        not isinstance(value, tuple)
        or not value
        or tuple(sorted(set(value))) != value
        or any(isinstance(item, bool) or not isinstance(item, int) or item <= 0 for item in value)
    ):
        _fail("adVNTR replay VNTR identifiers must be sorted unique positive integers")
    return value


@dataclass(frozen=True)
class ReplayCapture:
    """One exact flat capture file and its expected completed locus roster."""

    key: str
    filename: str
    sha256: str
    vntr_ids: tuple[int, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.key, str) or not self.key or self.key.strip() != self.key:
            _fail("adVNTR replay capture key must be nonempty trimmed text")
        if (
            not isinstance(self.filename, str)
            or not self.filename
            or Path(self.filename).name != self.filename
            or self.filename in {".", ".."}
        ):
            _fail("adVNTR replay capture filename must be one safe basename")
        _digest(self.sha256, "adVNTR replay capture")
        _positive_ids(self.vntr_ids)


@dataclass(frozen=True)
class ReplayPlan:
    """One exact upstream replay invocation and its semantic manifest/policy bytes."""

    argv_prefix: tuple[str, ...]
    argv: tuple[str, ...]
    capture_root: Path
    captures: tuple[ReplayCapture, ...]
    manifest_path: Path
    policy_path: Path
    output_directory: Path
    capture_policy: CapturePolicy
    caller_policy: CallerPolicyValues
    policy_json: bytes
    background_path: Path | None
    background_sha256: str | None


@dataclass(frozen=True)
class ReplayLocus:
    """One upstream locus result; failed capture audits remain explicitly unassessable."""

    key: str
    vntr_id: int
    assessable: bool
    raw_json: bytes
    sha256: str


@dataclass(frozen=True)
class ReplayResult:
    """Validated complete replay output with no locally recomputed decisions."""

    capabilities: AdvntrCapabilities
    output_sha256: str
    loci: tuple[ReplayLocus, ...]


def replay_manifest_document(captures: tuple[ReplayCapture, ...]) -> dict[str, object]:
    """Project the exact sorted capture-file and target roster expected upstream."""
    if not isinstance(captures, tuple) or not captures or any(not isinstance(item, ReplayCapture) for item in captures):
        _fail("adVNTR replay requires a nonempty ReplayCapture tuple")
    for item in captures:
        ReplayCapture(item.key, item.filename, item.sha256, item.vntr_ids)
    keys = tuple(item.key for item in captures)
    filenames = tuple(item.filename for item in captures)
    if keys != tuple(sorted(set(keys))) or len(filenames) != len(set(filenames)):
        _fail("adVNTR replay capture keys must be sorted unique and filenames distinct")
    return {
        "schema_version": "advntr-frameshift-replay-manifest-v1",
        "captures": [
            {"key": item.key, "filename": item.filename, "sha256": item.sha256, "vntr_ids": list(item.vntr_ids)}
            for item in captures
        ],
    }


def _absolute(path: object, label: str) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        _fail(f"{label} must be an absolute Path")
    return path


def _compact_json(document: object) -> bytes:
    try:
        return json.dumps(document, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode(
            "ascii"
        )
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise ValueError("adVNTR replay policy is not finite JSON-compatible content") from error


def build_replay_plan(
    argv_prefix: tuple[str, ...],
    *,
    capture_root: Path,
    captures: tuple[ReplayCapture, ...],
    manifest_path: Path,
    policy_path: Path,
    output_directory: Path,
    capture_policy: CapturePolicy,
    caller_policy: CallerPolicyValues,
    background_path: Path | None = None,
    background_sha256: str | None = None,
) -> ReplayPlan:
    """Build the exact installed replay argv without invoking or copying its evaluator."""
    if (
        not isinstance(argv_prefix, tuple)
        or not argv_prefix
        or any(not isinstance(token, str) or not token or "\x00" in token for token in argv_prefix)
    ):
        _fail("adVNTR executable prefix must be a nonempty tuple of command tokens")
    capture_root = _absolute(capture_root, "adVNTR replay capture root")
    manifest_path = _absolute(manifest_path, "adVNTR replay manifest")
    policy_path = _absolute(policy_path, "adVNTR replay policy")
    output_directory = _absolute(output_directory, "adVNTR replay output")
    if capture_root in manifest_path.parents:
        _fail("adVNTR replay manifest must remain outside its capture root")
    replay_manifest_document(captures)
    capture_policy_document(capture_policy)
    caller_policy_values_document(caller_policy)
    policy_document = replay_policy_document(capture_policy, caller_policy)
    mode = policy_document["caller_policy"]["mode"]  # type: ignore[index]
    if background_path is not None:
        background_path = _absolute(background_path, "adVNTR replay background")
    if (background_path is None) != (background_sha256 is None):
        _fail("adVNTR replay background path and digest must be present together")
    if mode == "exact" and background_path is None:
        _fail("adVNTR exact replay requires a background")
    if mode == "legacy" and background_path is not None:
        _fail("adVNTR legacy replay forbids a background")
    if background_sha256 is not None:
        _digest(background_sha256, "adVNTR replay background")
    argv = [
        *argv_prefix,
        "replay-frameshift",
        "--capture-root",
        str(capture_root),
        "--manifest",
        str(manifest_path),
        "--policy",
        str(policy_path),
    ]
    if background_path is not None:
        argv.extend(("--background", str(background_path)))
    argv.extend(("--output", str(output_directory)))
    return ReplayPlan(
        argv_prefix,
        tuple(argv),
        capture_root,
        captures,
        manifest_path,
        policy_path,
        output_directory,
        capture_policy,
        caller_policy,
        _compact_json(policy_document),
        background_path,
        background_sha256,
    )


def _require_plan(plan: ReplayPlan) -> ReplayPlan:
    if not isinstance(plan, ReplayPlan):
        _fail("adVNTR replay requires ReplayPlan")
    rebuilt = build_replay_plan(
        plan.argv_prefix,
        capture_root=plan.capture_root,
        captures=plan.captures,
        manifest_path=plan.manifest_path,
        policy_path=plan.policy_path,
        output_directory=plan.output_directory,
        capture_policy=plan.capture_policy,
        caller_policy=plan.caller_policy,
        background_path=plan.background_path,
        background_sha256=plan.background_sha256,
    )
    if rebuilt != plan:
        _fail("adVNTR replay plan differs from its typed argv or policy")
    return plan


def _mapping(value: object, expected: set[str], label: str) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != expected:
        _fail(f"{label} fields differ from the closed upstream contract")
    return value


def _producer(capabilities: AdvntrCapabilities) -> dict[str, object]:
    return {
        "package_version": capabilities.package_version,
        "build_id": capabilities.build_id,
        "source_revision": capabilities.source_revision,
    }


def _capture_record_hashes(plan: ReplayPlan) -> dict[tuple[str, int], str]:
    result: dict[tuple[str, int], str] = {}
    for capture in plan.captures:
        raw = read_regular_path(plan.capture_root / capture.filename)
        if hashlib.sha256(raw).hexdigest() != capture.sha256 or not raw.endswith(b"\n"):
            _fail("adVNTR replay capture bytes differ from the manifest")
        identifiers: list[int] = []
        for line in raw.splitlines():
            document = load_strict_json_object(line)
            locus = document.get("locus")
            vntr_id = locus.get("vntr_id") if isinstance(locus, dict) else None
            if isinstance(vntr_id, bool) or not isinstance(vntr_id, int) or vntr_id <= 0:
                _fail("adVNTR replay capture record lacks a positive VNTR identifier")
            identifiers.append(vntr_id)
            result[(capture.key, vntr_id)] = advntr_canonical_sha256(document)
        if tuple(sorted(identifiers)) != capture.vntr_ids or len(identifiers) != len(set(identifiers)):
            _fail("adVNTR replay capture record roster differs from the manifest")
    return result


def replay_locus_document(locus: ReplayLocus) -> dict[str, object]:
    """Return one fresh upstream result after rechecking its immutable bytes."""
    if not isinstance(locus, ReplayLocus) or hashlib.sha256(locus.raw_json).hexdigest() != locus.sha256:
        _fail("adVNTR replay locus differs from its immutable bytes")
    document = load_strict_json_object(locus.raw_json)
    audit = document.get("capture_audit")
    assessable = (
        isinstance(audit, dict)
        and not audit.get("attribution_outside_trials")
        and not audit.get("calibrated_policy_domain_errors")
    )
    if document.get("vntr_id") != locus.vntr_id or bool(assessable) != locus.assessable:
        _fail("adVNTR replay locus differs from its typed identity or audit status")
    return document


def _decode_locus(
    key: str,
    document: object,
    expected_vntr: int,
    expected_record_sha: str,
    expected_policy_sha: str,
    capabilities: AdvntrCapabilities,
    expects_background: bool,
) -> ReplayLocus:
    raw = _mapping(document, _LOCUS_RESULT_FIELDS, "adVNTR replay locus result")
    if (
        raw["schema_version"] != "advntr-frameshift-replay-result-v1"
        or raw["vntr_id"] != expected_vntr
        or raw["capture_record_sha256"] != expected_record_sha
        or raw["policy_sha256"] != expected_policy_sha
        or raw["capture_producer"] != _producer(capabilities)
        or raw["baseline_parity"] is not True
    ):
        _fail("adVNTR replay locus identity, producer, policy, or baseline parity differs")
    loaded_background_sha256 = _digest(
        raw["loaded_background_sha256"], "adVNTR loaded replay background", nullable=True
    )
    if (loaded_background_sha256 is not None) != expects_background:
        _fail("adVNTR loaded replay background binding differs from the requested caller mode")
    if not isinstance(raw["capture_assets"], dict):
        _fail("adVNTR replay capture assets must remain an object")
    for name in ("decision_visits", "calls", "warnings"):
        if not isinstance(raw[name], list):
            _fail(f"adVNTR replay {name} must remain a complete list")
    audit = _mapping(raw["capture_audit"], _AUDIT_FIELDS, "adVNTR replay capture audit")
    for name in _AUDIT_FIELDS:
        entries = audit[name]
        if not isinstance(entries, list) or any(not isinstance(item, str) for item in entries):
            _fail("adVNTR replay capture audit entries must be text lists")
    assessable = not audit["attribution_outside_trials"] and not audit["calibrated_policy_domain_errors"]
    encoded = _compact_json(raw)
    return ReplayLocus(key, expected_vntr, bool(assessable), encoded, hashlib.sha256(encoded).hexdigest())


def _decode_output(
    raw: bytes,
    plan: ReplayPlan,
    capabilities: AdvntrCapabilities,
    manifest_raw: bytes,
    policy_raw: bytes,
    record_hashes: dict[tuple[str, int], str],
) -> ReplayResult:
    document = _mapping(load_strict_json_object(raw), _OUTPUT_FIELDS, "adVNTR replay output")
    manifest_document = load_strict_json_object(manifest_raw)
    policy_document = load_strict_json_object(policy_raw)
    expected_policy_sha = advntr_canonical_sha256(policy_document)
    if (
        document["schema_version"] != "advntr-frameshift-replay-output-v1"
        or document["manifest_file_sha256"] != hashlib.sha256(manifest_raw).hexdigest()
        or document["manifest_sha256"] != advntr_canonical_sha256(manifest_document)
        or document["policy_file_sha256"] != hashlib.sha256(policy_raw).hexdigest()
        or document["policy_sha256"] != expected_policy_sha
        or document["background_file_sha256"] != plan.background_sha256
        or document["replay_producer"] != _producer(capabilities)
    ):
        _fail("adVNTR replay output manifest, policy, background, or producer binding differs")
    rows = document["results"]
    if not isinstance(rows, list):
        _fail("adVNTR replay results must be a list")
    expected = {capture.key: capture for capture in plan.captures}
    loci: list[ReplayLocus] = []
    seen: list[str] = []
    for value in rows:
        row = _mapping(value, _RESULT_FIELDS, "adVNTR replay result row")
        key = row["key"]
        if not isinstance(key, str) or key not in expected:
            _fail("adVNTR replay result key is undeclared")
        capture = expected[key]
        if row["capture_sha256"] != capture.sha256 or not isinstance(row["vntrs"], list):
            _fail("adVNTR replay result capture binding differs")
        identifiers: list[int] = []
        for value_vntr in row["vntrs"]:
            vntr = _mapping(value_vntr, _VNTR_FIELDS, "adVNTR replay VNTR row")
            vntr_id = vntr["vntr_id"]
            if not isinstance(vntr_id, int) or isinstance(vntr_id, bool) or (key, vntr_id) not in record_hashes:
                _fail("adVNTR replay VNTR row is undeclared")
            identifiers.append(vntr_id)
            loci.append(
                _decode_locus(
                    key,
                    vntr["result"],
                    vntr_id,
                    record_hashes[(key, vntr_id)],
                    expected_policy_sha,
                    capabilities,
                    plan.background_path is not None,
                )
            )
        if tuple(identifiers) != capture.vntr_ids:
            _fail("adVNTR replay result VNTR roster differs from its capture")
        seen.append(key)
    if tuple(seen) != tuple(expected):
        _fail("adVNTR replay result keys differ from the planned capture roster")
    return ReplayResult(capabilities, hashlib.sha256(raw).hexdigest(), tuple(loci))


def run_replay(
    plan: ReplayPlan,
    pin: AdvntrToolPin,
    *,
    runner: ProcessRunner = subprocess.run,
) -> ReplayResult:
    """Invoke packaged replay and retain audit failures as unassessable results."""
    plan = _require_plan(plan)
    capabilities = probe_advntr_capabilities(plan.argv_prefix, pin, runner=runner)
    if not plan.capture_root.is_dir() or os.path.lexists(plan.output_directory):
        _fail("adVNTR replay capture root must exist and output must be new")
    manifest_raw = read_regular_path(plan.manifest_path)
    policy_raw = read_regular_path(plan.policy_path)
    if load_strict_json_object(manifest_raw) != replay_manifest_document(plan.captures):
        _fail("adVNTR replay manifest content differs from the typed plan")
    if policy_raw.rstrip(b"\n") != plan.policy_json or load_strict_json_object(policy_raw) != load_strict_json_object(
        plan.policy_json
    ):
        _fail("adVNTR replay policy content differs from the typed plan")
    if plan.background_path is not None and plan.background_sha256 is not None:
        background_raw = read_regular_path(plan.background_path)
        if hashlib.sha256(background_raw).hexdigest() != plan.background_sha256:
            _fail("adVNTR replay background bytes differ from the typed plan")
    record_hashes = _capture_record_hashes(plan)
    completed = runner(plan.argv, capture_output=True, text=True, check=False)
    if not isinstance(completed, subprocess.CompletedProcess) or completed.returncode != 0:
        raise RuntimeError("adVNTR replay process failed")
    with SecureDirectoryReader.open(plan.output_directory, {"replay.json"}) as reader:
        raw = reader.read_file("replay.json")
    return _decode_output(raw, plan, capabilities, manifest_raw, policy_raw, record_hashes)
