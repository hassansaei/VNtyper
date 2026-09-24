"""Native adVNTR replay adapter for finite legacy cutoff grids."""

from __future__ import annotations

import hashlib
import logging
import subprocess
from collections.abc import Mapping
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
    capture_policy_for_caller,
    decode_capture_policy,
)
from vntyper.modules.advntr.advntr_capture import caller_policy_document as native_caller_policy_document
from vntyper.modules.advntr.advntr_replay import (
    ReplayCapture,
    ReplayLocus,
    build_replay_plan,
    replay_locus_document,
    replay_manifest_document,
    run_replay,
)
from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    caller_policy_values_document,
)
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object

logger = logging.getLogger(__name__)

_ADVNTR_PREFIX = "/components/advntr/calibrated_calling/"
_TUNABLE_POINTERS = {
    f"{_ADVNTR_PREFIX}cutoff",
    f"{_ADVNTR_PREFIX}minimum_read_support",
}
_PRODUCER_FIELDS = {"package_version", "build_id", "source_revision"}
_CALLER_FIELDS = {"schema_version", "mode", "cutoff", "minimum_read_support"}


@dataclass(frozen=True)
class AdvntrCutoffSample:
    """One sample-level native replay outcome with full opaque locus results."""

    key: str
    assessable: bool
    called_positive: bool | None
    loci: tuple[ReplayLocus, ...]


@dataclass(frozen=True)
class AdvntrCutoffPolicyResult:
    """One requested full policy projected onto a deduplicated adVNTR execution."""

    policy_id: str
    policy_sha256: str
    execution_id: str
    samples: tuple[AdvntrCutoffSample, ...]


@dataclass(frozen=True)
class AdvntrCutoffGridResult:
    """Complete native legacy replay grid and its installed producer identity."""

    output: Path
    baseline_policy_id: str
    capture_policy_sha256: str
    capabilities: AdvntrCapabilities
    policies: tuple[AdvntrCutoffPolicyResult, ...]
    sha256: str


@dataclass(frozen=True)
class _CaptureSnapshot:
    key: str
    raw: bytes
    sha256: str
    vntr_ids: tuple[int, ...]


@dataclass(frozen=True)
class _CaptureContext:
    captures: tuple[_CaptureSnapshot, ...]
    capture_policy: CapturePolicy
    native_caller_policy: dict[str, object]
    pin: AdvntrToolPin


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"{label} fields differ from the closed contract")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value or "\x00" in value:
        _fail(f"{label} must be nonempty trimmed text")
    return value


def _absolute(path: object, label: str) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        _fail(f"{label} must be an absolute Path")
    return path


def _capture_metadata(document: Mapping[str, object]) -> tuple[CapturePolicy, dict[str, object], AdvntrToolPin, int]:
    if (
        document.get("schema_version") != "advntr-frameshift-capture-v2"
        or document.get("completion") != "completed-vntr"
    ):
        _fail("adVNTR cutoff input must contain completed v2 capture records")
    producer = _object(document.get("producer"), _PRODUCER_FIELDS, "adVNTR cutoff capture producer")
    package_version = _text(producer["package_version"], "adVNTR capture package version")
    build_id = producer["build_id"]
    source_revision = producer["source_revision"]
    if not isinstance(build_id, str) or not isinstance(source_revision, str):
        _fail("adVNTR cutoff capture producer identity is incomplete")
    pin = AdvntrToolPin(package_version, build_id, source_revision)
    capture_policy = decode_capture_policy(document.get("capture_policy"))
    caller = dict(_object(document.get("caller_policy"), _CALLER_FIELDS, "adVNTR cutoff baseline caller policy"))
    if caller["schema_version"] != "advntr-frameshift-policy-v1":
        _fail("adVNTR cutoff baseline caller policy schema is unsupported")
    locus = document.get("locus")
    vntr_id = locus.get("vntr_id") if isinstance(locus, Mapping) else None
    if isinstance(vntr_id, bool) or not isinstance(vntr_id, int) or vntr_id <= 0:
        _fail("adVNTR cutoff capture locus must contain a positive VNTR identifier")
    return capture_policy, caller, pin, vntr_id


def _load_capture_context(capture_paths: Mapping[str, Path]) -> _CaptureContext:
    if not isinstance(capture_paths, Mapping) or not capture_paths:
        _fail("adVNTR cutoff captures must be a nonempty sample-to-Path mapping")
    snapshots: list[_CaptureSnapshot] = []
    policies: list[CapturePolicy] = []
    callers: list[dict[str, object]] = []
    pins: list[AdvntrToolPin] = []
    source_paths: set[Path] = set()
    for key, path in sorted(capture_paths.items()):
        sample_key = _text(key, "adVNTR cutoff sample key")
        source = _absolute(path, "adVNTR cutoff capture")
        normalized = source.resolve(strict=False)
        if normalized in source_paths:
            _fail("adVNTR cutoff capture paths must be distinct")
        source_paths.add(normalized)
        raw = read_regular_path(source)
        if not raw or not raw.endswith(b"\n") or any(not line for line in raw.splitlines()):
            _fail("adVNTR cutoff capture must be complete nonempty newline-terminated JSONL")
        digest = hashlib.sha256(raw).hexdigest()
        identifiers: list[int] = []
        for line in raw.splitlines():
            document = load_strict_json_object(line)
            if hashlib.sha256(line).hexdigest() != advntr_canonical_sha256(document):
                _fail("adVNTR cutoff capture records must use canonical upstream JSON")
            capture_policy, caller, pin, vntr_id = _capture_metadata(document)
            policies.append(capture_policy)
            callers.append(caller)
            pins.append(pin)
            identifiers.append(vntr_id)
        if tuple(sorted(set(identifiers))) != tuple(identifiers):
            _fail("adVNTR cutoff capture VNTR identifiers must be sorted and unique")
        snapshots.append(_CaptureSnapshot(sample_key, raw, digest, tuple(identifiers)))
    if any(policy != policies[0] for policy in policies[1:]):
        _fail("adVNTR cutoff captures must share one capture policy")
    if any(caller != callers[0] for caller in callers[1:]):
        _fail("adVNTR cutoff captures must share one baseline caller policy")
    if any(pin != pins[0] for pin in pins[1:]):
        _fail("adVNTR cutoff captures must share one producer identity")
    return _CaptureContext(tuple(snapshots), policies[0], callers[0], pins[0])


def _validate_policies(
    policies: Mapping[str, CallerPolicyValues], baseline_policy_id: str, context: _CaptureContext
) -> tuple[tuple[str, CallerPolicyValues], ...]:
    if not isinstance(policies, Mapping) or not policies:
        _fail("adVNTR cutoff policies must be a nonempty mapping")
    baseline_id = _text(baseline_policy_id, "adVNTR cutoff baseline policy ID")
    rows: list[tuple[str, CallerPolicyValues]] = []
    for policy_id, policy in sorted(policies.items()):
        name = _text(policy_id, "adVNTR cutoff policy ID")
        caller_policy_values_document(policy)
        if "advntr" not in policy.required_callers:
            _fail("adVNTR cutoff policies must require adVNTR")
        rows.append((name, policy))
    by_name = dict(rows)
    if baseline_id not in by_name:
        _fail("adVNTR cutoff baseline policy ID is undeclared")
    baseline = by_name[baseline_id]
    values = baseline.values
    if values[f"{_ADVNTR_PREFIX}mode"] != "legacy":
        _fail("adVNTR cutoff optimization currently requires a legacy baseline")
    if native_caller_policy_document(baseline) != context.native_caller_policy:
        _fail("adVNTR capture baseline caller policy differs from the declared baseline")
    if capture_policy_for_caller(context.capture_policy, baseline) != context.capture_policy:
        _fail("adVNTR capture policy differs from the declared baseline")
    fixed = set(ADVNTR_CALLER_POLICY_POINTERS) - _TUNABLE_POINTERS
    for _, policy in rows:
        if policy.values[f"{_ADVNTR_PREFIX}mode"] != "legacy" or any(
            policy.values[pointer] != values[pointer] for pointer in fixed
        ):
            _fail("adVNTR cutoff grid may vary only legacy cutoff and support")
        if capture_policy_for_caller(context.capture_policy, policy) != context.capture_policy:
            _fail("adVNTR cutoff candidate changes captured native semantics")
    return tuple(rows)


def advntr_signature(policy: CallerPolicyValues) -> tuple[object, ...]:
    """The adVNTR part of a complete policy, which decides its native execution.

    Two policies with equal signatures differ only in other callers' values, so the grid
    replays them once, and the optimize execution guard counts them once.

    Args:
        policy: A complete caller policy that includes adVNTR.

    Returns:
        The adVNTR policy values, in ``ADVNTR_CALLER_POLICY_POINTERS`` order.
    """
    return tuple(policy.values[pointer] for pointer in ADVNTR_CALLER_POLICY_POINTERS)


def _execution_id(policy: CallerPolicyValues) -> str:
    document = {pointer: policy.values[pointer] for pointer in ADVNTR_CALLER_POLICY_POINTERS}
    return f"advntr-{advntr_canonical_sha256(document)}"


def _write_private(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_bytes(raw)
    path.chmod(0o600)


def _sample_results(loci: tuple[ReplayLocus, ...], keys: tuple[str, ...]) -> tuple[AdvntrCutoffSample, ...]:
    by_key: dict[str, list[ReplayLocus]] = {key: [] for key in keys}
    for locus in loci:
        if locus.key not in by_key:
            _fail("adVNTR cutoff replay returned an undeclared sample")
        by_key[locus.key].append(locus)
    results = []
    for key in keys:
        sample_loci = tuple(by_key[key])
        if not sample_loci:
            _fail("adVNTR cutoff replay omitted a declared sample")
        assessable = all(locus.assessable for locus in sample_loci)
        called = None
        if assessable:
            called = any(bool(replay_locus_document(locus)["calls"]) for locus in sample_loci)
        results.append(AdvntrCutoffSample(key, assessable, called, sample_loci))
    return tuple(results)


def _result_document(
    baseline_policy_id: str,
    context: _CaptureContext,
    capabilities: AdvntrCapabilities,
    policies: tuple[AdvntrCutoffPolicyResult, ...],
) -> dict[str, object]:
    return {
        "schema_version": "calibration-advntr-cutoff-grid-result-v1",
        "baseline_policy_id": baseline_policy_id,
        "capture_policy_sha256": context.capture_policy.sha256,
        "tool_identity": advntr_capabilities_document(capabilities),
        "policies": [
            {
                "policy_id": row.policy_id,
                "policy_sha256": row.policy_sha256,
                "execution_id": row.execution_id,
                "samples": [
                    {
                        "key": sample.key,
                        "assessable": sample.assessable,
                        "called_positive": sample.called_positive,
                        "loci": [
                            {"vntr_id": locus.vntr_id, "assessable": locus.assessable, "sha256": locus.sha256}
                            for locus in sample.loci
                        ],
                    }
                    for sample in row.samples
                ],
            }
            for row in policies
        ],
    }


def _private_tree(root: Path) -> None:
    for path in root.rglob("*"):
        path.chmod(0o700 if path.is_dir() else 0o600)


def evaluate_advntr_cutoff_grid(
    capture_paths: Mapping[str, Path],
    policies: Mapping[str, CallerPolicyValues],
    *,
    baseline_policy_id: str,
    executable_path: Path,
    output: Path,
    runner: ProcessRunner = subprocess.run,
) -> AdvntrCutoffGridResult:
    """Replay each unique legacy adVNTR cutoff/support policy with the installed evaluator.

    Kestrel-only policy differences share one native replay execution. Audit failures are
    retained as unassessable and never converted into negative observations.

    Args:
        capture_paths: Sample key to complete upstream v2 capture JSONL.
        policies: Full caller policy values keyed by candidate ID, including baseline.
        baseline_policy_id: Policy whose native capture decisions must replay exactly.
        executable_path: Absolute packaged adVNTR executable.
        output: New private atomic output directory.
        runner: Subprocess seam for deterministic tests.

    Returns:
        Immutable policy/sample observations plus full opaque upstream locus results.

    Raises:
        ValueError: If inputs, capture provenance, policies, or upstream output differ.
        RuntimeError: If the native tool fails or atomic publication is unavailable.
    """
    executable = _absolute(executable_path, "adVNTR cutoff executable")
    destination = _absolute(output, "adVNTR cutoff output")
    if not callable(runner):
        _fail("adVNTR cutoff runner must be callable")
    context = _load_capture_context(capture_paths)
    policy_rows = _validate_policies(policies, baseline_policy_id, context)
    grouped: dict[tuple[object, ...], list[tuple[str, CallerPolicyValues]]] = {}
    for row in policy_rows:
        grouped.setdefault(advntr_signature(row[1]), []).append(row)
    built: list[AdvntrCutoffPolicyResult] = []
    observed_capabilities: AdvntrCapabilities | None = None
    document_holder: list[dict[str, object]] = []

    def produce(staging: Path) -> bool:
        nonlocal observed_capabilities
        capture_root = staging / "captures"
        capture_root.mkdir(mode=0o700)
        replay_captures = []
        for index, snapshot in enumerate(context.captures):
            filename = f"capture-{index:06d}-{snapshot.sha256}.jsonl"
            _write_private(capture_root / filename, snapshot.raw)
            replay_captures.append(ReplayCapture(snapshot.key, filename, snapshot.sha256, snapshot.vntr_ids))
        captures = tuple(replay_captures)
        manifest_path = staging / "manifest.json"
        _write_private(manifest_path, canonical_json_bytes(replay_manifest_document(captures)) + b"\n")
        outcomes: dict[tuple[object, ...], tuple[str, tuple[AdvntrCutoffSample, ...]]] = {}
        for signature, aliases in sorted(grouped.items(), key=lambda item: repr(item[0])):
            representative = aliases[0][1]
            execution_id = _execution_id(representative)
            execution_root = staging / "executions" / execution_id
            execution_root.mkdir(parents=True, mode=0o700)
            policy_path = execution_root / "policy.json"
            replay_output = execution_root / "replay"
            plan = build_replay_plan(
                (str(executable),),
                capture_root=capture_root,
                captures=captures,
                manifest_path=manifest_path,
                policy_path=policy_path,
                output_directory=replay_output,
                capture_policy=context.capture_policy,
                caller_policy=representative,
            )
            _write_private(policy_path, plan.policy_json + b"\n")
            replay = run_replay(plan, context.pin, runner=runner)
            if observed_capabilities is None:
                observed_capabilities = replay.capabilities
            elif replay.capabilities != observed_capabilities:
                _fail("adVNTR cutoff executions observed different tool identities")
            samples = _sample_results(replay.loci, tuple(snapshot.key for snapshot in context.captures))
            outcomes[signature] = (execution_id, samples)
            receipt = {
                "schema_version": "calibration-advntr-cutoff-execution-v1",
                "execution_id": execution_id,
                "policy_ids": [name for name, _ in aliases],
                "policy_sha256": advntr_canonical_sha256(load_strict_json_object(plan.policy_json)),
                "capture_policy_sha256": context.capture_policy.sha256,
                "argv": list(plan.argv),
                "tool_identity": advntr_capabilities_document(replay.capabilities),
                "output_sha256": replay.output_sha256,
            }
            _write_private(execution_root / "receipt.json", canonical_json_bytes(receipt) + b"\n")
        if observed_capabilities is None:
            _fail("adVNTR cutoff grid produced no native executions")
        for policy_id, policy in policy_rows:
            execution_id, samples = outcomes[advntr_signature(policy)]
            built.append(AdvntrCutoffPolicyResult(policy_id, policy.sha256, execution_id, samples))
        document = _result_document(baseline_policy_id, context, observed_capabilities, tuple(built))
        document_holder.append(document)
        _write_private(staging / "grid-result.json", canonical_json_bytes(document) + b"\n")
        _private_tree(staging)
        return True

    atomic_output(destination, produce)
    if observed_capabilities is None or not document_holder:
        raise RuntimeError("adVNTR cutoff grid completed without a result")
    document = document_holder[0]
    return AdvntrCutoffGridResult(
        destination,
        baseline_policy_id,
        context.capture_policy.sha256,
        observed_capabilities,
        tuple(built),
        advntr_canonical_sha256(document),
    )


__all__ = [
    "AdvntrCutoffGridResult",
    "AdvntrCutoffPolicyResult",
    "AdvntrCutoffSample",
    "advntr_signature",
    "evaluate_advntr_cutoff_grid",
]
