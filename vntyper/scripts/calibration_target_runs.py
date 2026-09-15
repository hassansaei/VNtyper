"""Local target-specific run commitments without opening scientific outcomes.

An artifact can have several predeclared policy arms. Their input identity stays
fixed while capture-policy identity distinguishes recapture from scalar replay.
Paths belong only to this local manifest, never to a portable model projection.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_exposure import require_digest
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_DIGESTS = ("policy_sha256", "input_sha256", "producer_sha256", "baseline_assets_sha256", "capture_policy_sha256")
_RUN_FIELDS = {"manifest_key", "execution_kind", "exit_code", "assets", *_DIGESTS}
_ASSET_FIELDS = {"path", "sha256", "size_bytes"}
_CALLER_ASSETS = {
    "input_alignment",
    "kestrel_capture",
    "kestrel_result",
    "advntr_capture",
    "advntr_result",
    "advntr_model",
    "advntr_background",
    "advntr_replay_manifest",
    "advntr_replay_policy",
    "advntr_replay_result",
}
_LENGTH_ASSETS = {"length_features", "input_alignment"}


@dataclass(frozen=True)
class TargetRunAsset:
    """Exact bytes expected at one local regular file."""

    path: Path
    sha256: str
    size_bytes: int


@dataclass(frozen=True)
class TargetRun:
    """One artifact/policy arm, whose process outcome is separate from genotype."""

    manifest_key: str
    policy_sha256: str
    input_sha256: str
    producer_sha256: str
    baseline_assets_sha256: str
    capture_policy_sha256: str
    execution_kind: str
    exit_code: int
    assets: Mapping[str, TargetRunAsset]
    vntr_ids: tuple[int, ...]


@dataclass(frozen=True)
class TargetRuns:
    """Complete local commitments; a role selector chooses what may be opened."""

    target: str
    runs: tuple[TargetRun, ...]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"target runs {label} fields differ from the closed contract")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip() or "\x00" in value:
        _fail(f"target runs {label} must be nonempty trimmed text")
    return value


def _integer(value: object, maximum: int, label: str) -> int:
    if type(value) is not int or not 0 <= value <= maximum:
        _fail(f"target runs {label} must be a bounded nonnegative integer")
    return value


def _asset(value: object) -> TargetRunAsset:
    raw = _object(value, _ASSET_FIELDS, "asset")
    text = _text(raw["path"], "asset path")
    path = Path(text)
    if not path.is_absolute() or ".." in path.parts or str(path) != text:
        _fail("target runs asset path must be an absolute normalized file path")
    return TargetRunAsset(
        path, require_digest(raw["sha256"], "asset digest"), _integer(raw["size_bytes"], 2**53 - 1, "size")
    )


def _run(value: object, target: str) -> TargetRun:
    fields = _RUN_FIELDS | ({"vntr_ids"} if target == "callers" else set())
    raw = _object(value, fields, "run")
    execution = _text(raw["execution_kind"], "execution kind")
    allowed = {"measurement"} if target == "length" else {"baseline-rerun", "scalar-replay", "recapture"}
    if execution not in allowed:
        _fail("target runs execution kind differs from its target")
    assets = raw["assets"]
    allowed_assets = _LENGTH_ASSETS if target == "length" else _CALLER_ASSETS
    if not isinstance(assets, Mapping) or not assets or set(assets) - allowed_assets:
        _fail("target runs asset roles are empty or differ from their target")
    decoded = {}
    for name, asset in assets.items():
        if not isinstance(name, str):
            _fail("target runs asset role must be text")
        decoded[name] = _asset(asset)
    ids = raw.get("vntr_ids", [])
    if not isinstance(ids, list) or any(type(value) is not int or not 1 <= value <= 2**53 - 1 for value in ids):
        _fail("target runs VNTR roster must contain positive integer identifiers")
    if ids != sorted(set(ids)):
        _fail("target runs VNTR roster must be sorted and unique")
    return TargetRun(
        manifest_key=_text(raw["manifest_key"], "artifact key"),
        policy_sha256=require_digest(raw["policy_sha256"], "policy digest"),
        input_sha256=require_digest(raw["input_sha256"], "input digest"),
        producer_sha256=require_digest(raw["producer_sha256"], "producer digest"),
        baseline_assets_sha256=require_digest(raw["baseline_assets_sha256"], "baseline assets digest"),
        capture_policy_sha256=require_digest(raw["capture_policy_sha256"], "capture policy digest"),
        execution_kind=execution,
        exit_code=_integer(raw["exit_code"], 255, "process exit code"),
        assets=MappingProxyType(decoded),
        vntr_ids=tuple(ids),
    )


def _run_document(run: TargetRun, target: str) -> dict[str, object]:
    return {
        "manifest_key": run.manifest_key,
        **{name: getattr(run, name) for name in _DIGESTS},
        "execution_kind": run.execution_kind,
        "exit_code": run.exit_code,
        "assets": {
            name: {"path": str(asset.path), "sha256": asset.sha256, "size_bytes": asset.size_bytes}
            for name, asset in run.assets.items()
        },
        **({"vntr_ids": list(run.vntr_ids)} if target == "callers" else {}),
    }


def decode_target_runs(value: object) -> TargetRuns:
    """Decode run metadata without reading paths, captures, truth or results.

    Args:
        value: Closed local ``calibration-runs-v2`` document.

    Returns:
        Immutable metadata retaining every artifact/policy arm.

    Raises:
        ValueError: For malformed, duplicate, unordered or inconsistent commitments.
    """
    raw = _object(value, {"schema_version", "target", "runs"}, "manifest")
    if raw["schema_version"] != "calibration-runs-v2" or raw["target"] not in ("callers", "length"):
        _fail("target runs schema or target is unsupported")
    target = _text(raw["target"], "target")
    values = raw["runs"]
    if not isinstance(values, list) or not values:
        _fail("target runs requires a nonempty run roster")
    runs = tuple(_run(row, target) for row in values)
    keys = [(run.manifest_key, run.policy_sha256) for run in runs]
    if keys != sorted(set(keys)):
        _fail("target runs artifact/policy pairs must be sorted and unique")
    inputs: dict[str, str] = {}
    for run in runs:
        if run.manifest_key in inputs and inputs[run.manifest_key] != run.input_sha256:
            _fail("target runs policy arms disagree about the artifact input bytes")
        inputs[run.manifest_key] = run.input_sha256
    document: dict[str, object] = {
        "schema_version": "calibration-runs-v2",
        "target": target,
        "runs": [_run_document(row, target) for row in runs],
    }
    return TargetRuns(target, runs, canonical_sha256(document))


def target_runs_document(manifest: TargetRuns) -> dict[str, object]:
    """Revalidate a typed local manifest before projecting its commitments.

    Args:
        manifest: Previously decoded manifest.

    Returns:
        JSON-compatible local document including paths.

    Raises:
        ValueError: If typed content or its canonical digest has changed.
    """
    if not isinstance(manifest, TargetRuns):
        _fail("target runs projection requires a typed manifest")
    document: dict[str, object] = {
        "schema_version": "calibration-runs-v2",
        "target": manifest.target,
        "runs": [_run_document(run, manifest.target) for run in manifest.runs],
    }
    if decode_target_runs(document).sha256 != manifest.sha256:
        _fail("target runs typed manifest digest differs")
    return document


def require_run_assets(target: str, run: TargetRun, required_callers: tuple[str, ...]) -> None:
    """Require completed evidence for the study's declared target and caller set.

    Args:
        target: The study target.
        run: One decoded run to extract.
        required_callers: Exact study caller set; empty for length.

    Raises:
        ValueError: For failed processes, missing evidence, or wrong-target assets.
    """
    _run(_run_document(run, target), target)
    if run.exit_code != 0:
        _fail("target run did not complete successfully; process failure is not a negative call")
    if target == "length":
        if required_callers or "length_features" not in run.assets:
            _fail("length extraction requires length features and no caller set")
        return
    if target != "callers" or required_callers not in (("kestrel",), ("advntr", "kestrel")):
        _fail("caller extraction requires the exact declared sorted caller set")
    required = {"kestrel_capture", "kestrel_result"}
    advntr_assets = {name for name in run.assets if name.startswith("advntr_")}
    if "advntr" in required_callers:
        required |= {"advntr_capture", "advntr_result", "advntr_model"}
        if run.execution_kind == "scalar-replay":
            required |= {"advntr_replay_manifest", "advntr_replay_policy", "advntr_replay_result"}
        if not run.vntr_ids:
            _fail("adVNTR extraction requires an explicit target roster")
    elif advntr_assets or run.vntr_ids:
        _fail("Kestrel-only extraction cannot acquire adVNTR evidence")
    if required - set(run.assets):
        _fail("target run is missing required caller evidence assets")


def select_target_runs(
    manifest: TargetRuns, artifact_keys: tuple[str, ...], policy_sha256s: tuple[str, ...]
) -> tuple[TargetRun, ...]:
    """Select the exact requested Cartesian roster using metadata alone.

    Args:
        manifest: Complete run commitments, possibly spanning several roles.
        artifact_keys: Authorized role artifact keys, sorted and unique.
        policy_sha256s: Predeclared policies to execute, sorted and unique.

    Returns:
        Exactly the requested artifact/policy runs, in canonical order.

    Raises:
        ValueError: For empty, duplicate, unordered or missing requests.
    """
    target_runs_document(manifest)
    if not artifact_keys or artifact_keys != tuple(sorted(set(artifact_keys))):
        _fail("target run role artifact keys must be nonempty sorted and unique")
    for key in artifact_keys:
        _text(key, "authorized artifact key")
    if not policy_sha256s or policy_sha256s != tuple(sorted(set(policy_sha256s))):
        _fail("target run role policies must be nonempty sorted and unique")
    for digest in policy_sha256s:
        require_digest(digest, "authorized policy digest")
    expected = {(key, policy) for key in artifact_keys for policy in policy_sha256s}
    selected = tuple(run for run in manifest.runs if (run.manifest_key, run.policy_sha256) in expected)
    if len(selected) != len(expected):
        _fail("target runs is missing an authorized artifact/policy arm")
    return selected
