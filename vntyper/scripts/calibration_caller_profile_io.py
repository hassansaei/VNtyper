"""Verified local research-profile loading for calibrated caller candidates."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn, cast

from vntyper.modules.advntr.advntr_calibration_policy import capture_policy_for_caller, decode_capture_policy
from vntyper.scripts.calibration_advntr_runtime_policy import (
    advntr_runtime_policy_document,
    build_advntr_runtime_policy,
    decode_advntr_runtime_policy,
    validate_advntr_runtime_policy,
)
from vntyper.scripts.calibration_artifact_io import load_object, verify_checksums
from vntyper.scripts.calibration_caller_artifacts import decode_caller_evaluation, decode_caller_role_evidence
from vntyper.scripts.calibration_caller_background_training import (
    TrainingBackground,
    decode_training_background_document,
)
from vntyper.scripts.calibration_caller_observations import _capture_records
from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile
from vntyper.scripts.calibration_caller_protocol import CallerProtocol
from vntyper.scripts.calibration_caller_roster import (
    CallerEligibleRoster,
    caller_eligible_roster_document,
    decode_caller_eligible_roster,
)
from vntyper.scripts.calibration_candidate import CandidateEnvelope, decode_candidate, validate_candidate_payload
from vntyper.scripts.calibration_payload import (
    CallerBundleDescriptor,
    PayloadManifest,
    caller_bundle_descriptor_document,
    decode_caller_bundle_descriptor,
    decode_payload_manifest,
    validate_payload_observations,
)
from vntyper.scripts.calibration_portable_background import validate_portable_background
from vntyper.scripts.calibration_role_source import RoleSource, decode_role_source
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.calibration_target_asset_io import read_target_asset, read_target_json
from vntyper.scripts.calibration_target_contract import CallerBaselinePlan, TargetStudy, decode_target_study
from vntyper.scripts.calibration_target_runs import TargetRun, TargetRuns, decode_target_runs, select_target_runs
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CallerResearchProfile:
    """Verified research-only caller candidate and selected full policy."""

    candidate: CandidateEnvelope
    study: TargetStudy
    selected_protocol_candidate_id: str
    payload_manifest: PayloadManifest
    payload_files: Mapping[str, bytes]
    decision_profile: ResolvedDecisionProfile
    descriptor: CallerBundleDescriptor
    projection_sha256: str


_RESEARCH_FILES = {
    "candidate.json",
    "checksums.json",
    "evaluation.json",
    "payload",
    "payload-manifest.json",
    "report.html",
    "runs.json",
    "selection-evidence.json",
    "selection-roster.json",
    "selection-source.json",
    "study.json",
    "training-source.json",
}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _asset(run: TargetRun, role: str) -> bytes:
    asset = run.assets.get(role)
    if asset is None:
        _fail(f"caller payload run lacks required {role} asset")
    return read_target_asset(asset)


def _checked_local_json(path: Path, label: str) -> object:
    """Read any canonical JSON root while retaining strict object-member checks."""
    try:
        raw = read_regular_path(path)
        value = load_strict_json_object(b'{"value":' + raw + b"}")["value"]
        if canonical_json_bytes(value) != raw:
            _fail(f"{label} must use canonical JSON bytes")
        return value
    except (OSError, UnicodeDecodeError, ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{label} is missing, unreadable, or invalid") from error


def build_caller_payload(
    study: TargetStudy,
    selected_policy,
    selection_source: RoleSource,
    selected_runs: tuple[TargetRun, ...],
    training_background: TrainingBackground | None = None,
) -> tuple[PayloadManifest, CallerBundleDescriptor, Mapping[str, bytes], ResolvedDecisionProfile]:
    protocol = cast(CallerProtocol, study.protocol)
    profile = build_caller_generated_profile(
        selected_policy,
        dataset_manifest_hash=selection_source.sha256,
        partition_manifest_hash=study.partitions.sha256,
        seed=protocol.seed,
        generator_version=study.baseline.producer.version,
    )
    files: dict[str, bytes] = {"decision-profile.json": profile.canonical_bytes}
    background_sha: str | None = None
    if "advntr" in selected_policy.required_callers:
        versions = study.baseline.producer.tool_versions
        if not versions.get("advntr") or not versions.get("advntr_build_id"):
            _fail("selected adVNTR runtime requires package and build identities")
        model_hashes = {run.assets["advntr_model"].sha256 for run in selected_runs}
        capture_policies = set()
        native_revisions: set[str] = set()
        for run in selected_runs:
            record = _capture_records(_asset(run, "advntr_capture"), run.vntr_ids)[0]
            producer, assets = record["producer"], record["assets"]
            if not isinstance(producer, Mapping) or not isinstance(assets, Mapping):
                _fail("selected adVNTR capture provenance fields differ")
            if (
                producer.get("package_version") != versions["advntr"]
                or producer.get("build_id") != versions["advntr_build_id"]
                or assets.get("model_sha256") != run.assets["advntr_model"].sha256
            ):
                _fail("selected adVNTR capture producer or model differs from the study and run")
            revision = producer.get("source_revision")
            if not isinstance(revision, str):
                _fail("selected adVNTR capture lacks its native source revision")
            native_revisions.add(revision)
            capture_policies.add(
                capture_policy_for_caller(decode_capture_policy(record["capture_policy"]), selected_policy)
            )
        if len(model_hashes) != 1 or len(capture_policies) != 1:
            _fail("selected adVNTR runs do not share frozen model and capture assets")
        raw_backgrounds = [
            canonical_json_bytes(validate_portable_background(read_target_json(run.assets["advntr_background"])))
            for run in selected_runs
            if "advntr_background" in run.assets
        ]
        if selected_policy.values["/components/advntr/calibrated_calling/mode"] == "exact" and (
            len(raw_backgrounds) != len(selected_runs)
        ):
            _fail("selected exact adVNTR runs require the frozen portable background")
        if raw_backgrounds:
            if len(set(raw_backgrounds)) != 1 or training_background is None:
                _fail("selected adVNTR portable backgrounds differ")
            if raw_backgrounds[0] != training_background.background_bytes:
                _fail("selected adVNTR background differs from training-only fitted bytes")
            files["background.json"] = training_background.background_bytes
            background_sha = training_background.background_sha256
        if len(native_revisions) != 1:
            _fail("selected adVNTR captures do not share one native source revision")
        runtime = build_advntr_runtime_policy(
            selected_policy,
            model_sha256=next(iter(model_hashes)),
            background_sha256=background_sha,
            capture_policy_sha256=next(iter(capture_policies)).sha256,
            advntr_revision=next(iter(native_revisions)),
        )
        files["advntr-policy.json"] = canonical_json_bytes(advntr_runtime_policy_document(runtime))
    descriptor = decode_caller_bundle_descriptor(
        {
            "schema_version": "caller-bundle-v2",
            "required_callers": list(selected_policy.required_callers),
            "components": {
                "decision-profile.json": hashlib.sha256(files["decision-profile.json"]).hexdigest(),
                "advntr-policy.json": None
                if "advntr-policy.json" not in files
                else hashlib.sha256(files["advntr-policy.json"]).hexdigest(),
                "background.json": background_sha,
            },
        }
    )
    files["caller-bundle.json"] = canonical_json_bytes(caller_bundle_descriptor_document(descriptor))
    ordered = MappingProxyType({name: files[name] for name in sorted(files)})
    manifest = decode_payload_manifest(
        [
            {"path": name, "size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
            for name, raw in ordered.items()
        ]
    )
    return manifest, descriptor, ordered, profile


def _require_profile_run_bindings(
    profile: CallerResearchProfile,
    runs: TargetRuns,
    keys: tuple[str, ...],
) -> None:
    """Require fixed-role native assets to equal the selected payload commitments."""
    policy = next(
        item.policy
        for item in cast(CallerProtocol, profile.study.protocol).candidates
        if item.candidate_id == profile.selected_protocol_candidate_id
    )
    if "advntr" not in policy.required_callers:
        return
    runtime = decode_advntr_runtime_policy(load_strict_json_object(profile.payload_files["advntr-policy.json"]))
    for run in select_target_runs(runs, keys, (policy.sha256,)):
        background = run.assets.get("advntr_background")
        if (
            run.assets["advntr_model"].sha256 != runtime.model_sha256
            or run.capture_policy_sha256 != runtime.capture_policy_sha256
            or (None if background is None else background.sha256) != runtime.background_sha256
        ):
            _fail("fixed caller run assets differ from the selected runtime bundle")


def load_caller_research_profile(profile_dir: Path) -> CallerResearchProfile:
    """Load a local caller fit bundle and recheck its selected policy and payload."""
    if not isinstance(profile_dir, Path) or not profile_dir.is_dir() or profile_dir.is_symlink():
        _fail("caller research profile must be a nonsymlink directory")
    inventory = {path.name for path in profile_dir.iterdir()}
    if inventory not in (_RESEARCH_FILES, {*_RESEARCH_FILES, "training-background"}):
        _fail("caller research profile inventory differs")
    verify_checksums(profile_dir)
    study = decode_target_study(load_object(profile_dir / "study.json", "caller research study"))
    runs = decode_target_runs(load_object(profile_dir / "runs.json", "caller research runs"))
    if study.target != "callers" or not isinstance(study.baseline, CallerBaselinePlan):
        _fail("caller research profile requires a caller study")
    protocol = cast(CallerProtocol, study.protocol)
    needs_background = any(
        "advntr" in policy.required_callers and policy.values["/components/advntr/calibrated_calling/mode"] == "exact"
        for policy in (protocol.baseline_policy, *(item.policy for item in protocol.candidates))
    )
    if needs_background != ("training-background" in inventory):
        _fail("caller research profile training-background inventory differs from its protocol")
    selection_source = decode_role_source(
        load_object(profile_dir / "selection-source.json", "caller selection source"),
        study=study,
        runs=runs,
        expected_role="policy-selection",
    )
    training_source = decode_role_source(
        load_object(profile_dir / "training-source.json", "caller training source"),
        study=study,
        runs=runs,
        expected_role="training",
    )
    roster = cast(CallerEligibleRoster, selection_source.roster)
    stored_roster = decode_caller_eligible_roster(
        _checked_local_json(profile_dir / "selection-roster.json", "caller selection roster")
    )
    if stored_roster != roster or caller_eligible_roster_document(stored_roster) != caller_eligible_roster_document(
        roster
    ):
        _fail("caller selection roster differs from its role source")
    evidence = decode_caller_role_evidence(load_object(profile_dir / "selection-evidence.json", "caller evidence"))
    evaluation = decode_caller_evaluation(
        load_object(profile_dir / "evaluation.json", "caller evaluation"),
        protocol=protocol,
        roster=roster,
        evidence=evidence,
    )
    selected = evaluation.selection.selected_candidate_id
    if evaluation.selection.status != "selected" or selected is None:
        _fail("caller research profile lacks a selected protocol candidate")
    payload_root = profile_dir / "payload"
    if not payload_root.is_dir() or payload_root.is_symlink():
        _fail("caller research payload must be a nonsymlink directory")
    manifest = decode_payload_manifest(
        _checked_local_json(profile_dir / "payload-manifest.json", "caller payload manifest")
    )
    names = {item.path for item in manifest.files}
    if {path.name for path in payload_root.iterdir()} != names:
        _fail("caller research payload inventory differs")
    payload_files = MappingProxyType({name: read_regular_path(payload_root / name) for name in sorted(names)})
    validate_payload_observations(
        manifest, {name: (len(raw), hashlib.sha256(raw).hexdigest()) for name, raw in payload_files.items()}
    )
    descriptor = decode_caller_bundle_descriptor(load_strict_json_object(payload_files["caller-bundle.json"]))
    candidate = decode_candidate(load_object(profile_dir / "candidate.json", "caller candidate"))
    validate_candidate_payload(candidate, manifest, expected_target="callers", caller_descriptor=descriptor)
    decision = parse_decision_profile(
        payload_files["decision-profile.json"],
        packaged_document=load_packaged_decision_profile().document,
        allow_caller_generated=True,
    )
    policy = next(item.policy for item in protocol.candidates if item.candidate_id == selected)
    training_background = None
    if needs_background:
        background_root = profile_dir / "training-background"
        if (
            not background_root.is_dir()
            or background_root.is_symlink()
            or {path.name for path in background_root.iterdir()}
            != {"capture", "native", "portable-background.json", "training-background.json"}
        ):
            _fail("caller training background inventory differs")
        for directory in (background_root / "capture", background_root / "native"):
            if not directory.is_dir() or directory.is_symlink():
                _fail("caller training background directories must be nonsymlink directories")
        background_bytes = read_regular_path(background_root / "portable-background.json")
        training_background = decode_training_background_document(
            load_object(background_root / "training-background.json", "caller training background"),
            background_bytes,
        )
        if (
            candidate.training_evidence_sha256 != training_background.training_evidence_sha256
            or training_background.study_sha256 != study.sha256
            or training_background.run_manifest_sha256 != runs.sha256
            or training_background.source_sha256 != training_source.sha256
            or payload_files.get("background.json") not in (None, background_bytes)
        ):
            _fail("caller candidate differs from its training-only background evidence")
        native_root = background_root / "native"
        if (
            not native_root.is_dir()
            or native_root.is_symlink()
            or {path.name for path in native_root.iterdir()} != set(training_background.fitter_artifact_sha256)
        ):
            _fail("caller training background native artifact inventory differs")
        for name, expected in training_background.fitter_artifact_sha256.items():
            if hashlib.sha256(read_regular_path(native_root / name)).hexdigest() != expected:
                _fail("caller training background native artifact digest differs")
    else:
        expected_training = canonical_sha256(
            {
                "schema_version": "caller-training-source-v1",
                "source_sha256": training_source.sha256,
                "runs_sha256": runs.sha256,
            }
        )
        if candidate.training_evidence_sha256 != expected_training:
            _fail("caller candidate differs from its training source evidence")
    selected_runs = select_target_runs(runs, selection_source.keys, (policy.sha256,))
    expected_manifest, expected_descriptor, expected_files, expected_profile = build_caller_payload(
        study, policy, selection_source, selected_runs, training_background
    )
    if (
        manifest != expected_manifest
        or descriptor != expected_descriptor
        or dict(payload_files) != dict(expected_files)
        or decision != expected_profile
    ):
        _fail("caller payload differs from its selected native run evidence")
    expected_profile = build_caller_generated_profile(
        policy,
        dataset_manifest_hash=selection_source.sha256,
        partition_manifest_hash=study.partitions.sha256,
        seed=protocol.seed,
        generator_version=study.baseline.producer.version,
    )
    if payload_files["decision-profile.json"] != expected_profile.canonical_bytes:
        _fail("caller generated profile differs from the selected full policy")
    background_raw = payload_files.get("background.json")
    if background_raw is not None:
        portable = validate_portable_background(load_strict_json_object(background_raw))
        if canonical_json_bytes(portable) != background_raw:
            _fail("caller portable background is not canonical")
    if "advntr" in policy.required_callers:
        runtime = decode_advntr_runtime_policy(load_strict_json_object(payload_files["advntr-policy.json"]))
        validate_advntr_runtime_policy(
            runtime,
            policy,
            background_raw_sha256=None if background_raw is None else hashlib.sha256(background_raw).hexdigest(),
        )
    if (
        candidate.study_sha256 != study.sha256
        or candidate.baseline_sha256 != study.baseline.policy.sha256
        or candidate.partition_sha256 != study.partitions.sha256
        or candidate.selection_evidence_sha256 != evaluation.sha256
        or descriptor.decision_profile_sha256 != decision.digest
        or decision.components["kestrel"] is None
        or policy.required_callers != candidate.applicability.required_callers
    ):
        _fail("caller research candidate, policy, profile, or study bindings differ")
    projection = canonical_sha256(
        {
            "schema_version": "caller-research-profile-projection-v1",
            "candidate_sha256": candidate.sha256,
            "study_sha256": study.sha256,
            "selected_protocol_candidate_id": selected,
            "payload_sha256": manifest.sha256,
            "selection_evidence_sha256": evaluation.sha256,
        }
    )
    return CallerResearchProfile(candidate, study, selected, manifest, payload_files, decision, descriptor, projection)
