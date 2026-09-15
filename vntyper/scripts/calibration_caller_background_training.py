"""Training-only staging and native fitting for calibrated adVNTR backgrounds."""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn, cast

from vntyper.modules.advntr.advntr_background import BackgroundFitResult, build_background_fit_plan, fit_background
from vntyper.modules.advntr.advntr_calibration_policy import AdvntrToolPin, advntr_canonical_sha256
from vntyper.scripts.calibration_artifact_io import write_json
from vntyper.scripts.calibration_caller_observations import _capture_records, decode_caller_truth
from vntyper.scripts.calibration_caller_protocol import CallerProtocol
from vntyper.scripts.calibration_caller_roster import CallerEligibleRoster
from vntyper.scripts.calibration_exposure import ExposureReceipt, exposure_receipt_document, require_digest
from vntyper.scripts.calibration_portable_background import project_portable_background
from vntyper.scripts.calibration_role_source import RoleSource, role_source_document
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.calibration_target_asset_io import read_target_asset, read_target_json
from vntyper.scripts.calibration_target_contract import CallerBaselinePlan, TargetStudy, target_study_document
from vntyper.scripts.calibration_target_runs import TargetRuns, select_target_runs, target_runs_document
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class TrainingBackground:
    """Portable background plus exact training-only native-fit provenance."""

    background_bytes: bytes
    background_sha256: str
    training_evidence_sha256: str
    study_sha256: str
    run_manifest_sha256: str
    source_sha256: str
    exposure_receipt_sha256: str
    negative_keys: tuple[str, ...]
    diagnostic_positive_keys: tuple[str, ...]
    excluded_unknown_keys: tuple[str, ...]
    diagnostic_policy_sha256: str
    fitter_artifact_sha256: MappingProxyType[str, str]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _payload(result: TrainingBackground) -> dict[str, object]:
    return {
        "schema_version": "calibration-caller-training-background-v1",
        "study_sha256": result.study_sha256,
        "run_manifest_sha256": result.run_manifest_sha256,
        "source_sha256": result.source_sha256,
        "exposure_receipt_sha256": result.exposure_receipt_sha256,
        "negative_keys": list(result.negative_keys),
        "diagnostic_positive_keys": list(result.diagnostic_positive_keys),
        "excluded_unknown_keys": list(result.excluded_unknown_keys),
        "background_sha256": result.background_sha256,
        "diagnostic_policy_sha256": result.diagnostic_policy_sha256,
        "fitter_artifact_sha256": dict(result.fitter_artifact_sha256),
        "training_evidence_sha256": result.training_evidence_sha256,
    }


def training_background_document(result: TrainingBackground) -> dict[str, object]:
    """Project a training-background receipt after revalidating bytes and digests."""
    if not isinstance(result, TrainingBackground):
        _fail("caller training background projection requires TrainingBackground")
    decoded = decode_training_background_document({**_payload(result), "sha256": result.sha256}, result.background_bytes)
    if decoded != result:
        _fail("caller training background receipt differs from its canonical content")
    return {**_payload(result), "sha256": result.sha256}


def decode_training_background_document(value: object, background_bytes: bytes) -> TrainingBackground:
    """Decode a training receipt against the exact portable background bytes."""
    fields = {
        "schema_version", "study_sha256", "run_manifest_sha256", "source_sha256", "exposure_receipt_sha256",
        "negative_keys", "diagnostic_positive_keys", "excluded_unknown_keys", "background_sha256",
        "diagnostic_policy_sha256", "fitter_artifact_sha256", "training_evidence_sha256", "sha256",
    }
    if not isinstance(value, Mapping) or set(value) != fields or value["schema_version"] != "calibration-caller-training-background-v1":
        _fail("caller training background fields or schema differ")
    if not isinstance(background_bytes, bytes) or hashlib.sha256(background_bytes).hexdigest() != value["background_sha256"]:
        _fail("caller training background bytes differ from their digest")
    key_fields = ("negative_keys", "diagnostic_positive_keys", "excluded_unknown_keys")
    keys: dict[str, tuple[str, ...]] = {}
    for field in key_fields:
        raw = value[field]
        if not isinstance(raw, list) or any(not isinstance(item, str) or not item for item in raw):
            _fail("caller training background membership must be text lists")
        if raw != sorted(set(raw)):
            _fail("caller training background membership must be sorted and unique")
        keys[field] = tuple(raw)
    key_sets = [set(keys[field]) for field in key_fields]
    if any(key_sets[left] & key_sets[right] for left in range(3) for right in range(left + 1, 3)) or not key_sets[0]:
        _fail("caller training background requires disjoint truth-negative controls")
    artifacts = value["fitter_artifact_sha256"]
    if not isinstance(artifacts, Mapping) or not artifacts:
        _fail("caller training background fitter artifacts must be a nonempty object")
    if any(not isinstance(name, str) or not name for name in artifacts):
        _fail("caller training background fitter artifact names are invalid")
    checked_artifacts = MappingProxyType({
        name: require_digest(artifacts[name], "fitter artifact") for name in sorted(artifacts)
    })
    result = TrainingBackground(
        background_bytes, require_digest(value["background_sha256"], "background"),
        require_digest(value["training_evidence_sha256"], "training evidence"),
        require_digest(value["study_sha256"], "study"), require_digest(value["run_manifest_sha256"], "runs"),
        require_digest(value["source_sha256"], "source"),
        require_digest(value["exposure_receipt_sha256"], "receipt"), keys["negative_keys"],
        keys["diagnostic_positive_keys"], keys["excluded_unknown_keys"],
        require_digest(value["diagnostic_policy_sha256"], "diagnostic policy"), checked_artifacts,
        require_digest(value["sha256"], "training background receipt"),
    )
    payload = _payload(result)
    expected_training = canonical_sha256({key: value for key, value in payload.items() if key != "training_evidence_sha256"})
    expected_sha = canonical_sha256(payload)
    if result.training_evidence_sha256 != expected_training or result.sha256 != expected_sha:
        _fail("caller training background receipt differs from its canonical content")
    return result


def _require_context(study: TargetStudy, runs: TargetRuns, source: RoleSource, receipt: ExposureReceipt) -> CallerProtocol:
    target_study_document(study)
    target_runs_document(runs)
    role_source_document(source, study=study, runs=runs)
    exposure_receipt_document(receipt)
    if study.target != "callers" or not isinstance(study.baseline, CallerBaselinePlan):
        _fail("caller background fitting requires a caller study")
    if source.role != "training" or source.study_sha256 != study.sha256 or source.run_manifest_sha256 != runs.sha256:
        _fail("caller background fitting requires the exact training source")
    expected = ("callers", "training", study.sha256, study.partitions.sha256, source.sha256, study.exposure_ledger_id)
    observed = (receipt.target, receipt.role, receipt.study_sha256, receipt.partition_sha256,
                receipt.evidence_sha256, receipt.exposure_ledger_id)
    if observed != expected:
        _fail("caller background fitting exposure receipt differs from the training source")
    membership = canonical_sha256(
        [{"namespace": name, "sha256": digest} for name, digest in source.identities]
    )
    if receipt.membership_sha256 != membership:
        _fail("caller background fitting receipt membership differs from the training source")
    protocol = cast(CallerProtocol, study.protocol)
    policies = (protocol.baseline_policy, *(item.policy for item in protocol.candidates))
    if not any(
        "advntr" in policy.required_callers
        and policy.values["/components/advntr/calibrated_calling/mode"] == "exact"
        for policy in policies
    ):
        _fail("caller background fitting requires a predeclared exact adVNTR candidate")
    return protocol


def _private_write(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_bytes(raw)
    path.chmod(0o600)


def _stage_training(
    root: Path, study: TargetStudy, source: RoleSource, runs: TargetRuns, protocol: CallerProtocol
) -> tuple[Path, Path, tuple[str, ...], tuple[str, ...], tuple[str, ...], AdvntrToolPin]:
    truth = decode_caller_truth(read_target_json(source.truth_asset), source.keys)
    selected = select_target_runs(runs, source.keys, (protocol.baseline_policy_sha256,))
    by_key = {run.manifest_key: run for run in selected}
    labels: list[dict[str, object]] = []
    negative: list[str] = []
    positive: list[str] = []
    unknown: list[str] = []
    producer: tuple[str, str, str] | None = None
    negative_groups: set[str] = set()
    for member in cast(CallerEligibleRoster, source.roster).members:
        row = truth.by_key[member.key]
        if row.genotype == "unknown":
            unknown.append(member.key)
            continue
        run = by_key[member.key]
        capture_raw = read_target_asset(run.assets["advntr_capture"])
        records = _capture_records(capture_raw, run.vntr_ids)
        first = records[0]
        raw_producer = first["producer"]
        raw_assets = first["assets"]
        raw_policy = first["capture_policy"]
        if not isinstance(raw_producer, Mapping) or not isinstance(raw_assets, Mapping):
            _fail("caller training capture provenance fields differ")
        current = (str(raw_producer.get("package_version")), str(raw_producer.get("build_id")),
                   str(raw_producer.get("source_revision")))
        if producer is None:
            producer = current
        if producer != current or any(record.get("producer") != raw_producer for record in records):
            _fail("caller training captures do not share one native producer")
        if (
            current[:2] != (
                study.baseline.producer.tool_versions.get("advntr"),
                study.baseline.producer.tool_versions.get("advntr_build_id"),
            )
            or raw_assets.get("model_sha256") != run.assets["advntr_model"].sha256
            or advntr_canonical_sha256(raw_policy) != run.capture_policy_sha256
        ):
            _fail("caller training capture producer, model, or policy differs from its run")
        if row.genotype == "negative" and member.group_key in negative_groups:
            _fail("caller background fitting requires one primary negative per independent group")
        if row.genotype == "negative":
            negative_groups.add(member.group_key)
        sample_id = f"member-{canonical_sha256(member.key)[:24]}"
        sink = root / "runs" / sample_id / "output" / "calibration.jsonl"
        _private_write(sink, capture_raw)
        labels.append(
            {"sample_id": sample_id, "truth": row.genotype == "positive", "partition": "training",
             "pair_id": member.group_key, "variant_class": row.genotype, "array_length": None}
        )
        (positive if row.genotype == "positive" else negative).append(member.key)
    if not negative:
        _fail("caller background fitting requires at least one truth-negative training control")
    labels_path = root / "labels.json"
    _private_write(labels_path, canonical_json_bytes({"samples": labels}))
    diagnostic = {
        "schema_version": "advntr-frameshift-policy-v1", "mode": "exact",
        "cutoff": 0.001, "minimum_read_support": 3,
    }
    diagnostic_path = root / "diagnostic-policy.json"
    _private_write(diagnostic_path, canonical_json_bytes(diagnostic))
    if producer is None:
        _fail("caller background fitting has no assessable training captures")
    return (
        labels_path, diagnostic_path, tuple(negative), tuple(positive), tuple(unknown),
        AdvntrToolPin(*producer),
    )


def _result(
    study: TargetStudy, runs: TargetRuns, source: RoleSource, receipt: ExposureReceipt,
    raw: bytes, native: BackgroundFitResult, negative: tuple[str, ...], positive: tuple[str, ...], unknown: tuple[str, ...],
) -> TrainingBackground:
    portable = raw
    background_sha = hashlib.sha256(portable).hexdigest()
    seed = TrainingBackground(
        portable, background_sha, "", study.sha256, runs.sha256, source.sha256, receipt.sha256,
        negative, positive, unknown, native.diagnostic_policy_sha256, native.artifact_sha256, "",
    )
    payload = _payload(seed)
    training_sha = canonical_sha256({key: value for key, value in payload.items() if key != "training_evidence_sha256"})
    seeded = replace(seed, training_evidence_sha256=training_sha)
    return replace(seeded, sha256=canonical_sha256(_payload(seeded)))


def fit_caller_training_background(
    study: TargetStudy,
    runs: TargetRuns,
    source: RoleSource,
    receipt: ExposureReceipt,
    *,
    argv_prefix: tuple[str, ...],
    output: Path,
) -> TrainingBackground:
    """Fit and portable-project one background using authorized training evidence only."""
    protocol = _require_context(study, runs, source, receipt)
    if not isinstance(output, Path) or not output.is_absolute() or os.path.lexists(output):
        _fail("caller background fit output must be a new absolute Path")
    try:
        output.mkdir(mode=0o700)
        capture_root = output / "capture"
        capture_root.mkdir(mode=0o700)
        labels, diagnostic, negative, positive, unknown, pin = _stage_training(
            capture_root, study, source, runs, protocol
        )
        plan = build_background_fit_plan(
            argv_prefix, capture_root=capture_root, labels_path=labels, diagnostic_policy_path=diagnostic,
            output_directory=output / "native", partition="training", profile="caller-training",
            folds=protocol.fold_count, insert_lengths=8, source_cohort=source.sha256,
            design="vntyper-caller-training-v1",
        )
        native = fit_background(plan, pin)
        native_raw = read_regular_path(plan.output_directory / "caller-training.background.json")
        portable = canonical_json_bytes(project_portable_background(load_strict_json_object(native_raw)))
        _private_write(output / "portable-background.json", portable)
        result = _result(study, runs, source, receipt, portable, native, negative, positive, unknown)
        write_json(output / "training-background.json", training_background_document(result))
        return result
    except BaseException:
        shutil.rmtree(output, ignore_errors=True)
        raise
