"""Caller target orchestration over exposed, hash-bound production evidence."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import NoReturn, cast

from vntyper.scripts.calibration_advntr_runtime_policy import decode_advntr_runtime_policy
from vntyper.scripts.calibration_artifact_io import load_object, write_checksums, write_json
from vntyper.scripts.calibration_caller_artifacts import (
    caller_evaluation_document,
    caller_role_evidence_document,
    decode_caller_role_evidence,
    render_caller_evaluation,
)
from vntyper.scripts.calibration_caller_background_training import (
    TrainingBackground,
    fit_caller_training_background,
)
from vntyper.scripts.calibration_caller_metrics import CallerObservation
from vntyper.scripts.calibration_caller_observations import (
    CallerTruth,
    decode_advntr_baseline_calls,
    decode_advntr_replay_calls,
    decode_caller_truth,
    native_caller_observation,
    replayed_caller_observation,
)
from vntyper.scripts.calibration_caller_profile_io import (
    CallerResearchProfile,
    load_caller_research_profile,
)
from vntyper.scripts.calibration_caller_profile_io import (
    build_caller_payload as _payload,
)
from vntyper.scripts.calibration_caller_protocol import CallerProtocol
from vntyper.scripts.calibration_caller_roster import CallerEligibleRoster, caller_eligible_roster_document
from vntyper.scripts.calibration_callers import (
    CallerEvaluationResult,
    CallerEvidenceRow,
    CallerPolicyEvidence,
    CallerReplayEquivalence,
    CallerRoleEvidence,
    EvaluationPhase,
    evaluate_caller_grid,
)
from vntyper.scripts.calibration_candidate import (
    candidate_applicability_document,
    candidate_document,
    candidate_producer_document,
    decode_candidate,
    validate_candidate_payload,
)
from vntyper.scripts.calibration_development_source import DevelopmentSource, decode_development_source
from vntyper.scripts.calibration_exposure import ExposureReceipt, exposure_receipt_document
from vntyper.scripts.calibration_exposure_io import record_exposure
from vntyper.scripts.calibration_kestrel_capture import decode_kestrel_capture
from vntyper.scripts.calibration_payload import payload_manifest_document
from vntyper.scripts.calibration_role_source import RoleSource, decode_role_source, role_source_document
from vntyper.scripts.calibration_target_asset_io import read_target_asset, read_target_json
from vntyper.scripts.calibration_target_contract import (
    CallerBaselinePlan,
    TargetStudy,
    decode_target_study,
    target_study_document,
)
from vntyper.scripts.calibration_target_runs import (
    TargetRun,
    TargetRuns,
    decode_target_runs,
    require_run_assets,
    select_target_runs,
    target_runs_document,
)
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _require_exposure(
    source: RoleSource | DevelopmentSource,
    study: TargetStudy,
    receipt: ExposureReceipt,
    partition_sha256: str,
) -> None:
    exposure_receipt_document(receipt)
    expected = {
        "target": "callers",
        "role": source.role if isinstance(source, RoleSource) else source.evidence_role,
        "study_sha256": study.sha256,
        "partition_sha256": partition_sha256,
        "evidence_sha256": source.sha256,
        "exposure_ledger_id": study.exposure_ledger_id,
        "membership_sha256": canonical_sha256(
            [{"namespace": name, "sha256": digest} for name, digest in source.identities]
        ),
    }
    if any(getattr(receipt, name) != value for name, value in expected.items()):
        _fail("caller extraction exposure receipt differs from the authorized source")


def _opened_truth(source: RoleSource | DevelopmentSource, opened: bytes | None) -> CallerTruth:
    if opened is None:
        value = read_target_json(source.truth_asset)
    else:
        if (
            not isinstance(opened, bytes)
            or len(opened) != source.truth_asset.size_bytes
            or hashlib.sha256(opened).hexdigest() != source.truth_asset.sha256
        ):
            _fail("opened caller truth differs from its sealed source commitment")
        value = load_strict_json_object(opened)
    return decode_caller_truth(value, source.keys)


def _source_digest(run: TargetRun, roles: tuple[str, ...]) -> str:
    return canonical_sha256(
        {"manifest_key": run.manifest_key, "policy_sha256": run.policy_sha256,
         "assets": {role: run.assets[role].sha256 for role in sorted(roles)}}
    )


def _asset(run: TargetRun, role: str) -> bytes:
    asset = run.assets.get(role)
    if asset is None:
        _fail(f"caller run lacks required {role} asset")
    return read_target_asset(asset)


def _advntr_baseline(run: TargetRun, required_callers: tuple[str, ...]) -> tuple[tuple[dict[str, object], ...], bool, bytes | None]:
    if "advntr" not in required_callers:
        return (), True, None
    capture = _asset(run, "advntr_capture")
    calls, assessable = decode_advntr_baseline_calls(capture, run.vntr_ids)
    return calls, assessable, capture


def _advntr_replay(run: TargetRun, capture_raw: bytes) -> tuple[tuple[dict[str, object], ...], bool]:
    background = run.assets.get("advntr_background")
    return decode_advntr_replay_calls(
        _asset(run, "advntr_replay_result"),
        manifest_raw=_asset(run, "advntr_replay_manifest"),
        policy_raw=_asset(run, "advntr_replay_policy"),
        capture_raw=capture_raw,
        expected_key=run.manifest_key,
        expected_vntr_ids=run.vntr_ids,
        expected_background_sha256=None if background is None else background.sha256,
    )


def _native_row(
    run: TargetRun,
    member,
    truth: CallerTruth,
    replay_row: CallerEvidenceRow,
    required_callers: tuple[str, ...],
) -> CallerEvidenceRow:
    roles: tuple[str, ...] = ("kestrel_capture", "kestrel_result")
    if "advntr" in required_callers:
        roles += ("advntr_capture", "advntr_result", "advntr_model")
    return native_caller_observation(
        member,
        truth,
        kestrel_tsv=_asset(run, "kestrel_result"),
        advntr_tsv=_asset(run, "advntr_result") if "advntr" in required_callers else None,
        disposition=replay_row.disposition,
        source_evidence_sha256=_source_digest(run, roles),
    )


def _policy_evidence(
    protocol: CallerProtocol,
    roster: CallerEligibleRoster,
    truth: CallerTruth,
    selected: tuple[TargetRun, ...],
    policy_ids: tuple[str, ...],
) -> tuple[tuple[CallerPolicyEvidence, ...], tuple[CallerReplayEquivalence, ...]]:
    available = {protocol.baseline_policy_sha256: protocol.baseline_policy}
    available.update({candidate.candidate_id: candidate.policy for candidate in protocol.candidates})
    policies = {identity: available[identity] for identity in policy_ids}
    by_arm = {(run.manifest_key, run.policy_sha256): run for run in selected}
    policy_rows: dict[str, list[CallerEvidenceRow]] = {identity: [] for identity in policies}
    baseline_replay: dict[str, dict[str, CallerObservation]] = {}
    capture_by_policy: dict[str, set[str]] = {identity: set() for identity in policies}
    kind_by_policy: dict[str, str] = {}
    for policy_id, policy in policies.items():
        for member in roster.members:
            run = by_arm[(member.key, policy_id)]
            require_run_assets("callers", run, policy.required_callers)
            if run.policy_sha256 != policy.sha256:
                _fail("caller run policy differs from the full declared policy")
            capture = decode_kestrel_capture(read_target_json(run.assets["kestrel_capture"]))
            if capture.provenance.capture_policy_sha256 != run.capture_policy_sha256:
                _fail("caller Kestrel capture policy differs from its run")
            baseline_calls, baseline_assessable, capture_raw = _advntr_baseline(run, policy.required_callers)
            source_roles = tuple(sorted(run.assets))
            source_sha = _source_digest(run, source_roles)
            if run.execution_kind == "scalar-replay":
                calls, assessable = _advntr_replay(run, cast(bytes, capture_raw)) if capture_raw is not None else ((), True)
                row = replayed_caller_observation(
                    member, truth, kestrel_capture=capture, policy=policy,
                    capture_policy_sha256=run.capture_policy_sha256, advntr_calls=calls,
                    advntr_assessable=assessable, source_evidence_sha256=source_sha,
                )
                baseline_row = replayed_caller_observation(
                    member, truth, kestrel_capture=capture, policy=protocol.baseline_policy,
                    capture_policy_sha256=run.capture_policy_sha256, advntr_calls=baseline_calls,
                    advntr_assessable=baseline_assessable, source_evidence_sha256=source_sha,
                )
                by_member = baseline_replay.setdefault(run.capture_policy_sha256, {})
                previous = by_member.setdefault(member.key, baseline_row.observation)
                if previous != baseline_row.observation:
                    _fail("caller baseline replay differs across scalar candidates")
            else:
                replay_row = replayed_caller_observation(
                    member, truth, kestrel_capture=capture, policy=policy,
                    capture_policy_sha256=run.capture_policy_sha256, advntr_calls=baseline_calls,
                    advntr_assessable=baseline_assessable, source_evidence_sha256=source_sha,
                )
                row = _native_row(run, member, truth, replay_row, policy.required_callers)
                if row.observation != replay_row.observation:
                    _fail("caller native production and capture replay observations differ")
            policy_rows[policy_id].append(row)
            capture_by_policy[policy_id].add(run.capture_policy_sha256)
            kind_by_policy.setdefault(policy_id, run.execution_kind)
            if kind_by_policy[policy_id] != run.execution_kind:
                _fail("caller policy arm mixes execution kinds")
    raw_policies: list[dict[str, object]] = []
    for policy_id in sorted(policies):
        captures = capture_by_policy[policy_id]
        if len(captures) != 1:
            _fail("caller policy arm must use one capture policy across its role")
        raw_policies.append(
            {
                "candidate_id": policy_id,
                "policy_sha256": policy_id,
                "execution_kind": kind_by_policy[policy_id],
                "capture_policy_sha256": next(iter(captures)),
                "rows": [_evidence_row_document(row) for row in policy_rows[policy_id]],
            }
        )
    draft = decode_caller_role_evidence(
        {"schema_version": "calibration-caller-role-evidence-v1", "phase": "policy-selection",
         "protocol_sha256": protocol.sha256, "eligible_roster_sha256": roster.sha256,
         "run_manifest_sha256": "0" * 64, "baseline_assets_sha256": "0" * 64,
         "policies": raw_policies, "replay_equivalence": []}
    )
    baseline = next(policy for policy in draft.policies if policy.candidate_id == protocol.baseline_policy_sha256)
    equivalence = tuple(
        CallerReplayEquivalence(
            capture_sha, protocol.baseline_policy_sha256, baseline.sha256,
            canonical_sha256([
                {"key": row.key, "group_key": row.group_key, "truth_positive": row.truth_positive,
                 "truth_variants": None if row.truth_variants is None else list(row.truth_variants),
                 "called_positive": row.called_positive, "called_variants": list(row.called_variants),
                 "tier_a_variants": list(row.tier_a_variants)} for row in rows
            ]), tuple(rows),
        )
        for capture_sha, by_member in sorted(baseline_replay.items())
        for rows in [tuple(by_member[member.key] for member in roster.members)]
    )
    return draft.policies, equivalence


def _observation_document(row) -> dict[str, object]:
    return {
        "key": row.key,
        "group_key": row.group_key,
        "truth_positive": row.truth_positive,
        "truth_variants": None if row.truth_variants is None else list(row.truth_variants),
        "called_positive": row.called_positive,
        "called_variants": list(row.called_variants),
        "tier_a_variants": list(row.tier_a_variants),
    }


def _evidence_row_document(row: CallerEvidenceRow) -> dict[str, object]:
    return {
        **_observation_document(row.observation),
        "disposition": row.disposition,
        "source_evidence_sha256": row.source_evidence_sha256,
    }


def _policy_document(policy: CallerPolicyEvidence) -> dict[str, object]:
    return {
        "candidate_id": policy.candidate_id,
        "policy_sha256": policy.policy_sha256,
        "execution_kind": policy.execution_kind,
        "capture_policy_sha256": policy.capture_policy_sha256,
        "rows": [_evidence_row_document(row) for row in policy.rows],
    }


def _equivalence_document(item: CallerReplayEquivalence) -> dict[str, object]:
    return {
        "capture_policy_sha256": item.capture_policy_sha256,
        "baseline_policy_sha256": item.baseline_policy_sha256,
        "baseline_rerun_sha256": item.baseline_rerun_sha256,
        "baseline_replay_sha256": item.baseline_replay_sha256,
        "baseline_replay_rows": [_observation_document(row) for row in item.baseline_replay_rows],
    }


def bind_caller_role_evidence(
    phase: EvaluationPhase,
    protocol_sha256: str,
    eligible_roster_sha256: str,
    run_manifest_sha256: str,
    baseline_assets_sha256: str,
    policies: tuple[CallerPolicyEvidence, ...],
    equivalence: tuple[CallerReplayEquivalence, ...],
) -> CallerRoleEvidence:
    """Bind decoded policy observations to their exact authorized source assets."""
    return decode_caller_role_evidence(
        {
            "schema_version": "calibration-caller-role-evidence-v1",
            "phase": phase,
            "protocol_sha256": protocol_sha256,
            "eligible_roster_sha256": eligible_roster_sha256,
            "run_manifest_sha256": run_manifest_sha256,
            "baseline_assets_sha256": baseline_assets_sha256,
            "policies": [_policy_document(policy) for policy in policies],
            "replay_equivalence": [_equivalence_document(item) for item in equivalence],
        }
    )


def load_caller_source_evidence(
    study: TargetStudy,
    runs: TargetRuns,
    source: RoleSource | DevelopmentSource,
    receipt: ExposureReceipt,
    *,
    phase: EvaluationPhase,
    fixed_candidate_id: str | None = None,
    opened_truth: bytes | None = None,
) -> CallerRoleEvidence:
    """Open authorized caller outcomes only after validating durable exposure."""
    if study.target != "callers" or not isinstance(study.baseline, CallerBaselinePlan):
        _fail("caller extraction requires a caller target study")
    protocol = cast(CallerProtocol, study.protocol)
    roster = source.roster
    if not isinstance(roster, CallerEligibleRoster) or runs.target != "callers":
        _fail("caller extraction requires a caller roster and run manifest")
    partition = study.partitions.sha256 if isinstance(source, RoleSource) else source.partition_sha256
    _require_exposure(source, study, receipt, partition)
    truth = _opened_truth(source, opened_truth)
    policy_ids = [protocol.baseline_policy_sha256]
    if fixed_candidate_id is None:
        policy_ids.extend(candidate.candidate_id for candidate in protocol.candidates)
    else:
        policy_ids.append(fixed_candidate_id)
    selected = select_target_runs(runs, source.keys, tuple(sorted(set(policy_ids))))
    producer_sha = canonical_sha256(candidate_producer_document(study.baseline.producer))
    if any(
        run.baseline_assets_sha256 != study.baseline.assets_sha256 or run.producer_sha256 != producer_sha
        for run in selected
    ):
        _fail("caller run producer or baseline assets differ from the study")
    ordered_policy_ids = tuple(sorted(set(policy_ids)))
    policies, equivalence = _policy_evidence(protocol, roster, truth, selected, ordered_policy_ids)
    allowed = set(policy_ids)
    policies = tuple(policy for policy in policies if policy.candidate_id in allowed)
    capture_ids = {policy.capture_policy_sha256 for policy in policies if policy.execution_kind == "scalar-replay"}
    equivalence = tuple(item for item in equivalence if item.capture_policy_sha256 in capture_ids)
    return bind_caller_role_evidence(
        phase,
        protocol.sha256,
        roster.sha256,
        runs.sha256,
        study.baseline.assets_sha256,
        policies,
        equivalence,
    )


def _source(root: Path, study: TargetStudy, runs: TargetRuns, role: str) -> RoleSource:
    return decode_role_source(
        load_object(root / "roles" / role / "source.json", f"caller {role} source"),
        study=study,
        runs=runs,
        expected_role=role,
    )


def _expose(ledger: Path, source: RoleSource, study: TargetStudy, evidence: Path, output: Path) -> ExposureReceipt:
    return record_exposure(
        ledger,
        expected_ledger_id=study.exposure_ledger_id,
        target="callers",
        role=source.role,
        study_sha256=study.sha256,
        partition_sha256=study.partitions.sha256,
        evidence_sha256=source.sha256,
        identities=[{"namespace": name, "sha256": digest} for name, digest in source.identities],
        forbidden_roots=(evidence, output),
    )


def _require_training_background_arms(
    protocol: CallerProtocol,
    runs: TargetRuns,
    keys: tuple[str, ...],
    training_background: TrainingBackground | None,
) -> None:
    policies = (protocol.baseline_policy, *(item.policy for item in protocol.candidates))
    for policy in policies:
        if (
            "advntr" not in policy.required_callers
            or policy.values["/components/advntr/calibrated_calling/mode"] != "exact"
        ):
            continue
        if training_background is None:
            _fail("exact caller policy arm lacks a training-only fitted background")
        for run in select_target_runs(runs, keys, (policy.sha256,)):
            background = run.assets.get("advntr_background")
            if background is None or background.sha256 != training_background.background_sha256:
                _fail("exact caller policy arm differs from the training-only fitted background")


def fit_caller_bundle(args: object, output: Path) -> bool:
    """Evaluate the finite caller grid and write one research-only selected bundle."""
    evidence_root = cast(Path, getattr(args, "evidence", None))
    ledger = cast(Path, getattr(args, "exposure_ledger", None))
    if not all(isinstance(path, Path) for path in (evidence_root, ledger, output)):
        _fail("caller fit requires Path evidence, ledger, and output")
    study = decode_target_study(load_object(evidence_root / "study.json", "caller study"))
    runs = decode_target_runs(load_object(evidence_root / "runs.json", "caller runs"))
    if study.target != "callers" or not isinstance(study.baseline, CallerBaselinePlan):
        _fail("caller fit requires a caller target study")
    protocol = cast(CallerProtocol, study.protocol)
    if getattr(args, "objective", None) != "caller-safety-v1":
        _fail("caller fit objective differs from the frozen caller protocol")
    training = _source(evidence_root, study, runs, "training")
    training_receipt = _expose(ledger, training, study, evidence_root, output)
    policies = (protocol.baseline_policy, *(item.policy for item in protocol.candidates))
    needs_background = any(
        "advntr" in policy.required_callers
        and policy.values["/components/advntr/calibrated_calling/mode"] == "exact"
        for policy in policies
    )
    training_background = None
    if needs_background:
        executable = getattr(args, "advntr_executable", None)
        if not isinstance(executable, Path):
            _fail("exact caller fit requires Path advntr_executable")
        training_background = fit_caller_training_background(
            study, runs, training, training_receipt,
            argv_prefix=(str(executable),), output=(output / "training-background").resolve(),
        )
        training_evidence_sha = training_background.training_evidence_sha256
    else:
        training_evidence_sha = canonical_sha256(
            {"schema_version": "caller-training-source-v1", "source_sha256": training.sha256,
             "runs_sha256": runs.sha256}
        )
    selection = _source(evidence_root, study, runs, "policy-selection")
    receipt = _expose(ledger, selection, study, evidence_root, output)
    _require_training_background_arms(protocol, runs, selection.keys, training_background)
    evidence = load_caller_source_evidence(study, runs, selection, receipt, phase="policy-selection")
    roster = cast(CallerEligibleRoster, selection.roster)
    result = evaluate_caller_grid(protocol, roster, evidence)
    write_json(output / "study.json", target_study_document(study))
    write_json(output / "runs.json", target_runs_document(runs))
    write_json(output / "training-source.json", role_source_document(training, study=study, runs=runs))
    write_json(output / "selection-source.json", role_source_document(selection, study=study, runs=runs))
    write_json(output / "selection-roster.json", caller_eligible_roster_document(roster))
    write_json(output / "selection-evidence.json", caller_role_evidence_document(evidence))
    write_json(output / "evaluation.json", caller_evaluation_document(result))
    (output / "report.html").write_text(
        render_caller_evaluation(result, protocol=protocol, roster=roster, evidence=evidence), encoding="utf-8"
    )
    if result.selection.status != "selected" or result.selection.selected_candidate_id is None:
        write_checksums(output)
        return False
    candidate_policy = next(
        candidate.policy for candidate in protocol.candidates
        if candidate.candidate_id == result.selection.selected_candidate_id
    )
    selected_runs = select_target_runs(runs, selection.keys, (candidate_policy.sha256,))
    manifest, descriptor, payload_files, _ = _payload(
        study, candidate_policy, selection, selected_runs, training_background
    )
    candidate_raw: dict[str, object] = {
        "schema_version": "calibration-candidate-v2", "target": "callers",
        "study_sha256": study.sha256, "baseline_sha256": study.baseline.policy.sha256,
        "partition_sha256": study.partitions.sha256, "training_evidence_sha256": training_evidence_sha,
        "selection_evidence_sha256": result.sha256, "payload_sha256": manifest.sha256,
        "applicability": candidate_applicability_document(study.applicability, target="callers"),
        "producer": candidate_producer_document(study.baseline.producer), "status": "research-only",
    }
    candidate_raw["candidate_id"] = canonical_sha256(candidate_raw)
    candidate = decode_candidate(candidate_raw)
    validate_candidate_payload(candidate, manifest, expected_target="callers", caller_descriptor=descriptor)
    payload_root = output / "payload"
    payload_root.mkdir()
    for name, raw in payload_files.items():
        (payload_root / name).write_bytes(raw)
    write_json(output / "candidate.json", candidate_document(candidate))
    write_json(output / "payload-manifest.json", payload_manifest_document(manifest))
    write_checksums(output)
    return True


def evaluate_fixed_caller_source(
    profile: CallerResearchProfile,
    runs: TargetRuns,
    source: RoleSource,
    receipt: ExposureReceipt,
    *,
    opened_truth: bytes | None = None,
) -> tuple[CallerEvaluationResult, str]:
    """Evaluate one fixed caller policy on an authorized validation or locked role."""
    if not isinstance(profile, CallerResearchProfile) or source.role not in {"validation", "locked-heldout"}:
        _fail("fixed caller evaluation requires a verified profile and confirmation role")
    _require_profile_run_bindings(profile, runs, source.keys)
    evidence = load_caller_source_evidence(
        profile.study, runs, source, receipt, phase=cast(EvaluationPhase, source.role),
        fixed_candidate_id=profile.selected_protocol_candidate_id, opened_truth=opened_truth,
    )
    roster = cast(CallerEligibleRoster, source.roster)
    protocol = cast(CallerProtocol, profile.study.protocol)
    result = evaluate_caller_grid(
        protocol, roster, evidence, fixed_candidate_id=profile.selected_protocol_candidate_id
    )
    return result, render_caller_evaluation(
        result, protocol=protocol, roster=roster, evidence=evidence
    )


def _require_profile_run_bindings(
    profile: CallerResearchProfile,
    runs: TargetRuns,
    keys: tuple[str, ...],
) -> None:
    """Require fixed-role native assets to equal the selected payload commitments."""
    policy = next(
        item.policy for item in cast(CallerProtocol, profile.study.protocol).candidates
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


def assess_caller_bundle(args: object, output: Path) -> bool:
    """Assess one fixed caller profile on explicitly previously examined evidence."""
    profile_root = cast(Path, getattr(args, "profile", None))
    intake_root = cast(Path, getattr(args, "intake", None))
    runs_path = cast(Path, getattr(args, "runs", None))
    ledger = cast(Path, getattr(args, "exposure_ledger", None))
    if not all(isinstance(path, Path) for path in (profile_root, intake_root, runs_path, ledger, output)):
        _fail("caller assessment requires Path profile, intake, runs, ledger, and output")
    profile = load_caller_research_profile(profile_root)
    runs = decode_target_runs(load_object(runs_path, "caller development runs"))
    if not intake_root.is_dir() or intake_root.is_symlink() or {path.name for path in intake_root.iterdir()} != {"source.json"}:
        _fail("caller development intake inventory differs")
    source = decode_development_source(
        load_object(intake_root / "source.json", "caller development source"), candidate=profile.candidate, runs=runs
    )
    _require_profile_run_bindings(profile, runs, source.keys)
    receipt = record_exposure(
        ledger, expected_ledger_id=profile.study.exposure_ledger_id, target="callers",
        role="development-assessment", study_sha256=profile.study.sha256,
        partition_sha256=source.partition_sha256, evidence_sha256=source.sha256,
        identities=[{"namespace": name, "sha256": digest} for name, digest in source.identities],
        forbidden_roots=(profile_root, intake_root, runs_path.parent, output),
    )
    evidence = load_caller_source_evidence(
        profile.study, runs, source, receipt, phase="development-assessment",
        fixed_candidate_id=profile.selected_protocol_candidate_id,
    )
    roster = cast(CallerEligibleRoster, source.roster)
    protocol = cast(CallerProtocol, profile.study.protocol)
    result = evaluate_caller_grid(
        protocol, roster, evidence, fixed_candidate_id=profile.selected_protocol_candidate_id
    )
    write_json(output / "assessment.json", caller_evaluation_document(result))
    (output / "report.html").write_text(
        render_caller_evaluation(result, protocol=protocol, roster=roster, evidence=evidence),
        encoding="utf-8",
    )
    write_checksums(output)
    return True
