"""Role-scoped length command orchestration over verified measurement artifacts."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import NoReturn, cast

from vntyper.scripts.calibration_artifact_io import load_object, write_checksums, write_json
from vntyper.scripts.calibration_candidate import (
    candidate_applicability_document,
    candidate_document,
    candidate_producer_document,
    decode_candidate,
)
from vntyper.scripts.calibration_development_source import (
    DevelopmentSource,
    decode_development_source,
    development_source_document,
)
from vntyper.scripts.calibration_exposure import ExposureReceipt, exposure_receipt_document
from vntyper.scripts.calibration_exposure_io import record_exposure
from vntyper.scripts.calibration_length import LengthTrainingRow, decode_length_training_metadata
from vntyper.scripts.calibration_length_artifacts import (
    build_length_training_artifact,
    length_training_artifact_document,
    validate_length_training_profile,
)
from vntyper.scripts.calibration_length_assessment import (
    assess_length_candidate,
    length_development_assessment_document,
)
from vntyper.scripts.calibration_length_evaluation import (
    EvaluationPhase,
    LengthEvaluationResult,
    LengthEvaluationRow,
    evaluate_length_hypotheses,
    length_evaluation_document,
)
from vntyper.scripts.calibration_length_evidence import (
    LengthRoleEvidence,
    bind_length_role_evidence,
    bind_length_training_evidence,
    length_role_evidence_document,
    length_training_evidence_document,
)
from vntyper.scripts.calibration_length_metrics import (
    LengthEligibleRoster,
    length_eligible_roster_document,
)
from vntyper.scripts.calibration_length_policy import (
    decode_length_source_truth,
    length_measurement_policy_sha256,
    length_protocol_for_roster,
)
from vntyper.scripts.calibration_length_profile import (
    LengthResearchProfile,
    build_length_payload,
    load_length_research_profile,
)
from vntyper.scripts.calibration_length_protocol import (
    LengthProtocol,
)
from vntyper.scripts.calibration_length_report import render_length_evaluation
from vntyper.scripts.calibration_payload import payload_manifest_document
from vntyper.scripts.calibration_role_source import RoleSource, decode_role_source, role_source_document
from vntyper.scripts.calibration_target_asset_io import read_target_json
from vntyper.scripts.calibration_target_contract import (
    LengthBaselinePlan,
    TargetStudy,
    decode_target_study,
    target_study_document,
)
from vntyper.scripts.calibration_target_runs import (
    TargetRuns,
    decode_target_runs,
    require_run_assets,
    select_target_runs,
)
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object
from vntyper.scripts.length_annotation import decode_length_annotation
from vntyper.scripts.length_estimation import EvidenceDomain
from vntyper.scripts.length_feature_artifact import decode_length_feature_artifact
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

logger = logging.getLogger(__name__)

_truth = decode_length_source_truth


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _require_exposure(source: RoleSource, study: TargetStudy, receipt: ExposureReceipt) -> None:
    exposure_receipt_document(receipt)
    expected = {
        "target": "length",
        "role": source.role,
        "study_sha256": study.sha256,
        "partition_sha256": study.partitions.sha256,
        "evidence_sha256": source.sha256,
        "exposure_ledger_id": study.exposure_ledger_id,
        "membership_sha256": canonical_sha256(
            [{"namespace": name, "sha256": digest} for name, digest in source.identities]
        ),
    }
    if any(getattr(receipt, name) != value for name, value in expected.items()):
        _fail("length extraction exposure receipt differs from the authorized source")


def load_length_source_rows(
    study: TargetStudy,
    runs: TargetRuns,
    source: RoleSource,
    receipt: ExposureReceipt,
    *,
    opened_truth: bytes | None = None,
) -> tuple[LengthTrainingRow | LengthEvaluationRow, ...]:
    """Open one authorized role's committed truth and measurement files.

    Args:
        study: Frozen length study.
        runs: Complete immutable run manifest.
        source: Metadata checked against the exact requested role.
        receipt: Actual durable exposure receipt, recorded before this call.

    Returns:
        Complete typed rows in the independently declared group order. Missing
        feature ratios remain in these rows for the engine's availability logic.

    Raises:
        ValueError: If authorization, process completion, provenance or bytes differ.
        OSError: If a committed file cannot be opened.
    """
    role_source_document(source, study=study, runs=runs)
    if not isinstance(study.baseline, LengthBaselinePlan) or not isinstance(source.roster, LengthEligibleRoster):
        _fail("length extraction received the wrong target")
    _require_exposure(source, study, receipt)
    policy = length_measurement_policy_sha256(study)
    selected = select_target_runs(runs, source.keys, (policy,))
    producer_sha = canonical_sha256(candidate_producer_document(study.baseline.producer))
    for run in selected:
        require_run_assets("length", run, ())
        if (
            run.capture_policy_sha256 != policy
            or run.baseline_assets_sha256 != study.baseline.sha256
            or run.producer_sha256 != producer_sha
        ):
            _fail("length run producer, capture policy or baseline commitment differs")
        alignment = run.assets.get("input_alignment")
        if alignment is not None and alignment.sha256 != run.input_sha256:
            _fail("length run alignment commitment differs from its input identity")
    if opened_truth is None:
        truth_raw = read_target_json(source.truth_asset)
    else:
        if (
            not isinstance(opened_truth, bytes)
            or len(opened_truth) != source.truth_asset.size_bytes
            or hashlib.sha256(opened_truth).hexdigest() != source.truth_asset.sha256
        ):
            _fail("opened length truth differs from its sealed source commitment")
        truth_raw = load_strict_json_object(opened_truth)
    truth = _truth(truth_raw, source.keys)
    measured = {}
    for run in selected:
        feature = decode_length_feature_artifact(read_target_json(run.assets["length_features"]))
        if (
            feature.manifest_key != run.manifest_key
            or feature.provenance.measurement_context.input_sha256 != run.input_sha256
            or feature.annotation_sha256 != study.baseline.annotation_sha256
            or feature.counting_policy_sha256 != study.baseline.counting_policy_sha256
        ):
            _fail("length feature body differs from its run input or measurement policy")
        measured[run.manifest_key] = feature
    rows: list[LengthTrainingRow | LengthEvaluationRow] = []
    for member in source.roster.members:
        domain = cast(EvidenceDomain, source.evidence_domains[member.key])
        if source.role == "training":
            rows.append(
                LengthTrainingRow(
                    member.key,
                    member.group_key,
                    "training",
                    measured[member.key],
                    TARGET_BOUNDARY_DEFINITION,
                    truth[member.key],
                    domain,
                )
            )
        else:
            rows.append(
                LengthEvaluationRow(
                    member.key,
                    member.group_key,
                    cast(EvaluationPhase, source.role),
                    measured[member.key],
                    TARGET_BOUNDARY_DEFINITION,
                    truth[member.key],
                    domain,
                )
            )
    return tuple(rows)


def _source(root: Path, study: TargetStudy, runs: TargetRuns, role: str) -> RoleSource:
    return decode_role_source(
        load_object(root / "roles" / role / "source.json", f"length {role} source"),
        study=study,
        runs=runs,
        expected_role=role,
    )


def _expose(ledger: Path, source: RoleSource, study: TargetStudy, evidence_root: Path, output: Path) -> ExposureReceipt:
    return record_exposure(
        ledger,
        expected_ledger_id=study.exposure_ledger_id,
        target="length",
        role=source.role,
        study_sha256=study.sha256,
        partition_sha256=study.partitions.sha256,
        evidence_sha256=source.sha256,
        identities=[{"namespace": namespace, "sha256": digest} for namespace, digest in source.identities],
        forbidden_roots=(evidence_root, output),
    )


def evaluate_fixed_length_source(
    profile: LengthResearchProfile,
    runs: TargetRuns,
    source: RoleSource,
    receipt: ExposureReceipt,
    *,
    opened_truth: bytes | None = None,
) -> tuple[LengthEvaluationResult, str]:
    """Evaluate the research profile's one fixed model on an authorized role.

    Args:
        profile: Verified local research fit and selected protocol candidate.
        runs: Exact target run commitments for the requested role.
        source: Validation or locked-heldout metadata authorized by the caller.
        receipt: Durable exposure receipt recorded before outcome access.
        opened_truth: Optional already consumed locked truth bytes.

    Returns:
        The immutable evaluation result and its deterministic offline report.

    Raises:
        ValueError: If profile, role, evidence, or fixed-model bindings differ.
    """
    if not isinstance(profile, LengthResearchProfile):
        _fail("fixed length evaluation requires a verified research profile")
    if source.role not in {"validation", "locked-heldout"}:
        _fail("fixed length evaluation requires validation or locked-heldout evidence")
    roster = cast(LengthEligibleRoster, source.roster)
    protocol = length_protocol_for_roster(cast(LengthProtocol, profile.study.protocol), roster)
    validate_length_training_profile(profile.training_profile, protocol)
    rows = cast(
        tuple[LengthEvaluationRow, ...],
        load_length_source_rows(
            profile.study,
            runs,
            source,
            receipt,
            opened_truth=opened_truth,
        ),
    )
    evidence = bind_length_role_evidence(
        rows,
        roster,
        protocol,
        study_sha256=profile.study.sha256,
        partition_sha256=profile.study.partitions.sha256,
        run_manifest_sha256=runs.sha256,
    )
    length_role_evidence_document(evidence)
    result = evaluate_length_hypotheses(
        profile.training_profile.artifact.outcomes,
        profile.training_profile.artifact.baseline,
        roster,
        protocol,
        rows,
        fixed_candidate_id=profile.selected_protocol_candidate_id,
    )
    return result, render_length_evaluation(
        result,
        roster,
        protocol,
        profile.training_profile.artifact.baseline,
    )


def _require_development_exposure(
    source: DevelopmentSource, profile: LengthResearchProfile, receipt: ExposureReceipt
) -> None:
    exposure_receipt_document(receipt)
    expected = {
        "target": "length",
        "role": "development-assessment",
        "study_sha256": profile.study.sha256,
        "partition_sha256": source.partition_sha256,
        "evidence_sha256": source.sha256,
        "exposure_ledger_id": profile.study.exposure_ledger_id,
        "membership_sha256": canonical_sha256(
            [{"namespace": name, "sha256": digest} for name, digest in source.identities]
        ),
    }
    if any(getattr(receipt, name) != value for name, value in expected.items()):
        _fail("length development exposure receipt differs from its source")


def load_length_development_evidence(
    profile: LengthResearchProfile,
    runs: TargetRuns,
    source: DevelopmentSource,
    receipt: ExposureReceipt,
) -> LengthRoleEvidence:
    """Open one previously examined development source after durable exposure.

    Args:
        profile: Verified research profile whose fixed model will be assessed.
        runs: Exact development measurement run commitments.
        source: Closed nonpromotable development source declaration.
        receipt: Durable receipt recorded before any truth or feature read.

    Returns:
        Complete typed development evidence ready for fixed-model assessment.

    Raises:
        ValueError: If source, exposure, run, feature, truth, or policy bindings differ.
    """
    development_source_document(source, candidate=profile.candidate, runs=runs)
    baseline = profile.study.baseline
    if (
        not isinstance(source.roster, LengthEligibleRoster)
        or not isinstance(baseline, LengthBaselinePlan)
        or runs.target != "length"
    ):
        _fail("length development assessment received the wrong target")
    _require_development_exposure(source, profile, receipt)
    policy = length_measurement_policy_sha256(profile.study)
    selected = select_target_runs(runs, source.keys, (policy,))
    producer_sha = canonical_sha256(candidate_producer_document(baseline.producer))
    measured = {}
    for run in selected:
        require_run_assets("length", run, ())
        if (
            run.capture_policy_sha256 != policy
            or run.baseline_assets_sha256 != baseline.sha256
            or run.producer_sha256 != producer_sha
        ):
            _fail("length development run differs from the profile measurement plan")
        alignment = run.assets.get("input_alignment")
        if alignment is not None and alignment.sha256 != run.input_sha256:
            _fail("length development run alignment differs from its input identity")
        feature = decode_length_feature_artifact(read_target_json(run.assets["length_features"]))
        if (
            feature.manifest_key != run.manifest_key
            or feature.provenance.measurement_context.input_sha256 != run.input_sha256
            or feature.annotation_sha256 != baseline.annotation_sha256
            or feature.counting_policy_sha256 != baseline.counting_policy_sha256
        ):
            _fail("length development feature differs from its run or measurement plan")
        measured[run.manifest_key] = feature
    truth = _truth(read_target_json(source.truth_asset), source.keys)
    rows = tuple(
        LengthEvaluationRow(
            member.key,
            member.group_key,
            "development-assessment",
            measured[member.key],
            TARGET_BOUNDARY_DEFINITION,
            truth[member.key],
            cast(EvidenceDomain, source.evidence_domains[member.key]),
        )
        for member in source.roster.members
    )
    protocol = length_protocol_for_roster(cast(LengthProtocol, profile.study.protocol), source.roster)
    return bind_length_role_evidence(
        rows,
        source.roster,
        protocol,
        study_sha256=profile.study.sha256,
        partition_sha256=source.partition_sha256,
        run_manifest_sha256=runs.sha256,
    )


def assess_length_bundle(args: object, output: Path) -> bool:
    """Assess one fixed research model on audited, previously examined evidence.

    Args:
        args: Namespace with profile, intake, runs, and external exposure ledger paths.
        output: Already staged output directory.

    Returns:
        True after the nonpromotable assessment and report are written.

    Raises:
        ValueError: If inputs, exposure, fixed model, or development evidence differ.
    """
    profile_root = cast(Path, getattr(args, "profile", None))
    intake_root = cast(Path, getattr(args, "intake", None))
    runs_path = cast(Path, getattr(args, "runs", None))
    ledger = cast(Path, getattr(args, "exposure_ledger", None))
    paths = (profile_root, intake_root, runs_path, ledger, output)
    if any(not isinstance(path, Path) for path in paths):
        _fail("length assessment requires Path profile, intake, runs, ledger, and output")
    profile = load_length_research_profile(profile_root)
    runs = decode_target_runs(load_object(runs_path, "length development runs"))
    if (
        not intake_root.is_dir()
        or intake_root.is_symlink()
        or {path.name for path in intake_root.iterdir()} != {"source.json"}
    ):
        _fail("length development intake inventory differs")
    source = decode_development_source(
        load_object(intake_root / "source.json", "length development source"),
        candidate=profile.candidate,
        runs=runs,
    )
    receipt = record_exposure(
        ledger,
        expected_ledger_id=profile.study.exposure_ledger_id,
        target="length",
        role="development-assessment",
        study_sha256=profile.study.sha256,
        partition_sha256=source.partition_sha256,
        evidence_sha256=source.sha256,
        identities=[{"namespace": name, "sha256": digest} for name, digest in source.identities],
        forbidden_roots=(profile_root, intake_root, runs_path.parent, output),
    )
    evidence = load_length_development_evidence(profile, runs, source, receipt)
    roster = cast(LengthEligibleRoster, source.roster)
    protocol = length_protocol_for_roster(cast(LengthProtocol, profile.study.protocol), roster)
    result = assess_length_candidate(
        profile.training_profile,
        evidence,
        roster,
        protocol,
        fixed_candidate_id=profile.selected_protocol_candidate_id,
    )
    write_json(output / "assessment.json", length_development_assessment_document(result.assessment))
    (output / "report.html").write_text(result.report_html, encoding="utf-8")
    write_checksums(output)
    return True


def fit_length_bundle(args: object, output: Path) -> bool:
    """Fit and select a research-only length candidate from exposed roles.

    Args:
        args: Namespace with evidence, exposure_ledger, and frozen objective.
        output: Already staged output directory.

    Returns:
        True when the policy-selection phase selects a fitted candidate.

    Raises:
        ValueError: If declarations, exposure, evidence, fit, or selection differ.
    """
    evidence_root = cast(Path, getattr(args, "evidence", None))
    ledger = cast(Path, getattr(args, "exposure_ledger", None))
    objective = getattr(args, "objective", None)
    if not isinstance(evidence_root, Path) or not isinstance(ledger, Path) or not isinstance(output, Path):
        _fail("length fit requires Path evidence, ledger, and output")
    study_raw = load_object(evidence_root / "study.json", "length study")
    study = decode_target_study(study_raw)
    runs = decode_target_runs(load_object(evidence_root / "runs.json", "length runs"))
    annotation_raw = load_object(evidence_root / "annotation.json", "length annotation")
    annotation = decode_length_annotation(annotation_raw)
    baseline = study.baseline
    protocol = study.protocol
    if (
        study.target != "length"
        or not isinstance(baseline, LengthBaselinePlan)
        or not isinstance(protocol, LengthProtocol)
    ):
        _fail("length fit requires a length target study")
    if objective != "length-total-v1":
        _fail("length fit objective differs from the frozen length protocol")
    if annotation.sha256 != baseline.annotation_sha256:
        _fail("length fit annotation differs from the study baseline plan")

    training_source = _source(evidence_root, study, runs, "training")
    training_receipt = _expose(ledger, training_source, study, evidence_root, output)
    training_rows = cast(
        tuple[LengthTrainingRow, ...], load_length_source_rows(study, runs, training_source, training_receipt)
    )
    training_roster = cast(LengthEligibleRoster, training_source.roster)
    training_evidence = bind_length_training_evidence(
        training_rows,
        training_roster,
        study_sha256=study.sha256,
        partition_sha256=study.partitions.sha256,
        run_manifest_sha256=runs.sha256,
    )
    metadata = decode_length_training_metadata(
        {
            "schema_version": "length-training-metadata-v1",
            "study_sha256": study.sha256,
            "training_evidence_sha256": training_evidence.sha256,
            "training_roster_sha256": training_roster.sha256,
            "annotation_sha256": baseline.annotation_sha256,
            "counting_policy_sha256": baseline.counting_policy_sha256,
            "applicability": candidate_applicability_document(study.applicability, target="length"),
            "qc": dict(protocol.qc),
            "producer": candidate_producer_document(baseline.producer),
            "maximum_condition_number": baseline.maximum_condition_number,
        },
        training_roster,
    )
    training_artifact = build_length_training_artifact(study, training_roster, metadata, training_evidence)

    selection_source = _source(evidence_root, study, runs, "policy-selection")
    selection_receipt = _expose(ledger, selection_source, study, evidence_root, output)
    selection_rows = cast(
        tuple[LengthEvaluationRow, ...], load_length_source_rows(study, runs, selection_source, selection_receipt)
    )
    selection_roster = cast(LengthEligibleRoster, selection_source.roster)
    selection_evidence = bind_length_role_evidence(
        selection_rows,
        selection_roster,
        protocol,
        study_sha256=study.sha256,
        partition_sha256=study.partitions.sha256,
        run_manifest_sha256=runs.sha256,
    )
    evaluation = evaluate_length_hypotheses(
        training_artifact.outcomes,
        training_artifact.baseline,
        selection_roster,
        protocol,
        selection_rows,
    )
    if evaluation.selection.status != "selected" or evaluation.selection.selected_candidate_id is None:
        _fail("length policy selection did not produce a research candidate")
    selected = next(
        outcome
        for outcome in training_artifact.outcomes
        if outcome.candidate_id == evaluation.selection.selected_candidate_id
    )
    if selected.model is None:
        _fail("length selected candidate has no fitted model")
    manifest, annotation_bytes, model_bytes = build_length_payload(selected.model, annotation)
    candidate_raw: dict[str, object] = {
        "schema_version": "calibration-candidate-v2",
        "target": "length",
        "study_sha256": study.sha256,
        "baseline_sha256": training_artifact.baseline.sha256,
        "partition_sha256": study.partitions.sha256,
        "training_evidence_sha256": training_evidence.sha256,
        "selection_evidence_sha256": evaluation.sha256,
        "payload_sha256": manifest.sha256,
        "applicability": candidate_applicability_document(selected.model.applicability, target="length"),
        "producer": candidate_producer_document(selected.model.producer),
        "status": "research-only",
    }
    candidate_raw["candidate_id"] = canonical_sha256(candidate_raw)
    candidate = decode_candidate(candidate_raw)
    payload_root = output / "payload"
    payload_root.mkdir()
    (payload_root / "length-annotation.json").write_bytes(annotation_bytes)
    (payload_root / "length-model.json").write_bytes(model_bytes)
    write_json(output / "candidate.json", candidate_document(candidate))
    write_json(output / "payload-manifest.json", payload_manifest_document(manifest))
    write_json(output / "study.json", target_study_document(study))
    write_json(output / "training-roster.json", length_eligible_roster_document(training_roster))
    write_json(output / "training-evidence.json", length_training_evidence_document(training_evidence))
    write_json(
        output / "training-artifact.json",
        length_training_artifact_document(training_artifact, study=study, training_roster=training_roster),
    )
    write_json(output / "selection-evidence.json", length_role_evidence_document(selection_evidence))
    write_json(output / "evaluation.json", length_evaluation_document(evaluation))
    (output / "report.html").write_text(
        render_length_evaluation(evaluation, selection_roster, protocol, training_artifact.baseline),
        encoding="utf-8",
    )
    write_checksums(output)
    return True
