"""Closed filesystem artifacts for the four calibration CLI operations."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from types import MappingProxyType
from typing import cast

from vntyper.scripts.calibration_artifact_io import (
    freeze_json,
    load_object,
    thaw_json,
    verify_checksums,
    write_checksums,
    write_json,
)
from vntyper.scripts.calibration_baseline import project_baseline
from vntyper.scripts.calibration_contract import CandidateMetrics, decode_evidence_manifest, decode_metrics
from vntyper.scripts.calibration_features import decode_feature_artifact, decode_label_artifact
from vntyper.scripts.calibration_manifest import connected_leakage_groups, decode_study_declaration
from vntyper.scripts.calibration_objective import (
    CandidateEvaluation,
    OutcomeObservation,
    calculate_metrics,
    count_free_parameters,
    select_candidate,
)
from vntyper.scripts.calibration_reporting import write_evaluation_artifacts
from vntyper.scripts.calibration_role_inputs import RoleInputs, build_role_inputs
from vntyper.scripts.calibration_run_extraction import decode_run_hashes
from vntyper.scripts.calibration_scalar_replay import replay_scalar_dominance
from vntyper.scripts.calibration_statistics import PairedObservation, paired_group_bootstrap
from vntyper.scripts.calibration_workflow import (
    ExtractedEvidence,
    evaluate_locked_candidate,
    extract_evidence,
    fit_candidate,
    validate_candidate,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import ResolvedDecisionProfile, resolve_decision_profile

_ROLES = ("training", "policy-selection", "validation", "locked-heldout")


def extract_artifact_bundle(truth_path: Path, study_path: Path, runs_path: Path, output: Path) -> None:
    """Validate raw study inputs and write role-separated immutable evidence."""
    truth = load_object(truth_path, "calibration truth")
    if set(truth) != {"schema_version", "labels"}:
        raise ValueError("calibration truth fields differ from the closed contract")
    if truth["schema_version"] != "calibration-truth-v1":
        raise ValueError("calibration truth schema version is unsupported")
    study_raw = load_object(study_path, "calibration study")
    study = decode_study_declaration(study_raw)
    labels = decode_label_artifact(truth["labels"])
    runs_raw = load_object(runs_path, "calibration runs")
    if set(runs_raw) != {"schema_version", "runs"} or runs_raw["schema_version"] != "calibration-runs-v1":
        raise ValueError("calibration runs fields or schema version differ")
    runs = runs_raw["runs"]
    if not isinstance(runs, Mapping) or any(not isinstance(key, str) for key in runs):
        raise ValueError("calibration runs must map manifest keys to immutable run declarations")
    extracted = extract_evidence(study, labels, runs)

    write_json(output / "study.json", study_raw)
    write_json(
        output / "groups.json",
        {
            "schema_version": "calibration-groups-v1",
            "connected_groups": dict(connected_leakage_groups(study.partitions)),
        },
    )
    feature_rows = [
        {
            "feature_key": row.feature_key,
            "manifest_key": row.manifest_key,
            "features": dict(row.features),
        }
        for row in extracted.features.rows
    ]
    label_rows = _rows(truth["labels"], "calibration labels")
    baseline = _mapping(thaw_json(extracted.baseline), "calibration baseline")
    members_by_role = {
        role: tuple(member.key for member in study.partitions.members if member.role == role) for role in _ROLES
    }
    for role, keys in members_by_role.items():
        role_inputs = build_role_inputs(keys, feature_rows, label_rows, baseline, extracted.run_hashes)
        role_root = output / "roles" / role
        role_root.mkdir(parents=True)
        write_json(role_root / "features.json", role_inputs.features)
        write_json(role_root / "labels.json", role_inputs.labels)
        write_json(role_root / "baseline.json", role_inputs.baseline)
        write_json(
            role_root / "run_hashes.json",
            {"schema_version": "calibration-run-hashes-v1", "run_hashes": role_inputs.run_hashes},
        )
        write_json(
            role_root / "manifest.json",
            _role_manifest(study, role, role_inputs, extracted.study_sha256),
        )
        if role == "locked-heldout":
            _write_locked_payload(role_root, study_raw, study, role_inputs)
        write_checksums(role_root)
    write_checksums(output)


def fit_artifact_bundle(evidence_path: Path, objective: str, output: Path) -> None:
    """Fit the frozen candidate grid from training and policy-selection only."""
    evidence = _load_roles(evidence_path, ("training", "policy-selection"))
    if objective != evidence.study.protocol.objective:
        raise ValueError("calibration fit objective differs from the snapshotted protocol")
    fitted = fit_candidate(evidence, objective=objective, evaluator=lambda profile: _evaluate(profile, evidence))
    output.mkdir(parents=True, exist_ok=True)
    (output / "decision_profile.json").write_bytes(fitted.profile.canonical_bytes)
    write_json(output / "metrics.json", _encode_metrics(fitted.evaluation.metrics))
    write_json(output / "evaluation.json", _encode_evaluation(fitted.evaluation))
    write_json(
        output / "fit_attestation.json",
        {
            "schema_version": "calibration-fit-attestation-v1",
            "profile_sha256": fitted.profile.digest,
            "evidence_sha256": evidence.dataset_sha256,
            "objective": objective,
            "baseline_reproduced": fitted.baseline_reproduced,
            "accessed_roles": list(fitted.accessed_roles),
        },
    )
    write_evaluation_artifacts(
        output,
        phase="fitted",
        profile=fitted.profile,
        evidence=evidence,
        evaluation=fitted.evaluation,
        accessed_roles=fitted.accessed_roles,
    )
    write_checksums(output)


def validate_artifact_bundle(profile_path: Path, evidence_path: Path, output: Path) -> None:
    """Evaluate one already-fixed generated profile on validation evidence only."""
    profile = resolve_decision_profile(profile_path)
    evidence = _load_roles(evidence_path, ("validation",))
    workflow = validate_candidate(profile, evidence)
    evaluation = _evaluate(profile, evidence)
    output.mkdir(parents=True, exist_ok=True)
    metrics_raw = _encode_metrics(evaluation.metrics)
    write_json(output / "metrics.json", metrics_raw)
    write_json(
        output / "attestation.json",
        _attestation(
            "validation",
            profile,
            evidence,
            canonical_sha256(metrics_raw),
            select_candidate((evaluation,), evidence.study.protocol) is not None,
        ),
    )
    write_json(
        output / "access.json",
        {"schema_version": "calibration-access-v1", "accessed_roles": list(workflow.accessed_roles)},
    )
    write_evaluation_artifacts(
        output,
        phase="validation",
        profile=profile,
        evidence=evidence,
        evaluation=evaluation,
        accessed_roles=workflow.accessed_roles,
    )
    write_checksums(output)


def evaluate_artifact_bundle(profile_path: Path, evidence_path: Path, output: Path) -> None:
    """Evaluate one fixed profile after one-use locked-payload precommit."""
    profile = resolve_decision_profile(profile_path)
    role_root = evidence_path / "roles" / "locked-heldout"
    manifest = decode_evidence_manifest(load_object(role_root / "evidence_manifest.json", "locked evidence manifest"))
    custody = evidence_path.parent / f".{evidence_path.name}.calibration-custody"

    def evaluator(raw: bytes) -> dict[str, object]:
        evidence = _decode_locked_payload(raw)
        if evidence.study.protocol.sha256 != manifest.protocol_sha256:
            raise ValueError("locked calibration protocol hash differs from its payload")
        if evidence.study.partitions.sha256 != manifest.partition_manifest_sha256:
            raise ValueError("locked calibration partition hash differs from its payload")
        if evidence.labels.sha256 != manifest.labels_sha256:
            raise ValueError("locked calibration labels hash differs from its payload")
        if canonical_sha256(thaw_json(evidence.baseline)) != manifest.baseline_sha256:
            raise ValueError("locked calibration baseline hash differs from its payload")
        evaluation = _evaluate(profile, evidence)
        write_evaluation_artifacts(
            output,
            phase="held-out",
            profile=profile,
            evidence=evidence,
            evaluation=evaluation,
            accessed_roles=("locked-heldout",),
        )
        return {
            "metrics": _encode_metrics(evaluation.metrics),
            "evidence_sha256": evidence.dataset_sha256,
            "accessed_roles": ["locked-heldout"],
            "passed": select_candidate((evaluation,), evidence.study.protocol) is not None,
        }

    locked = evaluate_locked_candidate(
        profile,
        role_root / "locked_payload.json",
        manifest,
        custody,
        evaluator=evaluator,
    )
    if not isinstance(locked.result, Mapping):
        raise ValueError("locked calibration evaluator returned an invalid result")
    metrics_raw = locked.result.get("metrics")
    if not isinstance(metrics_raw, Mapping):
        raise ValueError("locked calibration evaluator omitted metrics")
    metrics = decode_metrics(metrics_raw)
    evidence_sha256 = locked.result.get("evidence_sha256")
    if not isinstance(evidence_sha256, str):
        raise ValueError("locked calibration evaluator omitted its evidence hash")
    passed = locked.result.get("passed")
    if not isinstance(passed, bool):
        raise ValueError("locked calibration evaluator omitted its validation status")
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "metrics.json", dict(metrics_raw))
    attestation = _attestation_raw(
        "locked-heldout",
        profile.digest,
        manifest.protocol_sha256,
        evidence_sha256,
        metrics.sha256,
        passed,
    )
    write_json(output / "attestation.json", attestation)
    write_json(
        output / "custody_limitations.json",
        {
            "schema_version": "calibration-custody-limitations-v1",
            "local_custody_is_independent_proof": False,
            "statement": (
                "Local append-only custody guards are not proof of independent custody; "
                "closure requires a named external custodian."
            ),
        },
    )
    write_json(
        output / "access.json",
        {"schema_version": "calibration-access-v1", "accessed_roles": ["locked-heldout"]},
    )
    write_json(output / "custody_receipt.json", load_object(locked.receipt.path, "custody receipt"))
    write_checksums(output)


def _load_roles(root: Path, roles: tuple[str, ...]) -> ExtractedEvidence:
    verify_checksums(root)
    study_raw = load_object(root / "study.json", "calibration study")
    study = decode_study_declaration(study_raw)
    feature_rows: list[object] = []
    label_rows: list[object] = []
    expected_baseline_rows: list[object] = []
    observed_baseline_rows: list[object] = []
    run_hashes: dict[str, dict[str, str]] = {}
    for role in roles:
        if role not in _ROLES:
            raise ValueError(f"unsupported calibration evidence role: {role}")
        role_root = root / "roles" / role
        verify_checksums(role_root)
        features_document = load_object(role_root / "features.json", "features")
        labels_document = load_object(role_root / "labels.json", "labels")
        role_feature_rows = _rows(features_document, "features")
        role_label_rows = _rows(labels_document, "labels")
        baseline = load_object(role_root / "baseline.json", "baseline")
        expected = _mapping(baseline.get("expected"), "baseline expected")
        observed = _mapping(baseline.get("observed"), "baseline observed")
        role_expected_rows = _sequence(expected.get("rows"), "baseline expected rows")
        role_observed_rows = _sequence(observed.get("rows"), "baseline observed rows")
        run_document = load_object(role_root / "run_hashes.json", "run hashes")
        if set(run_document) != {"schema_version", "run_hashes"} or run_document["schema_version"] != (
            "calibration-run-hashes-v1"
        ):
            raise ValueError("calibration role run-hash fields or schema version differ")
        raw_hashes = _mapping(run_document.get("run_hashes"), "run hashes")
        role_run_hashes = decode_run_hashes(raw_hashes)
        inputs = RoleInputs(features_document, labels_document, baseline, role_run_hashes)
        observed_manifest = load_object(role_root / "manifest.json", "calibration role bundle manifest")
        if observed_manifest != _role_manifest(study, role, inputs, study.sha256):
            raise ValueError(f"calibration {role} role bundle does not match its study and artifacts")
        expected_keys = tuple(member.key for member in study.partitions.members if member.role == role)
        feature_keys = tuple(str(_mapping(row, "feature row")["manifest_key"]) for row in role_feature_rows)
        label_keys = tuple(str(_mapping(row, "label row")["manifest_key"]) for row in role_label_rows)
        if (
            feature_keys != expected_keys
            or label_keys != expected_keys
            or tuple(sorted(role_run_hashes)) != expected_keys
        ):
            raise ValueError(f"calibration {role} role bundle members do not match the partition declaration")
        feature_rows.extend(role_feature_rows)
        label_rows.extend(role_label_rows)
        expected_baseline_rows.extend(role_expected_rows)
        observed_baseline_rows.extend(role_observed_rows)
        run_hashes.update(role_run_hashes)
    feature_rows.sort(key=lambda row: str(_mapping(row, "feature row")["manifest_key"]))
    label_rows.sort(key=lambda row: str(_mapping(row, "label row")["manifest_key"]))
    expected_baseline_rows.sort(key=lambda row: str(_mapping(row, "baseline row")["manifest_key"]))
    observed_baseline_rows.sort(key=lambda row: str(_mapping(row, "baseline row")["manifest_key"]))
    features_raw = {"schema_version": "calibration-features-v1", "rows": feature_rows}
    labels_raw = {"schema_version": "calibration-labels-v1", "rows": label_rows}
    features = decode_feature_artifact(features_raw)
    labels = decode_label_artifact(labels_raw)
    baseline = project_baseline(
        expected_baseline_rows,
        observed_baseline_rows,
        {row.manifest_key: row for row in labels.rows},
    )
    dataset = canonical_sha256(
        {
            "study_sha256": study.sha256,
            "features_sha256": features.sha256,
            "labels_sha256": labels.sha256,
            "baseline_sha256": canonical_sha256(baseline),
            "run_artifact_sha256": run_hashes,
        }
    )
    return ExtractedEvidence(
        study,
        features,
        labels,
        cast(Mapping[str, object], freeze_json(baseline)),
        MappingProxyType({key: MappingProxyType(dict(value)) for key, value in sorted(run_hashes.items())}),
        study.sha256,
        dataset,
    )


def _evaluate(profile: ResolvedDecisionProfile, evidence: ExtractedEvidence) -> CandidateEvaluation:
    observations, baseline_observations = _observations(profile, evidence)
    protocol = evidence.study.protocol
    component = _mapping(profile.components["dominance"], "dominance component")
    summary = calculate_metrics(
        observations,
        profile_sha256=profile.digest,
        free_parameter_count=count_free_parameters(component),
        required_strata=protocol.required_strata,
    )
    groups = connected_leakage_groups(evidence.study.partitions)
    detection_pairs = []
    exact_pairs = []
    for candidate, baseline in zip(observations, baseline_observations, strict=True):
        if candidate.expected_identity is None:
            continue
        stratum = f"{candidate.assay_class}:{candidate.mutation_class}"
        group = groups[candidate.key]
        detection_pairs.append(
            PairedObservation(
                group,
                stratum,
                Fraction(baseline.selected_identity is not None),
                Fraction(candidate.selected_identity is not None),
            )
        )
        exact_pairs.append(
            PairedObservation(
                group,
                stratum,
                Fraction(baseline.selected_identity == baseline.expected_identity),
                Fraction(candidate.selected_identity == candidate.expected_identity),
            )
        )
    detection = paired_group_bootstrap(
        detection_pairs,
        required_strata=protocol.required_strata,
        iterations=protocol.bootstrap_iterations,
        seed=protocol.seed,
    )
    exact = paired_group_bootstrap(
        exact_pairs,
        required_strata=protocol.required_strata,
        iterations=protocol.bootstrap_iterations,
        seed=protocol.seed,
    )
    family_p = max(
        detection.one_sided_noninferiority_p_value,
        exact.one_sided_noninferiority_p_value,
    )
    return CandidateEvaluation(
        summary.metrics,
        detection.one_sided_lower,
        exact.one_sided_lower,
        summary.stratum_counts,
        family_p,
    )


def _observations(
    profile: ResolvedDecisionProfile, evidence: ExtractedEvidence
) -> tuple[tuple[OutcomeObservation, ...], tuple[OutcomeObservation, ...]]:
    component = _mapping(profile.components["dominance"], "dominance component")
    labels = {row.manifest_key: row for row in evidence.labels.rows}
    members = {member.key: member for member in evidence.study.partitions.members}
    rows = _sequence(_mapping(evidence.baseline["expected"], "baseline expected")["rows"], "baseline rows")
    if len(rows) != len(evidence.features.rows):
        raise ValueError("calibration baseline rows must align with role feature rows")
    candidate_rows = []
    baseline_rows = []
    for feature, baseline_raw in zip(evidence.features.rows, rows, strict=True):
        label = labels[feature.manifest_key]
        baseline = _mapping(baseline_raw, "baseline row")
        decision = replay_scalar_dominance(feature.features, component)
        selected = decision.selected_identity
        abstained = decision.abstention_reason is not None
        applicable = decision.applicable
        baseline_identity = baseline.get("canonical_identity")
        displayed, tier = _rendered_identity_projection(baseline, selected)
        assay = feature.features.get("assay_class")
        if not isinstance(assay, str) or not assay:
            raise ValueError("calibration feature row requires assay_class")
        baseline_tier_value = baseline.get("tier")
        baseline_tier = baseline_tier_value if isinstance(baseline_tier_value, str) else None
        candidate_rows.append(
            OutcomeObservation(
                feature.manifest_key,
                assay,
                label.mutation_class,
                label.expected_identity,
                label.expected_display_name,
                selected,
                displayed if isinstance(displayed, str) else None,
                tier if isinstance(tier, str) else None,
                abstained,
                applicable,
                True,
                baseline_tier,
            )
        )
        baseline_name = baseline.get("name")
        baseline_selected = baseline_identity if isinstance(baseline_identity, str) else None
        baseline_rows.append(
            OutcomeObservation(
                feature.manifest_key,
                assay,
                label.mutation_class,
                label.expected_identity,
                label.expected_display_name,
                baseline_selected,
                baseline_name if isinstance(baseline_name, str) else None,
                baseline_tier,
                baseline.get("abstention") is not None,
                True,
                True,
                baseline_tier,
            )
        )
        if members[feature.manifest_key].role == "training":
            candidate_rows.pop()
            baseline_rows.pop()
            continue
        if members[feature.manifest_key].role not in {"policy-selection", "validation", "locked-heldout"}:
            raise ValueError("calibration metrics cannot be calculated from training-role outcomes")
    return tuple(candidate_rows), tuple(baseline_rows)


def _rendered_identity_projection(
    baseline: Mapping[str, object], selected_identity: str | None
) -> tuple[str | None, str | None]:
    """Look up a selected identity in the closed non-feature rendering map."""
    projection = _mapping(baseline.get("identity_projection"), "baseline identity projection")
    for identity, raw in projection.items():
        if not isinstance(identity, str) or not identity:
            raise ValueError("calibration baseline identity projection keys must be non-empty strings")
        value = _mapping(raw, "baseline identity rendering")
        if set(value) != {"name", "tier"}:
            raise ValueError("calibration baseline identity rendering fields differ")
        if not isinstance(value["name"], str) or not value["name"]:
            raise ValueError("calibration baseline identity rendering name must be non-empty")
        if value["tier"] is not None and (not isinstance(value["tier"], str) or not value["tier"]):
            raise ValueError("calibration baseline identity rendering tier must be non-empty or null")
    if selected_identity is None:
        return None, None
    selected = projection.get(selected_identity)
    if selected is None:
        raise ValueError("calibration selected identity is absent from the closed baseline rendering projection")
    value = _mapping(selected, "selected baseline identity rendering")
    name = value["name"]
    tier = value["tier"]
    assert isinstance(name, str)
    assert tier is None or isinstance(tier, str)
    return name, tier


def _write_locked_payload(role_root: Path, study_raw: Mapping[str, object], study, inputs: RoleInputs) -> None:
    payload = {
        "schema_version": "calibration-locked-payload-v1",
        "study": dict(study_raw),
        "features": inputs.features,
        "labels": inputs.labels,
        "baseline": inputs.baseline,
        "run_hashes": inputs.run_hashes,
    }
    raw = canonical_json_bytes(payload)
    (role_root / "locked_payload.json").write_bytes(raw)
    evidence = {
        "schema_version": "calibration-evidence-v1",
        "role": "locked-heldout",
        "provenance": "external-custodian",
        "protocol_sha256": study.protocol.sha256,
        "partition_manifest_sha256": study.partitions.sha256,
        "features_sha256": hashlib.sha256(raw).hexdigest(),
        "labels_sha256": canonical_sha256(inputs.labels),
        "baseline_sha256": canonical_sha256(inputs.baseline),
    }
    write_json(role_root / "evidence_manifest.json", evidence)


def _decode_locked_payload(raw: bytes) -> ExtractedEvidence:
    payload = load_strict_json_object(raw)
    expected = {"schema_version", "study", "features", "labels", "baseline", "run_hashes"}
    if set(payload) != expected or payload["schema_version"] != "calibration-locked-payload-v1":
        raise ValueError("locked calibration payload fields or schema version differ")
    study = decode_study_declaration(payload["study"])
    features = decode_feature_artifact(payload["features"])
    labels = decode_label_artifact(payload["labels"])
    hashes = decode_run_hashes(_mapping(payload["run_hashes"], "run hashes"))
    baseline = dict(_mapping(payload["baseline"], "baseline"))
    dataset = canonical_sha256(
        {
            "study_sha256": study.sha256,
            "features_sha256": features.sha256,
            "labels_sha256": labels.sha256,
            "baseline_sha256": canonical_sha256(baseline),
            "run_artifact_sha256": hashes,
        }
    )
    return ExtractedEvidence(
        study,
        features,
        labels,
        cast(Mapping[str, object], freeze_json(baseline)),
        MappingProxyType({key: MappingProxyType(dict(value)) for key, value in hashes.items()}),
        study.sha256,
        dataset,
    )


def _role_manifest(study, role: str, inputs: RoleInputs, study_sha256: str) -> dict[str, object]:
    return {
        "schema_version": "calibration-role-bundle-v1",
        "role": role,
        "study_sha256": study_sha256,
        "protocol_sha256": study.protocol.sha256,
        "partition_manifest_sha256": study.partitions.sha256,
        "features_sha256": canonical_sha256(inputs.features),
        "labels_sha256": canonical_sha256(inputs.labels),
        "baseline_sha256": canonical_sha256(inputs.baseline),
        "run_hashes_sha256": canonical_sha256(inputs.run_hashes),
    }


def _encode_metrics(metrics: CandidateMetrics) -> dict[str, object]:
    return {
        "schema_version": "calibration-metrics-v1",
        "candidate_profile_sha256": metrics.candidate_profile_sha256,
        "wrong_tier_a_displayed_names": metrics.wrong_tier_a_displayed_names,
        "control_findings": metrics.control_findings,
        "wrong_displayed_names_all_tiers": metrics.wrong_displayed_names_all_tiers,
        "macro_exact_recovery": str(metrics.macro_exact_recovery),
        "binary_detection_sensitivity": str(metrics.binary_detection_sensitivity),
        "free_parameter_count": metrics.free_parameter_count,
        "abstention_fraction": str(metrics.abstention_fraction),
        "tier_a_reachable": metrics.tier_a_reachable,
        "applicability_matches": metrics.applicability_matches,
    }


def _encode_evaluation(evaluation: CandidateEvaluation) -> dict[str, object]:
    return {
        "schema_version": "calibration-candidate-evaluation-v1",
        "metrics_sha256": evaluation.metrics.sha256,
        "detection_lower_bound": str(evaluation.detection_lower_bound),
        "macro_exact_lower_bound": str(evaluation.macro_exact_lower_bound),
        "stratum_counts": list(evaluation.stratum_counts),
        "holm_adjusted_p_value": str(evaluation.holm_adjusted_p_value),
    }


def _attestation(
    role: str,
    profile: ResolvedDecisionProfile,
    evidence: ExtractedEvidence,
    metrics_sha256: str,
    passed: bool,
) -> dict[str, object]:
    return _attestation_raw(
        role,
        profile.digest,
        evidence.study.protocol.sha256,
        evidence.dataset_sha256,
        metrics_sha256,
        passed,
    )


def _attestation_raw(
    role: str,
    profile_sha256: str,
    protocol_sha256: str,
    evidence_sha256: str,
    metrics_sha256: str,
    passed: bool,
) -> dict[str, object]:
    return {
        "schema_version": "calibration-attestation-v1",
        "role": role,
        "status": "passed" if passed else "failed",
        "profile_sha256": profile_sha256,
        "protocol_sha256": protocol_sha256,
        "evidence_sha256": evidence_sha256,
        "metrics_sha256": metrics_sha256,
    }


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be an object")
    return value


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{label} must be a sequence")
    return value


def _rows(value: object, label: str) -> Sequence[object]:
    root = _mapping(value, label)
    return _sequence(root.get("rows"), f"{label} rows")


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"calibration {label} must be an integer")
    return value
