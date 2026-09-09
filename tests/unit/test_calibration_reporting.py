"""Tests for calibration reporting and candidate identity projection."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from typing import Literal

import pytest

from vntyper.scripts.calibration_features import (
    FeatureArtifact,
    FeatureRow,
    LabelArtifact,
    LabelRow,
)
from vntyper.scripts.calibration_manifest import StudyDeclaration, decode_study_declaration
from vntyper.scripts.calibration_objective import CandidateEvaluation, CandidateMetrics
from vntyper.scripts.calibration_reporting import _analyse, write_evaluation_artifacts
from vntyper.scripts.calibration_workflow import ExtractedEvidence
from vntyper.scripts.decision_profile import ResolvedDecisionProfile

pytestmark = pytest.mark.unit


def _groups(key: str) -> dict[str, list[str]]:
    return {
        namespace: [f"{namespace}:{key}"]
        for namespace in (
            "individual-family",
            "simulated-pair",
            "backbone-seed-lineage",
            "replicate-rerun",
            "depth-series-source",
            "batch",
            "repeat-context",
        )
    }


def _dummy_study() -> StudyDeclaration:
    members = [
        {
            "key": "sample_1",
            "role": "validation",
            "provenance": "development",
            "assay_class": "capture-short-read",
            "groups": _groups("sample_1"),
        },
        {
            "key": "sample_held",
            "role": "locked-heldout",
            "provenance": "external-custodian",
            "assay_class": "capture-short-read",
            "groups": _groups("sample_held"),
        },
        {
            "key": "sample_select",
            "role": "policy-selection",
            "provenance": "development",
            "assay_class": "capture-short-read",
            "groups": _groups("sample_select"),
        },
        {
            "key": "sample_train",
            "role": "training",
            "provenance": "development",
            "assay_class": "capture-short-read",
            "groups": _groups("sample_train"),
        },
    ]
    return decode_study_declaration(
        {
            "schema_version": "calibration-study-v1",
            "protocol": {
                "objective": "lexicographic-safety-v1",
                "bootstrap_iterations": 10000,
                "bootstrap_interval": "percentile",
                "multiplicity_method": "holm",
                "seed": 295,
                "maximum_free_parameters": 4,
                "minimum_stratum_count": 2,
                "maximum_abstention_fraction": 0.25,
                "assay_classes": ["capture-short-read"],
                "mutation_classes": ["duplication", "insertion"],
                "candidate_grid": {
                    "minimum_record_count_margin": [1, 2],
                    "minimum_record_share": [0.5, 0.75],
                    "minimum_record_share_margin": [0.0, 0.25],
                    "xd_veto": ["disabled", "missingness"],
                },
            },
            "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
        }
    )


def _build_evidence(
    *,
    candidate_identity: str | None,
    candidate_name: str,
    candidate_tier: str,
    truth_identity: str | None = "c.59dupC",
    truth_display_name: str | None = "59dupC",
    truth_status: Literal["mutated", "control"] = "mutated",
    tie: bool = False,
) -> tuple[ResolvedDecisionProfile, ExtractedEvidence]:
    dominance_component = {
        "enabled": True,
        "minimum_record_count_margin": 1,
        "minimum_record_share": 0.5,
        "minimum_record_share_margin": 0.1,
        "xd_veto": "disabled",
        "abstain_on_inadmissible_advntr": False,
    }
    profile = ResolvedDecisionProfile(
        profile_id="test_profile",
        profile_revision="1",
        profile_kind="generated",
        source="package",
        digest="a" * 64,
        canonical_bytes=b"{}",
        document={},
        components={"dominance": dominance_component},
    )
    feature_dict: dict[str, object] = {
        "assay_class": "capture-short-read",
        "haplotype_record_count_margin": 10,
        "haplotype_record_share": 0.8,
        "haplotype_record_share_margin": 0.5,
        "haplotype_record_tie": tie,
    }
    if candidate_identity is not None:
        feature_dict["canonical_identity"] = candidate_identity

    feature_row = FeatureRow(
        feature_key="feat_1",
        manifest_key="sample_1",
        features=feature_dict,
    )
    label_row = LabelRow(
        label_key="lbl_1",
        manifest_key="sample_1",
        truth_status=truth_status,
        expected_identity=truth_identity,
        expected_display_name=truth_display_name,
        mutation_class="duplication",
    )
    baseline_row = {
        "manifest_key": "sample_1",
        "name": "baseline_default_name",
        "tier": "C",
        "canonical_identity": "c.54_56delinsAT",
        "identity_projection": {
            "c.59dupC": {"name": candidate_name, "tier": candidate_tier},
            "c.54_56delinsAT": {"name": "frameshift +1, representation-limited", "tier": "A"},
        },
    }
    evidence = ExtractedEvidence(
        study=_dummy_study(),
        features=FeatureArtifact(rows=(feature_row,), sha256="f" * 64),
        labels=LabelArtifact(rows=(label_row,), sha256="l" * 64),
        baseline={
            "expected": {
                "aggregate": {"displayed": 0, "exact": 0, "wrong": 0, "control_findings": 0},
                "per_tier": {},
                "rows": (baseline_row,),
            }
        },
        run_hashes={},
        study_sha256="s" * 64,
        dataset_sha256="d" * 64,
        profile_dataset_sha256="p" * 64,
    )
    return profile, evidence


def _tier_metrics_map(analysis: dict[str, object]) -> dict[str, dict[str, object]]:
    raw = analysis["tier_metrics"]
    assert isinstance(raw, list)
    result: dict[str, dict[str, object]] = {}
    for item in raw:
        assert isinstance(item, dict)
        tier = item["tier"]
        assert isinstance(tier, str)
        result[tier] = item
    return result


def test_analyse_projects_candidate_identity_displayed_name_and_tier() -> None:
    """Reporting analyzes candidate's projected identity rather than baseline default row fields."""
    profile, evidence = _build_evidence(
        candidate_identity="c.59dupC",
        candidate_name="59dupC",
        candidate_tier="A",
    )
    analysis = _analyse(profile, evidence)
    tier_metrics = _tier_metrics_map(analysis)

    assert tier_metrics["A"]["displayed"] == 1
    assert tier_metrics["A"]["exact"] == 1
    assert tier_metrics["A"]["wrong"] == 0
    # Baseline had tier C, but candidate projected to tier A:
    assert tier_metrics["C"]["displayed"] == 0


def test_analyse_withholds_representation_limited_candidate_name() -> None:
    """When candidate identity projects to representation-limited, displayed count is not incremented."""
    profile, evidence = _build_evidence(
        candidate_identity="c.54_56delinsAT",
        candidate_name="frameshift +1, representation-limited",
        candidate_tier="A",
        truth_identity="c.54_56delinsAT",
        truth_display_name="54_56delinsAT",
    )
    analysis = _analyse(profile, evidence)
    tier_metrics = _tier_metrics_map(analysis)

    assert tier_metrics["A"]["displayed"] == 0
    assert tier_metrics["A"]["exact"] == 0
    assert tier_metrics["A"]["wrong"] == 0


def test_analyse_abstained_candidate_records_abstention() -> None:
    """When candidate dominance abstains, displayed count is zero and abstention recorded."""
    profile, evidence = _build_evidence(
        candidate_identity="c.59dupC",
        candidate_name="59dupC",
        candidate_tier="A",
        tie=True,
    )
    analysis = _analyse(profile, evidence)
    tier_metrics = _tier_metrics_map(analysis)

    assert tier_metrics["A"]["displayed"] == 0
    assert tier_metrics["A"]["exact"] == 0
    assert tier_metrics["A"]["wrong"] == 0
    abstentions = analysis["abstentions"]
    assert isinstance(abstentions, list)
    assert any(isinstance(a, dict) and a.get("reason") == "record-tie" for a in abstentions)


def test_write_evaluation_artifacts_validates_inputs(tmp_path: Path) -> None:
    """Unsupported phase or invalid types are rejected with ValueError."""
    profile, evidence = _build_evidence(
        candidate_identity="c.59dupC",
        candidate_name="59dupC",
        candidate_tier="A",
    )
    metrics = CandidateMetrics(
        candidate_profile_sha256="a" * 64,
        wrong_tier_a_displayed_names=0,
        control_findings=0,
        wrong_displayed_names_all_tiers=0,
        macro_exact_recovery=Fraction(1, 1),
        binary_detection_sensitivity=Fraction(1, 1),
        free_parameter_count=1,
        abstention_fraction=Fraction(0, 1),
        tier_a_reachable=True,
        applicability_matches=True,
        sha256="m" * 64,
    )
    evaluation = CandidateEvaluation(
        metrics=metrics,
        detection_lower_bound=Fraction(1, 1),
        macro_exact_lower_bound=Fraction(1, 1),
        stratum_counts=(1, 1),
        holm_adjusted_p_value=Fraction(1, 100),
    )

    with pytest.raises(ValueError, match="unsupported calibration reporting phase"):
        write_evaluation_artifacts(
            tmp_path,
            phase="invalid_phase",
            profile=profile,
            evidence=evidence,
            evaluation=evaluation,
            accessed_roles=("validation",),
        )
