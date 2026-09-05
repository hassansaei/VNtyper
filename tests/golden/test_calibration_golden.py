"""Development-only calibration replay and external-evidence blocker."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.golden import calibration_oracle
from tests.golden.calibration_oracle import (
    load_development_snapshot,
    materialize_development_fixture,
)
from tests.golden.identity_oracle import DisplayCounts
from vntyper.scripts.decision_profile import resolve_decision_profile
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.golden

REPO_ROOT = Path(__file__).parents[2]
PACKAGED_PROFILE_SHA256 = "0b13d07370491b3ea773e65144891cb30caebcae70b0ef98feb0f2c5ccd2f4a1"
PACKAGED_PROJECTION_SHA256 = "cfa5ec402a3a20096b76273c4347ff8b5975db942aa6dccf9f2d99474260236d"
FITTED_PROFILE_SHA256 = "b1a25f00e8cc828ca1117c56a63caaab39f71825d00472346ddf59f98f417b5e"
SOURCE_SHA256 = {
    "simulation/experiment1_dupC/ground_truth.csv": "007026c594f2182a4385aa49a9cd2892b5729c6e4c9b6441cb0e6e10b1458a73",
    "advntr/experiment1_dupC/pair_3000/mutated/pipeline_summary.json": (
        "542aad10253cd8fb1c8efc2cbabf73dddd3e112f5bc6c5cf93b6af6d6a29880d"
    ),
    "advntr/experiment1_dupC/pair_3000/mutated/kestrel/kestrel_result.tsv": (
        "b728d1321f03e95b9cebab892f84f59a2838dae4e6e6c1f652f22573e44a9654"
    ),
    "advntr/experiment1_dupC/pair_3001/mutated/pipeline_summary.json": (
        "4c5d6f843d3ec3bf883e89c8bf5279ba487dd88c5fd95ffb9b976e0d0eb67b67"
    ),
    "advntr/experiment1_dupC/pair_3001/mutated/kestrel/kestrel_result.tsv": (
        "753535b4ccefb75e47db043a68215beafb37b2e06f1c655893bfb689358fe7cc"
    ),
    "advntr/experiment1_dupC/pair_3003/mutated/pipeline_summary.json": (
        "5a0f8acadd8d04ca86d3cec1bd4607b0ee7788a10cc656710b51005ef1da03ee"
    ),
    "advntr/experiment1_dupC/pair_3003/mutated/kestrel/kestrel_result.tsv": (
        "e2f93d452e93ab5f324779e07b594d8edc01d1bc5916f4378743f3bccae65374"
    ),
    "advntr/experiment1_dupC/pair_3004/mutated/pipeline_summary.json": (
        "534b451037f3185dc3ecaf3a886fd26d838f6c1e360bd5412030c96498a54139"
    ),
    "advntr/experiment1_dupC/pair_3004/mutated/kestrel/kestrel_result.tsv": (
        "e929dad91812dda639ad825d497b34d1f6c204c83e53dd7fac5f8d7017820b46"
    ),
}


@pytest.fixture(scope="module")
def explicit_roots() -> tuple[Path, Path]:
    """Require both explicit corpus roots; absence is a golden failure, never a skip."""
    return calibration_oracle.require_explicit_roots(os.environ)


@pytest.fixture(scope="module")
def workflow_events() -> list[str]:
    """Record completed replay and CLI phases for the ordering contract."""
    return []


@pytest.fixture(scope="module")
def packaged_replay(explicit_roots: tuple[Path, Path], workflow_events: list[str]):
    """Complete the literal historical projection before any candidate workflow."""
    replay = load_development_snapshot(*explicit_roots)
    assert replay.mutated_samples == 200
    assert replay.control_samples == 200
    assert replay.public_identity_rows == 374
    assert replay.selected_locus_rows == 178
    assert replay.total == DisplayCounts(displayed=154, exact=136, wrong=18)
    assert replay.by_tier == {
        "A": DisplayCounts(displayed=53, exact=53, wrong=0),
        "B": DisplayCounts(displayed=101, exact=83, wrong=18),
        "C": DisplayCounts(displayed=0, exact=0, wrong=0),
    }
    assert replay.control_findings == 0
    profile = REPO_ROOT / "vntyper" / "profiles" / "decision_profile.json"
    projection = REPO_ROOT / "vntyper" / "profiles" / "decision_projection.json"
    assert hashlib.sha256(profile.read_bytes()).hexdigest() == PACKAGED_PROFILE_SHA256
    assert hashlib.sha256(projection.read_bytes()).hexdigest() == PACKAGED_PROJECTION_SHA256
    workflow_events.append("packaged-replay-complete")
    return replay


@pytest.fixture(scope="module")
def snapshot(packaged_replay):
    """Expose the completed replay with no per-test skip fallback."""
    return packaged_replay


@pytest.fixture(scope="module")
def completed_workflow(
    explicit_roots: tuple[Path, Path],
    packaged_replay,
    workflow_events: list[str],
    tmp_path_factory: pytest.TempPathFactory,
):
    """Materialize the independent bridge, then run the real extract and fit CLI."""
    assert packaged_replay.evidence_role == "previously-examined-development-simulation"
    work = tmp_path_factory.mktemp("calibration-golden")
    fixture = materialize_development_fixture(*explicit_roots, REPO_ROOT, work / "inputs")
    evidence = work / "evidence"
    candidate = work / "candidate"
    extract = _run_cli(
        fixture.truth_path.parent,
        "extract",
        "--truth",
        fixture.truth_path,
        "--partitions",
        fixture.study_path,
        "--runs",
        fixture.runs_path,
        "--output",
        evidence,
    )
    if extract.returncode != 0:
        raise AssertionError(f"calibration extract failed before fit: {extract.stderr}")
    workflow_events.append("extract-complete")
    fit = _run_cli(
        fixture.truth_path.parent,
        "fit",
        "--evidence",
        evidence,
        "--objective",
        "lexicographic-safety-v1",
        "--output",
        candidate,
    )
    if fit.returncode != 0:
        raise AssertionError(f"calibration fit failed: {fit.stderr}")
    workflow_events.append("fit-complete")
    return fixture, evidence, candidate, extract, fit


def test_missing_explicit_root_is_a_failure_not_a_skip() -> None:
    with pytest.raises(AssertionError, match="VNTYPER_SIM_ROOT"):
        calibration_oracle.require_explicit_roots({})
    with pytest.raises(AssertionError, match="VNTYPER_ADVNTR_ROOT"):
        calibration_oracle.require_explicit_roots({"VNTYPER_SIM_ROOT": str(REPO_ROOT)})


@pytest.mark.parametrize("incomplete_name", ["VNTYPER_SIM_ROOT", "VNTYPER_ADVNTR_ROOT"])
def test_present_but_incomplete_explicit_root_is_a_failure(
    explicit_roots: tuple[Path, Path], tmp_path: Path, incomplete_name: str
) -> None:
    sim_root = tmp_path / "simulation"
    advntr_root = tmp_path / "advntr"
    sim_root.mkdir()
    advntr_root.mkdir()

    with pytest.raises(AssertionError, match=rf"{incomplete_name}.*incomplete"):
        calibration_oracle.require_explicit_roots(
            {
                "VNTYPER_SIM_ROOT": str(sim_root if incomplete_name == "VNTYPER_SIM_ROOT" else explicit_roots[0]),
                "VNTYPER_ADVNTR_ROOT": str(
                    advntr_root if incomplete_name == "VNTYPER_ADVNTR_ROOT" else explicit_roots[1]
                ),
            }
        )


def test_calibration_oracle_has_no_production_import_path() -> None:
    scanned = calibration_oracle.assert_independent_import_closure(REPO_ROOT)

    assert Path(calibration_oracle.__file__).resolve() in scanned


def test_both_roots_complete_population_and_row_locus_counts_are_loaded(
    snapshot, explicit_roots: tuple[Path, Path]
) -> None:
    assert snapshot.sim_root == explicit_roots[0].resolve()
    assert snapshot.advntr_root == explicit_roots[1].resolve()
    assert snapshot.mutated_samples == 200
    assert snapshot.control_samples == 200
    assert snapshot.public_identity_rows == 374
    assert snapshot.selected_locus_rows == 178


def test_shipped_projection_is_reproduced_before_any_candidate_claim(snapshot) -> None:
    assert snapshot.total == DisplayCounts(displayed=154, exact=136, wrong=18)
    assert snapshot.by_tier == {
        "A": DisplayCounts(displayed=53, exact=53, wrong=0),
        "B": DisplayCounts(displayed=101, exact=83, wrong=18),
        "C": DisplayCounts(displayed=0, exact=0, wrong=0),
    }
    assert snapshot.control_findings == 0

    profile = REPO_ROOT / "vntyper" / "profiles" / "decision_profile.json"
    projection = REPO_ROOT / "vntyper" / "profiles" / "decision_projection.json"
    assert hashlib.sha256(profile.read_bytes()).hexdigest() == PACKAGED_PROFILE_SHA256
    assert hashlib.sha256(projection.read_bytes()).hexdigest() == PACKAGED_PROJECTION_SHA256


def test_examined_simulations_are_not_misrepresented_as_locked_heldout(snapshot) -> None:
    assert snapshot.evidence_role == "previously-examined-development-simulation"
    assert snapshot.eligible_for_independent_validation is False
    assert snapshot.eligible_for_locked_evaluate is False
    assert "previously examined development evidence" in snapshot.ineligibility_reason
    assert "neither independent external validation nor a custodian-locked held-out cohort" in (
        snapshot.ineligibility_reason
    )


def test_root_derived_schema_three_bridge_binds_both_sources_and_declares_legacy_limit(completed_workflow) -> None:
    fixture, _evidence, _candidate, _extract, _fit = completed_workflow

    assert {f"{binding.root_kind}/{binding.relative_path}": binding.sha256 for binding in fixture.source_bindings} == (
        SOURCE_SHA256
    )
    assert tuple(member.key for member in fixture.members) == (
        "held-pair-3004",
        "select-pair-3001",
        "train-pair-3000",
        "validate-pair-3003",
    )
    assert tuple(member.role for member in fixture.members) == (
        "locked-heldout",
        "policy-selection",
        "training",
        "validation",
    )
    assert tuple(
        (
            member.key,
            member.source_summary_schema,
            member.motif_context,
            member.pair_context,
            member.support,
            member.active_region_depth,
            member.depth_score,
            member.confidence,
            member.flag,
            member.tier,
        )
        for member in fixture.members
    ) == (
        ("held-pair-3004", 2, "H", "H-C", 291, 9718, 0.029944433010907594, "High_Precision*", "Not flagged", "B"),
        ("select-pair-3001", 2, "K", "K-6", 192, 12560, 0.015286624203821656, "High_Precision*", "Not flagged", "A"),
        ("train-pair-3000", 2, "B", "B-W", 643, 11781, 0.05457940752058399, "High_Precision*", "Not flagged", "A"),
        ("validate-pair-3003", 2, "K", "K-Q", 489, 11855, 0.041248418388865456, "High_Precision*", "Not flagged", "A"),
    )
    classification = _json(fixture.classification_path)
    assert classification == {
        "schema_version": "calibration-development-evidence-classification-v1",
        "classification": "previously-examined-development-simulation",
        "eligible_for_independent_validation": False,
        "eligible_for_locked_evaluate": False,
        "statement": (
            "This development fixture is neither an independent external cohort nor custodian-locked heldout evidence."
        ),
    }


def test_real_cli_extract_installs_exact_role_inventory_and_literal_replay(completed_workflow) -> None:
    fixture, evidence, _candidate, extract, _fit = completed_workflow

    assert extract.returncode == 0, extract.stderr
    assert sorted(path.name for path in evidence.iterdir()) == [
        "checksums.json",
        "groups.json",
        "profile_binding.json",
        "roles",
        "run_commitments.json",
        "study.json",
    ]
    locked_root = evidence / "roles" / "locked-heldout"
    assert sorted(path.name for path in locked_root.iterdir()) == [
        "checksums.json",
        "member_declaration.json",
        "run_commitments.json",
    ]
    calibration_oracle.verify_checksum_tree(locked_root)
    calibration_oracle.verify_checksum_tree(evidence)
    calibration_oracle.verify_exact_child_directories(
        evidence / "roles",
        ("locked-heldout", "policy-selection", "training", "validation"),
    )
    for role in ("training", "policy-selection", "validation"):
        assert sorted(path.name for path in (evidence / "roles" / role).iterdir()) == [
            "baseline.json",
            "checksums.json",
            "features.json",
            "labels.json",
            "manifest.json",
            "run_hashes.json",
        ]
        calibration_oracle.verify_checksum_tree(evidence / "roles" / role)

    expected = {member.key: member for member in fixture.members}
    for role, key in (
        ("policy-selection", "select-pair-3001"),
        ("training", "train-pair-3000"),
        ("validation", "validate-pair-3003"),
    ):
        features_document = _json(evidence / "roles" / role / "features.json")
        feature_rows = features_document["rows"]
        assert isinstance(feature_rows, list) and len(feature_rows) == 1
        features = feature_rows[0]
        assert isinstance(features, dict)
        baseline = _json(evidence / "roles" / role / "baseline.json")
        member = expected[key]
        assert features["manifest_key"] == key
        assert features["features"] == {
            "active_region_depth": member.active_region_depth,
            "advntr_evidence_disposition": "admissible",
            "alternate_kmer_path_depth": member.support,
            "assay_class": "development-schema2-bridge",
            "canonical_identity": "MUC1-X-60-coding-v1|60|59|-|C",
            "cooccurring_identity_count": 1,
            "context_local_representation_count": 1,
            "decision_policy": "legacy-selection-v1",
            "decision_profile_sha256": PACKAGED_PROFILE_SHA256,
            "depth_score": member.depth_score,
            "haplotype_record_count": 2,
            "haplotype_record_count_margin": 2,
            "haplotype_record_share": 1.0,
            "haplotype_record_share_margin": 1.0,
            "haplotype_record_tie": False,
            "motif_context": member.motif_context,
            "pair_context": member.pair_context,
            "reference_version": "hg19",
            "structural_gate": True,
            "tool_version": "2.0.26",
            "xd_availability_count": 2,
            "xd_availability_fraction": 1.0,
            "xd_interquartile_range": 4.0,
            "xd_median": 10.0,
        }
        expected_baseline = baseline["expected"]
        observed_baseline = baseline["observed"]
        assert isinstance(expected_baseline, dict) and expected_baseline == observed_baseline
        baseline_rows = expected_baseline["rows"]
        assert isinstance(baseline_rows, list) and len(baseline_rows) == 1
        row = baseline_rows[0]
        assert row == {
            "abstention": None,
            "canonical_identity": "MUC1-X-60-coding-v1|60|59|-|C",
            "confidence": member.confidence,
            "flag": member.flag,
            "identity_projection": {"MUC1-X-60-coding-v1|60|59|-|C": {"name": "59dupC", "tier": member.tier}},
            "manifest_key": key,
            "name": "59dupC",
            "order": 0,
            "support": member.support,
            "tie": False,
            "tier": member.tier,
        }


def test_locked_heldout_checksum_corruption_fails_independent_verification(completed_workflow, tmp_path: Path) -> None:
    """A changed locked declaration must fail the same independent checksum verifier."""
    _fixture, evidence, _candidate, _extract, _fit = completed_workflow
    corrupted = tmp_path / "locked-heldout"
    shutil.copytree(evidence / "roles" / "locked-heldout", corrupted)
    (corrupted / "member_declaration.json").write_text("{}\n", encoding="utf-8")

    with pytest.raises(AssertionError, match="calibration checksum differs"):
        calibration_oracle.verify_checksum_tree(corrupted)


def test_extra_evidence_role_directory_fails_independent_inventory(completed_workflow, tmp_path: Path) -> None:
    """A role unknown to the literal four-role study cannot hide in extract output."""
    _fixture, evidence, _candidate, _extract, _fit = completed_workflow
    corrupted = tmp_path / "roles"
    shutil.copytree(evidence / "roles", corrupted)
    (corrupted / "unreviewed-role").mkdir()

    with pytest.raises(AssertionError, match="calibration child-directory inventory differs"):
        calibration_oracle.verify_exact_child_directories(
            corrupted,
            ("locked-heldout", "policy-selection", "training", "validation"),
        )


def test_packaged_replay_completes_before_extract_and_fit(completed_workflow, workflow_events: list[str]) -> None:
    """Fixture dependencies must complete neutral replay before either CLI phase."""
    _fixture, _evidence, _candidate, extract, fit = completed_workflow
    assert (extract.returncode, fit.returncode) == (0, 0)
    assert workflow_events == ["packaged-replay-complete", "extract-complete", "fit-complete"]


def test_bridge_md5_is_explicitly_nonsecurity(
    explicit_roots: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The legacy summary checksum must opt out of security-sensitive MD5 use."""
    real_md5 = hashlib.md5

    def fips_md5(value: bytes, *, usedforsecurity: bool = True):
        if usedforsecurity:
            raise ValueError("MD5 is unavailable for security-sensitive use")
        return real_md5(value, usedforsecurity=usedforsecurity)

    monkeypatch.setattr(calibration_oracle.hashlib, "md5", fips_md5)
    fixture = materialize_development_fixture(*explicit_roots, REPO_ROOT, tmp_path / "inputs")

    assert len(fixture.members) == 4


def test_snapshot_constructor_binds_every_field_by_keyword(
    explicit_roots: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Eligibility booleans cannot shift position when the snapshot schema changes."""
    corpus = calibration_oracle.identity_oracle.load_golden_corpus(*explicit_roots)
    sentinel = object()

    def keyword_only_snapshot(
        *,
        sim_root: Path,
        advntr_root: Path,
        mutated_samples: int,
        control_samples: int,
        public_identity_rows: int,
        selected_locus_rows: int,
        total: DisplayCounts,
        by_tier: dict[str, DisplayCounts],
        control_findings: int,
        evidence_role: str,
        eligible_for_independent_validation: bool,
        eligible_for_locked_evaluate: bool,
        ineligibility_reason: str,
    ) -> object:
        assert sim_root == explicit_roots[0]
        assert advntr_root == explicit_roots[1]
        assert mutated_samples == control_samples == 200
        assert public_identity_rows == 374
        assert selected_locus_rows == 178
        assert total == DisplayCounts(154, 136, 18)
        assert by_tier["A"] == DisplayCounts(53, 53, 0)
        assert control_findings == 0
        assert evidence_role == "previously-examined-development-simulation"
        assert eligible_for_independent_validation is False
        assert eligible_for_locked_evaluate is False
        assert "previously examined development evidence" in ineligibility_reason
        return sentinel

    monkeypatch.setattr(calibration_oracle, "DevelopmentCalibrationSnapshot", keyword_only_snapshot)

    assert calibration_oracle.snapshot_from_corpus(corpus) is sentinel


def test_packaged_replay_precedes_fit_and_generated_profile_is_explicit_only(completed_workflow) -> None:
    _fixture, _evidence, candidate, _extract, fit = completed_workflow

    assert fit.returncode == 0, fit.stderr
    attestation = _json(candidate / "fit_attestation.json")
    profile = resolve_decision_profile(candidate / "decision_profile.json")
    run = resolve_run_configuration(candidate / "decision_profile.json")
    packaged = resolve_run_configuration(None)
    assert sorted(path.name for path in candidate.iterdir()) == [
        "abstentions.tsv",
        "checksums.json",
        "decision_profile.json",
        "evaluation.json",
        "fit_attestation.json",
        "grid.json",
        "intervals.json",
        "joint_surface.tsv",
        "metrics.json",
        "pr.tsv",
        "report.html",
        "report_data.json",
        "roc.tsv",
    ]
    calibration_oracle.verify_checksum_tree(candidate)
    assert attestation["baseline_reproduced"] is True
    assert attestation["accessed_roles"] == ["training", "policy-selection"]
    assert profile.digest == FITTED_PROFILE_SHA256
    assert attestation["profile_sha256"] == FITTED_PROFILE_SHA256
    assert profile.profile_kind == "generated"
    assert profile.source == "explicit-cli"
    assert run.decision_profile.digest == FITTED_PROFILE_SHA256
    assert packaged.decision_profile.digest == PACKAGED_PROFILE_SHA256
    assert packaged.decision_profile.profile_kind == "packaged"
    assert run.decision_profile.digest != packaged.decision_profile.digest
    assert profile.components["dominance"] == {
        "abstain_on_inadmissible_advntr": False,
        "enabled": True,
        "minimum_record_count_margin": 1,
        "minimum_record_share": 0.5,
        "minimum_record_share_margin": 0.0,
        "xd_veto": "disabled",
    }


def _run_cli(workdir: Path, operation: str, *arguments: object) -> subprocess.CompletedProcess[str]:
    argv = [
        sys.executable,
        "-m",
        "vntyper.cli",
        "calibrate",
        operation,
        *(str(argument) for argument in arguments),
    ]
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        value for value in (str(REPO_ROOT), environment.get("PYTHONPATH", "")) if value
    )
    return subprocess.run(argv, cwd=workdir, env=environment, capture_output=True, text=True, check=False)


def _json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value
