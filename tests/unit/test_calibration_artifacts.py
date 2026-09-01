"""Filesystem codecs and end-to-end calibration operation artifacts."""

import logging
import shutil
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.calibration_artifacts import (
    evaluate_artifact_bundle,
    extract_artifact_bundle,
    fit_artifact_bundle,
    validate_artifact_bundle,
)
from vntyper.scripts.calibration_contract import decode_attestation
from vntyper.scripts.calibration_statistics import BootstrapInterval
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.cli_calibrate import handle_calibrate
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


def test_extract_fit_validate_and_one_use_evaluate_round_trip(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    _run_cli(
        [
            "calibrate",
            "extract",
            "--truth",
            str(truth),
            "--partitions",
            str(partitions),
            "--runs",
            str(runs),
            "--output",
            str(evidence),
        ]
    )

    assert (evidence / "study.json").is_file()
    assert (evidence / "groups.json").is_file()
    assert (evidence / "roles" / "policy-selection" / "features.json").is_file()
    assert (evidence / "roles" / "locked-heldout" / "locked_payload.json").is_file()
    assert load_strict_json_object((evidence / "checksums.json").read_bytes())["schema_version"] == (
        "calibration-checksums-v1"
    )

    candidate = tmp_path / "candidate"
    _run_cli(
        [
            "calibrate",
            "fit",
            "--evidence",
            str(evidence),
            "--objective",
            "lexicographic-safety-v1",
            "--output",
            str(candidate),
        ]
    )
    profile = candidate / "decision_profile.json"
    assert profile.is_file()
    for name in (
        "report.html",
        "grid.json",
        "intervals.json",
        "roc.tsv",
        "pr.tsv",
        "joint_surface.tsv",
        "abstentions.tsv",
    ):
        assert (candidate / name).is_file(), name
    assert "optional minimum k-mer depth" in (candidate / "report.html").read_text(encoding="utf-8")
    assert load_strict_json_object((candidate / "fit_attestation.json").read_bytes())["accessed_roles"] == [
        "training",
        "policy-selection",
    ]

    validation = tmp_path / "validation"
    _run_cli(
        [
            "calibrate",
            "validate",
            "--profile",
            str(profile),
            "--evidence",
            str(evidence),
            "--output",
            str(validation),
        ]
    )
    validation_attestation = decode_attestation(load_strict_json_object((validation / "attestation.json").read_bytes()))
    assert validation_attestation.role == "validation"
    assert (validation / "report.html").is_file()
    assert (validation / "intervals.json").is_file()

    heldout = tmp_path / "heldout"
    _run_cli(
        [
            "calibrate",
            "evaluate",
            "--profile",
            str(profile),
            "--evidence",
            str(evidence),
            "--output",
            str(heldout),
        ]
    )
    heldout_attestation = load_strict_json_object((heldout / "attestation.json").read_bytes())
    assert decode_attestation(heldout_attestation).role == "locked-heldout"
    assert (heldout / "report.html").is_file()
    assert (heldout / "abstentions.tsv").is_file()
    limitations = load_strict_json_object((heldout / "custody_limitations.json").read_bytes())
    assert limitations["local_custody_is_independent_proof"] is False

    second = tmp_path / "heldout-again"
    with pytest.raises(ValueError, match="consumed|one-use"):
        _run_cli(
            [
                "calibrate",
                "evaluate",
                "--profile",
                str(profile),
                "--evidence",
                str(evidence),
                "--output",
                str(second),
            ]
        )
    assert not second.exists()


def test_fit_refuses_an_objective_different_from_the_snapshotted_protocol(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)

    with pytest.raises(ValueError, match="objective"):
        fit_artifact_bundle(evidence, "f1", tmp_path / "candidate")


def test_fit_verifies_only_authorized_role_payloads(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    locked_payload = evidence / "roles" / "locked-heldout" / "locked_payload.json"
    original_read = Path.read_bytes

    def guarded_read(path: Path) -> bytes:
        if path == locked_payload:
            raise AssertionError("fit opened locked held-out evidence")
        return original_read(path)

    candidate = tmp_path / "candidate"
    candidate.mkdir()
    with patch.object(Path, "read_bytes", guarded_read):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    assert (candidate / "decision_profile.json").is_file()


def test_fit_refuses_tampered_policy_selection_artifacts(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    features = evidence / "roles" / "policy-selection" / "features.json"
    features.write_bytes(features.read_bytes() + b" ")

    with pytest.raises(ValueError, match="checksum"):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", tmp_path / "candidate")


def test_validation_refuses_a_self_checksummed_role_directory_swap(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    validation_role = evidence / "roles" / "validation"
    shutil.rmtree(validation_role)
    shutil.copytree(evidence / "roles" / "policy-selection", validation_role)

    with pytest.raises(ValueError, match="role bundle|validation"):
        validate_artifact_bundle(
            candidate / "decision_profile.json",
            evidence,
            tmp_path / "validation-attestation",
        )


def test_locked_evaluation_refuses_manifest_protocol_mismatched_to_payload(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)
    manifest_path = evidence / "roles" / "locked-heldout" / "evidence_manifest.json"
    manifest = load_strict_json_object(manifest_path.read_bytes())
    manifest["protocol_sha256"] = "a" * 64
    manifest_path.write_bytes(canonical_json_bytes(manifest))

    with pytest.raises(ValueError, match="protocol"):
        evaluate_artifact_bundle(
            candidate / "decision_profile.json",
            evidence,
            tmp_path / "heldout-attestation",
        )


def test_fit_rejects_a_declared_but_unobserved_bootstrap_stratum(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    study = load_strict_json_object(partitions.read_bytes())
    protocol = study["protocol"]
    assert isinstance(protocol, dict)
    protocol["assay_classes"] = ["capture-short-read", "genome-short-read"]
    partitions.write_bytes(canonical_json_bytes(study))
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)

    with pytest.raises(ValueError, match="empty declared bootstrap strata.*genome-short-read:duplication"):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", tmp_path / "candidate")


def test_fit_does_not_infer_canonical_identity_correctness_from_display_name(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    truth_document = load_strict_json_object(truth.read_bytes())
    features = truth_document["features"]
    assert isinstance(features, dict)
    rows = features["rows"]
    assert isinstance(rows, list)
    for raw in rows:
        assert isinstance(raw, dict)
        feature_values = raw["features"]
        assert isinstance(feature_values, dict)
        feature_values["canonical_identity"] = f"wrong-{raw['manifest_key']}"
    truth.write_bytes(canonical_json_bytes(truth_document))
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)

    metrics = load_strict_json_object((candidate / "metrics.json").read_bytes())
    assert metrics["macro_exact_recovery"] == "0"
    assert metrics["wrong_displayed_names_all_tiers"] == 0


def test_candidate_p_value_is_the_conservative_maximum_of_both_required_endpoints(tmp_path: Path) -> None:
    truth, partitions, runs = _inputs(tmp_path)
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    extract_artifact_bundle(truth, partitions, runs, evidence)
    candidate = tmp_path / "candidate"
    detection = BootstrapInterval(Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 100), 10_000, 1)
    exact = BootstrapInterval(Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1, 25), 10_000, 1)

    with patch("vntyper.scripts.calibration_artifacts.paired_group_bootstrap", side_effect=(detection, exact)):
        fit_artifact_bundle(evidence, "lexicographic-safety-v1", candidate)

    evaluation = load_strict_json_object((candidate / "evaluation.json").read_bytes())
    assert evaluation["holm_adjusted_p_value"] == "1/25"


def _inputs(root: Path) -> tuple[Path, Path, Path]:
    keys = ("held", "select", "train", "validate")
    roles = ("locked-heldout", "policy-selection", "training", "validation")
    features = []
    labels = []
    rows = []
    members = []
    run_paths: dict[str, str] = {}
    for order, (key, role) in enumerate(zip(keys, roles, strict=True)):
        identity = f"identity-{key}"
        features.append(
            {
                "feature_key": f"f-{key}",
                "manifest_key": key,
                "features": {
                    "canonical_identity": identity,
                    "assay_class": "capture-short-read",
                    "haplotype_record_count": 2,
                    "haplotype_record_count_margin": 2,
                    "haplotype_record_share": 1.0,
                    "haplotype_record_share_margin": 1.0,
                    "haplotype_record_tie": False,
                    "xd_availability_fraction": 1.0,
                    "advntr_evidence_disposition": "admissible",
                },
            }
        )
        labels.append(
            {
                "label_key": f"l-{key}",
                "manifest_key": key,
                "truth_status": "mutated",
                "expected_identity": identity,
                "expected_display_name": "59dupC",
                "mutation_class": "duplication",
            }
        )
        rows.append(
            {
                "order": order,
                "name": "59dupC",
                "confidence": "High_Precision",
                "flag": "Not flagged",
                "tier": "A",
                "support": 10,
                "tie": False,
                "abstention": None,
            }
        )
        members.append(
            {
                "key": key,
                "role": role,
                "provenance": "external-custodian" if role == "locked-heldout" else "development",
                "groups": {
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
                },
            }
        )
        run = root / "runs" / key
        (run / "kestrel").mkdir(parents=True)
        (run / "pipeline_summary.json").write_bytes(canonical_json_bytes({"schema_version": 3}))
        (run / "kestrel" / "kestrel_pre_result.tsv").write_text("__Identity_Observation_Ordinal\n0\n", encoding="utf-8")
        (run / "kestrel" / "bam_identity_replay.v1.json").write_bytes(
            canonical_json_bytes({"schema_version": "bam-identity-replay-v1", "loci": []})
        )
        run_paths[key] = str(run)

    expected = {
        "aggregate": {"displayed": 4, "exact": 4, "wrong": 0, "control_findings": 0},
        "per_tier": {"A": {"displayed": 4, "exact": 4, "wrong": 0}},
        "rows": rows,
    }
    truth_value = {
        "schema_version": "calibration-truth-v1",
        "features": {"schema_version": "calibration-features-v1", "rows": features},
        "labels": {"schema_version": "calibration-labels-v1", "rows": labels},
        "baseline": {
            "schema_version": "calibration-baseline-replay-v1",
            "expected": expected,
            "observed": expected,
        },
    }
    study_value = {
        "schema_version": "calibration-study-v1",
        "protocol": {
            "objective": "lexicographic-safety-v1",
            "bootstrap_iterations": 10000,
            "bootstrap_interval": "percentile",
            "multiplicity_method": "holm",
            "seed": 295,
            "maximum_free_parameters": 4,
            "minimum_stratum_count": 1,
            "maximum_abstention_fraction": 0.25,
            "assay_classes": ["capture-short-read"],
            "mutation_classes": ["duplication"],
            "candidate_grid": {
                "minimum_record_count_margin": [1],
                "minimum_record_share": [0.5],
                "minimum_record_share_margin": [0.0],
                "xd_veto": ["disabled"],
            },
        },
        "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
    }
    runs_value = {"schema_version": "calibration-runs-v1", "runs": run_paths}
    truth_path = root / "truth.json"
    partitions_path = root / "partitions.json"
    runs_path = root / "runs.json"
    truth_path.write_bytes(canonical_json_bytes(truth_value))
    partitions_path.write_bytes(canonical_json_bytes(study_value))
    runs_path.write_bytes(canonical_json_bytes(runs_value))
    return truth_path, partitions_path, runs_path


def _run_cli(argv: list[str]) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    handle_calibrate(args, {}, parser, logging.INFO, None)
