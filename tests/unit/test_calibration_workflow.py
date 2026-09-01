"""In-memory calibration extraction, fitting, and validation workflow."""

from copy import deepcopy
from fractions import Fraction
from pathlib import Path

import pytest

from tests.unit.test_calibration_contract import synthetic_protocol
from vntyper.scripts.calibration_contract import CandidateMetrics
from vntyper.scripts.calibration_features import decode_feature_artifact, decode_label_artifact
from vntyper.scripts.calibration_manifest import decode_study_declaration
from vntyper.scripts.calibration_objective import CandidateEvaluation, count_free_parameters
from vntyper.scripts.calibration_statistics import PairedObservation, paired_group_bootstrap
from vntyper.scripts.calibration_workflow import extract_evidence, fit_candidate, validate_candidate

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


def _study():
    members = [
        {"key": "held", "role": "locked-heldout", "provenance": "external-custodian", "groups": _groups("held")},
        {"key": "select", "role": "policy-selection", "provenance": "development", "groups": _groups("select")},
        {"key": "train", "role": "training", "provenance": "development", "groups": _groups("train")},
        {"key": "validate", "role": "validation", "provenance": "development", "groups": _groups("validate")},
    ]
    return decode_study_declaration(
        {
            "schema_version": "calibration-study-v1",
            "protocol": synthetic_protocol(),
            "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
        }
    )


def _features():
    return decode_feature_artifact(
        {
            "schema_version": "calibration-features-v1",
            "rows": [
                {
                    "feature_key": f"f-{key}",
                    "manifest_key": key,
                    "features": {"assay_class": "capture-short-read", "haplotype_record_count": 2},
                }
                for key in ("held", "select", "train", "validate")
            ],
        }
    )


def _labels():
    return decode_label_artifact(
        {
            "schema_version": "calibration-labels-v1",
            "rows": [
                {
                    "label_key": f"l-{key}",
                    "manifest_key": key,
                    "truth_status": "mutated",
                    "expected_identity": f"identity-{key}",
                    "expected_display_name": "59dupC",
                    "mutation_class": "duplication",
                }
                for key in ("held", "select", "train", "validate")
            ],
        }
    )


def _runs(tmp_path: Path) -> dict[str, Path]:
    runs: dict[str, Path] = {}
    for key in ("held", "select", "train", "validate"):
        run = tmp_path / key
        (run / "kestrel").mkdir(parents=True, exist_ok=True)
        (run / "pipeline_summary.json").write_text('{"schema_version":3}\n', encoding="utf-8")
        (run / "kestrel" / "kestrel_pre_result.tsv").write_text("__Identity_Observation_Ordinal\n0\n", encoding="utf-8")
        (run / "kestrel" / "bam_identity_replay.v1.json").write_text(
            '{"schema_version":"bam-identity-replay-v1","loci":[]}\n', encoding="utf-8"
        )
        runs[key] = run
    return runs


def _baseline(*, reproduced: bool = True) -> dict[str, object]:
    expected = {
        "aggregate": {"displayed": 154, "exact": 136, "wrong": 18, "control_findings": 0},
        "per_tier": {"A": {"displayed": 53, "exact": 53, "wrong": 0}},
        "rows": [
            {
                "order": 0,
                "name": "59dupC",
                "confidence": "High_Precision",
                "flag": "Not flagged",
                "tier": "A",
                "support": 10,
                "tie": False,
                "abstention": None,
            }
        ],
    }
    observed = deepcopy(expected)
    if not reproduced:
        aggregate = observed["aggregate"]
        assert isinstance(aggregate, dict)
        aggregate["exact"] = 135
    return {"schema_version": "calibration-baseline-replay-v1", "expected": expected, "observed": observed}


def _evaluate(profile) -> CandidateEvaluation:
    free_parameters = count_free_parameters(profile.components["dominance"])
    metrics = CandidateMetrics(
        profile.digest,
        0,
        0,
        0,
        Fraction(3, 4),
        Fraction(4, 5),
        free_parameters,
        Fraction(1, 10),
        True,
        True,
        "a" * 64,
    )
    return CandidateEvaluation(metrics, Fraction(0), Fraction(0), (2,), Fraction(1, 1000))


def test_extract_snapshots_complete_study_features_labels_baseline_and_run_hashes(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _features(), _labels(), _runs(tmp_path), baseline=_baseline())

    assert evidence.study_sha256 == _study().sha256
    assert evidence.features.sha256 == _features().sha256
    assert evidence.labels.sha256 == _labels().sha256
    assert evidence.baseline["schema_version"] == "calibration-baseline-replay-v1"
    assert tuple(evidence.run_hashes) == ("held", "select", "train", "validate")


@pytest.mark.parametrize(
    "relative",
    [
        "pipeline_summary.json",
        "kestrel/kestrel_pre_result.tsv",
        "kestrel/bam_identity_replay.v1.json",
    ],
)
def test_extract_refuses_old_runs_missing_complete_replay_artifacts(tmp_path: Path, relative: str) -> None:
    runs = _runs(tmp_path)
    (runs["train"] / relative).unlink()

    with pytest.raises(ValueError, match="replay|artifact"):
        extract_evidence(_study(), _features(), _labels(), runs, baseline=_baseline())


def test_fit_reads_only_training_and_policy_selection_and_requires_baseline_replay(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _features(), _labels(), _runs(tmp_path), baseline=_baseline())
    candidate = fit_candidate(evidence, objective="lexicographic-safety-v1", evaluator=_evaluate)

    assert candidate.profile.profile_kind == "generated"
    assert candidate.accessed_roles == ("training", "policy-selection")
    assert candidate.baseline_reproduced
    assert candidate.evaluation.holm_adjusted_p_value == Fraction(3, 250)

    failed = extract_evidence(_study(), _features(), _labels(), _runs(tmp_path), baseline=_baseline(reproduced=False))
    with pytest.raises(ValueError, match="baseline"):
        fit_candidate(failed, objective="lexicographic-safety-v1", evaluator=_evaluate)


def test_fit_objective_must_match_snapshotted_protocol(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _features(), _labels(), _runs(tmp_path), baseline=_baseline())

    with pytest.raises(ValueError, match="objective"):
        fit_candidate(evidence, objective="f1", evaluator=_evaluate)


def test_fit_applies_holm_to_real_marginal_bootstrap_evidence_across_the_family(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _features(), _labels(), _runs(tmp_path), baseline=_baseline())
    rows = tuple(
        PairedObservation(
            f"group-{index:02d}",
            "capture-short-read:duplication",
            Fraction(index >= 13),
            Fraction(index < 13),
        )
        for index in range(20)
    )
    interval = paired_group_bootstrap(rows, iterations=10_000, seed=295)

    def marginal(profile) -> CandidateEvaluation:
        evaluation = _evaluate(profile)
        return CandidateEvaluation(
            evaluation.metrics,
            interval.one_sided_lower,
            interval.one_sided_lower,
            evaluation.stratum_counts,
            interval.one_sided_noninferiority_p_value,
        )

    assert interval.one_sided_noninferiority_p_value <= Fraction(1, 20)
    with pytest.raises(ValueError, match="no admissible candidate"):
        fit_candidate(evidence, objective="lexicographic-safety-v1", evaluator=marginal)


def test_validate_cannot_select_another_profile_or_open_heldout(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _features(), _labels(), _runs(tmp_path), baseline=_baseline())
    candidate = fit_candidate(evidence, objective="lexicographic-safety-v1", evaluator=_evaluate)
    attestation = validate_candidate(candidate.profile, evidence)

    assert attestation.role == "validation"
    assert attestation.profile_sha256 == candidate.profile.digest
    assert attestation.accessed_roles == ("validation",)
