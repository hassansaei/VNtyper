"""In-memory calibration extraction, fitting, and validation workflow."""

from fractions import Fraction
from pathlib import Path

import pytest

from tests.calibration_run_fixture import (
    ADVNTR_ARTIFACTS,
    CORE_ARTIFACTS,
    IDENTITY,
    OTHER_IDENTITY,
    refresh_run_hashes,
    write_schema_three_run,
)
from tests.unit.test_calibration_contract import synthetic_protocol
from vntyper.scripts.calibration_contract import CandidateMetrics
from vntyper.scripts.calibration_features import decode_label_artifact
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
        {
            "key": "held",
            "role": "locked-heldout",
            "provenance": "external-custodian",
            "assay_class": "capture-short-read",
            "groups": _groups("held"),
        },
        {
            "key": "select",
            "role": "policy-selection",
            "provenance": "development",
            "assay_class": "capture-short-read",
            "groups": _groups("select"),
        },
        {
            "key": "train",
            "role": "training",
            "provenance": "development",
            "assay_class": "capture-short-read",
            "groups": _groups("train"),
        },
        {
            "key": "validate",
            "role": "validation",
            "provenance": "development",
            "assay_class": "capture-short-read",
            "groups": _groups("validate"),
        },
    ]
    return decode_study_declaration(
        {
            "schema_version": "calibration-study-v1",
            "protocol": synthetic_protocol(),
            "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
        }
    )


def _labels(*, held_control: bool = False):
    return decode_label_artifact(
        {
            "schema_version": "calibration-labels-v1",
            "rows": [
                {
                    "label_key": f"l-{key}",
                    "manifest_key": key,
                    "truth_status": "control" if held_control and key == "held" else "mutated",
                    "expected_identity": None if held_control and key == "held" else f"identity-{key}",
                    "expected_display_name": None if held_control and key == "held" else "59dupC",
                    "mutation_class": "duplication",
                }
                for key in ("held", "select", "train", "validate")
            ],
        }
    )


def _runs(tmp_path: Path) -> dict[str, object]:
    runs: dict[str, object] = {}
    for key in ("held", "select", "train", "validate"):
        runs[key] = write_schema_three_run(tmp_path / key, key)
    return runs


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
    evidence = extract_evidence(_study(), _labels(), _runs(tmp_path))

    assert evidence.study_sha256 == _study().sha256
    assert evidence.features.rows[0].features["canonical_identity"] == ("MUC1-X-60-coding-v1|60|59|-|C")
    assert evidence.labels.sha256 == _labels().sha256
    assert evidence.baseline["schema_version"] == "calibration-baseline-replay-v1"
    assert tuple(evidence.run_hashes) == ("held", "select", "train", "validate")
    assert set(evidence.run_hashes["train"]) == set(CORE_ARTIFACTS)


@pytest.mark.parametrize(
    "relative",
    [
        "pipeline_summary.json",
        "provenance/decision_profile.json",
        "kestrel/kestrel_pre_result.tsv",
        "kestrel/bam_identity_replay.v1.json",
        "kestrel/kestrel_result.tsv",
    ],
)
def test_extract_refuses_old_runs_missing_complete_replay_artifacts(tmp_path: Path, relative: str) -> None:
    runs = _runs(tmp_path)
    train = runs["train"]
    assert isinstance(train, dict)
    (Path(str(train["root"])) / relative).unlink()

    with pytest.raises(ValueError, match="replay|artifact"):
        extract_evidence(_study(), _labels(), runs)


def test_extract_refuses_nearly_empty_schema_three_placeholders(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    train = runs["train"]
    assert isinstance(train, dict)
    root = Path(str(train["root"]))
    (root / "pipeline_summary.json").write_text('{"schema_version":3}\n', encoding="utf-8")
    (root / "kestrel" / "kestrel_pre_result.tsv").write_text("__Identity_Observation_Ordinal\n0\n", encoding="utf-8")
    refresh_run_hashes(train)

    with pytest.raises(ValueError, match="profile|identity|schema 3"):
        extract_evidence(_study(), _labels(), runs)


def test_extract_derives_feature_and_baseline_digests_from_verified_run_bytes(tmp_path: Path) -> None:
    first = extract_evidence(_study(), _labels(), _runs(tmp_path / "first"))
    changed_runs = _runs(tmp_path / "changed")
    changed_runs["train"] = write_schema_three_run(tmp_path / "replacement", "train", support=11)

    changed = extract_evidence(_study(), _labels(), changed_runs)

    assert changed.features.sha256 != first.features.sha256
    assert changed.dataset_sha256 != first.dataset_sha256


def test_extract_rejects_a_cross_member_run_root_swap(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    held = runs["held"]
    train = runs["train"]
    assert isinstance(held, dict) and isinstance(train, dict)
    held["root"], train["root"] = train["root"], held["root"]

    with pytest.raises(ValueError, match="SHA-256|artifact hash"):
        extract_evidence(_study(), _labels(), runs)


def test_extract_rejects_a_stale_final_result_even_when_its_declared_hash_is_updated(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    train = runs["train"]
    assert isinstance(train, dict)
    result = Path(str(train["root"])) / "kestrel" / "kestrel_result.tsv"
    result.write_text(result.read_text(encoding="utf-8").replace("59dupC", "60dupA"), encoding="utf-8")
    refresh_run_hashes(train)

    with pytest.raises(ValueError, match="summary.*final|final.*summary|parsed result"):
        extract_evidence(_study(), _labels(), runs)


@pytest.mark.parametrize("retain_pre_result", [False, True])
def test_extract_retains_a_no_finding_member_without_fabricating_a_baseline_call(
    tmp_path: Path, retain_pre_result: bool
) -> None:
    runs = _runs(tmp_path)
    runs["held"] = write_schema_three_run(
        tmp_path / "held-negative",
        "held",
        no_kestrel_finding=True,
        retain_pre_result=retain_pre_result,
    )
    evidence = extract_evidence(_study(), _labels(held_control=True), runs)

    held_feature = next(row for row in evidence.features.rows if row.manifest_key == "held")
    held_baseline = next(row for row in evidence.baseline["expected"]["rows"] if row["manifest_key"] == "held")
    assert held_feature.features["assay_class"] == "capture-short-read"
    assert held_feature.features.get("canonical_identity") == (IDENTITY if retain_pre_result else None)
    assert held_baseline == {
        "manifest_key": "held",
        "order": 0,
        "canonical_identity": None,
        "name": None,
        "confidence": None,
        "flag": None,
        "tier": None,
        "support": None,
        "tie": False,
        "abstention": None,
    }
    assert evidence.baseline["expected"]["aggregate"]["control_findings"] == 0


def test_extract_rejects_pre_result_diagnostics_that_disagree_with_identity_capture(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    train = runs["train"]
    assert isinstance(train, dict)
    pre_result = Path(str(train["root"])) / "kestrel" / "kestrel_pre_result.tsv"
    pre_result.write_text(
        pre_result.read_text(encoding="utf-8").replace("\tFalse\t", "\tTrue\t", 1),
        encoding="utf-8",
    )
    refresh_run_hashes(train)

    with pytest.raises(ValueError, match="diagnostics.*captured identity"):
        extract_evidence(_study(), _labels(), runs)


def test_extract_rejects_bam_identity_bindings_that_disagree_with_pre_result(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    train = runs["train"]
    assert isinstance(train, dict)
    replay = Path(str(train["root"])) / "kestrel" / "bam_identity_replay.v1.json"
    replay.write_bytes(replay.read_bytes().replace(IDENTITY.encode(), OTHER_IDENTITY.encode()))
    refresh_run_hashes(train)

    with pytest.raises(ValueError, match="BAM.*binding.*pre-result"):
        extract_evidence(_study(), _labels(), runs)


def test_extract_cross_binds_governed_advntr_evidence_and_support(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(tmp_path / "train-advntr", "train", with_advntr=True)

    evidence = extract_evidence(_study(), _labels(), runs)

    train = next(row for row in evidence.features.rows if row.manifest_key == "train")
    assert train.features["advntr_evidence_disposition"] == "admissible"
    assert train.features["advntr_sequencing_read_support"] == 7
    assert train.features["advntr_p_value"] == 0.001
    assert train.features["advntr_coverage"] == 42
    assert set(evidence.run_hashes["train"]) == {*CORE_ARTIFACTS, *ADVNTR_ARTIFACTS}

    train_run = runs["train"]
    assert isinstance(train_run, dict)
    snapshot = Path(str(train_run["root"])) / "provenance" / "advntr_artifact_evidence.json"
    snapshot.write_bytes(snapshot.read_bytes().replace(b'"schema_version"', b'"tampered_version"', 1))
    refresh_run_hashes(train_run)
    with pytest.raises(ValueError, match="adVNTR|artifact evidence|digest"):
        extract_evidence(_study(), _labels(), runs)


def test_fit_reads_only_training_and_policy_selection_and_requires_baseline_replay(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _labels(), _runs(tmp_path))
    candidate = fit_candidate(evidence, objective="lexicographic-safety-v1", evaluator=_evaluate)

    assert candidate.profile.profile_kind == "generated"
    assert candidate.accessed_roles == ("training", "policy-selection")
    assert candidate.baseline_reproduced
    assert candidate.evaluation.holm_adjusted_p_value == Fraction(3, 250)

    baseline = dict(evidence.baseline)
    baseline["observed"] = {"aggregate": {}, "per_tier": {}, "rows": []}
    failed = type(evidence)(
        evidence.study,
        evidence.features,
        evidence.labels,
        baseline,
        evidence.run_hashes,
        evidence.study_sha256,
        evidence.dataset_sha256,
    )
    with pytest.raises(ValueError, match="baseline"):
        fit_candidate(failed, objective="lexicographic-safety-v1", evaluator=_evaluate)


def test_fit_objective_must_match_snapshotted_protocol(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _labels(), _runs(tmp_path))

    with pytest.raises(ValueError, match="objective"):
        fit_candidate(evidence, objective="f1", evaluator=_evaluate)


def test_fit_applies_holm_to_real_marginal_bootstrap_evidence_across_the_family(tmp_path: Path) -> None:
    evidence = extract_evidence(_study(), _labels(), _runs(tmp_path))
    rows = tuple(
        PairedObservation(
            f"group-{index:02d}",
            "capture-short-read:duplication",
            Fraction(index >= 13),
            Fraction(index < 13),
        )
        for index in range(20)
    )
    interval = paired_group_bootstrap(
        rows,
        required_strata=("capture-short-read:duplication",),
        iterations=10_000,
        seed=295,
    )

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
    evidence = extract_evidence(_study(), _labels(), _runs(tmp_path))
    candidate = fit_candidate(evidence, objective="lexicographic-safety-v1", evaluator=_evaluate)
    attestation = validate_candidate(candidate.profile, evidence)

    assert attestation.role == "validation"
    assert attestation.profile_sha256 == candidate.profile.digest
    assert attestation.accessed_roles == ("validation",)
