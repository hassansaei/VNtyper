"""In-memory calibration extraction, fitting, and validation workflow."""

import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import cast

import pandas as pd
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
from vntyper.scripts import calibration_artifacts
from vntyper.scripts.calibration_contract import CandidateMetrics
from vntyper.scripts.calibration_features import decode_label_artifact
from vntyper.scripts.calibration_manifest import decode_study_declaration
from vntyper.scripts.calibration_objective import CandidateEvaluation, count_free_parameters
from vntyper.scripts.calibration_profiles import build_generated_profile
from vntyper.scripts.calibration_run_projection import with_complete_kestrel_candidate_projections
from vntyper.scripts.calibration_statistics import PairedObservation, paired_group_bootstrap
from vntyper.scripts.calibration_workflow import extract_evidence, fit_candidate, validate_candidate
from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.decision_profile import load_packaged_decision_profile, parse_decision_profile
from vntyper.scripts.molecular_identity import parse_molecular_identity
from vntyper.scripts.nomenclature_annotate import DominanceSeamOutcome, reconcile_caller_outputs
from vntyper.scripts.nomenclature_bam_evidence import BamIdentityEvidence, BamLocusEvidence
from vntyper.scripts.nomenclature_bam_replay import BamReplayArtifact, BamReplayLocus, write_bam_replay_artifact
from vntyper.scripts.profile_provenance import profile_summary_fields
from vntyper.scripts.run_configuration import resolve_run_configuration
from vntyper.version import __version__

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


def _study(*, protocol: dict[str, object] | None = None):
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
            "protocol": protocol or synthetic_protocol(),
            "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
        }
    )


def _labels(
    *,
    held_control: bool = False,
    expected_identity: str | None = None,
    expected_display_name: str = "59dupC",
):
    return decode_label_artifact(
        {
            "schema_version": "calibration-labels-v1",
            "rows": [
                {
                    "label_key": f"l-{key}",
                    "manifest_key": key,
                    "truth_status": "control" if held_control and key == "held" else "mutated",
                    "expected_identity": (
                        None if held_control and key == "held" else expected_identity or f"identity-{key}"
                    ),
                    "expected_display_name": None if held_control and key == "held" else expected_display_name,
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


def test_complete_candidate_rendering_projection_rejects_a_missing_shipped_map() -> None:
    with pytest.raises(ValueError, match="identity projection"):
        with_complete_kestrel_candidate_projections({}, ())


def _producer_shaped_ab_run(root: Path, key: str, profile) -> dict[str, object]:
    """Write complete A/B pre-selection and deduplicated replay before extraction."""
    run = write_schema_three_run(root, key)
    kestrel_root = root / "kestrel"
    pre_path = kestrel_root / "kestrel_pre_result.tsv"
    pre = pd.read_csv(pre_path, sep="\t", dtype=str, keep_default_na=False)
    selection = {
        "__Identity_Equivalent_Representation_Count": "1",
        "__Identity_Hypothesis_Count": "2",
        "__Identity_Group_Blocking_Gates": "[]",
        "__Identity_Group_Flags": "[]",
        "__Identity_Group_Context_Diverges": "false",
    }
    for column, value in selection.items():
        pre[column] = value
    pre["__Identity_Selected_Raw_Representation_Key"] = pre["__Identity_Raw_Representation_Key"]
    pre["__Identity_Selected_Observation_Ordinal"] = pre["__Identity_Observation_Ordinal"]
    pre["Nomenclature"] = "59dupC"
    second = pre.loc[0].copy()
    raw_key = json.dumps(
        {"source": "kestrel", "values": ["X-X", 62, "G", "GC"]},
        separators=(",", ":"),
        sort_keys=True,
    )
    second.update(
        {
            "POS": "62",
            "ALT": "GC",
            "Nomenclature": "58dupG",
            "Molecular_Identity": OTHER_IDENTITY,
            "__Identity_Raw_Representation_Key": raw_key,
            "__Identity_Molecular_Identity": OTHER_IDENTITY,
            "__Identity_Observation_Ordinal": "1",
            "__Identity_Selected_Raw_Representation_Key": raw_key,
            "__Identity_Selected_Observation_Ordinal": "1",
        }
    )
    pd.DataFrame([pre.loc[0], second]).to_csv(pre_path, sep="\t", index=False)

    identity_a = parse_molecular_identity(IDENTITY)
    identity_b = parse_molecular_identity(OTHER_IDENTITY)
    records = (
        BamIdentityEvidence((identity_a,), ((0,),), 7),
        BamIdentityEvidence((identity_b,), ((1,),), 8),
        BamIdentityEvidence((identity_b,), ((1,),), 12),
        BamIdentityEvidence((identity_b,), ((1,),), 16),
    )
    write_bam_replay_artifact(
        kestrel_root,
        BamReplayArtifact(
            (BamReplayLocus((0, 1), "observed", BamLocusEvidence(records, 4, {identity_a: 1, identity_b: 3})),)
        ),
    )

    (root / "provenance" / "decision_profile.json").write_bytes(profile.canonical_bytes)
    summary_path = root / "pipeline_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary.update(profile_summary_fields(profile))
    result = kestrel_root / "kestrel_result.tsv"
    result_frame = pd.read_csv(result, sep="\t", dtype=str, keep_default_na=False)
    result_frame["__Identity_Hypothesis_Count"] = "2"
    result_frame["Molecular_Identity_Status"] = "legacy-selected-among-multiple"
    result_frame["Identity_Hypothesis_Count"] = "2"
    result_frame.to_csv(result, sep="\t", index=False)
    summary["steps"][0]["md5sum"] = hashlib.md5(result.read_bytes()).hexdigest()
    summary["steps"][0]["parsed_result"]["data"] = result_frame.to_dict(orient="records")
    summary_path.write_bytes(canonical_json_bytes(summary))
    refresh_run_hashes(run)
    return run


def _explicit_custom_dominance_profile(*, enabled: bool):
    """Build a valid explicit-custom profile with the requested dominance state."""
    packaged = load_packaged_decision_profile()
    document = json.loads(packaged.canonical_bytes)
    document.update(
        {
            "profile_id": f"test-explicit-dominance-{str(enabled).lower()}",
            "profile_kind": "explicit-custom",
        }
    )
    document["inventory"]["/components/dominance/enabled"]["value"] = enabled
    return parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def _bind_profile_to_run(run: dict[str, object], profile) -> None:
    """Replace a fixture's snapshotted profile and refresh its closed hashes."""
    root = Path(str(run["root"]))
    (root / "provenance" / "decision_profile.json").write_bytes(profile.canonical_bytes)
    summary_path = root / "pipeline_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary.update(profile_summary_fields(profile))
    summary_path.write_bytes(canonical_json_bytes(summary))
    refresh_run_hashes(run)


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


def test_extract_rejects_a_complete_cross_member_declaration_swap(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["held"], runs["train"] = runs["train"], runs["held"]

    with pytest.raises(ValueError, match="member key|manifest key"):
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
    expected = evidence.baseline["expected"]
    assert isinstance(expected, Mapping)
    expected_rows = cast(Sequence[Mapping[str, object]], expected["rows"])
    held_baseline = next(row for row in expected_rows if row["manifest_key"] == "held")
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
        "identity_projection": {},
    }
    aggregate = expected["aggregate"]
    assert isinstance(aggregate, Mapping)
    assert aggregate["control_findings"] == 0


def test_extract_accepts_replay_for_only_the_retained_final_candidate(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(tmp_path / "two-pre-one-final", "train", extra_pre_result=True)

    evidence = extract_evidence(_study(), _labels(), runs)

    train = next(row for row in evidence.features.rows if row.manifest_key == "train")
    assert train.features["canonical_identity"] == IDENTITY
    assert train.features["cooccurring_identity_count"] == 2


def test_extract_generated_profile_with_disabled_dominance_uses_only_the_retained_candidate(tmp_path: Path) -> None:
    packaged = load_packaged_decision_profile()
    packaged_dominance = packaged.components["dominance"]
    assert isinstance(packaged_dominance, Mapping)
    dominance = {**packaged_dominance, "enabled": False}
    profile = build_generated_profile(
        dominance,
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=295,
        objective="lexicographic-safety-v1",
        generator_version=__version__,
    )
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(tmp_path / "generated-disabled", "train", extra_pre_result=True)
    train_run = runs["train"]
    assert isinstance(train_run, dict)
    _bind_profile_to_run(train_run, profile)

    evidence = extract_evidence(_study(), _labels(), runs)

    expected = evidence.baseline["expected"]
    assert isinstance(expected, Mapping)
    rows = cast(Sequence[Mapping[str, object]], expected["rows"])
    train = next(row for row in rows if row["manifest_key"] == "train")
    assert train["identity_projection"] == {IDENTITY: {"name": "59dupC", "tier": "A"}}


def test_extract_explicit_profile_with_enabled_dominance_uses_complete_candidate_universe(tmp_path: Path) -> None:
    profile = _explicit_custom_dominance_profile(enabled=True)
    runs = _runs(tmp_path)
    runs["train"] = _producer_shaped_ab_run(tmp_path / "explicit-enabled", "train", profile)

    evidence = extract_evidence(_study(), _labels(), runs)

    expected = evidence.baseline["expected"]
    assert isinstance(expected, Mapping)
    rows = cast(Sequence[Mapping[str, object]], expected["rows"])
    train = next(row for row in rows if row["manifest_key"] == "train")
    assert train["identity_projection"] == {
        IDENTITY: {"name": "59dupC", "tier": "A"},
        OTHER_IDENTITY: {"name": "58dupG", "tier": "A"},
    }


def test_extract_rejects_replay_for_a_stale_nonretained_pre_result_ordinal(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(tmp_path / "two-pre-one-final", "train", extra_pre_result=True)
    train = runs["train"]
    assert isinstance(train, dict)
    replay = Path(str(train["root"])) / "kestrel" / "bam_identity_replay.v1.json"
    stale = replay.read_bytes().replace(
        b'"candidate_observation_ordinals":[0]', b'"candidate_observation_ordinals":[1]'
    )
    replay.write_bytes(
        stale.replace(b'"candidate_observation_ordinals":[[0]]', b'"candidate_observation_ordinals":[[1]]')
    )
    refresh_run_hashes(train)

    with pytest.raises(ValueError, match="BAM replay.*ordinals|cover.*ordinals"):
        extract_evidence(_study(), _labels(), runs)


def test_extract_uses_advntr_when_kestrel_is_the_exact_negative_placeholder(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(
        tmp_path / "advntr-only",
        "train",
        with_advntr=True,
        no_kestrel_finding=True,
        advntr_identity=OTHER_IDENTITY,
        advntr_name="58dupG",
        selected_name="58dupG",
        reconciled_identity=OTHER_IDENTITY,
    )

    evidence = extract_evidence(_study(), _labels(), runs)

    expected = evidence.baseline["expected"]
    assert isinstance(expected, Mapping)
    expected_rows = cast(Sequence[Mapping[str, object]], expected["rows"])
    row = next(row for row in expected_rows if row["manifest_key"] == "train")
    assert row["canonical_identity"] == OTHER_IDENTITY
    assert row["name"] == "58dupG"
    assert row["support"] == 7


def test_extract_uses_the_reconciled_advntr_selected_whole_locus_verdict(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(
        tmp_path / "advntr-selected",
        "train",
        with_advntr=True,
        advntr_identity=OTHER_IDENTITY,
        advntr_name="58dupG",
        selected_name="58dupG",
        reconciled_identity=OTHER_IDENTITY,
    )

    evidence = extract_evidence(_study(), _labels(), runs)

    expected = evidence.baseline["expected"]
    assert isinstance(expected, Mapping)
    expected_rows = cast(Sequence[Mapping[str, object]], expected["rows"])
    row = next(row for row in expected_rows if row["manifest_key"] == "train")
    assert row["canonical_identity"] == OTHER_IDENTITY
    assert row["name"] == "58dupG"
    assert row["identity_projection"] == {
        IDENTITY: {"name": "59dupC", "tier": "A"},
        OTHER_IDENTITY: {"name": "58dupG", "tier": "A"},
    }


def test_extract_never_selects_identity_from_equal_display_names(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    runs["train"] = write_schema_three_run(
        tmp_path / "same-display-different-identity",
        "train",
        with_advntr=True,
        advntr_identity=OTHER_IDENTITY,
        advntr_name="59dupC",
        reconciled_identity=OTHER_IDENTITY,
    )

    evidence = extract_evidence(_study(), _labels(), runs)

    expected = evidence.baseline["expected"]
    assert isinstance(expected, Mapping)
    expected_rows = cast(Sequence[Mapping[str, object]], expected["rows"])
    row = next(row for row in expected_rows if row["manifest_key"] == "train")
    assert row["canonical_identity"] == OTHER_IDENTITY
    assert row["name"] == "59dupC"


def test_extract_binds_raw_profile_bytes_to_the_summary_digest(tmp_path: Path) -> None:
    runs = _runs(tmp_path)
    train = runs["train"]
    assert isinstance(train, dict)
    profile = Path(str(train["root"])) / "provenance" / "decision_profile.json"
    profile.write_bytes(profile.read_bytes() + b" \n")
    refresh_run_hashes(train)

    with pytest.raises(ValueError, match="profile.*SHA-256|SHA-256.*profile"):
        extract_evidence(_study(), _labels(), runs)


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
        evidence.profile_dataset_sha256,
    )
    with pytest.raises(ValueError, match="baseline"):
        fit_candidate(failed, objective="lexicographic-safety-v1", evaluator=_evaluate)


def test_fitted_profile_resolves_and_drives_production_whole_locus_reconciliation(tmp_path: Path) -> None:
    """Producer A/B evidence drives real fitting and then the production identity change."""
    protocol = synthetic_protocol()
    protocol["candidate_grid"] = {
        "minimum_record_count_margin": [2],
        "minimum_record_share": [0.75],
        "minimum_record_share_margin": [0.25],
        "xd_veto": ["disabled"],
    }
    study = _study(protocol=protocol)
    expected_component = {
        "enabled": True,
        "minimum_record_count_margin": 2,
        "minimum_record_share": 0.75,
        "minimum_record_share_margin": 0.25,
        "xd_veto": "disabled",
        "abstain_on_inadmissible_advntr": False,
    }
    producer_profile = build_generated_profile(
        expected_component,
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash=study.partitions.sha256,
        seed=study.protocol.seed,
        objective=study.protocol.objective,
        generator_version=__version__,
    )
    runs = {
        key: _producer_shaped_ab_run(tmp_path / "runs" / key, key, producer_profile)
        for key in ("held", "select", "train", "validate")
    }
    evidence = extract_evidence(
        study,
        _labels(expected_identity=OTHER_IDENTITY, expected_display_name="58dupG"),
        runs,
    )
    assert evidence.features.rows[0].features["canonical_identity"] == OTHER_IDENTITY
    assert evidence.features.rows[0].features["haplotype_record_count_margin"] == 2
    assert evidence.features.rows[0].features["haplotype_record_share"] == 0.75

    candidate = fit_candidate(
        evidence,
        objective="lexicographic-safety-v1",
        evaluator=lambda profile: calibration_artifacts._evaluate(profile, evidence),
    )
    fitted_dominance = candidate.profile.components["dominance"]
    assert isinstance(fitted_dominance, Mapping)
    assert dict(fitted_dominance) == expected_component
    profile_path = tmp_path / "fitted-profile.json"
    profile_path.write_bytes(candidate.profile.canonical_bytes)

    resolved = resolve_run_configuration(profile_path)
    train = runs["train"]
    assert isinstance(train, dict)
    root = Path(str(train["root"]))
    kestrel = root / "kestrel" / "kestrel_result.tsv"
    before = pd.read_csv(kestrel, sep="\t", dtype=str, keep_default_na=False)
    assert before.loc[0, "Nomenclature"] == "59dupC"

    assert resolved.dominance["enabled"] is True
    outcome = reconcile_caller_outputs(
        kestrel,
        None,
        resolved_component=resolved.nomenclature,
        dominance_component=resolved.dominance,
        custom_context_active=True,
    )

    written = pd.read_csv(kestrel, sep="\t", dtype=str, keep_default_na=False)
    assert isinstance(outcome, DominanceSeamOutcome)
    assert outcome.dominance_outcome == "selected"
    assert written.loc[0, "__Reconciled_Molecular_Identity"] == OTHER_IDENTITY
    assert written.loc[0, "Nomenclature"] == "58dupG"
    assert written.loc[0, "Nomenclature"] != before.loc[0, "Nomenclature"]
    assert written.loc[0, "__Dominance_Abstention_Reason"] == ""


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


def test_generated_profile_cannot_cross_studies_with_the_same_candidate_parameters(tmp_path: Path) -> None:
    study_a = _study()
    evidence_a = extract_evidence(study_a, _labels(), _runs(tmp_path / "a"))
    candidate = fit_candidate(evidence_a, objective="lexicographic-safety-v1", evaluator=_evaluate)
    protocol_b = synthetic_protocol()
    protocol_b["seed"] = 296
    evidence_b = extract_evidence(_study(protocol=protocol_b), _labels(), _runs(tmp_path / "b"))

    with pytest.raises(ValueError, match="study|protocol|dataset|seed|partition"):
        validate_candidate(candidate.profile, evidence_b)
