"""Length command adapters consume only committed, authorized role evidence."""

import hashlib
from dataclasses import replace
from importlib import import_module
from unittest.mock import patch

import pytest

from tests.unit.test_calibration_role_source import source_fixture
from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_calibration_target_runs import run_document
from tests.unit.test_length_estimation import features
from vntyper.scripts.calibration_candidate import candidate_producer_document
from vntyper.scripts.calibration_exposure import decode_exposure_receipt
from vntyper.scripts.calibration_role_source import decode_role_source
from vntyper.scripts.calibration_target_contract import decode_target_study
from vntyper.scripts.calibration_target_runs import decode_target_runs
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256
from vntyper.scripts.length_feature_artifact import encode_length_feature_artifact
from vntyper.scripts.length_model import TARGET_BOUNDARY_DEFINITION

pytestmark = pytest.mark.unit


def prepared_source(tmp_path):
    module = import_module("vntyper.scripts.calibration_length_controller")
    measured = features(manifest_key="artifact-0")
    raw_study = study_document()
    raw_study["baseline"]["annotation_sha256"] = measured.annotation_sha256
    raw_study["baseline"]["counting_policy_sha256"] = measured.counting_policy_sha256
    raw_study["applicability"]["counting_policy_sha256"] = measured.counting_policy_sha256
    study = decode_target_study(raw_study)
    feature_path = tmp_path / "features.json"
    feature_raw = canonical_json_bytes(encode_length_feature_artifact(measured))
    feature_path.write_bytes(feature_raw)
    raw_runs = run_document()
    raw_run = raw_runs["runs"][0]
    raw_run.update(
        manifest_key="artifact-0",
        input_sha256=measured.provenance.measurement_context.input_sha256,
        policy_sha256=module.length_measurement_policy_sha256(study),
        capture_policy_sha256=module.length_measurement_policy_sha256(study),
        baseline_assets_sha256=study.baseline.sha256,
        producer_sha256=canonical_sha256(candidate_producer_document(study.baseline.producer)),
    )
    raw_run["assets"]["length_features"] = {
        "path": str(feature_path),
        "sha256": hashlib.sha256(feature_raw).hexdigest(),
        "size_bytes": len(feature_raw),
    }
    runs = decode_target_runs(raw_runs)
    truth_path = tmp_path / "truth.json"
    truth_raw = canonical_json_bytes(
        {
            "schema_version": "calibration-length-truth-v1",
            "boundary_definition": TARGET_BOUNDARY_DEFINITION,
            "rows": [{"key": "artifact-0", "total_repeat_count": 110}],
        }
    )
    truth_path.write_bytes(truth_raw)
    _, _, raw_source = source_fixture()
    raw_source.update(study_sha256=study.sha256, run_manifest_sha256=runs.sha256)
    raw_source["truth_asset"] = {
        "path": str(truth_path),
        "sha256": hashlib.sha256(truth_raw).hexdigest(),
        "size_bytes": len(truth_raw),
    }
    raw_source["identities_by_key"]["artifact-0"][0]["sha256"] = raw_run["input_sha256"]
    source = decode_role_source(raw_source, study=study, runs=runs, expected_role="training")
    receipt = decode_exposure_receipt(
        {
            "schema_version": "calibration-target-exposure-receipt-v2",
            "target": "length",
            "role": "training",
            "study_sha256": study.sha256,
            "partition_sha256": study.partitions.sha256,
            "evidence_sha256": source.sha256,
            "membership_sha256": canonical_sha256(
                [{"namespace": name, "sha256": digest} for name, digest in source.identities]
            ),
            "exposure_ledger_id": study.exposure_ledger_id,
            "sequence": 1,
        }
    )
    return study, runs, source, receipt, measured


def test_role_rows_come_from_verified_feature_bodies_and_truth(tmp_path):
    module = import_module("vntyper.scripts.calibration_length_controller")
    study, runs, source, receipt, measured = prepared_source(tmp_path)
    rows = module.load_length_source_rows(study, runs, source, receipt)
    assert len(rows) == 1
    assert rows[0].key == "artifact-0"
    assert rows[0].features == measured
    assert rows[0].total_truth_repeat_count == 110
    assert rows[0].role == "training"


def test_wrong_exposure_receipt_precedes_all_asset_reads(tmp_path):
    module = import_module("vntyper.scripts.calibration_length_controller")
    study, runs, source, receipt, _ = prepared_source(tmp_path)
    with (
        patch.object(module, "read_target_json", side_effect=AssertionError("opened")) as reader,
        pytest.raises(ValueError),
    ):
        module.load_length_source_rows(study, runs, source, replace(receipt, evidence_sha256="0" * 64))
    reader.assert_not_called()

    receipt_raw = module.exposure_receipt_document(receipt)
    receipt_raw["evidence_sha256"] = "0" * 64
    independently_valid = decode_exposure_receipt(receipt_raw)
    with (
        patch.object(module, "read_target_json", side_effect=AssertionError("opened")) as reader,
        pytest.raises(ValueError, match="authorized source"),
    ):
        module.load_length_source_rows(study, runs, source, independently_valid)
    reader.assert_not_called()


def test_modified_truth_file_is_not_used(tmp_path):
    module = import_module("vntyper.scripts.calibration_length_controller")
    study, runs, source, receipt, _ = prepared_source(tmp_path)
    source.truth_asset.path.write_bytes(source.truth_asset.path.read_bytes().replace(b"110", b"999"))
    with pytest.raises(ValueError, match="digest"):
        module.load_length_source_rows(study, runs, source, receipt)


def test_opened_locked_truth_is_verified_and_never_reopened(tmp_path):
    module = import_module("vntyper.scripts.calibration_length_controller")
    study, runs, source, receipt, measured = prepared_source(tmp_path)
    opened = source.truth_asset.path.read_bytes()
    original_reader = module.read_target_json

    def guarded_reader(asset):
        if asset == source.truth_asset:
            raise AssertionError("reopened truth")
        return original_reader(asset)

    with patch.object(module, "read_target_json", side_effect=guarded_reader) as reader:
        rows = module.load_length_source_rows(study, runs, source, receipt, opened_truth=opened)
    reader.assert_called_once()  # The feature body is still a separately committed asset.
    assert rows[0].features == measured
    with pytest.raises(ValueError, match="opened length truth"):
        module.load_length_source_rows(study, runs, source, receipt, opened_truth=opened + b" ")


@pytest.mark.parametrize(
    "value,keys",
    [
        ({}, ("k",)),
        (
            {"schema_version": "calibration-length-truth-v1", "boundary_definition": "wrong", "rows": []},
            (),
        ),
        (
            {
                "schema_version": "calibration-length-truth-v1",
                "boundary_definition": TARGET_BOUNDARY_DEFINITION,
                "rows": {},
            },
            (),
        ),
        (
            {
                "schema_version": "calibration-length-truth-v1",
                "boundary_definition": TARGET_BOUNDARY_DEFINITION,
                "rows": [{"key": "k", "total_repeat_count": True}],
            },
            ("k",),
        ),
        (
            {
                "schema_version": "calibration-length-truth-v1",
                "boundary_definition": TARGET_BOUNDARY_DEFINITION,
                "rows": [{"key": "k", "total_repeat_count": 1.5}],
            },
            ("k",),
        ),
        (
            {
                "schema_version": "calibration-length-truth-v1",
                "boundary_definition": TARGET_BOUNDARY_DEFINITION,
                "rows": [{"key": "k", "total_repeat_count": 2}],
            },
            ("different",),
        ),
    ],
)
def test_length_truth_decoder_rejects_malformed_nonintegral_or_wrong_roster(value, keys):
    module = import_module("vntyper.scripts.calibration_length_controller")
    with pytest.raises(ValueError):
        module._truth(value, keys)


def fit_fixture(tmp_path):
    from argparse import Namespace

    from tests.unit.test_calibration_length import metadata_document, roster_for, training_row
    from tests.unit.test_calibration_length_evaluation import evaluation_protocol
    from tests.unit.test_calibration_manifest import _manifest, _member
    from tests.unit.test_length_features import _annotation, _annotation_raw, _context
    from vntyper.scripts.calibration_artifact_io import write_json
    from vntyper.scripts.calibration_exposure_io import initialize_exposure_ledger
    from vntyper.scripts.calibration_length_metrics import decode_length_eligible_roster
    from vntyper.scripts.calibration_length_protocol import length_protocol_document
    from vntyper.scripts.calibration_manifest import connected_leakage_groups
    from vntyper.scripts.length_features import DepthPosition, extract_length_features

    module = import_module("vntyper.scripts.calibration_length_controller")
    root = tmp_path / "evidence"
    root.mkdir()
    ledger = tmp_path / "exposure.jsonl"
    ledger_id = initialize_exposure_ledger(ledger)
    values = {"t1": 1, "t2": 2, "t3": 3, "s1": 1, "s2": 3, "v1": 1, "v2": 3, "l1": 1, "l2": 3}
    annotations = _annotation()
    measured = {}
    for key, value in values.items():
        context = _context(annotations, manifest_key=key, input_sha256=canonical_sha256(key))
        depth = tuple(
            DepthPosition(annotations.contig, position, count, tuple(f"pair-{i}" for i in range(100)))
            for position, count in enumerate((100, 100 * value, 100 * value, 100, 100, 100))
        )
        measured[key] = extract_length_features(depth, annotations, context)
    raw_study = study_document()
    metadata_rows = tuple(training_row(key, values[key], 10 + 50 * values[key]) for key in ("t1", "t2", "t3"))
    metadata = metadata_document(metadata_rows, roster_for(metadata_rows))
    raw_study["applicability"] = metadata["applicability"]
    raw_study["baseline"].update(
        {
            key: metadata[key]
            for key in ("annotation_sha256", "counting_policy_sha256", "producer", "maximum_condition_number")
        }
    )
    roles = {
        key: (
            "training"
            if key.startswith("t")
            else "policy-selection"
            if key.startswith("s")
            else "validation"
            if key.startswith("v")
            else "locked-heldout"
        )
        for key in values
    }
    members = [_member(key, role, assay_class=measured["t1"].assay_class) for key, role in sorted(roles.items())]
    for member in members:
        if member["role"] == "locked-heldout":
            member["provenance"] = "external-custodian"
    raw_study["partitions"] = _manifest(*members)
    raw_study["exposure_ledger_id"] = ledger_id
    study = decode_target_study(raw_study)
    groups = connected_leakage_groups(study.partitions)
    rosters = {
        role: [
            {"key": key, "group_key": groups[key], "strata": ["all"]}
            for key in sorted(values, key=lambda key: groups[key])
            if roles[key] == role
        ]
        for role in ("training", "policy-selection", "validation", "locked-heldout")
    }
    from tests.unit.test_calibration_length_evaluation import evaluation_rows

    protocol, _ = evaluation_protocol(evaluation_rows(), ("affine-a", "affine-A"))
    raw_protocol = length_protocol_document(protocol)
    raw_protocol["eligible_roster_sha256"] = decode_length_eligible_roster(rosters["policy-selection"]).sha256
    raw_study["protocol"] = raw_protocol
    study = decode_target_study(raw_study)
    write_json(root / "study.json", raw_study)
    write_json(root / "annotation.json", _annotation_raw())
    raw_runs = {"schema_version": "calibration-runs-v2", "target": "length", "runs": []}
    for key in sorted(values):
        path = root / f"{key}-features.json"
        write_json(path, encode_length_feature_artifact(measured[key]))
        row = run_document()["runs"][0]
        row.update(
            manifest_key=key,
            policy_sha256=module.length_measurement_policy_sha256(study),
            capture_policy_sha256=module.length_measurement_policy_sha256(study),
            input_sha256=canonical_sha256(key),
            baseline_assets_sha256=study.baseline.sha256,
            producer_sha256=canonical_sha256(candidate_producer_document(study.baseline.producer)),
        )
        row["assets"]["length_features"] = {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
        }
        raw_runs["runs"].append(row)
    runs = decode_target_runs(raw_runs)
    write_json(root / "runs.json", raw_runs)
    for role, roster in rosters.items():
        keys = sorted(member["key"] for member in roster)
        truth_path = root / "roles" / role / "truth.json"
        write_json(
            truth_path,
            {
                "schema_version": "calibration-length-truth-v1",
                "boundary_definition": TARGET_BOUNDARY_DEFINITION,
                "rows": [{"key": key, "total_repeat_count": 10 + 50 * values[key]} for key in keys],
            },
        )
        source = {
            "schema_version": "calibration-role-source-v2",
            "study_sha256": study.sha256,
            "role": role,
            "roster": roster,
            "excluded": [],
            "evidence_domains": dict.fromkeys(keys, "synthetic"),
            "identities_by_key": {
                key: [
                    {"namespace": "physical-readset", "sha256": canonical_sha256(key)},
                    {"namespace": "specimen", "sha256": canonical_sha256("specimen-" + key)},
                ]
                for key in keys
            },
            "truth_asset": {
                "path": str(truth_path),
                "sha256": hashlib.sha256(truth_path.read_bytes()).hexdigest(),
                "size_bytes": truth_path.stat().st_size,
            },
            "run_manifest_sha256": runs.sha256,
        }
        write_json(root / "roles" / role / "source.json", source)
    return Namespace(evidence=root, exposure_ledger=ledger, objective="length-total-v1"), study


def test_fit_bundle_executes_real_fit_selection_and_emits_research_candidate(tmp_path):
    from vntyper.scripts.calibration_artifact_io import load_object
    from vntyper.scripts.calibration_candidate import decode_candidate

    module = import_module("vntyper.scripts.calibration_length_controller")
    args, study = fit_fixture(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    assert module.fit_length_bundle(args, output) is True
    candidate = decode_candidate(load_object(output / "candidate.json", "candidate"))
    assert candidate.study_sha256 == study.sha256
    assert candidate.status == "research-only"
    model = load_object(output / "payload" / "length-model.json", "model")
    assert model["coefficients"] == pytest.approx([50])
    assert model["intercept"] == pytest.approx(10)
    assert (output / "report.html").is_file()
    assert not (output / "validation.json").exists()
    history = args.exposure_ledger.read_text()
    assert '"role":"training"' in history
    assert '"role":"policy-selection"' in history
    assert '"role":"validation"' not in history
    profile = module.load_length_research_profile(output)
    assert profile.candidate == candidate
    assert profile.selected_protocol_candidate_id == "affine-a"
    assert profile.training_profile.artifact.baseline.mean_total_repeat_count == pytest.approx(110)
    assert set(profile.payload_files) == {"length-annotation.json", "length-model.json"}
    assert profile.payload_files["length-model.json"] == (output / "payload" / "length-model.json").read_bytes()
    assert len(profile.projection_sha256) == 64


def test_fit_rejects_wrong_objective_before_exposing_any_role(tmp_path):
    module = import_module("vntyper.scripts.calibration_length_controller")
    args, _ = fit_fixture(tmp_path)
    args.objective = "caller-safety-v1"
    output = tmp_path / "output"
    output.mkdir()
    before = args.exposure_ledger.read_bytes()
    with pytest.raises(ValueError, match="objective"):
        module.fit_length_bundle(args, output)
    assert args.exposure_ledger.read_bytes() == before
    assert not tuple(output.iterdir())


@pytest.mark.parametrize("mode", ["extra-file", "payload-bytes", "payload-inventory", "wrong-path-type"])
def test_research_profile_rejects_inventory_or_payload_tampering(tmp_path, mode):
    module = import_module("vntyper.scripts.calibration_length_controller")
    args, _ = fit_fixture(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    assert module.fit_length_bundle(args, output) is True
    if mode == "extra-file":
        (output / "unexpected.json").write_text("{}")
    elif mode == "payload-bytes":
        path = output / "payload" / "length-model.json"
        path.write_bytes(path.read_bytes() + b" ")
    elif mode == "payload-inventory":
        (output / "payload" / "unexpected.json").write_text("{}")
    profile_path = output if mode != "wrong-path-type" else str(output)
    with pytest.raises(ValueError, match="inventory|payload|directory"):
        module.load_length_research_profile(profile_path)


def test_fixed_profile_evaluates_one_validation_source_without_fit_or_selection(tmp_path):
    from vntyper.scripts.calibration_artifact_io import load_object
    from vntyper.scripts.calibration_exposure_io import record_exposure
    from vntyper.scripts.calibration_role_source import decode_role_source

    module = import_module("vntyper.scripts.calibration_length_controller")
    args, study = fit_fixture(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    assert module.fit_length_bundle(args, output) is True
    profile = module.load_length_research_profile(output)
    runs = decode_target_runs(load_object(args.evidence / "runs.json", "runs"))
    source = decode_role_source(
        load_object(args.evidence / "roles" / "validation" / "source.json", "source"),
        study=study,
        runs=runs,
        expected_role="validation",
    )
    receipt = record_exposure(
        args.exposure_ledger,
        expected_ledger_id=study.exposure_ledger_id,
        target="length",
        role="validation",
        study_sha256=study.sha256,
        partition_sha256=study.partitions.sha256,
        evidence_sha256=source.sha256,
        identities=[{"namespace": name, "sha256": digest} for name, digest in source.identities],
    )
    with patch("vntyper.scripts.calibration_length.fit_length_hypotheses", side_effect=AssertionError("refit")):
        result, report = module.evaluate_fixed_length_source(profile, runs, source, receipt)
    assert result.phase == "validation"
    assert result.selection.status == "not-applicable"
    assert sum(candidate.status == "evaluated" for candidate in result.candidates) == 1
    assert "Total diploid repeat-count calibration: validation" in report


def test_assess_bundle_records_development_exposure_and_never_refits_or_selects(tmp_path):
    from argparse import Namespace

    from vntyper.scripts.calibration_artifact_io import load_object, write_json
    from vntyper.scripts.calibration_development_source import decode_development_source
    from vntyper.scripts.calibration_exposure_io import record_exposure
    from vntyper.scripts.calibration_length_assessment import decode_length_development_assessment

    module = import_module("vntyper.scripts.calibration_length_controller")
    args, _ = fit_fixture(tmp_path)
    profile_root = tmp_path / "profile"
    profile_root.mkdir()
    assert module.fit_length_bundle(args, profile_root) is True
    profile = module.load_length_research_profile(profile_root)
    runs = decode_target_runs(load_object(args.evidence / "runs.json", "runs"))
    role_raw = load_object(args.evidence / "roles" / "validation" / "source.json", "source")
    source_raw = {
        "schema_version": "calibration-development-source-v1",
        "target": "length",
        "evidence_role": "development-assessment",
        "promotion_eligible": False,
        "candidate_sha256": profile.candidate.sha256,
        "study_sha256": profile.study.sha256,
        "run_manifest_sha256": runs.sha256,
        "partition_sha256": canonical_sha256({"development": role_raw["roster"]}),
        "roster": role_raw["roster"],
        "identities_by_key": role_raw["identities_by_key"],
        "evidence_domains": role_raw["evidence_domains"],
        "previously_examined": {row["key"]: True for row in role_raw["roster"]},
        "truth_asset": role_raw["truth_asset"],
    }
    development = tmp_path / "development"
    development.mkdir()
    write_json(development / "source.json", source_raw)
    source = decode_development_source(source_raw, candidate=profile.candidate, runs=runs)
    assessment_root = tmp_path / "assessment"
    assessment_root.mkdir()
    assess_args = Namespace(
        profile=profile_root,
        intake=development,
        runs=args.evidence / "runs.json",
        exposure_ledger=args.exposure_ledger,
    )
    with (
        patch("vntyper.scripts.calibration_length.fit_length_hypotheses", side_effect=AssertionError("refit")),
        patch.object(module, "evaluate_length_hypotheses", side_effect=AssertionError("selection")),
    ):
        assert module.assess_length_bundle(assess_args, assessment_root) is True
    receipt = record_exposure(
        args.exposure_ledger,
        expected_ledger_id=profile.study.exposure_ledger_id,
        target="length",
        role="development-assessment",
        study_sha256=profile.study.sha256,
        partition_sha256=source.partition_sha256,
        evidence_sha256=source.sha256,
        identities=[{"namespace": name, "sha256": digest} for name, digest in source.identities],
    )
    evidence = module.load_length_development_evidence(profile, runs, source, receipt)
    protocol = module.length_protocol_for_roster(profile.study.protocol, source.roster)
    decoded = decode_length_development_assessment(
        load_object(assessment_root / "assessment.json", "assessment"),
        profile=profile.training_profile,
        evidence=evidence,
        roster=source.roster,
        protocol=protocol,
        fixed_candidate_id=profile.selected_protocol_candidate_id,
    )
    assert decoded.assessment.promotion_eligible is False
    assert decoded.assessment.selection_status == "not-applicable"
    assert (assessment_root / "report.html").read_text() == decoded.report_html
    assert {path.name for path in assessment_root.iterdir()} == {
        "assessment.json",
        "checksums.json",
        "report.html",
    }
    history = args.exposure_ledger.read_text()
    assert '"role":"development-assessment"' in history
    assert not (assessment_root / "candidate.json").exists()
