"""Bound length training artifacts and development-only assessment."""

import hashlib
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from importlib import import_module

import pytest

from tests.unit.test_calibration_length import metadata_document, protocol_for, roster_for, training_row
from tests.unit.test_calibration_length_evaluation import evaluation_protocol, evaluation_rows
from tests.unit.test_calibration_target_contract import study_document
from tests.unit.test_length_estimation import features
from vntyper.scripts import calibration_length_evidence as evidence_contract

pytestmark = pytest.mark.unit


def _study_and_metadata(rows, roster, *kinds):
    target = import_module("vntyper.scripts.calibration_target_contract")
    training = import_module("vntyper.scripts.calibration_length")
    protocol_contract = import_module("vntyper.scripts.calibration_length_protocol")
    protocol = protocol_for(*kinds)
    raw_metadata = metadata_document(rows, roster)
    raw_study = study_document()
    raw_study["protocol"] = protocol_contract.length_protocol_document(protocol)
    raw_study["applicability"] = deepcopy(raw_metadata["applicability"])
    for member in raw_study["partitions"]["members"]:
        member["assay_class"] = rows[0].features.assay_class
    raw_study["baseline"].update(
        annotation_sha256=raw_metadata["annotation_sha256"],
        counting_policy_sha256=raw_metadata["counting_policy_sha256"],
        maximum_condition_number=raw_metadata["maximum_condition_number"],
        producer=deepcopy(raw_metadata["producer"]),
    )
    study = target.decode_target_study(raw_study)
    evidence = evidence_contract.bind_length_training_evidence(
        rows,
        roster,
        study_sha256=study.sha256,
        partition_sha256="6" * 64,
        run_manifest_sha256="7" * 64,
    )
    raw_metadata["study_sha256"] = study.sha256
    raw_metadata["training_evidence_sha256"] = evidence.sha256
    metadata = training.decode_length_training_metadata(raw_metadata, roster)
    return study, metadata, evidence


def _training_artifact(*kinds):
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    rows = tuple(training_row(str(index), value, 10 + 50 * value) for index, value in enumerate((1, 2, 3), 1))
    roster = roster_for(rows)
    study, metadata, evidence = _study_and_metadata(rows, roster, *kinds)
    return artifacts.build_length_training_artifact(study, roster, metadata, evidence), study, roster


def _training_profile(*kinds):
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    artifact, study, roster = _training_artifact(*kinds)
    return artifacts.bind_length_training_profile(artifact, study, roster)


def test_training_artifact_runs_fit_and_binds_the_exact_study_baseline_plan():
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    artifact, study, roster = _training_artifact("affine-A", "physical-A")

    assert artifact.study_sha256 == study.sha256
    assert artifact.protocol_sha256 == study.protocol.sha256
    assert artifact.training_roster_sha256 == roster.sha256
    assert artifact.baseline.study_sha256 == study.sha256
    assert [item.status for item in artifact.outcomes] == ["fitted", "ineligible"]
    assert (
        artifacts.decode_length_training_artifact(
            artifacts.length_training_artifact_document(artifact, study=study, training_roster=roster),
            study=study,
            training_roster=roster,
        )
        == artifact
    )
    with pytest.raises(FrozenInstanceError):
        artifact.study_sha256 = "0" * 64


def test_training_evidence_digest_is_derived_from_exact_opened_rows_and_required_by_metadata():
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    rows = tuple(training_row(str(index), value, 10 + 50 * value) for index, value in enumerate((1, 2, 3), 1))
    roster = roster_for(rows)
    study, _, evidence = _study_and_metadata(rows, roster, "affine-A")
    document = evidence_contract.length_training_evidence_document(evidence)
    assert document["rows"][0]["features_sha256"] == rows[0].features.sha256
    assert evidence_contract.decode_length_training_evidence(document, rows=rows, roster=roster) == evidence
    with pytest.raises(ValueError, match="partition_sha256"):
        evidence_contract.length_training_evidence_document(replace(evidence, partition_sha256="bad"))

    changed_rows = (rows[0], rows[1], replace(rows[2], total_truth_repeat_count=161))
    with pytest.raises(ValueError, match="supplied typed rows"):
        evidence_contract.decode_length_training_evidence(document, rows=changed_rows, roster=roster)
    raw_metadata = metadata_document(
        rows,
        roster,
        study_sha256=study.sha256,
        training_evidence_sha256="0" * 64,
    )
    changed_metadata = import_module("vntyper.scripts.calibration_length").decode_length_training_metadata(
        raw_metadata, roster
    )
    with pytest.raises(ValueError, match="study, roster, or metadata"):
        artifacts.build_length_training_artifact(
            study,
            roster,
            changed_metadata,
            evidence,
        )


@pytest.mark.parametrize(
    "change,match",
    [
        ("study", "study"),
        ("annotation", "annotation"),
        ("counting", "counting"),
        ("condition", "condition"),
        ("producer", "producer"),
        ("applicability", "applicability"),
        ("qc", "QC"),
    ],
)
def test_training_artifact_rejects_metadata_that_drifted_from_the_study_plan(change, match):
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    rows = tuple(training_row(str(index), value, 10 + 50 * value) for index, value in enumerate((1, 2, 3), 1))
    roster = roster_for(rows)
    study, _, evidence = _study_and_metadata(rows, roster, "affine-A")
    raw = metadata_document(
        rows,
        roster,
        study_sha256=study.sha256,
        training_evidence_sha256=evidence.sha256,
    )
    if change == "study":
        raw["study_sha256"] = "0" * 64
    elif change == "annotation":
        raw["annotation_sha256"] = "0" * 64
    elif change == "counting":
        raw["counting_policy_sha256"] = "0" * 64
        raw["applicability"]["counting_policy_sha256"] = "0" * 64
    elif change == "condition":
        raw["maximum_condition_number"] += 1
    elif change == "producer":
        raw["producer"]["version"] = "changed"
    elif change == "applicability":
        raw["applicability"]["domain"] = "external"
    else:
        raw["qc"]["minimum_denominator_mean_depth"] = 11
    metadata = import_module("vntyper.scripts.calibration_length").decode_length_training_metadata(raw, roster)
    with pytest.raises(ValueError, match=match):
        artifacts.build_length_training_artifact(study, roster, metadata, evidence)


def test_training_artifact_decoder_rejects_unknown_fields_hash_and_nested_model_drift():
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    artifact, study, roster = _training_artifact("affine-A")
    document = artifacts.length_training_artifact_document(artifact, study=study, training_roster=roster)
    for changed in (
        {**document, "unknown": None},
        {**document, "sha256": "0" * 64},
    ):
        with pytest.raises(ValueError):
            artifacts.decode_length_training_artifact(changed, study=study, training_roster=roster)
    nested = deepcopy(document)
    nested["outcomes"][0]["model"]["intercept"] += 1
    with pytest.raises(ValueError):
        artifacts.decode_length_training_artifact(nested, study=study, training_roster=roster)


def test_training_projector_and_role_specific_profile_reject_forged_typed_content():
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    artifact, study, roster = _training_artifact("affine-A")
    with pytest.raises(ValueError, match="typed artifact"):
        artifacts.length_training_artifact_document({}, study=study, training_roster=roster)
    with pytest.raises(ValueError, match="study, protocol, or roster"):
        artifacts.length_training_artifact_document(
            replace(artifact, protocol_sha256="0" * 64), study=study, training_roster=roster
        )
    profile = artifacts.bind_length_training_profile(artifact, study, roster)
    assert profile.artifact is artifact
    with pytest.raises(ValueError, match="contextual typed profile"):
        artifacts.validate_length_training_profile({}, study.protocol)
    with pytest.raises(ValueError, match="canonical content"):
        artifacts.validate_length_training_profile(
            replace(profile, artifact=replace(artifact, sha256="0" * 64)), study.protocol
        )


@pytest.mark.parametrize(
    "change,match",
    [
        ({"key": " bad"}, "row key"),
        ({"group_key": ""}, "group key"),
        ({"truth_boundary_definition": "other"}, "boundary"),
        ({"total_truth_repeat_count": True}, "truth"),
        ({"evidence_domain": "guessed"}, "domain"),
    ],
)
def test_role_evidence_rejects_malformed_typed_rows(change, match):
    _, rows, protocol, roster = _development_evidence()
    malformed = (replace(rows[0], **change), rows[1])
    with pytest.raises(ValueError, match=match):
        evidence_contract.bind_length_role_evidence(
            malformed,
            roster,
            protocol,
            study_sha256="a" * 64,
            partition_sha256="b" * 64,
            run_manifest_sha256="c" * 64,
        )


def test_role_evidence_projection_and_decoder_reject_mutability_schema_and_digest_drift():
    evidence, rows, protocol, roster = _development_evidence()
    document = evidence_contract.length_role_evidence_document(evidence)
    with pytest.raises(ValueError, match="immutable typed evidence"):
        evidence_contract.length_role_evidence_document(replace(evidence, rows=list(rows)))
    with pytest.raises(ValueError, match="digest"):
        evidence_contract.length_role_evidence_document(replace(evidence, sha256="0" * 64))
    with pytest.raises(ValueError, match="phase"):
        evidence_contract.length_role_evidence_document(replace(evidence, phase=[]))
    for changed in (
        {**document, "schema_version": "calibration-length-role-evidence-v0"},
        {**document, "phase": "training"},
        {**document, "rows": [dict(document["rows"][0], unknown=None), document["rows"][1]]},
    ):
        with pytest.raises(ValueError):
            evidence_contract.decode_length_role_evidence(changed, rows=rows, roster=roster, protocol=protocol)


def _development_evidence():
    rows = evaluation_rows(phase="development-assessment")
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    evidence = evidence_contract.bind_length_role_evidence(
        rows,
        roster,
        protocol,
        study_sha256="a" * 64,
        partition_sha256="b" * 64,
        run_manifest_sha256="c" * 64,
    )
    return evidence, rows, protocol, roster


def test_role_evidence_is_closed_and_rederived_from_opened_typed_rows():
    evidence, rows, protocol, roster = _development_evidence()
    document = evidence_contract.length_role_evidence_document(evidence)
    assert document["rows"][0]["features_sha256"] == rows[0].features.sha256
    assert "features" not in document["rows"][0]
    assert (
        evidence_contract.decode_length_role_evidence(document, rows=rows, roster=roster, protocol=protocol) == evidence
    )
    with pytest.raises(ValueError, match="typed rows"):
        evidence_contract.decode_length_role_evidence(
            document, rows=tuple(reversed(rows)), roster=roster, protocol=protocol
        )
    changed = deepcopy(document)
    changed["rows"][0]["total_truth_repeat_count"] += 1
    with pytest.raises(ValueError, match="typed rows"):
        evidence_contract.decode_length_role_evidence(changed, rows=rows, roster=roster, protocol=protocol)


def test_development_assessment_evaluates_one_fixed_model_without_fit_selection_or_authority(monkeypatch):
    artifacts = import_module("vntyper.scripts.calibration_length_artifacts")
    assessment_contract = import_module("vntyper.scripts.calibration_length_assessment")
    evidence, _, protocol, roster = _development_evidence()
    profile = _training_profile("affine-A")
    artifact = profile.artifact
    monkeypatch.setattr(artifacts, "fit_length_hypotheses", lambda *_a, **_k: pytest.fail("assessment refit"))

    output = assessment_contract.assess_length_candidate(
        profile, evidence, roster, protocol, fixed_candidate_id="affine-a"
    )

    assessment = output.assessment
    assert assessment.evidence_role == "development-assessment"
    assert assessment.promotion_eligible is False
    assert assessment.selection_status == "not-applicable"
    assert assessment.candidate_id == "affine-a"
    assert assessment.model_sha256 == artifact.outcomes[0].model.sha256
    assert assessment.input_evidence_sha256 == evidence.sha256
    assert assessment.report_sha256 == hashlib.sha256(output.report_html.encode()).hexdigest()
    assert "Validation passed" not in output.report_html
    document = assessment_contract.length_development_assessment_document(assessment)
    assert document["promotion_eligible"] is False
    assert "authority" not in document and "receipt" not in document
    assert (
        assessment_contract.decode_length_development_assessment(
            document,
            profile=profile,
            evidence=evidence,
            roster=roster,
            protocol=protocol,
            fixed_candidate_id="affine-a",
        ).assessment
        == assessment
    )


def test_development_assessment_records_exploratory_out_of_domain_values_separately():
    assessment_contract = import_module("vntyper.scripts.calibration_length_assessment")
    evidence, rows, protocol, roster = _development_evidence()
    profile = _training_profile("affine-A")
    changed_rows = tuple(replace(row, evidence_domain="external") for row in rows)
    changed_evidence = evidence_contract.bind_length_role_evidence(
        changed_rows,
        roster,
        protocol,
        study_sha256=evidence.study_sha256,
        partition_sha256=evidence.partition_sha256,
        run_manifest_sha256=evidence.run_manifest_sha256,
    )

    output = assessment_contract.assess_length_candidate(
        profile, changed_evidence, roster, protocol, fixed_candidate_id="affine-a"
    )

    candidate = next(item for item in output.assessment.evaluation.candidates if item.candidate_id == "affine-a")
    assert all(item.prediction is None for item in candidate.predictions)
    assert all(item.availability_reasons == ("unsupported_evidence_domain",) for item in candidate.predictions)
    assert [item.prediction for item in output.assessment.exploratory] == pytest.approx([60, 160])
    assert all(item.reasons == () for item in output.assessment.exploratory)

    missing = replace(
        changed_rows[0],
        features=features(depths=(100, 200, 200, 0, 0, 100), manifest_key=changed_rows[0].key),
    )
    missing_evidence = evidence_contract.bind_length_role_evidence(
        (missing, changed_rows[1]),
        roster,
        protocol,
        study_sha256=evidence.study_sha256,
        partition_sha256=evidence.partition_sha256,
        run_manifest_sha256=evidence.run_manifest_sha256,
    )
    unavailable = assessment_contract.assess_length_candidate(
        profile, missing_evidence, roster, protocol, fixed_candidate_id="affine-a"
    )
    assert unavailable.assessment.exploratory[0].prediction is None
    assert unavailable.assessment.exploratory[0].reasons == ("missing_A",)


def test_assessment_projection_and_decoder_fail_closed_on_forged_local_evidence():
    assessment_contract = import_module("vntyper.scripts.calibration_length_assessment")
    evidence, _, protocol, roster = _development_evidence()
    profile = _training_profile("affine-A")
    output = assessment_contract.assess_length_candidate(
        profile, evidence, roster, protocol, fixed_candidate_id="affine-a"
    )
    assessment = output.assessment
    for forged in (
        replace(assessment, promotion_eligible=True),
        replace(assessment, selection_status="selected"),
        replace(assessment, candidate_id="other"),
        replace(assessment, exploratory=tuple(reversed(assessment.exploratory))),
        replace(assessment, sha256="0" * 64),
    ):
        with pytest.raises(ValueError):
            assessment_contract.length_development_assessment_document(forged)
    invalid_prediction = replace(assessment.exploratory[0], prediction=float("inf"), reasons=())
    with pytest.raises(ValueError, match="prediction"):
        assessment_contract.length_development_assessment_document(
            replace(assessment, exploratory=(invalid_prediction, assessment.exploratory[1]))
        )
    document = assessment_contract.length_development_assessment_document(assessment)
    changed = deepcopy(document)
    changed["report_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="recomputed evidence"):
        assessment_contract.decode_length_development_assessment(
            changed,
            profile=profile,
            evidence=evidence,
            roster=roster,
            protocol=protocol,
            fixed_candidate_id="affine-a",
        )
    with pytest.raises(ValueError, match="must be an object"):
        assessment_contract.decode_length_development_assessment(
            None,
            profile=profile,
            evidence=evidence,
            roster=roster,
            protocol=protocol,
            fixed_candidate_id="affine-a",
        )


@pytest.mark.parametrize("phase", ["policy-selection", "validation", "locked-heldout"])
def test_assessment_rejects_promotion_roles_and_wrong_artifact_bindings(phase):
    assessment_contract = import_module("vntyper.scripts.calibration_length_assessment")
    rows = evaluation_rows(phase=phase)
    protocol, roster = evaluation_protocol(rows, ("affine-a", "affine-A"))
    evidence = evidence_contract.bind_length_role_evidence(
        rows,
        roster,
        protocol,
        study_sha256="a" * 64,
        partition_sha256="b" * 64,
        run_manifest_sha256="c" * 64,
    )
    profile = _training_profile("affine-A")
    with pytest.raises(ValueError, match="development-assessment"):
        assessment_contract.assess_length_candidate(profile, evidence, roster, protocol, fixed_candidate_id="affine-a")
