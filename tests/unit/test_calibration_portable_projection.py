"""Aggregate-free portable approval projection and consistency checks."""

from dataclasses import replace

import pytest

from vntyper.scripts.calibration_candidate import decode_candidate
from vntyper.scripts.calibration_target_attestation import (
    decode_target_locked_attestation,
    decode_target_validation_attestation,
)
from vntyper.scripts.calibration_target_authority import decode_target_custodian_authority
from vntyper.scripts.calibration_target_completion import decode_target_completion
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def _candidate(*, target: str = "length", payload_sha256: str = "f" * 64):
    applicability = {
        "domain": "synthetic",
        "assemblies": ["GRCh38"],
        "assay_classes": ["capture"],
        "input_scopes": ["full"],
        "preprocessing_ids": ["alignment-v1"],
    }
    if target == "callers":
        applicability["required_callers"] = ["kestrel"]
    else:
        applicability.update(
            {
                "aligner_name": "synthetic-aligner",
                "aligner_version": "1.0",
                "aligner_arguments_sha256": "1" * 64,
                "primary_secondary_marking": "primary-only",
                "counting_policy_sha256": "2" * 64,
            }
        )
    raw = {
        "schema_version": "calibration-candidate-v2",
        "target": target,
        "study_sha256": "2" * 64,
        "baseline_sha256": "5" * 64,
        "partition_sha256": "4" * 64,
        "training_evidence_sha256": "d" * 64,
        "selection_evidence_sha256": "e" * 64,
        "payload_sha256": payload_sha256,
        "applicability": applicability,
        "producer": {
            "name": "synthetic-calibration",
            "version": "1.0",
            "source_revision": "a" * 40,
            "tool_versions": {"depth-tool": "1.0"},
            "feature_schema_sha256": "0" * 64,
        },
        "status": "research-only",
    }
    raw["candidate_id"] = canonical_sha256(raw)
    return decode_candidate(raw)


def _artifacts(*, target: str = "length", validation_status: str = "passed", locked_status: str = "passed"):
    candidate = _candidate(target=target)
    validation_raw = {
        "schema_version": "calibration-target-validation-attestation-v2",
        "target": target,
        "role": "validation",
        "status": validation_status,
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": candidate.study_sha256,
        "protocol_sha256": "3" * 64,
        "partition_sha256": candidate.partition_sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "evidence_sha256": "6" * 64,
        "run_manifest_sha256": "7" * 64,
        "metrics_sha256": "8" * 64,
        "exposure_ledger_id": "9" * 64,
        "exposure_receipt_sha256": "a" * 64,
    }
    validation = decode_target_validation_attestation(validation_raw)
    authority_raw = {
        "schema_version": "calibration-target-custodian-authority-v2",
        "authority_kind": "external-custodian",
        "custodian_name": "Independent Example Repository",
        "attestation_id": "IER-332-0001",
        "status": "authorized",
        "role": "locked-heldout",
        "target": target,
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": candidate.study_sha256,
        "protocol_sha256": validation.protocol_sha256,
        "partition_sha256": candidate.partition_sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "validation_evidence_sha256": validation.evidence_sha256,
        "locked_heldout_evidence_sha256": "b" * 64,
        "validation_run_manifest_sha256": validation.run_manifest_sha256,
        "locked_heldout_run_manifest_sha256": "c" * 64,
        "locked_payload_sha256": "d" * 64,
        "validation_attestation_sha256": validation.sha256,
        "exposure_ledger_id": validation.exposure_ledger_id,
        "validation_exposure_receipt_sha256": validation.exposure_receipt_sha256,
    }
    authority = decode_target_custodian_authority(authority_raw)
    locked_raw = {
        "schema_version": "calibration-target-locked-attestation-v2",
        "target": target,
        "role": "locked-heldout",
        "status": locked_status,
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": candidate.study_sha256,
        "protocol_sha256": validation.protocol_sha256,
        "partition_sha256": candidate.partition_sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "evidence_sha256": authority.locked_heldout_evidence_sha256,
        "run_manifest_sha256": authority.locked_heldout_run_manifest_sha256,
        "metrics_sha256": "e" * 64,
        "exposure_ledger_id": validation.exposure_ledger_id,
        "exposure_receipt_sha256": "f" * 64,
        "validation_attestation_sha256": validation.sha256,
        "custodian_authority_sha256": authority.sha256,
        "custody_consumption_receipt_sha256": "0" * 64,
    }
    locked = decode_target_locked_attestation(locked_raw)
    completion_raw = {
        "schema_version": "calibration-target-completion-v2",
        "status": "completed",
        "target": target,
        "candidate_sha256": candidate.sha256,
        "candidate_id": candidate.candidate_id,
        "study_sha256": candidate.study_sha256,
        "protocol_sha256": validation.protocol_sha256,
        "partition_sha256": candidate.partition_sha256,
        "baseline_sha256": candidate.baseline_sha256,
        "validation_attestation_sha256": validation.sha256,
        "locked_heldout_attestation_sha256": locked.sha256,
        "custodian_authority_sha256": authority.sha256,
        "locked_heldout_evidence_sha256": locked.evidence_sha256,
        "locked_heldout_run_manifest_sha256": locked.run_manifest_sha256,
        "exposure_ledger_id": validation.exposure_ledger_id,
        "validation_exposure_receipt_sha256": validation.exposure_receipt_sha256,
        "locked_heldout_exposure_receipt_sha256": locked.exposure_receipt_sha256,
        "custody_consumption_receipt_sha256": locked.custody_consumption_receipt_sha256,
    }
    completion = decode_target_completion(completion_raw)
    return candidate, validation, locked, authority, completion


@pytest.mark.parametrize("target", ["callers", "length"])
def test_portable_approval_is_aggregate_free_and_binds_runtime_metadata(target: str) -> None:
    from vntyper.scripts.calibration_portable_projection import (
        build_portable_approval,
        portable_approval_document,
        validate_portable_approval_candidate,
    )

    candidate, validation, locked, authority, completion = _artifacts(target=target)
    approval = build_portable_approval(candidate, validation, locked, authority, completion)
    document = portable_approval_document(approval)

    assert approval.disposition == "passed-validation-and-locked-heldout"
    assert approval.exposure_ledger_id == validation.exposure_ledger_id
    assert document["candidate_sha256"] == candidate.sha256
    assert document["custodian_authority_sha256"] == authority.sha256
    assert not ({"metrics", "rows", "members", "custodian_name", "attestation_id"} & set(document))
    validate_portable_approval_candidate(approval, candidate)


@pytest.mark.parametrize("component", ["candidate", "validation", "locked", "authority", "completion"])
def test_portable_builder_rejects_each_inconsistent_component(component: str) -> None:
    from vntyper.scripts.calibration_portable_projection import build_portable_approval

    artifacts = list(_artifacts())
    index = ("candidate", "validation", "locked", "authority", "completion").index(component)
    artifacts[index] = replace(artifacts[index], candidate_id="f" * 64)
    with pytest.raises(ValueError):
        build_portable_approval(*artifacts)


def test_portable_builder_rejects_canonical_but_cross_bound_artifacts() -> None:
    from vntyper.scripts.calibration_portable_projection import build_portable_approval
    from vntyper.scripts.calibration_target_attestation import target_validation_attestation_document

    candidate, validation, locked, authority, completion = _artifacts()
    document = target_validation_attestation_document(validation)
    document["baseline_sha256"] = "f" * 64
    different_validation = decode_target_validation_attestation(document)
    with pytest.raises(ValueError, match="bindings"):
        build_portable_approval(candidate, different_validation, locked, authority, completion)


def test_portable_builder_rejects_canonical_protocol_and_lineage_mismatches() -> None:
    from vntyper.scripts.calibration_portable_projection import build_portable_approval
    from vntyper.scripts.calibration_target_attestation import (
        target_locked_attestation_document,
        target_validation_attestation_document,
    )
    from vntyper.scripts.calibration_target_authority import target_custodian_authority_document

    candidate, validation, locked, authority, completion = _artifacts()
    locked_document = target_locked_attestation_document(locked)
    locked_document["protocol_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="protocol"):
        build_portable_approval(
            candidate,
            validation,
            decode_target_locked_attestation(locked_document),
            authority,
            completion,
        )

    authority_document = target_custodian_authority_document(authority)
    authority_document["validation_evidence_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="lineage"):
        build_portable_approval(
            candidate,
            validation,
            locked,
            decode_target_custodian_authority(authority_document),
            completion,
        )

    validation_document = target_validation_attestation_document(validation)
    validation_document["status"] = "failed"
    with pytest.raises(ValueError, match="passed"):
        build_portable_approval(
            candidate,
            decode_target_validation_attestation(validation_document),
            locked,
            authority,
            completion,
        )


@pytest.mark.parametrize("status_kind", ["validation", "locked"])
def test_portable_builder_requires_both_passed_attestations(status_kind: str) -> None:
    from vntyper.scripts.calibration_portable_projection import build_portable_approval

    artifacts = _artifacts(
        validation_status="failed" if status_kind == "validation" else "passed",
        locked_status="failed" if status_kind == "locked" else "passed",
    )
    with pytest.raises(ValueError, match="passed"):
        build_portable_approval(*artifacts)


def test_portable_decoder_is_closed_and_runtime_rejects_wrong_candidate() -> None:
    from vntyper.scripts.calibration_portable_projection import (
        build_portable_approval,
        decode_portable_approval,
        portable_approval_document,
        validate_portable_approval_candidate,
    )

    artifacts = _artifacts()
    approval = build_portable_approval(*artifacts)
    document = portable_approval_document(approval)
    assert decode_portable_approval(document) == approval
    document["metrics"] = {}
    with pytest.raises(ValueError, match="fields"):
        decode_portable_approval(document)
    with pytest.raises(ValueError, match="candidate"):
        validate_portable_approval_candidate(approval, _candidate(payload_sha256="0" * 64))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("schema_version", "calibration-portable-approval-v1"),
        ("disposition", "authorized"),
        ("target", "dominance"),
        ("payload_sha256", "F" * 64),
    ],
)
def test_portable_decoder_rejects_invalid_schema_disposition_target_and_digest(field: str, value: object) -> None:
    from vntyper.scripts.calibration_portable_projection import (
        build_portable_approval,
        decode_portable_approval,
        portable_approval_document,
    )

    document = portable_approval_document(build_portable_approval(*_artifacts()))
    document[field] = value
    with pytest.raises(ValueError):
        decode_portable_approval(document)


@pytest.mark.parametrize("component", range(5))
def test_portable_builder_requires_decoded_typed_components(component: int) -> None:
    from vntyper.scripts.calibration_portable_projection import build_portable_approval

    artifacts = list(_artifacts())
    artifacts[component] = {}
    with pytest.raises(ValueError):
        build_portable_approval(*artifacts)  # type: ignore[arg-type]


def test_portable_projection_revalidates_typed_identity() -> None:
    from vntyper.scripts.calibration_portable_projection import build_portable_approval, portable_approval_document

    approval = build_portable_approval(*_artifacts())
    with pytest.raises(ValueError, match="canonical"):
        portable_approval_document(replace(approval, sha256="f" * 64))
    with pytest.raises(ValueError, match="typed"):
        portable_approval_document({})  # type: ignore[arg-type]
