"""Closed study/profile binding commitments for calibration candidates."""

from dataclasses import replace

import pytest

from tests.unit.test_calibration_contract import synthetic_protocol
from vntyper.scripts.calibration_locked_artifacts import build_locked_declaration_documents, decode_locked_payload
from vntyper.scripts.calibration_manifest import decode_study_declaration
from vntyper.scripts.calibration_study_binding import build_study_binding, decode_study_binding, validate_study_binding

pytestmark = pytest.mark.unit


def _study():
    members = []
    for key, role in (
        ("held", "locked-heldout"),
        ("select", "policy-selection"),
        ("train", "training"),
        ("validate", "validation"),
    ):
        members.append(
            {
                "key": key,
                "role": role,
                "provenance": "external-custodian" if role == "locked-heldout" else "development",
                "assay_class": "capture-short-read",
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
    return decode_study_declaration(
        {
            "schema_version": "calibration-study-v1",
            "protocol": synthetic_protocol(),
            "partitions": {"schema_version": "calibration-partitions-v1", "members": members},
        }
    )


def test_study_binding_round_trips_and_binds_role_and_run_commitments() -> None:
    study = _study()
    manifests = {"training": {"role": "training"}, "policy-selection": {"role": "policy-selection"}}
    run_commitments = {"train": {"summary.json": "a" * 64}, "select": {"summary.json": "b" * 64}}

    binding = build_study_binding(study, manifests, run_commitments)
    decoded = decode_study_binding(binding.document)
    validate_study_binding(decoded, study, manifests, run_commitments)

    assert decoded.dataset_manifest_sha256 == binding.dataset_manifest_sha256
    assert decoded.study_sha256 == study.sha256
    assert decoded.protocol_sha256 == study.protocol.sha256
    assert decoded.partition_manifest_sha256 == study.partitions.sha256


def test_study_binding_rejects_cross_study_and_commitment_tampering() -> None:
    study = _study()
    manifests = {"training": {"role": "training"}, "policy-selection": {"role": "policy-selection"}}
    runs = {"train": {"summary.json": "a" * 64}, "select": {"summary.json": "b" * 64}}
    binding = build_study_binding(study, manifests, runs)

    with pytest.raises(ValueError, match="study|binding"):
        validate_study_binding(replace(binding, study_sha256="f" * 64), study, manifests, runs)
    with pytest.raises(ValueError, match="role manifest|binding"):
        validate_study_binding(binding, study, {**manifests, "training": {"role": "changed"}}, runs)
    with pytest.raises(ValueError, match="run commitment|binding"):
        validate_study_binding(binding, study, manifests, {**runs, "train": {"summary.json": "c" * 64}})


def test_study_binding_decoder_rejects_unknown_fields_and_bad_dataset_digest() -> None:
    study = _study()
    binding = build_study_binding(
        study,
        {"training": {"role": "training"}, "policy-selection": {"role": "policy-selection"}},
        {"train": {"summary.json": "a" * 64}, "select": {"summary.json": "b" * 64}},
    )
    unknown = dict(binding.document)
    unknown["extra"] = True
    with pytest.raises(ValueError, match="fields"):
        decode_study_binding(unknown)
    wrong = dict(binding.document)
    wrong["dataset_manifest_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="dataset"):
        decode_study_binding(wrong)


def test_locked_declaration_is_value_free_and_rejects_wrong_members() -> None:
    study = _study()
    runs = {"held": {"root": "/custodian/run", "artifacts": {"summary.json": "a" * 64}}}
    declaration, commitments = build_locked_declaration_documents(study, ("held",), runs)

    assert declaration["members"][0]["key"] == "held"  # type: ignore[index]
    assert "provenance" not in declaration["members"][0]  # type: ignore[operator,index]
    assert commitments["runs"] == runs
    with pytest.raises(ValueError, match="keys"):
        build_locked_declaration_documents(study, ("validate",), runs)


def test_locked_payload_decoder_rejects_non_payload_documents() -> None:
    with pytest.raises(ValueError, match="payload fields"):
        decode_locked_payload(b"{}\n")
