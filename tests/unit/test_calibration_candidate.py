"""Research candidates bind their target, provenance and complete payload."""

from copy import deepcopy
from importlib import import_module

import pytest

from vntyper.scripts.calibration_payload import decode_payload_manifest
from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def candidate_document(target="length"):
    applicability = {
        "domain": "synthetic",
        "assemblies": ["GRCh38"],
        "assay_classes": ["capture"],
        "input_scopes": ["full"],
        "preprocessing_ids": ["alignment-v1"],
    }
    if target == "callers":
        applicability["required_callers"] = ["advntr", "kestrel"]
    raw = {
        "schema_version": "calibration-candidate-v2",
        "target": target,
        "study_sha256": "a" * 64,
        "baseline_sha256": "b" * 64,
        "partition_sha256": "c" * 64,
        "training_evidence_sha256": "d" * 64,
        "selection_evidence_sha256": "e" * 64,
        "payload_sha256": "f" * 64,
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
    return resign(raw)


def resign(raw):
    raw["candidate_id"] = canonical_sha256({key: value for key, value in raw.items() if key != "candidate_id"})
    return raw


@pytest.mark.parametrize("target", ["callers", "length"])
def test_candidate_roundtrip_and_nested_immutability(target):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document(target)
    expected = deepcopy(raw)
    decoded = candidates.decode_candidate(raw)
    assert decoded.target == target
    assert decoded.candidate_id == raw["candidate_id"]
    assert decoded.sha256 == canonical_sha256(raw)
    assert decoded.applicability.assemblies == ("GRCh38",)
    assert decoded.producer.tool_versions["depth-tool"] == "1.0"
    raw["applicability"]["assemblies"].append("GRCh37")
    raw["producer"]["tool_versions"]["depth-tool"] = "2.0"
    assert candidates.candidate_document(decoded) == expected
    with pytest.raises(TypeError):
        decoded.producer.tool_versions["depth-tool"] = "2.0"


def test_candidate_identity_changes_with_every_bound_input():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    baseline = candidates.decode_candidate(candidate_document())
    for key in (
        "study_sha256",
        "baseline_sha256",
        "partition_sha256",
        "training_evidence_sha256",
        "selection_evidence_sha256",
        "payload_sha256",
    ):
        changed = candidate_document()
        changed[key] = "1" * 64
        decoded = candidates.decode_candidate(resign(changed))
        assert decoded.candidate_id != baseline.candidate_id
        assert decoded.sha256 != baseline.sha256


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", "calibration-candidate-v1"),
        ("target", "dominance"),
        ("status", "validated"),
        ("candidate_id", "0" * 64),
        ("payload_sha256", "f" * 63),
        ("baseline_sha256", "Z" * 64),
    ],
)
def test_candidate_rejects_unknown_version_target_status_and_stale_identity(field, value):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    raw[field] = value
    with pytest.raises(ValueError, match=field):
        candidates.decode_candidate(raw)


@pytest.mark.parametrize("where", ["root", "applicability", "producer"])
@pytest.mark.parametrize("change", ["missing", "extra"])
def test_candidate_is_closed_at_each_object(where, change):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    obj = raw if where == "root" else raw[where]
    if change == "extra":
        obj["extra"] = True
    else:
        del obj[next(iter(obj))]
    with pytest.raises(ValueError, match="fields"):
        candidates.decode_candidate(raw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("domain", "production"),
        ("assemblies", []),
        ("assemblies", ["GRCh38", "GRCh37"]),
        ("assemblies", ["GRCh38", "GRCh38"]),
        ("assay_classes", [""]),
        ("input_scopes", ["unknown"]),
        ("preprocessing_ids", [True]),
    ],
)
def test_applicability_rejects_ambiguous_or_empty_domains(field, value):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    raw["applicability"][field] = value
    with pytest.raises(ValueError, match=field):
        candidates.decode_candidate(raw)


def test_length_candidate_cannot_claim_caller_composition():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    raw["applicability"]["required_callers"] = ["kestrel"]
    with pytest.raises(ValueError, match="fields"):
        candidates.decode_candidate(raw)


@pytest.mark.parametrize("callers", [[], ["advntr"], ["unknown"], ["kestrel", "advntr"], ["kestrel", "kestrel"]])
def test_caller_set_is_sorted_and_includes_pipeline_kestrel(callers):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document("callers")
    raw["applicability"]["required_callers"] = callers
    with pytest.raises(ValueError, match="required_callers"):
        candidates.decode_candidate(raw)


def test_kestrel_only_candidate_is_valid_and_identified_build_is_optional_for_research():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document("callers")
    raw["applicability"]["required_callers"] = ["kestrel"]
    raw["producer"]["source_revision"] = None
    decoded = candidates.decode_candidate(resign(raw))
    assert decoded.applicability.required_callers == ("kestrel",)
    assert decoded.producer.source_revision is None
    assert decoded.status == "research-only"


@pytest.mark.parametrize(
    "field,value",
    [
        ("name", ""),
        ("version", None),
        ("source_revision", "abc"),
        ("source_revision", "A" * 40),
        ("tool_versions", {}),
        ("tool_versions", {"tool": ""}),
        ("tool_versions", {1: "v1"}),
        ("feature_schema_sha256", "bad"),
    ],
)
def test_producer_provenance_is_strict(field, value):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    raw["producer"][field] = value
    with pytest.raises(ValueError, match=field):
        candidates.decode_candidate(raw)


def test_payload_binding_checks_target_and_entire_manifest():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    manifest = decode_payload_manifest([{"path": "length-model.json", "size_bytes": 5, "sha256": "a" * 64}])
    raw = candidate_document()
    raw["payload_sha256"] = manifest.sha256
    decoded = candidates.decode_candidate(resign(raw))
    assert candidates.validate_candidate_payload(decoded, manifest, expected_target="length") is None
    with pytest.raises(ValueError, match="target"):
        candidates.validate_candidate_payload(decoded, manifest, expected_target="callers")
    with pytest.raises(ValueError, match="payload"):
        candidates.validate_candidate_payload(
            candidates.decode_candidate(candidate_document()), manifest, expected_target="length"
        )


@pytest.mark.parametrize("raw", [None, [], "candidate"])
def test_candidate_requires_an_object(raw):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    with pytest.raises(ValueError, match="fields"):
        candidates.decode_candidate(raw)


def test_projection_and_binding_refuse_untyped_objects():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    with pytest.raises(ValueError, match="CandidateEnvelope"):
        candidates.candidate_document(candidate_document())
    with pytest.raises(ValueError, match="CandidateEnvelope"):
        candidates.validate_candidate_payload({}, None, expected_target="length")
    with pytest.raises(ValueError, match="PayloadManifest"):
        candidates.validate_candidate_payload(
            candidates.decode_candidate(candidate_document()), None, expected_target="length"
        )
