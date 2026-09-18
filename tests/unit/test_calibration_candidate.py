"""Research candidates bind their target, provenance and complete payload."""

from copy import deepcopy
from dataclasses import replace
from importlib import import_module

import pytest

from vntyper.scripts.calibration_payload import (
    CALLER_BUNDLE_DESCRIPTOR_PATH,
    decode_caller_bundle_descriptor,
    decode_payload_manifest,
    payload_manifest_document,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256

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


def test_public_applicability_and_producer_helpers_preserve_candidate_contracts():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document("length")
    applicability = candidates.decode_candidate_applicability(raw["applicability"], target="length")
    producer = candidates.decode_candidate_producer(raw["producer"])
    assert candidates.candidate_applicability_document(applicability, target="length") == raw["applicability"]
    assert candidates.candidate_producer_document(producer) == raw["producer"]

    mutable = replace(producer, tool_versions=dict(producer.tool_versions))
    with pytest.raises(ValueError, match="immutable"):
        candidates.candidate_producer_document(mutable)
    with pytest.raises(ValueError, match="applicability"):
        candidates.candidate_applicability_document({}, target="length")
    with pytest.raises(ValueError, match="target"):
        candidates.decode_candidate_applicability(raw["applicability"], target=[])  # type: ignore[arg-type]


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


@pytest.mark.parametrize(
    "field,value",
    [
        ("aligner_name", ""),
        ("aligner_version", None),
        ("aligner_arguments_sha256", "bad"),
        ("primary_secondary_marking", True),
        ("counting_policy_sha256", "F" * 64),
    ],
)
def test_length_applicability_requires_explicit_aligner_and_counting_identity(field, value):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    raw["applicability"][field] = value
    with pytest.raises(ValueError, match=field):
        candidates.decode_candidate(raw)


@pytest.mark.parametrize(
    "field",
    [
        "aligner_name",
        "aligner_version",
        "aligner_arguments_sha256",
        "primary_secondary_marking",
        "counting_policy_sha256",
    ],
)
def test_length_applicability_requires_every_measurement_identity_field(field):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document()
    del raw["applicability"][field]
    with pytest.raises(ValueError, match="fields"):
        candidates.decode_candidate(raw)


def test_caller_candidate_cannot_claim_length_measurement_identity():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document("callers")
    raw["applicability"]["aligner_name"] = "synthetic-aligner"
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


def caller_bundle_payload(callers=("advntr", "kestrel"), *, include_background=False):
    descriptor_raw = {
        "schema_version": "caller-bundle-v2",
        "required_callers": list(callers),
        "components": {
            "decision-profile.json": "a" * 64,
            "advntr-policy.json": "b" * 64 if "advntr" in callers else None,
            "background.json": "c" * 64 if include_background else None,
        },
    }
    descriptor = decode_caller_bundle_descriptor(descriptor_raw)
    rows = [
        {
            "path": CALLER_BUNDLE_DESCRIPTOR_PATH,
            "size_bytes": len(canonical_json_bytes(descriptor_raw)),
            "sha256": descriptor.sha256,
        },
        {"path": "decision-profile.json", "size_bytes": 102, "sha256": "a" * 64},
    ]
    if "advntr" in callers:
        rows.append({"path": "advntr-policy.json", "size_bytes": 103, "sha256": "b" * 64})
    if include_background:
        rows.append({"path": "background.json", "size_bytes": 104, "sha256": "c" * 64})
    return descriptor, decode_payload_manifest(sorted(rows, key=lambda row: row["path"]))


def bound_candidate(target, manifest, *, callers=("advntr", "kestrel")):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    raw = candidate_document(target)
    raw["payload_sha256"] = manifest.sha256
    if target == "callers":
        raw["applicability"]["required_callers"] = list(callers)
    return candidates.decode_candidate(resign(raw))


@pytest.mark.parametrize(
    "callers,include_background",
    [(("kestrel",), False), (("advntr", "kestrel"), False), (("advntr", "kestrel"), True)],
)
def test_caller_payload_requires_matching_closed_descriptor_and_conditional_files(callers, include_background):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    descriptor, manifest = caller_bundle_payload(callers, include_background=include_background)
    candidate = bound_candidate("callers", manifest, callers=callers)

    assert (
        candidates.validate_candidate_payload(
            candidate, manifest, expected_target="callers", caller_descriptor=descriptor
        )
        is None
    )


def test_caller_payload_rejects_missing_or_mismatched_descriptor():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    descriptor, manifest = caller_bundle_payload()
    candidate = bound_candidate("callers", manifest)
    with pytest.raises(ValueError, match="descriptor"):
        candidates.validate_candidate_payload(candidate, manifest, expected_target="callers")

    kestrel_descriptor, _ = caller_bundle_payload(("kestrel",))
    with pytest.raises(ValueError, match="required_callers"):
        candidates.validate_candidate_payload(
            candidate, manifest, expected_target="callers", caller_descriptor=kestrel_descriptor
        )

    forged = replace(descriptor, decision_profile_sha256="d" * 64)
    with pytest.raises(ValueError, match="canonical"):
        candidates.validate_candidate_payload(candidate, manifest, expected_target="callers", caller_descriptor=forged)


@pytest.mark.parametrize("path", ["extra.json", "advntr-policy.json"])
def test_caller_payload_rejects_files_outside_descriptor_layout(path):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    descriptor, manifest = caller_bundle_payload(("kestrel",))
    rows = payload_manifest_document(manifest)
    rows.append({"path": path, "size_bytes": 1, "sha256": "d" * 64})
    changed = decode_payload_manifest(sorted(rows, key=lambda row: row["path"]))
    candidate = bound_candidate("callers", changed, callers=("kestrel",))
    with pytest.raises(ValueError, match="file set"):
        candidates.validate_candidate_payload(
            candidate, changed, expected_target="callers", caller_descriptor=descriptor
        )


@pytest.mark.parametrize("path", [CALLER_BUNDLE_DESCRIPTOR_PATH, "advntr-policy.json"])
def test_caller_payload_rejects_a_missing_required_file(path):
    candidates = import_module("vntyper.scripts.calibration_candidate")
    descriptor, manifest = caller_bundle_payload()
    rows = [row for row in payload_manifest_document(manifest) if row["path"] != path]
    changed = decode_payload_manifest(rows)
    candidate = bound_candidate("callers", changed)
    with pytest.raises(ValueError, match="file set"):
        candidates.validate_candidate_payload(
            candidate, changed, expected_target="callers", caller_descriptor=descriptor
        )


def test_caller_payload_rejects_component_or_descriptor_digest_mismatch():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    descriptor, manifest = caller_bundle_payload()
    for path in ("decision-profile.json", CALLER_BUNDLE_DESCRIPTOR_PATH):
        rows = payload_manifest_document(manifest)
        for row in rows:
            if row["path"] == path:
                row["sha256"] = "d" * 64
        changed = decode_payload_manifest(rows)
        candidate = bound_candidate("callers", changed)
        with pytest.raises(ValueError, match="digest"):
            candidates.validate_candidate_payload(
                candidate, changed, expected_target="callers", caller_descriptor=descriptor
            )

    rows = payload_manifest_document(manifest)
    next(row for row in rows if row["path"] == CALLER_BUNDLE_DESCRIPTOR_PATH)["size_bytes"] += 1
    changed = decode_payload_manifest(rows)
    candidate = bound_candidate("callers", changed)
    with pytest.raises(ValueError, match="size"):
        candidates.validate_candidate_payload(
            candidate, changed, expected_target="callers", caller_descriptor=descriptor
        )


def test_length_payload_rejects_caller_descriptor():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    descriptor = decode_caller_bundle_descriptor(
        {
            "schema_version": "caller-bundle-v2",
            "required_callers": ["kestrel"],
            "components": {
                "decision-profile.json": "a" * 64,
                "advntr-policy.json": None,
                "background.json": None,
            },
        }
    )
    manifest = decode_payload_manifest([{"path": "length-model.json", "size_bytes": 5, "sha256": "a" * 64}])
    candidate = bound_candidate("length", manifest)
    with pytest.raises(ValueError, match="length.*descriptor"):
        candidates.validate_candidate_payload(
            candidate, manifest, expected_target="length", caller_descriptor=descriptor
        )


def test_public_boundaries_revalidate_candidate_and_payload_content():
    candidates = import_module("vntyper.scripts.calibration_candidate")
    manifest = decode_payload_manifest([{"path": "length-model.json", "size_bytes": 5, "sha256": "a" * 64}])
    candidate = bound_candidate("length", manifest)

    forged_candidate = replace(candidate, candidate_id="0" * 64, sha256="0" * 64)
    with pytest.raises(ValueError, match="candidate_id"):
        candidates.candidate_document(forged_candidate)
    with pytest.raises(ValueError, match="candidate_id"):
        candidates.validate_candidate_payload(forged_candidate, manifest, expected_target="length")

    forged_manifest = replace(manifest, sha256="f" * 64)
    with pytest.raises(ValueError, match="canonical"):
        candidates.validate_candidate_payload(candidate, forged_manifest, expected_target="length")

    mutable_producer = replace(candidate.producer, tool_versions={"depth-tool": "1.0"})
    with pytest.raises(ValueError, match="decoded immutable"):
        candidates.candidate_document(replace(candidate, producer=mutable_producer))


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
