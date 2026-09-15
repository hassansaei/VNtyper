"""Strict non-executable total-length research model contracts."""

from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from importlib import import_module

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def model_document(feature="A"):
    return {
        "schema_version": "length-model-v1",
        "target": {
            "name": "total_diploid_repeat_count",
            "boundary_definition": "complete-core-plus-invariant-units-v1",
        },
        "unit": "repeat_units",
        "model_kind": f"affine-{feature}",
        "feature_order": [feature],
        "intercept": 10.0,
        "coefficients": [50.0],
        "annotation_sha256": "a" * 64,
        "counting_policy_sha256": "b" * 64,
        "applicability": {
            "domain": "synthetic",
            "assemblies": ["synthetic-build-v1"],
            "assay_classes": ["synthetic-short-read"],
            "input_scopes": ["regional"],
            "preprocessing_ids": ["synthetic-preprocessing-v1"],
            "aligner_name": "synthetic-aligner",
            "aligner_version": "1.0",
            "aligner_arguments_sha256": "c" * 64,
            "primary_secondary_marking": "primary-only",
            "counting_policy_sha256": "b" * 64,
        },
        "qc": {
            "minimum_denominator_mean_depth": 10.0,
            "minimum_denominator_covered_fraction": 0.9,
            "minimum_denominator_supporting_fragments": 100,
            "fragment_evidence_kind": "read-pair-identity-qc-proxy",
        },
        "feature_bounds": {feature: {"minimum": 0.5, "maximum": 3.0}},
        "study_sha256": "d" * 64,
        "training_evidence_sha256": "e" * 64,
        "producer": {
            "name": "synthetic-length-fit",
            "version": "1.0",
            "source_revision": "f" * 40,
            "tool_versions": {"numpy": "1.26"},
            "feature_schema_sha256": "0" * 64,
        },
    }


@pytest.mark.parametrize("feature", ["A", "F"])
def test_affine_model_roundtrip_is_hash_bound_deeply_immutable_and_research_only(feature):
    m = import_module("vntyper.scripts.length_model")
    raw = model_document(feature)
    expected = deepcopy(raw)
    model = m.decode_length_model(raw)
    assert model.sha256 == canonical_sha256(expected)
    assert model.feature_order == (feature,)
    assert model.coefficients == (50.0,)
    assert model.target_name == "total_diploid_repeat_count"
    assert m.encode_length_model(model) == expected
    assert not hasattr(model, "promotion_eligible")
    raw["coefficients"][0] = 99
    assert model.coefficients == (50.0,)
    with pytest.raises(TypeError):
        model.feature_bounds[feature] = model.feature_bounds[feature]
    with pytest.raises(FrozenInstanceError):
        model.intercept = 0


@pytest.mark.parametrize("where", ["root", "target", "applicability", "qc", "bounds", "producer"])
@pytest.mark.parametrize("change", ["missing", "extra"])
def test_every_model_object_is_closed(where, change):
    m = import_module("vntyper.scripts.length_model")
    raw = model_document()
    targets = {
        "root": raw,
        "target": raw["target"],
        "applicability": raw["applicability"],
        "qc": raw["qc"],
        "bounds": raw["feature_bounds"]["A"],
        "producer": raw["producer"],
    }
    target = targets[where]
    if change == "missing":
        del target[next(iter(target))]
    else:
        target["extra"] = True
    with pytest.raises(ValueError, match="fields"):
        m.decode_length_model(raw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", "length-model-v0"),
        ("target", {"name": "allele_length", "boundary_definition": "complete-core-plus-invariant-units-v1"}),
        ("target", {"name": "total_diploid_repeat_count", "boundary_definition": "ambiguous"}),
        ("unit", "bp"),
        ("model_kind", "affine-A+F"),
        ("model_kind", []),
        ("model_kind", {}),
        ("feature_order", ["F"]),
        ("intercept", True),
        ("intercept", float("nan")),
        ("intercept", 10**1000),
        ("coefficients", ["50 * A"]),
        ("coefficients", [float("inf")]),
        ("coefficients", [1, 2]),
        ("annotation_sha256", "A" * 64),
        ("study_sha256", "bad"),
    ],
)
def test_model_rejects_wrong_identity_expressions_booleans_and_nonfinite_values(field, value):
    m = import_module("vntyper.scripts.length_model")
    raw = model_document()
    raw[field] = value
    with pytest.raises(ValueError):
        m.decode_length_model(raw)


@pytest.mark.parametrize("kind", ["physical-A", "physical-F"])
def test_physical_models_fail_until_geometry_conversion_evidence_is_implemented(kind):
    m = import_module("vntyper.scripts.length_model")
    raw = model_document(kind[-1])
    raw["model_kind"] = kind
    with pytest.raises(ValueError, match="physical_model_geometry_evidence_unsupported"):
        m.decode_length_model(raw)


@pytest.mark.parametrize(
    "section,field,value",
    [
        ("qc", "minimum_denominator_mean_depth", 0),
        ("qc", "minimum_denominator_covered_fraction", 1.1),
        ("qc", "minimum_denominator_supporting_fragments", True),
        ("qc", "fragment_evidence_kind", "read-count"),
        ("feature_bounds", "minimum", True),
        ("feature_bounds", "minimum", float("-inf")),
        ("feature_bounds", "maximum", float("inf")),
        ("feature_bounds", "maximum", 0.5),
    ],
)
def test_qc_and_feature_bounds_are_strict(section, field, value):
    m = import_module("vntyper.scripts.length_model")
    raw = model_document()
    target = raw["qc"] if section == "qc" else raw["feature_bounds"]["A"]
    target[field] = value
    with pytest.raises(ValueError):
        m.decode_length_model(raw)


def test_model_requires_one_matching_bound_and_counting_policy_identity():
    m = import_module("vntyper.scripts.length_model")
    raw = model_document()
    raw["feature_bounds"] = {"F": {"minimum": 0.5, "maximum": 3.0}}
    with pytest.raises(ValueError, match="feature_bounds"):
        m.decode_length_model(raw)
    raw = model_document()
    raw["applicability"]["counting_policy_sha256"] = "9" * 64
    with pytest.raises(ValueError, match="counting_policy_sha256"):
        m.decode_length_model(raw)


@pytest.mark.parametrize("value", [None, b"pickle", "A * 50 + 10", [model_document()]])
def test_model_decoder_accepts_only_a_json_object(value):
    m = import_module("vntyper.scripts.length_model")
    with pytest.raises(ValueError, match="fields"):
        m.decode_length_model(value)


def test_model_encoder_rejects_forged_or_mutable_typed_values():
    m = import_module("vntyper.scripts.length_model")
    model = m.decode_length_model(model_document())
    with pytest.raises(ValueError, match="digest"):
        m.encode_length_model(replace(model, sha256="0" * 64))
    with pytest.raises(ValueError, match="immutable"):
        m.encode_length_model(replace(model, feature_bounds=dict(model.feature_bounds)))
