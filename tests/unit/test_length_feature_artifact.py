"""Closed serialized features retain complete, validated measurement provenance."""

import json
from copy import deepcopy
from dataclasses import replace
from importlib import import_module

import pytest

from tests.unit.test_length_estimation import features
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_feature_provenance import encode_length_feature_provenance
from vntyper.scripts.length_features import encode_length_features

pytestmark = pytest.mark.unit


def _module():
    return import_module("vntyper.scripts.length_feature_artifact")


def _document(measured=None):
    measured = features() if measured is None else measured
    return {
        "schema_version": "length-feature-artifact-v1",
        "features": encode_length_features(measured),
        "provenance": encode_length_feature_provenance(measured.provenance),
    }


@pytest.mark.parametrize(
    "options",
    [{}, {"missing_region": "CORE"}, {"depths": (0, 0, 0, 0, 0, 0)}, {"fragment_evidence": False}],
)
def test_roundtrip_preserves_values_unavailability_and_exact_hashes(options):
    measured = features(**options)
    encoded = _module().encode_length_feature_artifact(measured)
    assert encoded == _document(measured)
    decoded = _module().decode_length_feature_artifact(json.loads(json.dumps(encoded)))
    assert decoded == measured
    assert decoded.sha256 == canonical_sha256(encoded["features"])
    assert decoded.provenance_sha256 == canonical_sha256(encoded["provenance"])
    assert decoded.provenance == measured.provenance


@pytest.mark.parametrize("level", ["bundle", "features", "row", "regions", "region"])
@pytest.mark.parametrize("change", ["extra", "missing", "nonobject"])
def test_every_document_layer_is_closed(level, change):
    raw = _document()
    container, key = {
        "bundle": ({"root": raw}, "root"),
        "features": (raw, "features"),
        "row": (raw["features"]["rows"], 0),
        "regions": (raw["features"]["rows"][0], "regions"),
        "region": (raw["features"]["rows"][0]["regions"], "CORE"),
    }[level]
    if change == "extra":
        container[key]["undeclared"] = 1
    elif change == "missing":
        del container[key][next(iter(container[key]))]
    else:
        container[key] = []
    if level == "bundle":
        raw = container[key]
    with pytest.raises(ValueError):
        _module().decode_length_feature_artifact(raw)


@pytest.mark.parametrize("rows", [[], [None], [None, None], {}, ()])
def test_bundle_requires_exactly_one_feature_row(rows):
    raw = _document()
    raw["features"]["rows"] = rows
    with pytest.raises(ValueError):
        _module().decode_length_feature_artifact(raw)


@pytest.mark.parametrize("section", ["bundle", "features"])
def test_unsupported_schema_refused(section):
    raw = _document()
    (raw if section == "bundle" else raw["features"])["schema_version"] = "future"
    with pytest.raises(ValueError):
        _module().decode_length_feature_artifact(raw)


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), "2", 3])
def test_ratios_are_finite_numeric_and_match_regions(value):
    raw = _document()
    raw["features"]["rows"][0]["A"] = value
    with pytest.raises(ValueError):
        _module().decode_length_feature_artifact(raw)


def test_provenance_cannot_be_swapped_between_measurements():
    raw = _document()
    raw["provenance"] = encode_length_feature_provenance(features(manifest_key="other").provenance)
    with pytest.raises(ValueError, match="provenance"):
        _module().decode_length_feature_artifact(raw)


def test_derived_mean_and_denominator_qc_must_remain_consistent():
    for field, value in (("mean_depth", 199), ("supporting_fragment_count", 1), ("length_bp", True)):
        raw = _document()
        raw["features"]["rows"][0]["regions"]["CORE"][field] = value
        with pytest.raises(ValueError):
            _module().decode_length_feature_artifact(raw)


@pytest.mark.parametrize(
    "field,value", [("status", 1), ("status", "measured-later"), ("reasons", ()), ("reasons", [1])]
)
def test_status_and_reason_types_are_closed(field, value):
    raw = _document()
    raw["features"]["rows"][0][field] = value
    with pytest.raises(ValueError):
        _module().decode_length_feature_artifact(raw)


def test_output_is_independent_of_mutable_input_and_encoder_checks_typed_forgery():
    raw = _document()
    original = deepcopy(raw)
    decoded = _module().decode_length_feature_artifact(raw)
    raw["features"]["rows"][0]["regions"]["CORE"]["mean_depth"] = 1
    assert _module().encode_length_feature_artifact(decoded) == original
    with pytest.raises(ValueError):
        _module().encode_length_feature_artifact(replace(decoded, sha256="0" * 64))
