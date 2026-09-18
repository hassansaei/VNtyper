"""Finite predeclared length hypotheses and reference simulation gate contracts."""

from copy import deepcopy
from importlib import import_module

import pytest

from vntyper.scripts.canonical_json import canonical_sha256

pytestmark = pytest.mark.unit


def protocol_document():
    return {
        "schema_version": "length-protocol-v1",
        "objective": "length-total-v1",
        "seed": 7,
        "fold_count": 3,
        "grouping_rule": "connected-leakage-groups-v1",
        "representative_preprocessing_priority": ["primary", "secondary"],
        "maximum_free_parameters": 2,
        "maximum_candidate_count": 4,
        "candidate_grid": [
            {"candidate_id": "affine-f", "model_kind": "affine-F"},
            {"candidate_id": "physical-f", "model_kind": "physical-F"},
        ],
        "baseline": "training-mean-v1",
        "tie_margin_repeat_units": 1,
        "required_strata": ["long", "nominal"],
        "declared_exclusions": [],
        "uncertainty": {
            "bootstrap_iterations": 10000,
            "bootstrap_interval": "percentile",
            "confidence": 0.95,
            "binomial_bound": "one-sided-exact",
        },
        "multiplicity": {"mandatory": "intersection-union", "exploratory": "holm"},
        "acceptance": {
            "minimum_independent_count": 60,
            "maximum_mae": 10,
            "minimum_relative_mae_improvement": 0.2,
            "paired_error_difference_upper_limit": 0,
            "tolerance_absolute": 10,
            "tolerance_relative": 0.1,
            "minimum_tolerance_lower_bound": 0.9,
            "minimum_availability_lower_bound": 0.95,
        },
        "qc": {
            "minimum_denominator_mean_depth": 10,
            "minimum_denominator_covered_fraction": 0.9,
            "minimum_denominator_supporting_fragments": 100,
        },
        "feature_bounds_rule": "training-min-max-expand-10-percent-v1",
        "truth_distribution_sha256": "a" * 64,
        "eligible_roster_sha256": "b" * 64,
    }


def test_protocol_is_hash_bound_and_preserves_priority_order():
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    result = m.decode_length_protocol(raw)
    assert result.candidates[0].model_kind == "affine-F"
    assert result.candidates[0].free_parameters == 2
    assert result.candidates[1].free_parameters == 0
    assert result.representative_preprocessing_priority == ("primary", "secondary")
    assert result.qc_sha256 == canonical_sha256(raw["qc"])
    assert m.length_protocol_document(result) == raw
    raw["candidate_grid"][0]["model_kind"] = "physical-A"
    assert result.candidates[0].model_kind == "affine-F"
    changed = protocol_document()
    changed["representative_preprocessing_priority"].reverse()
    assert m.decode_length_protocol(changed).sha256 != result.sha256


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", "calibration-protocol-v1"),
        ("objective", "caller-safety-v1"),
        ("seed", True),
        ("seed", -1),
        ("seed", 1.1),
        ("fold_count", 1),
        ("fold_count", False),
        ("grouping_rule", "per-read"),
        ("representative_preprocessing_priority", []),
        ("representative_preprocessing_priority", ["primary", "primary"]),
        ("maximum_free_parameters", 1),
        ("maximum_free_parameters", 3),
        ("maximum_candidate_count", 1),
        ("maximum_candidate_count", 5),
        ("baseline", "validation-mean"),
        ("tie_margin_repeat_units", 2),
        ("required_strata", []),
        ("required_strata", ["nominal", "long"]),
        ("declared_exclusions", [" "]),
        ("feature_bounds_rule", "validation-min-max"),
        ("truth_distribution_sha256", "A" * 64),
        ("eligible_roster_sha256", "B" * 64),
    ],
)
def test_protocol_rejects_invalid_or_unplanned_rules(field, value):
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    raw[field] = value
    with pytest.raises(ValueError):
        m.decode_length_protocol(raw)


@pytest.mark.parametrize(
    "grid",
    [
        [],
        [{"candidate_id": "joint", "model_kind": "affine-A+F"}],
        [{"candidate_id": "constant", "model_kind": "constant"}],
        [{"candidate_id": "affine-f", "model_kind": "affine-F", "coefficients": [1]}],
        [{"candidate_id": "same", "model_kind": "affine-A"}, {"candidate_id": "same", "model_kind": "affine-F"}],
        [{"candidate_id": "z", "model_kind": "affine-A"}, {"candidate_id": "a", "model_kind": "affine-F"}],
        [{"candidate_id": "a", "model_kind": "affine-A"}, {"candidate_id": "b", "model_kind": "affine-A"}],
    ],
)
def test_candidate_search_is_explicit_unique_finite_and_cannot_promote_baseline_or_joint(grid):
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    raw["candidate_grid"] = grid
    with pytest.raises(ValueError):
        m.decode_length_protocol(raw)


@pytest.mark.parametrize(
    "section,key,value",
    [
        ("uncertainty", "bootstrap_iterations", 1000),
        ("uncertainty", "bootstrap_interval", "normal"),
        ("uncertainty", "confidence", 0.9),
        ("uncertainty", "binomial_bound", "two-sided"),
        ("multiplicity", "mandatory", "any-gate"),
        ("multiplicity", "exploratory", "none"),
        ("acceptance", "minimum_independent_count", 0),
        ("acceptance", "maximum_mae", float("nan")),
        ("acceptance", "minimum_relative_mae_improvement", -0.1),
        ("acceptance", "paired_error_difference_upper_limit", 0.1),
        ("acceptance", "paired_error_difference_upper_limit", True),
        ("acceptance", "paired_error_difference_upper_limit", float("-inf")),
        ("acceptance", "tolerance_absolute", 0),
        ("acceptance", "tolerance_relative", 1.1),
        ("acceptance", "minimum_tolerance_lower_bound", True),
        ("acceptance", "minimum_availability_lower_bound", float("inf")),
        ("qc", "minimum_denominator_mean_depth", 0),
        ("qc", "minimum_denominator_covered_fraction", 1.1),
        ("qc", "minimum_denominator_supporting_fragments", True),
    ],
)
def test_nested_acceptance_and_uncertainty_are_validated_before_outcomes(section, key, value):
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    raw[section][key] = value
    with pytest.raises(ValueError):
        m.decode_length_protocol(raw)


@pytest.mark.parametrize("section", [None, "acceptance", "qc", "uncertainty", "multiplicity"])
def test_every_protocol_object_is_closed(section):
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    target = raw if section is None else raw[section]
    target["unknown"] = 1
    with pytest.raises(ValueError):
        m.decode_length_protocol(raw)


def test_predeclared_external_limits_are_hash_bound_without_becoming_reference_simulation_limits():
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    changed = deepcopy(raw)
    changed["acceptance"]["maximum_mae"] = 5
    assert m.decode_length_protocol(changed).sha256 != m.decode_length_protocol(raw).sha256


def test_stricter_signed_paired_error_limit_is_valid_and_hash_bound():
    m = import_module("vntyper.scripts.calibration_length_protocol")
    raw = protocol_document()
    raw["acceptance"]["paired_error_difference_upper_limit"] = -0.5
    result = m.decode_length_protocol(raw)
    assert result.acceptance["paired_error_difference_upper_limit"] == -0.5
    assert m.length_protocol_document(result) == raw


def test_protocol_document_refuses_forged_typed_digest():
    from dataclasses import replace

    m = import_module("vntyper.scripts.calibration_length_protocol")
    result = m.decode_length_protocol(protocol_document())
    with pytest.raises(ValueError, match="digest"):
        m.length_protocol_document(replace(result, sha256="0" * 64))


def test_protocol_document_refuses_mutable_typed_collections():
    from dataclasses import replace

    m = import_module("vntyper.scripts.calibration_length_protocol")
    result = m.decode_length_protocol(protocol_document())
    with pytest.raises(ValueError, match="immutable"):
        m.length_protocol_document(replace(result, acceptance=dict(result.acceptance)))
    with pytest.raises(ValueError, match="immutable"):
        m.length_protocol_document(replace(result, qc=dict(result.qc)))
