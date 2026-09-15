"""Paired calibration truth and explicit length conversion contracts."""

from __future__ import annotations

import copy
import json
from collections.abc import Callable
from dataclasses import FrozenInstanceError
from importlib.resources import files
from typing import Any

import pytest

from tests.unit.test_calibration_intake_contract import synthetic_intake
from vntyper.scripts.calibration_intake_contract import TruthRecord, decode_intake
from vntyper.scripts.calibration_truth import (
    converted_length_target,
    decode_conversion_registry,
    length_target,
)

pytestmark = pytest.mark.unit


def _truth(
    allele_1: int | float | None,
    allele_2: int | float | None,
    *,
    unit: str = "repeat-count",
    boundary_definition: str = "target-boundary-v1",
    measurement: str = "exact",
    conversion_id: str | None = None,
    status: str = "confirmed",
    repeat_unit_bp: int = 60,
) -> TruthRecord:
    raw = synthetic_intake()
    raw_truth = raw["truth"][0]  # type: ignore[index]
    raw_truth["status"] = status
    raw_truth["length"] = {
        "allele_1": allele_1,
        "allele_2": allele_2,
        "unit": unit,
        "repeat_unit_bp": repeat_unit_bp,
        "boundary_definition": boundary_definition,
        "measurement": measurement,
        "lower_bound": 100 if measurement == "interval" else None,
        "upper_bound": 160 if measurement in {"interval", "censored"} else None,
        "conversion_id": conversion_id,
    }
    return decode_intake(raw).truth[0]


def _registry() -> dict[str, Any]:
    return {
        "schema_version": "calibration-length-conversions-v1",
        "conversions": [
            {
                "conversion_id": "bp-source-to-target-v1",
                "source_unit": "bp",
                "source_boundary_definition": "bp-boundary-v1",
                "target_boundary_definition": "target-boundary-v1",
                "repeat_unit_bp": 60,
                "allele_1_nonrepeat_bp": 30,
                "allele_2_nonrepeat_bp": 90,
            },
            {
                "conversion_id": "legacy-count-to-target-v1",
                "source_unit": "repeat-count",
                "source_boundary_definition": "legacy-count-boundary-v1",
                "target_boundary_definition": "target-boundary-v1",
                "repeat_unit_bp": 60,
                "allele_1_offset": 1,
                "allele_2_offset": 2,
            },
        ],
    }


def _additive_registry() -> dict[str, Any]:
    return {
        "schema_version": "calibration-length-conversions-v2",
        "conversions": [
            {
                "conversion_id": "canonical-only-plus-nine-terminals-v1",
                "source_unit": "repeat-count",
                "source_boundary_definition": "canonical-variable-only-units-v1",
                "target_boundary_definition": "complete-core-plus-invariant-units-v1",
                "repeat_unit_bp": 60,
                "allele_1_adjustment": 9,
                "allele_2_adjustment": 9,
            }
        ],
    }


@pytest.mark.parametrize(
    ("alleles", "expected"),
    [
        ((40, 90), 130.0),
        ((90, 40), 130.0),
        ((100, 180), 280.0),
        ((180, 100), 280.0),
    ],
)
def test_length_target_sums_exact_repeat_count_pairs_without_ordering(
    alleles: tuple[int, int], expected: float
) -> None:
    truth = _truth(*alleles)

    assert length_target(truth) == expected


@pytest.mark.parametrize("alleles", [(40, None), (None, 90), (None, None)])
def test_length_target_returns_none_for_an_incomplete_exact_pair(alleles: tuple[int | None, int | None]) -> None:
    assert length_target(_truth(*alleles)) is None


@pytest.mark.parametrize(
    ("measurement", "status"),
    [("interval", "confirmed"), ("censored", "confirmed"), ("exact", "disputed")],
)
def test_length_target_preserves_but_excludes_nonexact_or_unconfirmed_truth(measurement: str, status: str) -> None:
    truth = _truth(
        None if measurement != "exact" else 40,
        None if measurement != "exact" else 90,
        measurement=measurement,
        status=status,
    )

    assert truth.length is not None
    assert truth.length.measurement == measurement
    assert truth.status == status
    assert length_target(truth) is None


def test_length_target_rejects_a_unit_or_boundary_that_requires_conversion() -> None:
    with pytest.raises(ValueError, match="conversion"):
        length_target(
            _truth(2430, 5490, unit="bp", boundary_definition="bp-boundary-v1", conversion_id="bp-source-to-target-v1")
        )
    with pytest.raises(ValueError, match="conversion"):
        length_target(
            _truth(41, 92, boundary_definition="legacy-count-boundary-v1", conversion_id="legacy-count-to-target-v1")
        )


def test_bp_conversion_subtracts_each_declared_boundary_and_divides_by_the_declared_unit() -> None:
    registry = decode_conversion_registry(_registry())
    truth = _truth(
        2430,
        5490,
        unit="bp",
        boundary_definition="bp-boundary-v1",
        conversion_id="bp-source-to-target-v1",
    )

    assert converted_length_target(truth, registry, target_boundary_definition="target-boundary-v1") == 130.0
    assert len(registry.sha256) == 64
    with pytest.raises(FrozenInstanceError):
        registry.conversions["bp-source-to-target-v1"].repeat_unit_bp = 61  # type: ignore[misc]
    with pytest.raises(TypeError):
        registry.conversions["new"] = registry.conversions["bp-source-to-target-v1"]  # type: ignore[index]


def test_repeat_count_boundary_conversion_applies_exact_paired_offsets() -> None:
    registry = decode_conversion_registry(_registry())
    truth = _truth(
        41,
        92,
        boundary_definition="legacy-count-boundary-v1",
        conversion_id="legacy-count-to-target-v1",
    )

    assert converted_length_target(truth, registry, target_boundary_definition="target-boundary-v1") == 130.0


def test_additive_repeat_count_conversion_can_include_nine_terminal_units_per_allele() -> None:
    registry = decode_conversion_registry(_additive_registry())
    truth = _truth(
        40,
        90,
        boundary_definition="canonical-variable-only-units-v1",
        conversion_id="canonical-only-plus-nine-terminals-v1",
    )

    assert registry.schema_version == "calibration-length-conversions-v2"
    assert (
        converted_length_target(
            truth,
            registry,
            target_boundary_definition="complete-core-plus-invariant-units-v1",
        )
        == 148.0
    )


def test_packaged_terminal_count_conversion_is_the_reviewed_additive_contract() -> None:
    path = files("vntyper").joinpath("data/length/grch38-count-conversions-v2.json")
    raw = json.loads(path.read_text(encoding="utf-8"))

    assert raw == _additive_registry()
    assert decode_conversion_registry(raw).schema_version == "calibration-length-conversions-v2"


@pytest.mark.parametrize("adjustment", [True, 0.5, float("inf")])
def test_additive_repeat_count_conversion_requires_finite_integral_adjustments(adjustment: object) -> None:
    raw = _additive_registry()
    raw["conversions"][0]["allele_1_adjustment"] = adjustment

    with pytest.raises(ValueError, match="adjustment"):
        decode_conversion_registry(raw)


def test_repeat_count_boundary_conversion_binds_the_repeat_unit() -> None:
    registry = decode_conversion_registry(_registry())
    truth = _truth(
        41,
        92,
        boundary_definition="legacy-count-boundary-v1",
        conversion_id="legacy-count-to-target-v1",
        repeat_unit_bp=61,
    )

    with pytest.raises(ValueError, match="repeat unit"):
        converted_length_target(truth, registry, target_boundary_definition="target-boundary-v1")


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda registry: registry.update(extra=True), "fields differ"),
        (lambda registry: registry["conversions"][0].update(extra=True), "fields differ"),
        (lambda registry: registry["conversions"][0].update(repeat_unit_bp=True), "repeat unit"),
        (lambda registry: registry["conversions"][0].update(allele_1_nonrepeat_bp=float("inf")), "nonrepeat"),
        (lambda registry: registry["conversions"][0].update(conversion_id=""), "conversion id"),
    ],
)
def test_conversion_registry_rejects_unknown_bool_nonfinite_and_empty_fields(
    mutate: Callable[[dict[str, Any]], None], message: str
) -> None:
    raw = _registry()
    mutate(raw)

    with pytest.raises(ValueError, match=message):
        decode_conversion_registry(raw)


def test_conversion_registry_rejects_duplicate_ids_and_unknown_aliases() -> None:
    duplicate = _registry()
    duplicate["conversions"].append(copy.deepcopy(duplicate["conversions"][0]))
    aliased = _registry()
    aliased["conversions"][0]["source_unit"] = "base-pairs"  # type: ignore[index]

    with pytest.raises(ValueError, match="conversion ids.*unique"):
        decode_conversion_registry(duplicate)
    with pytest.raises(ValueError, match="source unit"):
        decode_conversion_registry(aliased)


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ({"schema_version": "calibration-length-conversions-v3", "conversions": []}, "schema version"),
        ({"schema_version": [], "conversions": []}, "schema version"),
        ({"schema_version": "calibration-length-conversions-v1", "conversions": []}, "non-empty list"),
        ({"schema_version": "calibration-length-conversions-v1", "conversions": [42]}, "must be an object"),
    ],
)
def test_conversion_registry_rejects_wrong_schema_empty_registry_and_nonobjects(
    value: dict[str, object], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        decode_conversion_registry(value)


@pytest.mark.parametrize(
    ("truth_kwargs", "message"),
    [
        ({"allele_1": 2431, "allele_2": 5490}, "integral"),
        ({"allele_1": 30, "allele_2": 5490}, "positive"),
        ({"allele_1": True, "allele_2": 5490}, "allele 1"),
    ],
)
def test_bp_conversion_refuses_rounding_nonpositive_and_bool_measurements(
    truth_kwargs: dict[str, Any], message: str
) -> None:
    registry = decode_conversion_registry(_registry())
    allele_1 = truth_kwargs["allele_1"]
    allele_2 = truth_kwargs["allele_2"]

    if isinstance(allele_1, bool):
        with pytest.raises(ValueError, match=message):
            _truth(
                allele_1,
                allele_2,
                unit="bp",
                boundary_definition="bp-boundary-v1",
                conversion_id="bp-source-to-target-v1",
            )
        return
    truth = _truth(
        allele_1,
        allele_2,
        unit="bp",
        boundary_definition="bp-boundary-v1",
        conversion_id="bp-source-to-target-v1",
    )
    with pytest.raises(ValueError, match=message):
        converted_length_target(truth, registry, target_boundary_definition="target-boundary-v1")


def test_conversion_requires_exact_registry_and_target_bindings() -> None:
    registry = decode_conversion_registry(_registry())
    truth = _truth(
        2430,
        5490,
        unit="bp",
        boundary_definition="bp-boundary-v1",
        conversion_id="bp-source-to-target-v1",
    )
    wrong_repeat_unit = _truth(
        2430,
        5490,
        unit="bp",
        boundary_definition="bp-boundary-v1",
        conversion_id="bp-source-to-target-v1",
        repeat_unit_bp=61,
    )

    with pytest.raises(ValueError, match="target boundary"):
        converted_length_target(truth, registry, target_boundary_definition="another-target")
    with pytest.raises(ValueError, match="repeat unit"):
        converted_length_target(wrong_repeat_unit, registry, target_boundary_definition="target-boundary-v1")


def test_conversion_requires_known_id_and_exact_source_bindings() -> None:
    registry = decode_conversion_registry(_registry())
    unknown = _truth(
        2430,
        5490,
        unit="bp",
        boundary_definition="bp-boundary-v1",
        conversion_id="unknown-conversion",
    )
    wrong_source_unit = _truth(
        40,
        90,
        boundary_definition="bp-boundary-v1",
        conversion_id="bp-source-to-target-v1",
    )
    wrong_source_boundary = _truth(
        2430,
        5490,
        unit="bp",
        boundary_definition="other-bp-boundary",
        conversion_id="bp-source-to-target-v1",
    )

    for truth, message in (
        (unknown, "unknown.*conversion id"),
        (wrong_source_unit, "source unit"),
        (wrong_source_boundary, "source boundary"),
    ):
        with pytest.raises(ValueError, match=message):
            converted_length_target(truth, registry, target_boundary_definition="target-boundary-v1")


def test_conversion_passes_through_matching_truth_and_excludes_incomplete_pairs() -> None:
    registry = decode_conversion_registry(_registry())
    direct = _truth(40, 90)
    incomplete = _truth(
        2430,
        None,
        unit="bp",
        boundary_definition="bp-boundary-v1",
        conversion_id="bp-source-to-target-v1",
    )

    assert converted_length_target(direct, registry, target_boundary_definition="target-boundary-v1") == 130.0
    assert converted_length_target(incomplete, registry, target_boundary_definition="target-boundary-v1") is None


def test_conversion_rejects_an_unregistered_boundary_and_invalid_api_values() -> None:
    registry = decode_conversion_registry(_registry())
    mismatched = _truth(40, 90, boundary_definition="legacy-count-boundary-v1")

    with pytest.raises(ValueError, match="boundary mismatch"):
        converted_length_target(mismatched, registry, target_boundary_definition="target-boundary-v1")
    with pytest.raises(ValueError, match="ConversionRegistry"):
        converted_length_target(mismatched, object(), target_boundary_definition="target-boundary-v1")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="target boundary"):
        converted_length_target(mismatched, registry, target_boundary_definition="")


def test_exact_repeat_count_truth_must_be_integral_and_bp_truth_requires_conversion() -> None:
    with pytest.raises(ValueError, match="integral"):
        _truth(40.5, 90)
    with pytest.raises(ValueError, match="conversion id"):
        _truth(2430, 5490, unit="bp", boundary_definition="bp-boundary-v1")


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"unit": "copies"}, "length unit"),
        ({"measurement": "approximate"}, "length measurement"),
        ({"allele_1": 0}, "allele 1"),
        ({"allele_2": float("nan")}, "allele 2"),
        ({"measurement": "exact", "lower_bound": 1}, "forbids interval bounds"),
        ({"measurement": "interval", "allele_1": 40, "lower_bound": 1}, "forbids exact allele"),
        ({"measurement": "interval", "allele_1": None, "allele_2": None}, "requires a directional bound"),
        (
            {"measurement": "interval", "allele_1": None, "allele_2": None, "lower_bound": 200, "upper_bound": 100},
            "must not exceed",
        ),
        ({"repeat_unit_bp": True}, "repeat unit"),
    ],
)
def test_length_truth_rejects_invalid_numeric_and_measurement_combinations(
    updates: dict[str, object], message: str
) -> None:
    raw = synthetic_intake()
    raw_length: dict[str, object] = {
        "allele_1": 40,
        "allele_2": 90,
        "unit": "repeat-count",
        "repeat_unit_bp": 60,
        "boundary_definition": "target-boundary-v1",
        "measurement": "exact",
        "lower_bound": None,
        "upper_bound": None,
        "conversion_id": None,
    }
    raw_length.update(updates)
    raw["truth"][0]["length"] = raw_length

    with pytest.raises(ValueError, match=message):
        decode_intake(raw)
