"""Strict synthetic contracts for VNTR length annotations."""

from __future__ import annotations

import pytest

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_annotation import (
    Interval,
    decode_length_annotation,
    encode_length_annotation,
    one_based_closed_to_zero_based_half_open,
)

pytestmark = pytest.mark.unit


def _annotation(**changes: object) -> dict[str, object]:
    value: dict[str, object] = {
        "schema_version": "length-annotation-v1",
        "assembly": "synthetic-build-v1",
        "contig": "synthetic-contig",
        "accepted_contig_aliases": ["alias-contig"],
        "reference_fasta_sha256": "a" * 64,
        "coordinate_system": "zero-based-half-open",
        "boundary_definition": "complete-core-plus-invariant-units-v1",
        "repeat_unit_bp": 2,
        "regions": {
            "CORE": [{"start": 10, "end": 14}],
            "INVARIANT": [{"start": 16, "end": 18}],
            "ARRAY": {"start": 10, "end": 19},
            "LEFT_FLANK": {"start": 7, "end": 10},
            "RIGHT_FLANK": {"start": 19, "end": 22},
        },
        "array_boundary_geometry": {"array_only_bp": 3, "target_only_bp": 0},
        "target_boundary_conversion_sha256": None,
        "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
        "annotation_provenance": "generated-synthetic-test-contract",
        "annotation_version": "synthetic-1",
    }
    value.update(changes)
    return value


def test_one_based_closed_coordinates_convert_without_an_off_by_one_error() -> None:
    assert one_based_closed_to_zero_based_half_open(11, 14) == Interval(start=10, end=14)


@pytest.mark.parametrize(("start", "end"), [(0, 1), (2, 1), (True, 2), (1, False)])
def test_coordinate_conversion_rejects_invalid_one_based_boundaries(start: object, end: object) -> None:
    with pytest.raises(ValueError, match="one-based closed"):
        one_based_closed_to_zero_based_half_open(start, end)


def test_decode_preserves_independent_array_geometry_and_has_a_deterministic_digest() -> None:
    raw = _annotation()

    annotation = decode_length_annotation(raw)

    assert annotation.core == (Interval(10, 14),)
    assert annotation.invariant == (Interval(16, 18),)
    assert annotation.array == Interval(10, 19)
    assert annotation.array.length_bp == 9
    assert annotation.array_only_bp == 3
    assert annotation.target_only_bp == 0
    assert annotation.array.length_bp % annotation.repeat_unit_bp == 1
    assert annotation.reference_core_repeat_count == 2
    assert annotation.reference_invariant_repeat_count == 1
    assert annotation.physical_a_compatible is False
    assert annotation.physical_f_compatible is False
    assert encode_length_annotation(annotation) == raw
    assert annotation.sha256 == canonical_sha256(raw)


def test_core_and_invariant_repeat_units_must_be_disjoint_and_invariant_fixed_width() -> None:
    for core, invariant, message in (
        ([{"start": 10, "end": 14}], [{"start": 16, "end": 19}], "complete repeat units"),
        ([{"start": 10, "end": 14}], [{"start": 13, "end": 17}], "disjoint"),
        (
            [{"start": 10, "end": 14}, {"start": 13, "end": 15}],
            [{"start": 16, "end": 18}],
            "disjoint",
        ),
    ):
        raw = _annotation()
        assert isinstance(raw["regions"], dict)
        regions = dict(raw["regions"])
        regions["CORE"] = core
        regions["INVARIANT"] = invariant
        raw["regions"] = regions
        with pytest.raises(ValueError, match=message):
            decode_length_annotation(raw)


def test_complete_core_region_with_indel_bases_does_not_invent_an_integral_reference_count() -> None:
    raw = _annotation()
    assert isinstance(raw["regions"], dict)
    raw["regions"] = {
        **raw["regions"],
        "CORE": [{"start": 10, "end": 15}],
        "ARRAY": {"start": 10, "end": 20},
        "RIGHT_FLANK": {"start": 20, "end": 23},
    }
    raw["array_boundary_geometry"] = {"array_only_bp": 3, "target_only_bp": 0}

    annotation = decode_length_annotation(raw)

    assert annotation.core == (Interval(10, 15),)
    assert annotation.reference_core_repeat_count is None


def test_physical_a_rejects_a_core_without_fixed_width_reference_count() -> None:
    raw = _annotation()
    assert isinstance(raw["regions"], dict)
    raw["regions"] = {
        **raw["regions"],
        "CORE": [{"start": 10, "end": 15}],
        "ARRAY": {"start": 10, "end": 20},
        "RIGHT_FLANK": {"start": 20, "end": 23},
    }
    raw["array_boundary_geometry"] = {"array_only_bp": 3, "target_only_bp": 0}
    raw["target_boundary_conversion_sha256"] = "b" * 64
    raw["physical_hypothesis_compatibility"] = {"physical_A": True, "physical_F": False}

    with pytest.raises(ValueError, match="fixed-width CORE reference count"):
        decode_length_annotation(raw)


@pytest.mark.parametrize(
    ("region", "interval"),
    [
        ("LEFT_FLANK", {"start": 9, "end": 11}),
        ("RIGHT_FLANK", {"start": 18, "end": 21}),
        ("RIGHT_FLANK", {"start": 8, "end": 10}),
    ],
)
def test_flanks_are_unique_nonoverlapping_bases_on_their_declared_side(region: str, interval: dict[str, int]) -> None:
    raw = _annotation()
    assert isinstance(raw["regions"], dict)
    regions = dict(raw["regions"])
    regions[region] = interval
    raw["regions"] = regions

    with pytest.raises(ValueError, match="flank"):
        decode_length_annotation(raw)


def test_boundary_geometry_is_checked_against_the_independent_array_and_target() -> None:
    raw = _annotation()
    raw["array_boundary_geometry"] = {"array_only_bp": 2, "target_only_bp": 0}

    with pytest.raises(ValueError, match="boundary geometry"):
        decode_length_annotation(raw)


def test_boundary_geometry_records_target_bases_outside_the_array() -> None:
    raw = _annotation()
    assert isinstance(raw["regions"], dict)
    raw["regions"] = {
        **raw["regions"],
        "INVARIANT": [{"start": 20, "end": 22}],
        "RIGHT_FLANK": {"start": 22, "end": 25},
    }
    raw["array_boundary_geometry"] = {"array_only_bp": 5, "target_only_bp": 2}

    annotation = decode_length_annotation(raw)

    assert (annotation.array_only_bp, annotation.target_only_bp) == (5, 2)


def test_missing_core_and_invariant_are_explicit_without_fabricating_a_geometry() -> None:
    raw = _annotation()
    assert isinstance(raw["regions"], dict)
    regions = dict(raw["regions"])
    regions["CORE"] = None
    regions["INVARIANT"] = None
    raw["regions"] = regions
    raw["array_boundary_geometry"] = {"array_only_bp": None, "target_only_bp": None}

    annotation = decode_length_annotation(raw)

    assert annotation.core is None
    assert annotation.invariant is None
    assert annotation.reference_core_repeat_count is None
    assert annotation.reference_invariant_repeat_count is None
    assert annotation.array == Interval(10, 19)


def test_physical_hypotheses_require_an_explicit_target_boundary_conversion() -> None:
    raw = _annotation()
    raw["physical_hypothesis_compatibility"] = {"physical_A": True, "physical_F": False}

    with pytest.raises(ValueError, match="conversion"):
        decode_length_annotation(raw)

    raw["target_boundary_conversion_sha256"] = "b" * 64
    annotation = decode_length_annotation(raw)
    assert annotation.physical_a_compatible is True


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"schema_version": "length-annotation-v2"}, "schema version"),
        ({"coordinate_system": "one-based-closed"}, "coordinate system"),
        ({"reference_fasta_sha256": "not-a-digest"}, "digest"),
        ({"repeat_unit_bp": True}, "repeat unit"),
        ({"repeat_unit_bp": 0}, "repeat unit"),
        ({"accepted_contig_aliases": ["alias-contig", "alias-contig"]}, "aliases"),
        ({"accepted_contig_aliases": ["synthetic-contig"]}, "aliases"),
    ],
)
def test_annotation_rejects_wrong_versions_types_digests_and_aliases(change: dict[str, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        decode_length_annotation(_annotation(**change))


def test_motif_names_cannot_stand_in_for_reviewed_repeat_class_boundaries() -> None:
    raw = _annotation()
    raw["motif_name"] = "ambiguous-end-repeat-name"

    with pytest.raises(ValueError, match="fields differ"):
        decode_length_annotation(raw)


def test_annotation_values_are_deeply_immutable() -> None:
    annotation = decode_length_annotation(_annotation())

    with pytest.raises(AttributeError):
        annotation.assembly = "changed"  # type: ignore[misc]
    with pytest.raises(TypeError):
        annotation.accepted_contig_aliases[0] = "changed"  # type: ignore[index]
