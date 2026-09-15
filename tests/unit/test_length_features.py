"""Pure synthetic arithmetic and failure contracts for length features."""

from __future__ import annotations

import math

import pytest

from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation
from vntyper.scripts.length_features import (
    DepthPosition,
    LengthFeatureProvenance,
    encode_length_features,
    extract_length_features,
)

pytestmark = pytest.mark.unit


def _annotation_raw(**changes: object) -> dict[str, object]:
    value: dict[str, object] = {
        "schema_version": "length-annotation-v1",
        "assembly": "synthetic-build-v1",
        "contig": "synthetic-contig",
        "accepted_contig_aliases": ["alias-contig"],
        "reference_fasta_sha256": "a" * 64,
        "coordinate_system": "zero-based-half-open",
        "boundary_definition": "complete-core-plus-invariant-units-v1",
        "repeat_unit_bp": 1,
        "regions": {
            "CORE": [{"start": 1, "end": 3}],
            "INVARIANT": [{"start": 3, "end": 5}],
            "ARRAY": {"start": 1, "end": 5},
            "LEFT_FLANK": {"start": 0, "end": 1},
            "RIGHT_FLANK": {"start": 5, "end": 6},
        },
        "array_boundary_geometry": {"array_only_bp": 0, "target_only_bp": 0},
        "target_boundary_conversion_sha256": None,
        "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
        "annotation_provenance": "generated-synthetic-test-contract",
        "annotation_version": "synthetic-1",
    }
    value.update(changes)
    return value


def _annotation(**changes: object) -> LengthAnnotation:
    return decode_length_annotation(_annotation_raw(**changes))


def _provenance(annotation: LengthAnnotation, **changes: object) -> LengthFeatureProvenance:
    values: dict[str, object] = {
        "manifest_key": "synthetic-member-1",
        "assembly": annotation.assembly,
        "assay_class": "synthetic-short-read",
        "input_scope": "regional",
        "contig": annotation.contig,
        "reference_fasta_sha256": annotation.reference_fasta_sha256,
        "annotation_sha256": annotation.sha256,
        "counting_policy_sha256": "b" * 64,
        "provenance_sha256": "c" * 64,
        "aligner_name": "synthetic-aligner",
        "aligner_version": "1.0",
        "aligner_arguments_sha256": "d" * 64,
        "primary_secondary_marking": "primary-only",
        "preprocessing_id": "synthetic-preprocessing-v1",
    }
    values.update(changes)
    return LengthFeatureProvenance(**values)  # type: ignore[arg-type]


def _depths(
    values: tuple[float, ...] = (5, 0, 20, 5, 5, 5),
    *,
    contig: str = "synthetic-contig",
    fragments: bool = False,
) -> tuple[DepthPosition, ...]:
    return tuple(
        DepthPosition(
            contig=contig,
            position_zero_based=position,
            depth=depth,
            supporting_fragment_ids=(f"fragment-{position}",) if fragments else None,
        )
        for position, depth in enumerate(values)
    )


def test_zero_depths_contribute_to_lengths_means_and_covered_fractions() -> None:
    annotation = _annotation()

    features = extract_length_features(_depths(), annotation, _provenance(annotation))

    core = features.regions["CORE"]
    invariant = features.regions["INVARIANT"]
    array = features.regions["ARRAY"]
    assert core is not None and invariant is not None and array is not None
    assert (core.length_bp, core.depth_sum, core.mean_depth, core.covered_fraction) == (2, 20.0, 10.0, 0.5)
    assert (invariant.length_bp, invariant.depth_sum, invariant.mean_depth, invariant.covered_fraction) == (
        2,
        10.0,
        5.0,
        1.0,
    )
    assert (array.length_bp, array.depth_sum, array.mean_depth, array.covered_fraction) == (4, 30.0, 7.5, 0.75)
    assert features.a == 2.0
    assert features.f == 1.5
    assert features.status == "measured"
    assert features.reasons == ()
    assert all(region is not None and region.supporting_fragment_count is None for region in features.regions.values())


def test_multiplying_every_depth_preserves_both_ratios() -> None:
    annotation = _annotation()
    provenance = _provenance(annotation)

    baseline = extract_length_features(_depths(), annotation, provenance)
    scaled = extract_length_features(_depths(tuple(value * 3 for value in (5, 0, 20, 5, 5, 5))), annotation, provenance)

    assert (baseline.a, baseline.f) == (2.0, 1.5)
    assert (scaled.a, scaled.f) == (2.0, 1.5)


def test_an_expected_zero_covered_position_changes_the_region_denominator() -> None:
    raw = _annotation_raw()
    raw["regions"] = {
        "CORE": [{"start": 1, "end": 4}],
        "INVARIANT": [{"start": 4, "end": 6}],
        "ARRAY": {"start": 1, "end": 6},
        "LEFT_FLANK": {"start": 0, "end": 1},
        "RIGHT_FLANK": {"start": 6, "end": 7},
    }
    annotation = decode_length_annotation(raw)
    depths = _depths((5, 0, 20, 0, 5, 5, 5))

    features = extract_length_features(depths, annotation, _provenance(annotation))

    core = features.regions["CORE"]
    assert core is not None
    assert (core.length_bp, core.depth_sum, core.mean_depth, core.covered_fraction) == (3, 20.0, 20.0 / 3.0, 1.0 / 3.0)


def test_exact_region_union_does_not_invent_positions_in_annotation_gaps() -> None:
    raw = _annotation_raw()
    raw["regions"] = {
        "CORE": [{"start": 1, "end": 2}],
        "INVARIANT": [{"start": 4, "end": 5}],
        "ARRAY": {"start": 1, "end": 3},
        "LEFT_FLANK": {"start": 0, "end": 1},
        "RIGHT_FLANK": {"start": 5, "end": 6},
    }
    raw["array_boundary_geometry"] = {"array_only_bp": 1, "target_only_bp": 1}
    annotation = decode_length_annotation(raw)
    depths = tuple(DepthPosition("synthetic-contig", position, 5) for position in (0, 1, 2, 4, 5))

    features = extract_length_features(depths, annotation, _provenance(annotation))

    assert features.regions["ARRAY"] is not None
    assert features.regions["ARRAY"].length_bp == 2


def test_depth_input_order_does_not_change_serialized_features_or_digest() -> None:
    annotation = _annotation()
    provenance = _provenance(annotation)
    ordered = extract_length_features(_depths(), annotation, provenance)

    reversed_input = extract_length_features(tuple(reversed(_depths())), annotation, provenance)

    assert encode_length_features(reversed_input) == encode_length_features(ordered)
    assert reversed_input.sha256 == ordered.sha256


def test_contig_alias_is_accepted_but_mixed_alias_spellings_are_rejected() -> None:
    annotation = _annotation()
    aliased = extract_length_features(
        _depths(contig="alias-contig"),
        annotation,
        _provenance(annotation, contig="alias-contig"),
    )
    assert aliased.a == 2.0

    mixed = list(_depths())
    mixed[-1] = DepthPosition("alias-contig", 5, 5)
    with pytest.raises(ValueError, match="mixed contigs"):
        extract_length_features(tuple(mixed), annotation, _provenance(annotation))


@pytest.mark.parametrize(
    ("depths", "message"),
    [
        ((*_depths()[:-1],), "missing"),
        (((*_depths(), DepthPosition("synthetic-contig", 5, 5))), "duplicate"),
        (((*_depths(), DepthPosition("synthetic-contig", 6, 5))), "out-of-range"),
    ],
)
def test_every_position_in_the_exact_declared_region_union_occurs_once(
    depths: tuple[DepthPosition, ...], message: str
) -> None:
    annotation = _annotation()
    with pytest.raises(ValueError, match=message):
        extract_length_features(depths, annotation, _provenance(annotation))


@pytest.mark.parametrize("depth", [-1, math.inf, math.nan, True, "5"])
def test_depth_position_rejects_negative_nonfinite_boolean_and_nonnumeric_depth(depth: object) -> None:
    with pytest.raises(ValueError, match="depth"):
        DepthPosition("synthetic-contig", 0, depth)  # type: ignore[arg-type]


@pytest.mark.parametrize("position", [-1, 0.5, True])
def test_depth_position_requires_an_explicit_zero_based_integer_position(position: object) -> None:
    with pytest.raises(ValueError, match="zero-based"):
        DepthPosition("synthetic-contig", position, 1)  # type: ignore[arg-type]


def test_depth_only_input_never_invents_fragment_support_counts() -> None:
    annotation = _annotation()
    features = extract_length_features(_depths(), annotation, _provenance(annotation))

    assert {name: region.supporting_fragment_count for name, region in features.regions.items() if region} == {
        "CORE": None,
        "INVARIANT": None,
        "ARRAY": None,
        "LEFT_FLANK": None,
        "RIGHT_FLANK": None,
    }


def test_fragment_support_is_deduplicated_per_region_without_changing_base_depth() -> None:
    annotation = _annotation()
    depths = tuple(
        DepthPosition(
            "synthetic-contig",
            position,
            depth,
            ("shared-pair", f"position-{position}"),
        )
        for position, depth in enumerate((5, 0, 20, 5, 5, 5))
    )

    features = extract_length_features(depths, annotation, _provenance(annotation))

    core = features.regions["CORE"]
    array = features.regions["ARRAY"]
    assert core is not None and array is not None
    assert (core.depth_sum, core.supporting_fragment_count) == (20.0, 3)
    assert (array.depth_sum, array.supporting_fragment_count) == (30.0, 5)


def test_partial_fragment_metadata_is_rejected_instead_of_undercounting_support() -> None:
    annotation = _annotation()
    depths = list(_depths(fragments=True))
    depths[0] = DepthPosition("synthetic-contig", 0, 5, None)

    with pytest.raises(ValueError, match="fragment metadata"):
        extract_length_features(tuple(depths), annotation, _provenance(annotation))


def test_duplicate_fragment_identity_at_one_position_is_rejected() -> None:
    with pytest.raises(ValueError, match="fragment IDs"):
        DepthPosition("synthetic-contig", 0, 5, ("same-pair", "same-pair"))


def test_partial_annotation_keeps_f_when_core_and_invariant_are_unavailable() -> None:
    raw = _annotation_raw()
    assert isinstance(raw["regions"], dict)
    regions = dict(raw["regions"])
    regions["CORE"] = None
    regions["INVARIANT"] = None
    raw["regions"] = regions
    raw["array_boundary_geometry"] = {"array_only_bp": None, "target_only_bp": None}
    annotation = decode_length_annotation(raw)

    features = extract_length_features(_depths(), annotation, _provenance(annotation))

    assert features.regions["CORE"] is None
    assert features.regions["INVARIANT"] is None
    assert features.a is None
    assert features.f == 1.5
    assert features.status == "partial"
    assert features.reasons == ("missing_core_annotation", "missing_invariant_annotation")


def test_zero_denominators_produce_null_ratios_and_stable_reasons() -> None:
    annotation = _annotation()

    features = extract_length_features(_depths((0, 0, 20, 0, 0, 0)), annotation, _provenance(annotation))

    assert features.a is None
    assert features.f is None
    assert features.status == "partial"
    assert features.reasons == ("zero_invariant_mean_depth", "zero_combined_flank_mean_depth")


def test_completely_unavailable_annotation_emits_an_unavailable_empty_row() -> None:
    raw = _annotation_raw()
    raw["regions"] = dict.fromkeys(("CORE", "INVARIANT", "ARRAY", "LEFT_FLANK", "RIGHT_FLANK"))
    raw["array_boundary_geometry"] = {"array_only_bp": None, "target_only_bp": None}
    annotation = decode_length_annotation(raw)

    features = extract_length_features((), annotation, _provenance(annotation))

    assert all(region is None for region in features.regions.values())
    assert (features.a, features.f, features.status) == (None, None, "unavailable")


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"assembly": "wrong-build"}, "assembly mismatch"),
        ({"contig": "wrong-contig"}, "contig mismatch"),
        ({"reference_fasta_sha256": "d" * 64}, "reference digest mismatch"),
        ({"annotation_sha256": "d" * 64}, "annotation digest mismatch"),
    ],
)
def test_provenance_must_match_annotation_build_contig_and_digests(change: dict[str, object], message: str) -> None:
    annotation = _annotation()
    with pytest.raises(ValueError, match=message):
        extract_length_features(_depths(), annotation, _provenance(annotation, **change))


@pytest.mark.parametrize(
    "change",
    [
        {"annotation_sha256": True},
        {"counting_policy_sha256": "short"},
        {"provenance_sha256": "C" * 64},
        {"manifest_key": ""},
    ],
)
def test_provenance_rejects_boolean_or_malformed_hashes_and_empty_identifiers(change: dict[str, object]) -> None:
    annotation = _annotation()
    with pytest.raises(ValueError):
        _provenance(annotation, **change)


def test_provenance_explicitly_types_alignment_marking_preprocessing_and_scope() -> None:
    annotation = _annotation()
    provenance = _provenance(annotation)

    assert (
        provenance.aligner_name,
        provenance.aligner_version,
        provenance.aligner_arguments_sha256,
        provenance.primary_secondary_marking,
        provenance.preprocessing_id,
        provenance.input_scope,
    ) == (
        "synthetic-aligner",
        "1.0",
        "d" * 64,
        "primary-only",
        "synthetic-preprocessing-v1",
        "regional",
    )

    with pytest.raises(ValueError, match="input scope"):
        _provenance(annotation, input_scope="unknown")


def test_feature_encoding_has_the_exact_closed_schema_and_no_fragment_identifiers() -> None:
    annotation = _annotation()
    features = extract_length_features(_depths(fragments=True), annotation, _provenance(annotation))

    encoded = encode_length_features(features)

    assert set(encoded) == {"schema_version", "rows"}
    assert encoded["schema_version"] == "length-features-v1"
    rows = encoded["rows"]
    assert isinstance(rows, list) and len(rows) == 1
    row = rows[0]
    assert set(row) == {
        "manifest_key",
        "assembly",
        "assay_class",
        "input_scope",
        "annotation_sha256",
        "counting_policy_sha256",
        "provenance_sha256",
        "regions",
        "A",
        "F",
        "status",
        "reasons",
    }
    assert set(row["regions"]) == {"CORE", "INVARIANT", "ARRAY", "LEFT_FLANK", "RIGHT_FLANK"}
    assert set(row["regions"]["CORE"]) == {
        "length_bp",
        "depth_sum",
        "mean_depth",
        "covered_fraction",
        "supporting_fragment_count",
    }
    assert "fragment-1" not in repr(encoded)
    assert features.sha256 == canonical_sha256(encoded)
    assert encode_length_features(features) == encoded


def test_feature_values_and_nested_region_mapping_are_immutable() -> None:
    annotation = _annotation()
    features = extract_length_features(_depths(), annotation, _provenance(annotation))

    with pytest.raises(AttributeError):
        features.a = 10.0  # type: ignore[misc]
    with pytest.raises(TypeError):
        features.regions["CORE"] = None  # type: ignore[index]
