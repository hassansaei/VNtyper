"""Contracts for complete, canonical VNtyper decision profiles."""

from __future__ import annotations

import copy
import hashlib
import re
from pathlib import Path

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile_schema import (
    ValidationClass,
    component_projection,
    flatten_decision_projection,
    validate_complete_inventory,
)

pytestmark = pytest.mark.unit

PACKAGE_ROOT = Path(__file__).parents[2] / "vntyper"
PROFILE_PATH = PACKAGE_ROOT / "profiles" / "decision_profile.json"
DIGEST_PATH = PACKAGE_ROOT / "profiles" / "decision_profile.sha256"
PROJECTION_PATH = PACKAGE_ROOT / "profiles" / "decision_projection.json"


def _packaged_profile() -> dict[str, object]:
    return load_strict_json_object(PROFILE_PATH.read_bytes())


def _field(profile: dict[str, object], pointer: str) -> dict[str, object]:
    inventory = profile["inventory"]
    assert isinstance(inventory, dict)
    value = inventory[pointer]
    assert isinstance(value, dict)
    return value


def test_packaged_profile_is_canonical_and_matches_sidecar() -> None:
    profile = _packaged_profile()
    raw = PROFILE_PATH.read_bytes()

    assert raw == canonical_json_bytes(profile)
    assert raw.endswith(b"\n") and not raw.endswith(b"\n\n")
    digest_bytes = DIGEST_PATH.read_bytes()
    assert re.fullmatch(rb"[0-9a-f]{64}\n", digest_bytes)
    digest = digest_bytes[:-1].decode("ascii")
    assert digest == hashlib.sha256(raw).hexdigest() == canonical_sha256(profile)
    validate_complete_inventory(profile)


def test_inventory_is_an_exact_flattening_of_the_checked_in_projection() -> None:
    profile = _packaged_profile()
    projection = load_strict_json_object(PROJECTION_PATH.read_bytes())
    flattened = flatten_decision_projection(projection)
    inventory = profile["inventory"]
    assert isinstance(inventory, dict)

    assert set(inventory) == set(flattened)
    assert {pointer: field["value"] for pointer, field in inventory.items()} == flattened
    assert set(projection) == {"advntr", "cross_match", "dominance", "kestrel", "nomenclature", "shark"}
    assert component_projection(profile, "shark") == {}
    for component in projection:
        assert component_projection(profile, component) == projection[component]


def test_every_inventory_leaf_has_exactly_one_validation_class() -> None:
    profile = _packaged_profile()
    inventory = profile["inventory"]
    assert isinstance(inventory, dict)

    observed = {field["class"] for field in inventory.values()}
    assert observed == {member.value for member in ValidationClass}
    assert all(isinstance(field["class"], str) for field in inventory.values())


@pytest.mark.parametrize(
    ("pointer", "value", "unit", "comparator", "inclusive"),
    [
        (
            "/components/kestrel/confidence_assignment/depth_score_thresholds/low",
            0.00469,
            "depth-score-ratio",
            "gte",
            True,
        ),
        ("/components/kestrel/alt_filtering/gg_depth_score_threshold", 0.00469, "depth-score-ratio", "gte", True),
        (
            "/components/kestrel/confidence_assignment/depth_score_thresholds/high",
            0.00515,
            "depth-score-ratio",
            "gte",
            True,
        ),
        (
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/low",
            20,
            "alternate-kmer-path-depth",
            "lte",
            True,
        ),
        (
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low",
            21,
            "alternate-kmer-path-depth",
            "gte",
            True,
        ),
        (
            "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high",
            100,
            "alternate-kmer-path-depth",
            "gte-and-upper-exclusive",
            False,
        ),
        (
            "/components/kestrel/confidence_assignment/var_active_region_threshold",
            200,
            "active-region-kmer-depth",
            "lte",
            True,
        ),
        ("/components/nomenclature/thresholds/bam_flank", 8, "base-pairs-per-side", "eq", True),
        (
            "/components/nomenclature/thresholds/bam_thin_haplotype_record_support",
            3,
            "resolved-haplotype-records",
            "lt",
            False,
        ),
        (
            "/components/nomenclature/identity_reconciliation/kestrel_min_alternate_kmer_path_depth",
            5,
            "alternate-kmer-path-depth",
            "gte",
            True,
        ),
        (
            "/components/nomenclature/identity_reconciliation/advntr_min_sequencing_read_support",
            5,
            "sequencing-reads",
            "gte",
            True,
        ),
    ],
)
def test_fixed_safety_boundaries_retain_value_unit_and_semantics(
    pointer: str, value: object, unit: str, comparator: str, inclusive: bool
) -> None:
    field = _field(_packaged_profile(), pointer)

    assert field == {
        "class": ValidationClass.FIXED_SAFETY.value,
        "value": value,
        "unit": unit,
        "comparator": comparator,
        "inclusive": inclusive,
    }


def test_low_reporting_floor_and_independent_gg_gate_are_linked() -> None:
    profile = _packaged_profile()
    low = _field(profile, "/components/kestrel/confidence_assignment/depth_score_thresholds/low")
    gg = _field(profile, "/components/kestrel/alt_filtering/gg_depth_score_threshold")

    assert low["value"] == gg["value"] == 0.00469
    assert low is not gg


def test_every_numeric_leaf_declares_unit_comparator_and_inclusivity() -> None:
    inventory = _packaged_profile()["inventory"]
    assert isinstance(inventory, dict)

    for pointer, field in inventory.items():
        value = field["value"]
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            assert set(field) == {"class", "value", "unit", "comparator", "inclusive"}, pointer


def test_rfc8785_hash_ignores_object_order_but_preserves_array_order() -> None:
    profile = _packaged_profile()
    reversed_root = dict(reversed(list(profile.items())))
    assert canonical_sha256(reversed_root) == canonical_sha256(profile)

    changed = copy.deepcopy(profile)
    inventory = changed["inventory"]
    assert isinstance(inventory, dict)
    ordered = sorted(pointer for pointer in inventory if "/selection/sort_order/" in pointer)
    first = inventory[ordered[0]]["value"]
    inventory[ordered[0]]["value"] = inventory[ordered[-1]]["value"]
    inventory[ordered[-1]]["value"] = first
    assert canonical_sha256(changed) != canonical_sha256(profile)


def test_issue_268_model_guard_is_not_a_profile_field() -> None:
    packaged = _packaged_profile()
    custom = copy.deepcopy(packaged)
    custom["profile_kind"] = "explicit-custom"
    custom["profile_id"] = "unit-test-issue-268"
    inventory = custom["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/advntr/model/required_version"] = {
        "class": "fixed-safety",
        "value": "2.0.4",
    }

    with pytest.raises(ValueError, match="inventory fields differ"):
        validate_complete_inventory(custom, packaged_profile=packaged)


def test_fixed_safety_and_generated_explicit_fields_cannot_drift() -> None:
    packaged = _packaged_profile()
    explicit = copy.deepcopy(packaged)
    explicit["profile_kind"] = "explicit-custom"
    explicit["profile_id"] = "unit-test-explicit"
    _field(explicit, "/components/kestrel/confidence_assignment/alt_depth_thresholds/low")["value"] = 19
    with pytest.raises(ValueError, match="fixed-safety field"):
        validate_complete_inventory(explicit, packaged_profile=packaged)

    generated = copy.deepcopy(packaged)
    generated["profile_kind"] = "generated"
    generated["profile_id"] = "unit-test-generated"
    generated["generated_metadata"] = {
        "packaged_base_hash": canonical_sha256(packaged),
        "generator_name": "unit-test",
        "generator_version": "1",
        "objective": "lexicographic-safety-v1",
        "dataset_manifest_hash": "a" * 64,
        "partition_manifest_hash": "b" * 64,
        "seed": 295,
    }
    generated_inventory = generated["inventory"]
    assert isinstance(generated_inventory, dict)
    explicit_pointer = next(
        pointer
        for pointer, field in generated_inventory.items()
        if field["class"] == ValidationClass.EXPLICIT_CUSTOM.value and isinstance(field["value"], str)
    )
    generated_inventory[explicit_pointer]["value"] += "-changed"
    with pytest.raises(ValueError, match="generated profile must copy explicit-custom field"):
        validate_complete_inventory(generated, packaged_profile=packaged)


def _different_json_value(value: object) -> object:
    if isinstance(value, bool):
        return not value
    if isinstance(value, (int, float)):
        return value + 1
    if isinstance(value, str):
        return f"{value}-mutated"
    raise AssertionError(f"fixed-safety mutation fixture has unsupported value: {value!r}")


def _fixed_fields() -> dict[str, dict[str, object]]:
    inventory = _packaged_profile()["inventory"]
    assert isinstance(inventory, dict)
    return {
        pointer: field
        for pointer, field in inventory.items()
        if isinstance(field, dict) and field["class"] == ValidationClass.FIXED_SAFETY.value
    }


_FIXED_FIELDS = _fixed_fields()
_FIXED_NUMERIC_SEMANTICS = [
    (pointer, key)
    for pointer, field in _FIXED_FIELDS.items()
    if isinstance(field["value"], (int, float)) and not isinstance(field["value"], bool)
    for key in ("unit", "comparator", "inclusive")
]


@pytest.mark.parametrize("pointer", sorted(_FIXED_FIELDS))
def test_every_fixed_safety_value_has_a_failing_mutation(pointer: str) -> None:
    packaged = _packaged_profile()
    explicit = copy.deepcopy(packaged)
    explicit["profile_kind"] = "explicit-custom"
    explicit["profile_id"] = "unit-test-fixed-value-mutation"
    field = _field(explicit, pointer)
    field["value"] = _different_json_value(field["value"])

    with pytest.raises(ValueError, match="immutable fixed-safety field differs|critical fixed-safety field differs"):
        validate_complete_inventory(explicit, packaged_profile=packaged)


@pytest.mark.parametrize(("pointer", "key"), _FIXED_NUMERIC_SEMANTICS)
def test_every_fixed_numeric_unit_comparator_and_inclusivity_has_a_failing_mutation(pointer: str, key: str) -> None:
    packaged = _packaged_profile()
    explicit = copy.deepcopy(packaged)
    explicit["profile_kind"] = "explicit-custom"
    explicit["profile_id"] = "unit-test-fixed-semantics-mutation"
    field = _field(explicit, pointer)
    field[key] = not field[key] if key == "inclusive" else f"{field[key]}-mutated"

    with pytest.raises(ValueError, match="semantics differ|critical fixed-safety field differs"):
        validate_complete_inventory(explicit, packaged_profile=packaged)


@pytest.mark.parametrize(
    "pointer",
    [
        "/components/advntr/model/path",
        "/components/advntr/model/sha256",
        "/components/advntr/model/schema_version",
        "/components/advntr/model/vid",
        "/components/advntr/model/genomic_interval",
        "/components/advntr/model/window_bp",
        "/components/advntr/model/n_segments",
        "/components/advntr/model/n_distinct_segments",
        "/components/advntr/model/max_segment_len",
        "/components/advntr/model/minimum_compatible_advntr_version",
        "/components/advntr/model/installed_advntr_version",
    ],
)
def test_every_issue_268_model_and_version_guard_is_rejected_as_an_unknown_key(pointer: str) -> None:
    packaged = _packaged_profile()
    explicit = copy.deepcopy(packaged)
    explicit["profile_kind"] = "explicit-custom"
    explicit["profile_id"] = "unit-test-excluded-issue-268-guard"
    inventory = explicit["inventory"]
    assert isinstance(inventory, dict)
    inventory[pointer] = {"class": "fixed-safety", "value": "forbidden"}

    with pytest.raises(ValueError, match="inventory fields differ"):
        validate_complete_inventory(explicit, packaged_profile=packaged)


def test_generated_metadata_is_closed_and_mutable_grid_is_bounded() -> None:
    packaged = _packaged_profile()
    generated = copy.deepcopy(packaged)
    generated["profile_kind"] = "generated"
    generated["profile_id"] = "unit-test-generated"
    generated["generated_metadata"] = {
        "packaged_base_hash": canonical_sha256(packaged),
        "generator_name": "unit-test",
        "generator_version": "1",
        "objective": "lexicographic-safety-v1",
        "dataset_manifest_hash": "a" * 64,
        "partition_manifest_hash": "b" * 64,
        "seed": 295,
    }
    _field(generated, "/components/dominance/minimum_record_count_margin")["value"] = 2
    _field(generated, "/components/dominance/minimum_record_share")["value"] = 0.75
    _field(generated, "/components/dominance/minimum_record_share_margin")["value"] = 0.25
    _field(generated, "/components/dominance/xd_veto")["value"] = "missingness"
    validate_complete_inventory(generated, packaged_profile=packaged)

    generated_metadata = generated["generated_metadata"]
    assert isinstance(generated_metadata, dict)
    generated_metadata["unexpected"] = True
    with pytest.raises(ValueError, match="generated_metadata fields differ"):
        validate_complete_inventory(generated, packaged_profile=packaged)
