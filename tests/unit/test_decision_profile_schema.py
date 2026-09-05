"""Contracts for complete, canonical VNtyper decision profiles."""

from __future__ import annotations

import copy
import hashlib
import re
from collections import Counter
from collections.abc import Mapping
from pathlib import Path

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import load_packaged_decision_profile, parse_decision_profile
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
EXPECTED_PACKAGED_PROFILE_SHA256 = "0b13d07370491b3ea773e65144891cb30caebcae70b0ef98feb0f2c5ccd2f4a1"


def _packaged_profile() -> dict[str, object]:
    return load_strict_json_object(PROFILE_PATH.read_bytes())


def _field(profile: dict[str, object], pointer: str) -> dict[str, object]:
    inventory = profile["inventory"]
    assert isinstance(inventory, dict)
    value = inventory[pointer]
    assert isinstance(value, dict)
    return value


def _resolved_pointer(components: Mapping[str, object], pointer: str) -> object:
    """Resolve one inventory pointer against immutable runtime components."""
    tokens = pointer.removeprefix("/").split("/")
    assert tokens.pop(0) == "components"
    value: object = components
    for token in tokens:
        if isinstance(value, Mapping):
            value = value[token]
        else:
            assert isinstance(value, tuple)
            value = value[int(token)]
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
    assert digest == EXPECTED_PACKAGED_PROFILE_SHA256
    validate_complete_inventory(profile)


def test_governed_advntr_evidence_and_cross_match_disposition_are_independently_fixed() -> None:
    profile = _packaged_profile()
    inventory = profile["inventory"]
    assert isinstance(inventory, dict)

    governed_prefix = "/components/advntr/flagging_rules/Polymorphic_Call/"
    governed = {pointer for pointer in inventory if pointer.startswith(governed_prefix)}
    assert governed
    assert all(inventory[pointer]["class"] == ValidationClass.FIXED_SAFETY.value for pointer in governed)
    assert _field(profile, "/components/cross_match/required_advntr_evidence_disposition") == {
        "class": ValidationClass.FIXED_SAFETY.value,
        "value": "admissible",
    }


def test_polymorphic_rule_states_are_exactly_the_governed_artifact_states() -> None:
    projection = load_strict_json_object(PROJECTION_PATH.read_bytes())
    advntr = projection["advntr"]
    assert isinstance(advntr, dict)
    evidence = advntr["artifact_evidence"]
    rules = advntr["flagging_rules"]
    assert isinstance(evidence, dict) and isinstance(rules, dict)

    governed = evidence["active_states"]
    polymorphic = rules["Polymorphic_Call"]
    assert polymorphic["all"][0]["right"]["literal"] == governed


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
    assert Counter(field["class"] for field in inventory.values()) == {
        ValidationClass.FIXED_SAFETY.value: 173,
        ValidationClass.EXPLICIT_CUSTOM.value: 69,
        ValidationClass.GENERATED_MUTABLE.value: 6,
    }


def test_one_complete_explicit_profile_changes_every_mutable_stage_argument() -> None:
    """Every field advertised as operator-editable reaches a non-default runtime value."""
    packaged = _packaged_profile()
    explicit = copy.deepcopy(packaged)
    explicit["profile_kind"] = "explicit-custom"
    explicit["profile_id"] = "unit-test-all-explicit-fields"
    inventory = explicit["inventory"]
    assert isinstance(inventory, dict)
    mutable_pointers = {
        pointer
        for pointer, field in inventory.items()
        if field["class"] in {ValidationClass.EXPLICIT_CUSTOM.value, ValidationClass.GENERATED_MUTABLE.value}
    }
    column_replacements = {
        "NumberOfSupportingReads": "MeanCoverage",
        "RU": "Variant",
        "Kestrel_Allele_Change": "Kestrel_REF",
        "Advntr_Allele_Change": "Advntr_REF",
        "Kestrel_Variant_Type": "Kestrel_ALT",
        "Advntr_Variant_Type": "Advntr_ALT",
        "REF": "ALT",
        "ALT": "REF",
        "Depth_Score": "haplo_count",
        "Motif": "Variant",
        "confidence_priority": "is_unflagged",
        "is_unflagged": "Depth_Score",
        "haplo_count": "POS",
        "POS": "confidence_priority",
    }
    operator_replacements = {"lt": "eq", "eq": "casefold_eq", "casefold_eq": "eq"}
    dominance_replacements = {
        "/components/dominance/enabled": True,
        "/components/dominance/minimum_record_count_margin": 2,
        "/components/dominance/minimum_record_share": 0.75,
        "/components/dominance/minimum_record_share_margin": 0.25,
        "/components/dominance/xd_veto": "missingness",
        "/components/dominance/abstain_on_inadmissible_advntr": True,
    }
    changed: set[str] = set()
    for pointer in sorted(mutable_pointers):
        field = inventory[pointer]
        value = field["value"]
        if pointer in dominance_replacements:
            field["value"] = dominance_replacements[pointer]
        elif pointer.endswith("/operator"):
            field["value"] = operator_replacements[value]
        elif "/duplicate_flagging/sort_by/" in pointer and pointer.endswith("/column"):
            field["value"] = "POS"
        elif pointer.endswith("/column"):
            field["value"] = column_replacements[value]
        elif "/duplicate_flagging/group_by/" in pointer:
            field["value"] = {"REF": "ALT", "ALT": "REF"}[value]
        elif pointer.endswith("/artifact_flags/0"):
            field["value"] = "Low_Depth_Conserved_Motifs"
        elif isinstance(value, bool):
            field["value"] = not value
        elif isinstance(value, (int, float)):
            field["value"] = value + 10
        elif isinstance(value, list):
            field["value"] = ["CUSTOM"]
        elif isinstance(value, str):
            field["value"] = f"{value}-custom"
        else:
            raise AssertionError(f"no valid differential fixture for {pointer}: {value!r}")
        if field["value"] != value:
            changed.add(pointer)

    assert changed == mutable_pointers
    validate_complete_inventory(explicit, packaged_profile=packaged)
    packaged_resolved = load_packaged_decision_profile()
    explicit_resolved = parse_decision_profile(
        canonical_json_bytes(explicit),
        packaged_document=packaged_resolved.document,
    )
    for pointer in mutable_pointers:
        assert _resolved_pointer(explicit_resolved.components, pointer) != _resolved_pointer(
            packaged_resolved.components, pointer
        )


def test_packaged_inventory_rejects_an_invalid_fixed_output_format() -> None:
    packaged = _packaged_profile()
    _field(packaged, "/components/advntr/settings/output_format")["value"] = "bam"

    with pytest.raises(ValueError, match="output_format"):
        validate_complete_inventory(packaged)


@pytest.mark.parametrize(
    ("pointer", "value", "unit", "comparator", "inclusive"),
    [
        (
            "/components/kestrel/confidence_assignment/reporting_floor",
            0.00469,
            "depth-score-ratio",
            "gte",
            True,
        ),
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


def test_reporting_floor_and_independent_gg_gate_are_linked() -> None:
    profile = _packaged_profile()
    floor = _field(profile, "/components/kestrel/confidence_assignment/reporting_floor")
    gg = _field(profile, "/components/kestrel/alt_filtering/gg_depth_score_threshold")

    assert floor["value"] == gg["value"] == 0.00469
    assert floor is not gg


def test_revision_1_profile_lacking_reporting_floor_fails_closed_with_inventory_fields_differ() -> None:
    """A revision-1 custom profile lacking reporting_floor must fail closed (#311)."""
    packaged = _packaged_profile()
    custom = copy.deepcopy(packaged)
    custom["profile_kind"] = "explicit-custom"
    custom["profile_id"] = "unit-test-revision-1-lacks-floor"
    inventory = custom["inventory"]
    assert isinstance(inventory, dict)
    inventory.pop("/components/kestrel/confidence_assignment/reporting_floor", None)

    with pytest.raises(ValueError, match="inventory fields differ"):
        validate_complete_inventory(custom, packaged_profile=packaged)


def test_divergent_reporting_floor_and_gg_depth_threshold_raises_in_critical_fields() -> None:
    """Critical field validation must reject divergent reporting_floor and gg_depth_score_threshold (#311)."""
    packaged = _packaged_profile()
    custom = copy.deepcopy(packaged)
    custom["profile_kind"] = "explicit-custom"
    custom["profile_id"] = "unit-test-divergent-floor-gg"
    inventory = custom["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/kestrel/confidence_assignment/reporting_floor"]["value"] = 0.00500

    with pytest.raises(
        ValueError,
        match="independent GG depth-score minimum must equal the reporting floor|critical fixed-safety field differs|immutable fixed-safety field differs",
    ):
        validate_complete_inventory(custom, packaged_profile=packaged)


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
    if isinstance(value, dict):
        return {"unexpected": True}
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
    _field(generated, "/components/dominance/enabled")["value"] = True
    _field(generated, "/components/dominance/minimum_record_count_margin")["value"] = 2
    _field(generated, "/components/dominance/minimum_record_share")["value"] = 0.75
    _field(generated, "/components/dominance/minimum_record_share_margin")["value"] = 0.25
    _field(generated, "/components/dominance/xd_veto")["value"] = "missingness"
    _field(generated, "/components/dominance/abstain_on_inadmissible_advntr")["value"] = True
    generated_inventory = generated["inventory"]
    assert isinstance(generated_inventory, dict)
    assert all(
        field["value"] != _field(packaged, pointer)["value"]
        for pointer, field in generated_inventory.items()
        if field["class"] == ValidationClass.GENERATED_MUTABLE.value
    )
    validate_complete_inventory(generated, packaged_profile=packaged)

    generated_metadata = generated["generated_metadata"]
    assert isinstance(generated_metadata, dict)
    generated_metadata["unexpected"] = True
    with pytest.raises(ValueError, match="generated_metadata fields differ"):
        validate_complete_inventory(generated, packaged_profile=packaged)
