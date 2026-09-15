"""Tests for versioned caller-generated decision profiles."""

from __future__ import annotations

from dataclasses import replace
from typing import cast

import pytest

from tests.unit.test_calibration_caller_policy import policy_document
from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    KESTREL_CALLER_POLICY_POINTERS,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_profiles import build_generated_profile
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.decision_profile import load_packaged_decision_profile, parse_decision_profile

pytestmark = pytest.mark.unit


def _policy(*, advntr: bool):
    document = policy_document(include_advntr=advntr)
    values = document["values"]
    assert isinstance(values, dict)
    values["/components/kestrel/confidence_assignment/reporting_floor"] = 0.003
    values["/components/kestrel/alt_filtering/gg_depth_score_threshold"] = 0.004
    if advntr:
        values["/components/advntr/calibrated_calling/rare_unit_fraction"] = None
    return decode_caller_policy_values(document)


def _build(*, advntr: bool = False):
    from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile

    return build_caller_generated_profile(
        _policy(advntr=advntr),
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=295,
        generator_version="unit-test",
    )


def test_caller_profile_v2_enumerates_exact_kestrel_capability_and_preserves_v1() -> None:
    policy = _policy(advntr=False)
    generated = _build()
    document = load_strict_json_object(generated.canonical_bytes)
    metadata = document["generated_metadata"]
    inventory = document["inventory"]
    assert isinstance(metadata, dict) and isinstance(inventory, dict)

    assert document["schema_version"] == 2
    assert metadata["generation_target"] == "callers"
    assert metadata["generated_pointers"] == list(KESTREL_CALLER_POLICY_POINTERS)
    assert metadata["required_callers"] == ["kestrel"]
    assert not any(pointer.startswith("/components/advntr/calibrated_calling/") for pointer in inventory)
    for pointer in KESTREL_CALLER_POLICY_POINTERS:
        assert inventory[pointer]["class"] == "generated-mutable"
        assert inventory[pointer]["value"] == policy.values[pointer]

    packaged = load_packaged_decision_profile()
    legacy = build_generated_profile(
        cast(dict[str, object], packaged.components["dominance"]),
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=295,
        objective="lexicographic-safety-v1",
        generator_version="unit-test",
        packaged_profile=packaged,
    )
    assert load_strict_json_object(legacy.canonical_bytes)["schema_version"] == 1


def test_advntr_profile_adds_only_the_closed_conditional_subtree_with_nullable_numeric_semantics() -> None:
    generated = _build(advntr=True)
    document = load_strict_json_object(generated.canonical_bytes)
    metadata = document["generated_metadata"]
    inventory = document["inventory"]
    assert isinstance(metadata, dict) and isinstance(inventory, dict)

    assert metadata["generated_pointers"] == sorted((*KESTREL_CALLER_POLICY_POINTERS, *ADVNTR_CALLER_POLICY_POINTERS))
    assert metadata["required_callers"] == ["advntr", "kestrel"]
    assert {pointer for pointer in inventory if pointer.startswith("/components/advntr/calibrated_calling/")} == set(
        ADVNTR_CALLER_POLICY_POINTERS
    )
    rare = inventory["/components/advntr/calibrated_calling/rare_unit_fraction"]
    assert rare == {
        "class": "generated-mutable",
        "value": None,
        "unit": "eligible-repeat-unit-fraction",
        "comparator": "gte-when-enabled",
        "inclusive": True,
    }


def test_caller_profile_parser_rejects_unlisted_changes_and_tampered_policy_binding() -> None:
    packaged = load_packaged_decision_profile()
    generated = _build(advntr=True)
    document = load_strict_json_object(generated.canonical_bytes)
    inventory = document["inventory"]
    metadata = document["generated_metadata"]
    assert isinstance(inventory, dict) and isinstance(metadata, dict)

    fixed = "/components/nomenclature/thresholds/bam_flank"
    inventory[fixed]["value"] = 9
    with pytest.raises(ValueError, match="fixed-safety|non-caller"):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)

    document = load_strict_json_object(generated.canonical_bytes)
    metadata = document["generated_metadata"]
    assert isinstance(metadata, dict)
    metadata["caller_policy_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="caller policy"):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def test_caller_profile_public_validation_rejects_omitted_capability_pointer() -> None:
    packaged = load_packaged_decision_profile()
    generated = _build()
    document = load_strict_json_object(generated.canonical_bytes)
    metadata = document["generated_metadata"]
    assert isinstance(metadata, dict)
    metadata["generated_pointers"] = list(KESTREL_CALLER_POLICY_POINTERS[:-1])

    with pytest.raises(ValueError, match="generated_pointers"):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def test_caller_profile_builder_revalidates_directly_forged_typed_policy() -> None:
    from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile

    forged = replace(_policy(advntr=False), sha256="f" * 64)
    with pytest.raises(ValueError, match="canonical content|digest"):
        build_caller_generated_profile(
            forged,
            dataset_manifest_hash="a" * 64,
            partition_manifest_hash="b" * 64,
            seed=295,
            generator_version="unit-test",
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("objective", "other", "objective"),
        ("generation_target", "length", "generation_target"),
        ("required_callers", ["kestrel", "advntr"], "required_callers"),
    ],
)
def test_caller_profile_parser_rejects_noncanonical_generation_metadata(
    field: str, value: object, message: str
) -> None:
    packaged = load_packaged_decision_profile()
    document = load_strict_json_object(_build(advntr=True).canonical_bytes)
    metadata = document["generated_metadata"]
    assert isinstance(metadata, dict)
    metadata[field] = value

    with pytest.raises(ValueError, match=message):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def test_caller_profile_parser_rejects_changed_caller_field_semantics() -> None:
    packaged = load_packaged_decision_profile()
    document = load_strict_json_object(_build().canonical_bytes)
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory[KESTREL_CALLER_POLICY_POINTERS[0]]["unit"] = "other"

    with pytest.raises(ValueError, match="semantics"):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)
