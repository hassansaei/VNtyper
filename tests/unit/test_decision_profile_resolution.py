"""Strict one-profile resolution and frozen run-context contracts."""

from __future__ import annotations

import copy
import json
from collections.abc import Mapping
from pathlib import Path

import pytest

from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.decision_profile import (
    load_packaged_decision_profile,
    parse_decision_profile,
    resolve_decision_profile,
)
from vntyper.scripts.run_configuration import resolve_compatibility_component, resolve_run_configuration

pytestmark = pytest.mark.unit

PROFILE_PATH = Path(__file__).parents[2] / "vntyper" / "profiles" / "decision_profile.json"


def _custom_document() -> dict[str, object]:
    document = load_strict_json_object(PROFILE_PATH.read_bytes())
    document["profile_id"] = "unit-test-explicit"
    document["profile_revision"] = "test-1"
    document["profile_kind"] = "explicit-custom"
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    pointer = "/components/kestrel/duplicate_flagging/flag_name"
    inventory[pointer]["value"] = "Unit_Test_Duplicate"
    return document


def test_no_argument_resolves_verified_packaged_profile() -> None:
    resolved = resolve_decision_profile()

    assert resolved.source == "package"
    assert resolved.profile_kind == "packaged"
    assert resolved.canonical_bytes == PROFILE_PATH.read_bytes()
    assert len(resolved.digest) == 64
    assert resolved.components["shark"] == {}


def test_explicit_profile_is_complete_canonicalized_and_never_overlaid(tmp_path: Path) -> None:
    document = _custom_document()
    path = tmp_path / "custom.json"
    path.write_text(json.dumps(document, indent=2), encoding="utf-8")

    resolved = resolve_decision_profile(path)

    assert resolved.source == "explicit-cli"
    assert resolved.profile_kind == "explicit-custom"
    assert resolved.canonical_bytes == canonical_json_bytes(document)
    kestrel = resolved.components["kestrel"]
    assert isinstance(kestrel, Mapping)
    duplicate_flagging = kestrel["duplicate_flagging"]
    assert isinstance(duplicate_flagging, Mapping)
    assert duplicate_flagging["flag_name"] == "Unit_Test_Duplicate"


@pytest.mark.parametrize(
    ("raw", "message"),
    [
        ('{"schema_version":1,"schema_version":1}', "duplicate JSON object key"),
        ('{"value":NaN}', "non-finite JSON constant"),
        ("[]", "top-level JSON value must be an object"),
    ],
)
def test_parser_rejects_ambiguous_or_nonstandard_json(raw: str, message: str) -> None:
    packaged = load_packaged_decision_profile()

    with pytest.raises(ValueError, match=message):
        parse_decision_profile(raw, packaged_document=packaged.document)


@pytest.mark.parametrize(
    ("pointer", "value", "message"),
    [
        (
            "/components/kestrel/flagging_rules/False_Positive_4bp_Insertion/all/0/operator",
            "definitely-not-an-operator",
            "operator",
        ),
        ("/components/kestrel/motif_filtering/position_threshold", -500, "position_threshold"),
        ("/components/kestrel/selection/confidence_priority/High_Precision", -1, "confidence_priority"),
        ("/components/kestrel/selection/sort_order/0/column", "not-a-column", "sort_order"),
        ("/components/kestrel/artifact_flags/0", "Missing_Flag", "artifact_flags"),
        ("/components/kestrel/duplicate_flagging/enabled", "yes", "enabled"),
        ("/components/kestrel/duplicate_flagging/sort_by/0/column", "missing", "sort_by"),
        ("/components/kestrel/selection/unflagged_value", "", "unflagged_value"),
    ],
)
def test_resolution_rejects_semantically_invalid_mutable_values(pointer: str, value: object, message: str) -> None:
    packaged = load_packaged_decision_profile()
    document = _custom_document()
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory[pointer]["value"] = value

    with pytest.raises(ValueError, match=message):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def test_resolution_rejects_a_governed_rule_evidence_mismatch() -> None:
    packaged = load_packaged_decision_profile()
    document = _custom_document()
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    pointer = "/components/advntr/flagging_rules/Polymorphic_Call/all/0/right/literal/0"
    inventory[pointer]["value"] = "NOT_A_GOVERNED_STATE"

    with pytest.raises(ValueError, match="fixed-safety|governed|Polymorphic"):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def test_run_resolution_rejects_non_neutral_generated_dominance_until_it_has_a_consumer(tmp_path: Path) -> None:
    packaged = load_packaged_decision_profile()
    document = json.loads(packaged.canonical_bytes)
    document["profile_id"] = "unit-test-generated-active"
    document["profile_revision"] = "test-1"
    document["profile_kind"] = "generated"
    document["generated_metadata"] = {
        "packaged_base_hash": packaged.digest,
        "generator_name": "unit-test",
        "generator_version": "1",
        "objective": "lexicographic-safety-v1",
        "dataset_manifest_hash": "a" * 64,
        "partition_manifest_hash": "b" * 64,
        "seed": 295,
    }
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/dominance/enabled"]["value"] = True
    path = tmp_path / "generated.json"
    path.write_bytes(canonical_json_bytes(document))

    with pytest.raises(ValueError, match="dominance.*not active|neutral"):
        resolve_run_configuration(path)


@pytest.mark.parametrize("kind", ["packaged", "unsupported"])
def test_explicit_path_cannot_select_packaged_or_unknown_kind(tmp_path: Path, kind: str) -> None:
    document = _custom_document()
    document["profile_kind"] = kind
    path = tmp_path / "bad-kind.json"
    path.write_bytes(canonical_json_bytes(document))

    with pytest.raises(ValueError, match="profile_kind|profile kind|explicit CLI"):
        resolve_decision_profile(path)


def test_missing_unknown_and_immutable_fields_fail_without_fallback(tmp_path: Path) -> None:
    packaged = load_packaged_decision_profile()
    for _label, mutate, message in (
        (
            "missing",
            lambda inventory: inventory.pop("/components/cross_match/required_advntr_evidence_disposition"),
            "inventory fields differ",
        ),
        (
            "unknown",
            lambda inventory: inventory.__setitem__(
                "/components/cross_match/unknown", {"class": "explicit-custom", "value": True}
            ),
            "inventory fields differ",
        ),
        (
            "immutable",
            lambda inventory: inventory["/components/nomenclature/thresholds/bam_flank"].__setitem__("value", 9),
            "fixed-safety field",
        ),
    ):
        document = _custom_document()
        inventory = document["inventory"]
        assert isinstance(inventory, dict)
        mutate(inventory)
        with pytest.raises(ValueError, match=message, check=lambda error: bool(str(error))):
            parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


@pytest.mark.parametrize("name", ["missing.json", "directory"])
def test_unreadable_explicit_path_fails_instead_of_falling_back(tmp_path: Path, name: str) -> None:
    path = tmp_path / name
    if name == "directory":
        path.mkdir()

    with pytest.raises(ValueError, match="cannot read explicit decision profile"):
        resolve_decision_profile(path)


def test_resolved_components_are_recursively_immutable() -> None:
    run = resolve_run_configuration()

    with pytest.raises(TypeError):
        run.kestrel["new"] = True  # type: ignore[index]
    confidence = run.kestrel["confidence_assignment"]
    assert isinstance(confidence, dict | tuple) is False
    with pytest.raises(TypeError):
        confidence["new"] = True  # type: ignore[index]


def test_resolved_profile_document_is_recursively_immutable() -> None:
    profile = resolve_decision_profile()

    with pytest.raises(TypeError):
        profile.document["profile_revision"] = "mutated"  # type: ignore[index]
    inventory = profile.document["inventory"]
    assert isinstance(inventory, Mapping)
    field = inventory["/components/kestrel/duplicate_flagging/flag_name"]
    assert isinstance(field, Mapping)
    with pytest.raises(TypeError):
        field["value"] = "mutated"  # type: ignore[index]


def test_run_configuration_freezes_excluded_stage_runtime_configuration() -> None:
    run = resolve_run_configuration()

    assert run.kestrel_runtime["kestrel_settings"]["java_memory"] == "12g"  # type: ignore[index]
    assert run.advntr_runtime["settings"] == {"additional_commands": "", "threads": None}
    assert run.shark_runtime["shark_settings"]["muc1_region_fasta_hg19"]  # type: ignore[index]
    with pytest.raises(TypeError):
        run.advntr_runtime["settings"]["threads"] = 8  # type: ignore[index]


def test_custom_context_cannot_implicitly_load_a_packaged_component() -> None:
    with pytest.raises(ValueError, match="custom Kestrel run context requires an explicit resolved component"):
        resolve_compatibility_component("kestrel", None, custom_context_active=True)


def test_legacy_context_can_resolve_the_packaged_component() -> None:
    resolved = resolve_compatibility_component("kestrel", None, custom_context_active=False)

    assert resolved == resolve_run_configuration().kestrel


def test_generated_profile_requires_exact_closed_metadata() -> None:
    packaged = load_packaged_decision_profile()
    document = load_strict_json_object(PROFILE_PATH.read_bytes())
    document["profile_id"] = "unit-test-generated"
    document["profile_revision"] = "test-1"
    document["profile_kind"] = "generated"
    document["generated_metadata"] = {
        "packaged_base_hash": packaged.digest,
        "generator_name": "unit-test",
        "generator_version": "1",
        "objective": "lexicographic-safety-v1",
        "dataset_manifest_hash": "a" * 64,
        "partition_manifest_hash": "b" * 64,
        "seed": 295,
    }
    baseline = copy.deepcopy(document)
    parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)

    document["generated_metadata"]["extra"] = True
    with pytest.raises(ValueError, match="generated_metadata fields differ"):
        parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)
    assert baseline != document
