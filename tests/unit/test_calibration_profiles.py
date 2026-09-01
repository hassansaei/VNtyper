"""Generated profiles may change only dominance/abstention leaves."""

import hashlib
from collections.abc import Mapping
from pathlib import Path

import pytest

from vntyper.scripts.calibration_profiles import build_generated_profile, validate_generated_allowlist
from vntyper.scripts.decision_profile import load_packaged_decision_profile, resolve_decision_profile

pytestmark = pytest.mark.unit

_MUTABLE_POINTERS = {
    "/components/dominance/enabled",
    "/components/dominance/minimum_record_count_margin",
    "/components/dominance/minimum_record_share",
    "/components/dominance/minimum_record_share_margin",
    "/components/dominance/xd_veto",
    "/components/dominance/abstain_on_inadmissible_advntr",
}


def _active_component() -> dict[str, object]:
    return {
        "enabled": True,
        "minimum_record_count_margin": 2,
        "minimum_record_share": 0.75,
        "minimum_record_share_margin": 0.25,
        "xd_veto": "missingness",
        "abstain_on_inadmissible_advntr": True,
    }


def _build():
    return build_generated_profile(
        _active_component(),
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=295,
        objective="lexicographic-safety-v1",
        generator_version="2.0.20",
    )


def test_generation_changes_only_six_allowlisted_dominance_leaves_and_metadata() -> None:
    packaged = load_packaged_decision_profile()
    generated = _build()
    packaged_inventory = packaged.document["inventory"]
    generated_inventory = generated.document["inventory"]
    assert isinstance(packaged_inventory, Mapping) and isinstance(generated_inventory, Mapping)

    changed = {
        pointer
        for pointer in packaged_inventory
        if packaged_inventory[pointer]["value"] != generated_inventory[pointer]["value"]  # type: ignore[index]
    }

    assert changed == _MUTABLE_POINTERS
    assert generated.profile_kind == "generated"
    assert generated.source == "explicit-cli"
    assert generated.components["dominance"] == _active_component()
    validate_generated_allowlist(generated, packaged)


def test_generated_profile_records_complete_provenance_and_reloads_explicitly(tmp_path: Path) -> None:
    generated = _build()
    destination = tmp_path / "candidate.json"
    destination.write_bytes(generated.canonical_bytes)

    reloaded = resolve_decision_profile(destination)
    metadata = reloaded.document["generated_metadata"]

    assert reloaded.digest == generated.digest
    assert metadata["packaged_base_hash"] == load_packaged_decision_profile().digest  # type: ignore[index]
    assert metadata["generator_name"] == "vntyper-calibrate"  # type: ignore[index]
    assert metadata["generator_version"] == "2.0.20"  # type: ignore[index]
    assert metadata["objective"] == "lexicographic-safety-v1"  # type: ignore[index]
    assert metadata["dataset_manifest_hash"] == "a" * 64  # type: ignore[index]
    assert metadata["partition_manifest_hash"] == "b" * 64  # type: ignore[index]
    assert metadata["seed"] == 295  # type: ignore[index]
    assert not str(destination).startswith(str(Path("vntyper/profiles")))


def test_generation_does_not_mutate_packaged_profile_bytes_or_files() -> None:
    packaged_path = Path("vntyper/profiles/decision_profile.json")
    sidecar_path = Path("vntyper/profiles/decision_profile.sha256")
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (packaged_path, sidecar_path)}

    _build()

    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before} == before
    assert resolve_decision_profile().profile_kind == "packaged"
    assert not resolve_decision_profile().components["dominance"]["enabled"]  # type: ignore[index]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("enabled", "yes"),
        ("minimum_record_count_margin", -1),
        ("minimum_record_share", 1.1),
        ("minimum_record_share_margin", True),
        ("xd_veto", "winner"),
        ("abstain_on_inadmissible_advntr", 1),
    ],
)
def test_generated_dominance_values_remain_inside_closed_schema(field: str, value: object) -> None:
    component = _active_component()
    component[field] = value

    with pytest.raises(ValueError):
        build_generated_profile(
            component,
            dataset_manifest_hash="a" * 64,
            partition_manifest_hash="b" * 64,
            seed=295,
            objective="lexicographic-safety-v1",
            generator_version="2.0.20",
        )


def test_fixed_safety_and_excluded_fields_cannot_be_smuggled_into_generation() -> None:
    component = _active_component()
    component["bam_flank"] = 9

    with pytest.raises(ValueError, match="fields|allowlist"):
        build_generated_profile(
            component,
            dataset_manifest_hash="a" * 64,
            partition_manifest_hash="b" * 64,
            seed=295,
            objective="lexicographic-safety-v1",
            generator_version="2.0.20",
        )
