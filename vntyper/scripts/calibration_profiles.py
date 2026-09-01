"""Generated opt-in decision profiles restricted to dominance/abstention leaves."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from vntyper.scripts.calibration_objective import count_free_parameters
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)

GENERATED_DOMINANCE_POINTERS = {
    "/components/dominance/enabled": "enabled",
    "/components/dominance/minimum_record_count_margin": "minimum_record_count_margin",
    "/components/dominance/minimum_record_share": "minimum_record_share",
    "/components/dominance/minimum_record_share_margin": "minimum_record_share_margin",
    "/components/dominance/xd_veto": "xd_veto",
    "/components/dominance/abstain_on_inadmissible_advntr": "abstain_on_inadmissible_advntr",
}


def build_generated_profile(
    dominance_component: Mapping[str, object],
    *,
    dataset_manifest_hash: str,
    partition_manifest_hash: str,
    seed: int,
    objective: str,
    generator_version: str,
    packaged_profile: ResolvedDecisionProfile | None = None,
) -> ResolvedDecisionProfile:
    """Build one complete explicit generated profile without mutating package data."""
    count_free_parameters(dominance_component)
    if objective != "lexicographic-safety-v1":
        raise ValueError("generated profile objective must be lexicographic-safety-v1")
    if not isinstance(generator_version, str) or not generator_version:
        raise ValueError("generated profile generator version must be a non-empty string")
    packaged = packaged_profile or load_packaged_decision_profile()
    if not isinstance(packaged, ResolvedDecisionProfile) or packaged.profile_kind != "packaged":
        raise ValueError("generated profile base must be the verified packaged profile")
    document = _mutable_packaged_document(packaged)
    inventory = cast(dict[str, dict[str, object]], document["inventory"])
    for pointer, component_key in GENERATED_DOMINANCE_POINTERS.items():
        inventory[pointer]["value"] = dominance_component[component_key]
    identity_payload = {
        "packaged_base_hash": packaged.digest,
        "dominance": dict(dominance_component),
        "objective": objective,
        "dataset_manifest_hash": dataset_manifest_hash,
        "partition_manifest_hash": partition_manifest_hash,
        "seed": seed,
        "generator_version": generator_version,
    }
    document.update(
        {
            "profile_id": f"vntyper-generated-{canonical_sha256(identity_payload)[:16]}",
            "profile_revision": "1",
            "profile_kind": "generated",
            "generated_metadata": {
                "packaged_base_hash": packaged.digest,
                "generator_name": "vntyper-calibrate",
                "generator_version": generator_version,
                "objective": objective,
                "dataset_manifest_hash": dataset_manifest_hash,
                "partition_manifest_hash": partition_manifest_hash,
                "seed": seed,
            },
        }
    )
    generated = parse_decision_profile(
        canonical_json_bytes(document),
        packaged_document=packaged.document,
    )
    validate_generated_allowlist(generated, packaged)
    return generated


def validate_generated_allowlist(
    generated: ResolvedDecisionProfile,
    packaged: ResolvedDecisionProfile,
) -> None:
    """Prove all non-generated inventory fields remain byte-for-value identical."""
    if not isinstance(generated, ResolvedDecisionProfile) or generated.profile_kind != "generated":
        raise ValueError("generated allowlist requires a generated ResolvedDecisionProfile")
    if not isinstance(packaged, ResolvedDecisionProfile) or packaged.profile_kind != "packaged":
        raise ValueError("generated allowlist requires a packaged ResolvedDecisionProfile base")
    raw_generated = generated.document.get("inventory")
    raw_packaged = packaged.document.get("inventory")
    if not isinstance(raw_generated, Mapping) or not isinstance(raw_packaged, Mapping):
        raise ValueError("generated allowlist requires complete decision inventories")
    if set(raw_generated) != set(raw_packaged):
        raise ValueError("generated profile inventory fields differ from the packaged allowlist")
    for pointer in raw_packaged:
        generated_field = raw_generated[pointer]
        packaged_field = raw_packaged[pointer]
        if pointer not in GENERATED_DOMINANCE_POINTERS:
            if generated_field != packaged_field:
                raise ValueError(f"generated profile changed a non-allowlisted field: {pointer}")
            continue
        if not isinstance(generated_field, Mapping) or not isinstance(packaged_field, Mapping):
            raise ValueError(f"generated profile field is malformed: {pointer}")
        generated_metadata = {key: value for key, value in generated_field.items() if key != "value"}
        packaged_metadata = {key: value for key, value in packaged_field.items() if key != "value"}
        if generated_metadata != packaged_metadata:
            raise ValueError(f"generated profile changed field semantics outside its value: {pointer}")
    metadata = generated.document.get("generated_metadata")
    if not isinstance(metadata, Mapping) or metadata.get("packaged_base_hash") != packaged.digest:
        raise ValueError("generated profile metadata does not bind the packaged base hash")


def _mutable_packaged_document(packaged: ResolvedDecisionProfile) -> dict[str, object]:
    return load_strict_json_object(packaged.canonical_bytes)
