"""Build complete versioned decision profiles from calibrated caller policies."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from vntyper.scripts.calibration_caller_policy import (
    ADVNTR_CALLER_POLICY_POINTERS,
    KESTREL_CALLER_POLICY_POINTERS,
    CallerPolicyValues,
    caller_policy_values_document,
    decode_caller_policy_values,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)
from vntyper.scripts.decision_profile_schema import CALLER_ADVNTR_FIELD_METADATA

_OBJECTIVE = "caller-safety-v1"


def _require_policy(policy: CallerPolicyValues) -> CallerPolicyValues:
    if not isinstance(policy, CallerPolicyValues):
        raise ValueError("caller-generated profile requires CallerPolicyValues")
    decoded = decode_caller_policy_values(caller_policy_values_document(policy))
    if decoded != policy:
        raise ValueError("caller-generated profile policy differs from its canonical content or digest")
    return decoded


def _caller_field(pointer: str, value: object) -> dict[str, object]:
    unit, comparator, inclusive = CALLER_ADVNTR_FIELD_METADATA[pointer]
    field: dict[str, object] = {"class": "generated-mutable", "value": value}
    if unit is not None:
        field.update({"unit": unit, "comparator": comparator, "inclusive": inclusive})
    return field


def build_caller_generated_profile(
    policy: CallerPolicyValues,
    *,
    dataset_manifest_hash: str,
    partition_manifest_hash: str,
    seed: int,
    generator_version: str,
    packaged_profile: ResolvedDecisionProfile | None = None,
) -> ResolvedDecisionProfile:
    """Build a complete schema-v2 profile for one validated caller policy.

    Args:
        policy: Exact Kestrel and conditional adVNTR calibrated values.
        dataset_manifest_hash: Canonical calibration dataset commitment.
        partition_manifest_hash: Canonical role-partition commitment.
        seed: Registered non-negative calibration seed.
        generator_version: Non-empty generator identity.
        packaged_profile: Optional verified packaged baseline.

    Returns:
        A validated complete generated profile.

    Raises:
        ValueError: If any typed input, hash, policy, or packaged baseline is invalid.
    """
    validated_policy = _require_policy(policy)
    if not isinstance(generator_version, str) or not generator_version:
        raise ValueError("caller-generated profile generator version must be a non-empty string")
    packaged = packaged_profile or load_packaged_decision_profile()
    if not isinstance(packaged, ResolvedDecisionProfile) or packaged.profile_kind != "packaged":
        raise ValueError("caller-generated profile base must be the verified packaged profile")
    document = load_strict_json_object(packaged.canonical_bytes)
    inventory = cast(dict[str, dict[str, object]], document["inventory"])
    for pointer in KESTREL_CALLER_POLICY_POINTERS:
        inventory[pointer]["class"] = "generated-mutable"
        inventory[pointer]["value"] = validated_policy.values[pointer]
    if "advntr" in validated_policy.required_callers:
        for pointer in ADVNTR_CALLER_POLICY_POINTERS:
            inventory[pointer] = _caller_field(pointer, validated_policy.values[pointer])
    generated_pointers = sorted(validated_policy.values)
    identity_payload: Mapping[str, object] = {
        "caller_policy_sha256": validated_policy.sha256,
        "dataset_manifest_hash": dataset_manifest_hash,
        "generator_version": generator_version,
        "packaged_base_hash": packaged.digest,
        "partition_manifest_hash": partition_manifest_hash,
        "seed": seed,
    }
    document.update(
        {
            "schema_version": 2,
            "profile_id": f"vntyper-caller-generated-{canonical_sha256(identity_payload)[:16]}",
            "profile_revision": "caller-v2",
            "profile_kind": "generated",
            "generated_metadata": {
                "packaged_base_hash": packaged.digest,
                "generator_name": "vntyper-calibrate",
                "generator_version": generator_version,
                "objective": _OBJECTIVE,
                "dataset_manifest_hash": dataset_manifest_hash,
                "partition_manifest_hash": partition_manifest_hash,
                "seed": seed,
                "generation_target": "callers",
                "caller_policy_sha256": validated_policy.sha256,
                "generated_pointers": generated_pointers,
                "required_callers": list(validated_policy.required_callers),
            },
        }
    )
    return parse_decision_profile(
        canonical_json_bytes(document),
        packaged_document=packaged.document,
        allow_caller_generated=True,
    )
