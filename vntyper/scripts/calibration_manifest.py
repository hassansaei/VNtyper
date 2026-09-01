"""Role-bound partition manifests and transitive leakage groups."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.calibration_contract import (
    CalibrationProtocol,
    CalibrationRole,
    EvidenceProvenance,
    decode_protocol,
)
from vntyper.scripts.canonical_json import canonical_sha256

CalibrationOperation = Literal["fit", "validate", "evaluate"]

GROUP_NAMESPACES = (
    "individual-family",
    "simulated-pair",
    "backbone-seed-lineage",
    "replicate-rerun",
    "depth-series-source",
    "batch",
    "repeat-context",
)
_ROLES = frozenset({"training", "policy-selection", "validation", "locked-heldout"})
_PROVENANCE = frozenset({"development", "external-custodian"})
_OPERATION_ROLES: dict[str, tuple[CalibrationRole, ...]] = {
    "fit": ("training", "policy-selection"),
    "validate": ("validation",),
    "evaluate": ("locked-heldout",),
}


@dataclass(frozen=True)
class PartitionMember:
    """One immutable study member and all predeclared leakage keys."""

    key: str
    role: CalibrationRole
    provenance: EvidenceProvenance
    assay_class: str
    groups: Mapping[str, tuple[str, ...]]


@dataclass(frozen=True)
class PartitionManifest:
    """Complete role-bound manifest with no cross-role connected component."""

    members: tuple[PartitionMember, ...]
    sha256: str


@dataclass(frozen=True)
class StudyDeclaration:
    """Complete canonical PARTITIONS declaration with all four roles."""

    protocol: CalibrationProtocol
    partitions: PartitionManifest
    sha256: str


def decode_study_declaration(value: object) -> StudyDeclaration:
    """Decode the complete protocol and four immutable split roles."""
    root = _exact_object(value, {"schema_version", "protocol", "partitions"}, "calibration study declaration")
    if root["schema_version"] != "calibration-study-v1":
        raise ValueError("calibration study schema version must be calibration-study-v1")
    protocol = decode_protocol(root["protocol"])
    partitions = decode_partition_manifest(root["partitions"])
    roles = {member.role for member in partitions.members}
    if roles != set(_ROLES):
        raise ValueError("calibration study declaration must contain all four partition roles")
    undeclared_assays = sorted({member.assay_class for member in partitions.members} - set(protocol.assay_classes))
    if undeclared_assays:
        raise ValueError(f"calibration partition assay classes are not declared by the protocol: {undeclared_assays}")
    return StudyDeclaration(protocol, partitions, canonical_sha256(root))


def decode_partition_manifest(value: object) -> PartitionManifest:
    """Decode a strict manifest and reject direct or transitive split leakage."""
    root = _exact_object(value, {"schema_version", "members"}, "calibration partition manifest")
    if root["schema_version"] != "calibration-partitions-v1":
        raise ValueError("calibration partition manifest schema version must be calibration-partitions-v1")
    raw_members = root["members"]
    if not isinstance(raw_members, list) or not raw_members:
        raise ValueError("calibration partition manifest members must be a non-empty list")
    members = tuple(_decode_member(raw) for raw in raw_members)
    keys = tuple(member.key for member in members)
    if keys != tuple(sorted(keys)) or len(keys) != len(set(keys)):
        raise ValueError("calibration partition member keys must be unique and increasing")
    manifest = PartitionManifest(members, canonical_sha256(root))
    groups = connected_leakage_groups(manifest)
    roles_by_group: dict[str, set[CalibrationRole]] = {}
    for member in members:
        roles_by_group.setdefault(groups[member.key], set()).add(member.role)
    crossing = {group: roles for group, roles in roles_by_group.items() if len(roles) > 1}
    if crossing:
        raise ValueError(f"transitive calibration leakage groups cross roles: {crossing}")
    return manifest


def connected_leakage_groups(manifest: PartitionManifest) -> Mapping[str, str]:
    """Return a deterministic component digest for every manifest member."""
    if not isinstance(manifest, PartitionManifest):
        raise ValueError("calibration partition manifest must be a PartitionManifest")
    parent = {member.key: member.key for member in manifest.members}

    def find(key: str) -> str:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            first, second = sorted((left_root, right_root))
            parent[second] = first

    owners: dict[tuple[str, str], str] = {}
    for member in manifest.members:
        for namespace, values in member.groups.items():
            for value in values:
                group_key = (namespace, value)
                previous = owners.setdefault(group_key, member.key)
                union(member.key, previous)
    components: dict[str, list[str]] = {}
    for key in parent:
        components.setdefault(find(key), []).append(key)
    digests = {key: canonical_sha256(sorted(components[find(key)])) for key in sorted(parent)}
    return MappingProxyType(digests)


def require_operation_roles(
    manifest: PartitionManifest,
    operation: CalibrationOperation,
    *,
    requested_roles: Sequence[str] | None = None,
) -> tuple[CalibrationRole, ...]:
    """Authorize only the exact partition roles needed by one operation."""
    if not isinstance(manifest, PartitionManifest):
        raise ValueError("calibration partition manifest must be a PartitionManifest")
    expected = _OPERATION_ROLES.get(operation)
    if expected is None:
        raise ValueError(f"unsupported calibration operation: {operation}")
    present = {member.role for member in manifest.members}
    if not set(expected) <= present:
        raise ValueError(f"calibration {operation} requires roles {expected}")
    if requested_roles is not None and tuple(requested_roles) != expected:
        raise ValueError(f"calibration {operation} may access only roles {expected}")
    return expected


def _decode_member(value: object) -> PartitionMember:
    raw = _exact_object(
        value,
        {"key", "role", "provenance", "assay_class", "groups"},
        "calibration partition member",
    )
    key = raw["key"]
    if not isinstance(key, str) or not key:
        raise ValueError("calibration partition member key must be a non-empty string")
    role = raw["role"]
    if role not in _ROLES:
        raise ValueError(f"unsupported calibration partition role: {role!r}")
    provenance = raw["provenance"]
    if provenance not in _PROVENANCE:
        raise ValueError(f"unsupported calibration partition provenance: {provenance!r}")
    if role == "locked-heldout" and provenance != "external-custodian":
        raise ValueError("locked held-out partition members require external custodian provenance")
    assay_class = raw["assay_class"]
    if not isinstance(assay_class, str) or not assay_class:
        raise ValueError("calibration partition assay class must be a non-empty string")
    raw_groups = _exact_object(raw["groups"], set(GROUP_NAMESPACES), "calibration member groups")
    groups = {namespace: _group_values(raw_groups[namespace], namespace) for namespace in GROUP_NAMESPACES}
    return PartitionMember(
        key,
        cast(CalibrationRole, role),
        cast(EvidenceProvenance, provenance),
        assay_class,
        MappingProxyType(groups),
    )


def _group_values(value: object, namespace: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not value or any(not isinstance(item, str) or not item for item in value):
        raise ValueError(f"calibration group {namespace} must be a non-empty string list")
    values = tuple(value)
    if values != tuple(sorted(values)) or len(values) != len(set(values)):
        raise ValueError(f"calibration group {namespace} values must be unique and increasing")
    return values


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value
