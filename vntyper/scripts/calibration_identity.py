"""Pure duplicate-evidence and biological-linkage audit before role extraction."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_intake_contract import (
    InputArtifact,
    IntakeDeclaration,
    TruthRecord,
    decode_intake,
    encode_intake,
)
from vntyper.scripts.calibration_manifest import PartitionMember, connected_leakage_groups, decode_partition_manifest
from vntyper.scripts.calibration_read_fingerprints import LogicalReadFingerprint
from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class ArtifactFingerprint:
    """Physical file hashes and logical read identities for one declared artifact."""

    artifact_key: str
    byte_sha256: str
    mate_byte_sha256: str | None
    logical: LogicalReadFingerprint


@dataclass(frozen=True)
class IdentityAudit:
    """Immutable grouping and representative decisions without opening truth files."""

    intake_sha256: str
    preprocessing_priority: tuple[str, ...]
    fingerprints: Mapping[str, ArtifactFingerprint]
    specimen_groups: Mapping[str, str]
    specimen_linkage: Mapping[str, Mapping[str, tuple[str, ...]]]
    execution_representatives: Mapping[str, str]
    primary_artifact_by_specimen: Mapping[str, str]
    primary_artifact_by_group: Mapping[str, str]
    quarantined_specimens: Mapping[str, tuple[str, ...]]
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _validated_intake(declaration: IntakeDeclaration) -> IntakeDeclaration:
    decoded = decode_intake(encode_intake(declaration))
    if decoded.sha256 != declaration.sha256:
        _fail("identity audit intake digest differs from its canonical content")
    return decoded


def _digest(value: object) -> None:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail("identity audit fingerprints require lowercase SHA256 digests")


def _validate_fingerprints(
    artifacts: tuple[InputArtifact, ...], fingerprints: Mapping[str, ArtifactFingerprint]
) -> dict[str, ArtifactFingerprint]:
    if not isinstance(fingerprints, Mapping) or set(fingerprints) != {artifact.key for artifact in artifacts}:
        _fail("identity audit fingerprints must match the exact artifact roster")
    validated: dict[str, ArtifactFingerprint] = {}
    for artifact in artifacts:
        fingerprint = fingerprints[artifact.key]
        if not isinstance(fingerprint, ArtifactFingerprint) or fingerprint.artifact_key != artifact.key:
            _fail("identity audit fingerprint key differs from its artifact")
        _digest(fingerprint.byte_sha256)
        if artifact.format == "FASTQ_PAIR":
            _digest(fingerprint.mate_byte_sha256)
        elif fingerprint.mate_byte_sha256 is not None:
            _fail("single-file artifact fingerprint forbids a mate file hash")
        if artifact.expected_sha256 is not None and artifact.expected_sha256 != fingerprint.byte_sha256:
            _fail("identity audit observed input differs from its expected file hash")
        logical = fingerprint.logical
        if not isinstance(logical, LogicalReadFingerprint):
            _fail("identity audit requires typed logical read fingerprints")
        for value in (logical.alignment_sha256, logical.named_sequence_sha256, logical.unnamed_sequence_sha256):
            _digest(value)
        if (
            isinstance(logical.primary_record_count, bool)
            or not isinstance(logical.primary_record_count, int)
            or logical.primary_record_count <= 0
        ):
            _fail("identity audit requires a positive primary record count")
        if (
            not isinstance(logical.reasons, tuple)
            or any(
                reason not in {"hard_clipped_sequence", "missing_sequence", "missing_qualities"}
                for reason in logical.reasons
            )
            or logical.reasons != tuple(sorted(set(logical.reasons)))
        ):
            _fail("identity audit logical reconstruction reasons are invalid")
        if not isinstance(logical.sequence_identity_reliable, bool) or logical.sequence_identity_reliable != (
            not logical.reasons
        ):
            _fail("identity audit logical reliability disagrees with reconstruction defects")
        validated[artifact.key] = fingerprint
    return validated


def _truth_conflicts(rows: Sequence[TruthRecord]) -> set[str]:
    reasons: set[str] = set()
    if any(row.status == "disputed" for row in rows):
        reasons.add("disputed_truth")
    confirmed = [row for row in rows if row.status == "confirmed"]
    if len({row.genotype for row in confirmed if row.genotype != "unknown"}) > 1:
        reasons.add("genotype_truth_conflict")
    if len({row.variants for row in confirmed if row.variants}) > 1:
        reasons.add("variant_truth_conflict")
    lengths = {
        (
            length.unit,
            length.repeat_unit_bp,
            length.boundary_definition,
            length.conversion_id,
            tuple(sorted((length.allele_1, length.allele_2))),
        )
        for row in confirmed
        if (length := row.length) is not None
        and length.measurement == "exact"
        and length.allele_1 is not None
        and length.allele_2 is not None
    }
    if len(lengths) > 1:
        reasons.add("length_truth_conflict")
    return reasons


def _policy_key(artifact: InputArtifact) -> tuple[str, ...]:
    return (artifact.format, artifact.assembly, artifact.assay_class, artifact.input_scope, artifact.preprocessing_id)


def _truth_quarantine(declaration: IntakeDeclaration) -> dict[str, set[str]]:
    # Explicit individual linkage joins labels; family linkage must not join truth.
    identities = {
        specimen.key: ("individual", specimen.individual_key)
        if specimen.individual_key is not None
        else ("specimen", specimen.key)
        for specimen in declaration.specimens
    }
    rows: dict[tuple[str, str], list[TruthRecord]] = {}
    for truth in declaration.truth:
        rows.setdefault(identities[truth.specimen_key], []).append(truth)
    conflicts = {identity: _truth_conflicts(values) for identity, values in rows.items()}
    return {key: set(conflicts.get(identity, set())) for key, identity in identities.items()}


def _logical_links(
    declaration: IntakeDeclaration,
    fingerprints: Mapping[str, ArtifactFingerprint],
    groups: dict[str, dict[str, set[str]]],
    quarantine: dict[str, set[str]],
) -> None:
    specimens = {row.key: row for row in declaration.specimens}
    artifacts = {row.key: row for row in declaration.artifacts}
    for kind in ("named_sequence_sha256", "unnamed_sequence_sha256"):
        matches: dict[str, list[str]] = {}
        for key, fingerprint in fingerprints.items():
            matches.setdefault(getattr(fingerprint.logical, kind), []).append(key)
        for digest, keys in matches.items():
            if len({fingerprints[key].logical.primary_record_count for key in keys}) > 1:
                _fail("identical logical read digest has inconsistent occurrence counts")
            specimen_keys = {artifacts[key].specimen_key for key in keys}
            for key in specimen_keys:
                groups[key]["replicate-rerun"].add(f"logical:{kind}:{digest}")
            if len(specimen_keys) < 2:
                continue
            individuals = {specimens[key].individual_key for key in specimen_keys}
            same_declared_individual = len(individuals) == 1 and None not in individuals
            if same_declared_individual:
                continue
            named_count = len({fingerprints[key].logical.named_sequence_sha256 for key in keys})
            reason = "duplicate_input_identity_requires_adjudication"
            if kind == "unnamed_sequence_sha256":
                if named_count == 1:
                    continue
                reason = "sequence_only_collision_requires_adjudication"
            for key in specimen_keys:
                quarantine[key].add(reason)


def _audit_document(audit: IdentityAudit) -> dict[str, object]:
    return {
        "schema_version": "calibration-identity-audit-v1",
        "intake_sha256": audit.intake_sha256,
        "preprocessing_priority": list(audit.preprocessing_priority),
        "fingerprints": {
            key: {
                "byte_sha256": item.byte_sha256,
                "mate_byte_sha256": item.mate_byte_sha256,
                "alignment_sha256": item.logical.alignment_sha256,
                "named_sequence_sha256": item.logical.named_sequence_sha256,
                "unnamed_sequence_sha256": item.logical.unnamed_sequence_sha256,
                "primary_record_count": item.logical.primary_record_count,
                "sequence_identity_reliable": item.logical.sequence_identity_reliable,
                "reasons": list(item.logical.reasons),
            }
            for key, item in audit.fingerprints.items()
        },
        "specimen_groups": dict(audit.specimen_groups),
        "specimen_linkage": {
            key: {name: list(values) for name, values in groups.items()}
            for key, groups in audit.specimen_linkage.items()
        },
        "execution_representatives": dict(audit.execution_representatives),
        "primary_artifact_by_specimen": dict(audit.primary_artifact_by_specimen),
        "primary_artifact_by_group": dict(audit.primary_artifact_by_group),
        "quarantined_specimens": {key: list(values) for key, values in audit.quarantined_specimens.items()},
    }


def resolve_identities(
    declaration: IntakeDeclaration,
    fingerprints: Mapping[str, ArtifactFingerprint],
    *,
    preprocessing_priority: Sequence[str],
) -> IdentityAudit:
    """Audit evidence duplicates and transitive biological linkage before extraction.

    Args:
        declaration: Strict normalized intake with explicit biological identities.
        fingerprints: Complete physical/logical identities supplied by input readers.
        preprocessing_priority: Predeclared preferred processing order, then artifact key.

    Returns:
        Immutable local audit. Read matches share an evidence group but never
        merge specimen identities. Truth/identity conflicts are quarantined for
        explicit adjudication; consumers must not fit those unresolved records.

    Raises:
        ValueError: If identities are inconsistent, fingerprints are incomplete,
            linkage crosses roles, or confirmatory identity is not assessable.
    """
    declaration = _validated_intake(declaration)
    if (
        not isinstance(preprocessing_priority, (tuple, list))
        or not preprocessing_priority
        or any(not isinstance(value, str) or not value or value != value.strip() for value in preprocessing_priority)
        or len(set(preprocessing_priority)) != len(preprocessing_priority)
    ):
        _fail("identity audit preprocessing priority must contain unique non-empty values")
    priority = tuple(preprocessing_priority)
    rank = {value: index for index, value in enumerate(priority)}
    if any(artifact.preprocessing_id not in rank for artifact in declaration.artifacts):
        _fail("identity audit artifact preprocessing is absent from the frozen priority")
    fps = _validate_fingerprints(declaration.artifacts, fingerprints)
    groups = {
        assignment.specimen_key: {name: set(values) for name, values in assignment.groups.items()}
        for assignment in declaration.assignments
    }
    quarantine = _truth_quarantine(declaration)
    for specimen in declaration.specimens:
        biological = groups[specimen.key]["individual-family"]
        biological.add("specimen:" + canonical_sha256(specimen.key))
        for namespace, value in (("individual", specimen.individual_key), ("family", specimen.family_key)):
            if value is not None:
                biological.add(namespace + ":" + canonical_sha256(value))
        if specimen.identity_status == "unresolved":
            quarantine[specimen.key].add("unresolved_biological_identity")
    by_specimen: dict[str, list[InputArtifact]] = {}
    physical: dict[tuple[object, ...], list[InputArtifact]] = {}
    same_bytes: dict[tuple[str, str | None], LogicalReadFingerprint] = {}
    for artifact in declaration.artifacts:
        fingerprint = fps[artifact.key]
        by_specimen.setdefault(artifact.specimen_key, []).append(artifact)
        groups[artifact.specimen_key]["replicate-rerun"].add(
            "artifact-replicate:" + canonical_sha256(artifact.replicate_group)
        )
        if not fingerprint.logical.sequence_identity_reliable:
            quarantine[artifact.specimen_key].add("unreliable_sequence_identity")
        byte_key = (fingerprint.byte_sha256, fingerprint.mate_byte_sha256)
        previous = same_bytes.setdefault(byte_key, fingerprint.logical)
        if previous != fingerprint.logical:
            _fail("byte-identical artifacts have inconsistent logical fingerprints")
        physical.setdefault((*byte_key, *_policy_key(artifact)), []).append(artifact)
    _logical_links(declaration, fps, groups, quarantine)
    identity_manifest = decode_partition_manifest(
        {
            "schema_version": "calibration-partitions-v1",
            "members": [
                {
                    "key": assignment.specimen_key,
                    "role": assignment.role,
                    "provenance": assignment.provenance,
                    "assay_class": "identity-audit-only",
                    "groups": {name: sorted(values) for name, values in groups[assignment.specimen_key].items()},
                }
                for assignment in declaration.assignments
            ],
        }
    )
    connected = connected_leakage_groups(identity_manifest)
    for assignment in declaration.assignments:
        if assignment.role in {"validation", "locked-heldout"} and quarantine[assignment.specimen_key]:
            _fail("confirmatory identity audit requires adjudication of all identity and truth conflicts")
    execution = {
        artifact.key: min(item.key for item in members) for members in physical.values() for artifact in members
    }
    primary = {
        key: min(artifacts, key=lambda item: (rank[item.preprocessing_id], item.key)).key
        for key, artifacts in by_specimen.items()
    }
    artifact_by_key = {artifact.key: artifact for artifact in declaration.artifacts}
    group_primary: dict[str, str] = {}
    for key, artifact_key in sorted(
        primary.items(), key=lambda item: (rank[artifact_by_key[item[1]].preprocessing_id], item[1])
    ):
        group_primary.setdefault(connected[key], artifact_key)
    linkage = MappingProxyType(
        {
            key: MappingProxyType({name: tuple(sorted(values)) for name, values in namespaces.items()})
            for key, namespaces in groups.items()
        }
    )
    audit = IdentityAudit(
        declaration.sha256,
        priority,
        MappingProxyType(fps),
        connected,
        linkage,
        MappingProxyType(dict(sorted(execution.items()))),
        MappingProxyType(dict(sorted(primary.items()))),
        MappingProxyType(dict(sorted(group_primary.items()))),
        MappingProxyType({key: tuple(sorted(values)) for key, values in sorted(quarantine.items()) if values}),
        "",
    )
    return replace(audit, sha256=canonical_sha256(_audit_document(audit)))


def build_partition_members(declaration: IntakeDeclaration, audit: IdentityAudit) -> tuple[PartitionMember, ...]:
    """Project audited input artifacts to role members while keeping artifact keys.

    Args:
        declaration: The exact intake bound by the audit.
        audit: Complete identity audit with frozen representative priority.

    Returns:
        Input-artifact members, including quarantined development artifacts for
        audit visibility. Locked membership-only specimens have no input member;
        their assay/asset metadata remains the separate custodian's responsibility.
        Callers must use quarantine reasons when constructing eligible views.

    Raises:
        ValueError: If the audit is forged, stale, or inconsistent with the intake.
    """
    if not isinstance(audit, IdentityAudit):
        _fail("identity audit must be an IdentityAudit")
    declaration = _validated_intake(declaration)
    expected = resolve_identities(declaration, audit.fingerprints, preprocessing_priority=audit.preprocessing_priority)
    if expected != audit:
        _fail("identity audit differs from its canonical content or intake")
    assignments = {row.specimen_key: row for row in declaration.assignments}
    return tuple(
        PartitionMember(
            artifact.key,
            assignments[artifact.specimen_key].role,
            assignments[artifact.specimen_key].provenance,
            artifact.assay_class,
            expected.specimen_linkage[artifact.specimen_key],
        )
        for artifact in declaration.artifacts
    )
