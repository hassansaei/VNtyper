"""Atomic producer for normalized local calibration intake audit bundles."""

from __future__ import annotations

import hashlib
import logging
import os
import re
import stat
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_identity import (
    ArtifactFingerprint,
    IdentityAudit,
    build_partition_members,
    identity_audit_document,
    resolve_identities,
)
from vntyper.scripts.calibration_intake_contract import IntakeDeclaration, decode_intake, encode_intake
from vntyper.scripts.calibration_manifest import GROUP_NAMESPACES, PartitionMember, decode_partition_manifest
from vntyper.scripts.calibration_read_io import fingerprint_input_artifact
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)

_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_PAYLOAD_NAMES = ("dedup_audit.json", "normalized.json", "partitions.json")
_BUNDLE_NAMES = frozenset((*_PAYLOAD_NAMES, "provenance.json"))


@dataclass(frozen=True)
class PinnedCramReference:
    """One explicit local reference binding for a declared CRAM assembly."""

    path: Path
    sha256: str


@dataclass(frozen=True)
class PreparedIntakeBundle:
    """Installed local intake artifacts and the validated objects they encode."""

    declaration: IntakeDeclaration
    audit: IdentityAudit
    partition_members: tuple[PartitionMember, ...]
    output: Path
    provenance_sha256: str


@dataclass(frozen=True)
class _SourceSnapshot:
    descriptor: int
    identity: tuple[int, int, int, int, int]
    raw: bytes


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _file_identity(metadata: os.stat_result) -> tuple[int, int, int, int, int]:
    return metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns


def _read_descriptor(descriptor: int) -> bytes:
    chunks: list[bytes] = []
    while True:
        chunk = os.read(descriptor, 1024 * 1024)
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)


def _open_source(path: Path) -> _SourceSnapshot:
    if not isinstance(path, Path) or not _NOFOLLOW:
        _fail("calibration intake declaration requires a Path and O_NOFOLLOW support")
    local_path = path if path.is_absolute() else Path.cwd() / path
    try:
        descriptor = os.open(local_path, os.O_RDONLY | os.O_NONBLOCK | _CLOEXEC | _NOFOLLOW)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            _fail("calibration intake declaration must be a regular non-symlink file")
        raw = _read_descriptor(descriptor)
        after = os.fstat(descriptor)
        if _file_identity(before) != _file_identity(after):
            raise RuntimeError("calibration intake declaration changed while it was read")
        return _SourceSnapshot(descriptor, _file_identity(after), raw)
    except OSError as error:
        if "descriptor" in locals():
            os.close(descriptor)
        raise ValueError("calibration intake declaration is unreadable or is a symlink") from error
    except BaseException:
        if "descriptor" in locals():
            os.close(descriptor)
        raise


def _assert_source_unchanged(path: Path, snapshot: _SourceSnapshot) -> None:
    local_path = path if path.is_absolute() else Path.cwd() / path
    try:
        descriptor_metadata = os.fstat(snapshot.descriptor)
        path_metadata = os.stat(local_path, follow_symlinks=False)
    except OSError:
        raise RuntimeError("calibration intake declaration changed during bundle preparation") from None
    if (
        not stat.S_ISREG(path_metadata.st_mode)
        or _file_identity(descriptor_metadata) != snapshot.identity
        or _file_identity(path_metadata) != snapshot.identity
    ):
        raise RuntimeError("calibration intake declaration changed during bundle preparation")


def _validate_priority(priority: tuple[str, ...]) -> tuple[str, ...]:
    if (
        not isinstance(priority, tuple)
        or not priority
        or any(not isinstance(value, str) or not value or value.strip() != value for value in priority)
        or len(priority) != len(set(priority))
    ):
        _fail("calibration intake preprocessing priority must contain unique non-empty strings")
    return priority


def _validate_local_paths(declaration: IntakeDeclaration) -> None:
    if not declaration.artifacts:
        _fail("calibration intake producer has no_local_artifacts")
    for artifact in declaration.artifacts:
        if not Path(artifact.path).is_absolute() or (
            artifact.mate_path is not None and not Path(artifact.mate_path).is_absolute()
        ):
            _fail("calibration intake artifact and mate paths must be absolute")


def _validate_references(
    declaration: IntakeDeclaration, values: Mapping[str, PinnedCramReference]
) -> Mapping[str, PinnedCramReference]:
    required = {artifact.assembly for artifact in declaration.artifacts if artifact.format == "CRAM"}
    if not isinstance(values, Mapping) or set(values) != required:
        _fail("calibration intake CRAM reference assemblies must match the exact required set")
    references: dict[str, PinnedCramReference] = {}
    for assembly in sorted(required):
        reference = values[assembly]
        if (
            not isinstance(reference, PinnedCramReference)
            or not isinstance(reference.path, Path)
            or not reference.path.is_absolute()
            or not isinstance(reference.sha256, str)
            or _SHA256.fullmatch(reference.sha256) is None
        ):
            _fail("calibration intake CRAM reference binding is invalid")
        references[assembly] = reference
    return MappingProxyType(references)


def _partition_document(members: tuple[PartitionMember, ...]) -> dict[str, object]:
    return {
        "schema_version": "calibration-partitions-v1",
        "members": [
            {
                "key": member.key,
                "role": member.role,
                "provenance": member.provenance,
                "assay_class": member.assay_class,
                "groups": {name: list(member.groups[name]) for name in GROUP_NAMESPACES},
            }
            for member in members
        ],
    }


def _write_private_json(path: Path, value: object) -> bytes:
    raw = canonical_json_bytes(value)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | _CLOEXEC | _NOFOLLOW, 0o600)
    try:
        os.fchmod(descriptor, 0o600)
        view = memoryview(raw)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("calibration intake artifact write made no progress")
            view = view[written:]
    finally:
        os.close(descriptor)
    return raw


def _file_row(name: str, raw: bytes) -> dict[str, object]:
    return {"path": name, "size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def _write_bundle(
    staging: Path,
    declaration: IntakeDeclaration,
    audit: IdentityAudit,
    partition_document: dict[str, object],
) -> str:
    payloads = {
        "normalized.json": canonical_json_bytes(encode_intake(declaration)),
        "dedup_audit.json": canonical_json_bytes(identity_audit_document(declaration, audit)),
        "partitions.json": canonical_json_bytes(partition_document),
    }
    for name in _PAYLOAD_NAMES:
        _write_private_json(staging / name, load_strict_json_object(payloads[name]))
    provenance: dict[str, object] = {
        "schema_version": "calibration-intake-provenance-v1",
        "intake_sha256": declaration.sha256,
        "identity_audit_sha256": audit.sha256,
        "partition_sha256": canonical_sha256(partition_document),
        "files": [_file_row(name, payloads[name]) for name in _PAYLOAD_NAMES],
    }
    provenance_raw = _write_private_json(staging / "provenance.json", provenance)
    _verify_bundle(staging, provenance)
    return hashlib.sha256(provenance_raw).hexdigest()


def _verify_bundle(staging: Path, provenance: dict[str, object]) -> None:
    if {path.name for path in staging.iterdir()} != _BUNDLE_NAMES:
        _fail("calibration intake bundle file inventory differs from its exact manifest")
    rows = provenance["files"]
    if not isinstance(rows, list) or [row.get("path") for row in rows if isinstance(row, dict)] != list(_PAYLOAD_NAMES):
        _fail("calibration intake provenance file manifest is invalid")
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"path", "size_bytes", "sha256"}:
            _fail("calibration intake provenance file manifest is invalid")
        path = staging / str(row["path"])
        metadata = path.lstat()
        if not stat.S_ISREG(metadata.st_mode) or stat.S_IMODE(metadata.st_mode) != 0o600:
            _fail("calibration intake bundle artifacts must be private regular files")
        raw = path.read_bytes()
        if row != _file_row(path.name, raw):
            _fail("calibration intake bundle artifact differs from its manifest")
    provenance_path = staging / "provenance.json"
    if not stat.S_ISREG(provenance_path.lstat().st_mode) or stat.S_IMODE(provenance_path.stat().st_mode) != 0o600:
        _fail("calibration intake bundle provenance must be a private regular file")
    if load_strict_json_object(provenance_path.read_bytes()) != provenance:
        _fail("calibration intake bundle provenance differs after writing")


def prepare_intake_bundle(
    declaration_path: Path,
    output: Path,
    *,
    preprocessing_priority: tuple[str, ...],
    cram_references: Mapping[str, PinnedCramReference],
    temporary_parent: Path | None = None,
) -> PreparedIntakeBundle:
    """Validate local intake, audit input identities, and atomically install its bundle.

    Artifact and mate paths must be absolute. This producer preserves the exact
    canonical declaration and its digest rather than rewriting paths. CRAM references
    are supplied separately and keyed by the exact declared assembly.

    Args:
        declaration_path: Strict calibration-intake-v1 JSON document.
        output: New local bundle directory.
        preprocessing_priority: Frozen unique preprocessing preference order.
        cram_references: Exact CRAM assembly-to-local-reference bindings.
        temporary_parent: Optional existing directory for fingerprint sort files.

    Returns:
        Validated in-memory audit result bound to the installed artifact directory.

    Raises:
        ValueError: If arguments, intake, reference bindings, identities, or output differ.
        RuntimeError: If source evidence changes or atomic publication is unavailable.
    """
    if not isinstance(output, Path):
        _fail("calibration intake output must be a Path")
    priority = _validate_priority(preprocessing_priority)
    if temporary_parent is not None and (
        not isinstance(temporary_parent, Path) or not temporary_parent.is_dir() or temporary_parent.is_symlink()
    ):
        _fail("calibration intake temporary parent must be an existing non-symlink directory")
    prepared: list[PreparedIntakeBundle] = []

    def produce(staging: Path) -> bool:
        snapshot = _open_source(declaration_path)
        try:
            declaration = decode_intake(load_strict_json_object(snapshot.raw))
            _validate_local_paths(declaration)
            references = _validate_references(declaration, cram_references)
            fingerprints: dict[str, ArtifactFingerprint] = {}
            for artifact in declaration.artifacts:
                reference = references.get(artifact.assembly) if artifact.format == "CRAM" else None
                fingerprints[artifact.key] = fingerprint_input_artifact(
                    artifact,
                    reference_path=None if reference is None else reference.path,
                    reference_sha256=None if reference is None else reference.sha256,
                    temporary_parent=temporary_parent,
                )
            audit = resolve_identities(declaration, fingerprints, preprocessing_priority=priority)
            members = build_partition_members(declaration, audit)
            partition_document = _partition_document(members)
            partition = decode_partition_manifest(partition_document)
            provenance_sha256 = _write_bundle(staging, declaration, audit, partition_document)
            _assert_source_unchanged(declaration_path, snapshot)
            prepared.append(PreparedIntakeBundle(declaration, audit, partition.members, output, provenance_sha256))
        finally:
            os.close(snapshot.descriptor)
        return True

    atomic_output(output, produce)
    if len(prepared) != 1:
        raise RuntimeError("calibration intake producer completed without a result")
    return prepared[0]
