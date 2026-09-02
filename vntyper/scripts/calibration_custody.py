"""Append-only local safeguards for one-use locked calibration evidence.

These records make accidental local reuse fail closed. They complement, but do
not replace, independent custodian controls and external evidence attestation.
"""

from __future__ import annotations

import hashlib
import os
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object


@dataclass(frozen=True)
class Precommit:
    """Durable candidate/evidence binding written before evidence access."""

    path: Path
    profile_sha256: str
    protocol_sha256: str
    evidence_sha256: str


@dataclass(frozen=True)
class ConsumptionReceipt:
    """Append-only proof that one evidence payload has been consumed locally."""

    path: Path
    profile_sha256: str
    protocol_sha256: str
    evidence_sha256: str


@dataclass(frozen=True)
class OpenedLockedPayload:
    """Locked bytes paired with the receipt created before they were read."""

    payload: bytes
    receipt: ConsumptionReceipt


def write_precommit(
    custody_dir: Path,
    profile_sha256: str,
    protocol_sha256: str,
    evidence_sha256: str,
) -> Precommit:
    """Durably bind a fixed candidate to locked evidence before opening it."""
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    protocol = _digest(protocol_sha256, "protocol")
    evidence = _digest(evidence_sha256, "evidence")
    path = root / "precommits" / f"{profile}.{evidence}.json"
    payload = {
        "schema_version": "calibration-precommit-v1",
        "profile_sha256": profile,
        "protocol_sha256": protocol,
        "evidence_sha256": evidence,
    }
    _exclusive_write(
        path,
        canonical_json_bytes(payload),
        "calibration precommit already exists; the one-use evidence may already be consumed",
    )
    return Precommit(path, profile, protocol, evidence)


def record_consumption(
    custody_dir: Path,
    profile_sha256: str,
    evidence_sha256: str,
    protocol_sha256: str,
) -> ConsumptionReceipt:
    """Consume an evidence hash exactly once, before its payload is exposed."""
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    evidence = _digest(evidence_sha256, "evidence")
    protocol = _digest(protocol_sha256, "protocol")
    path = root / "consumed" / f"{evidence}.json"
    payload = {
        "schema_version": "calibration-consumption-v1",
        "profile_sha256": profile,
        "protocol_sha256": protocol,
        "evidence_sha256": evidence,
    }
    _exclusive_write(
        path, canonical_json_bytes(payload), "locked evidence is already consumed under the one-use policy"
    )
    return ConsumptionReceipt(path, profile, protocol, evidence)


def open_locked_payload(payload_path: Path, precommit: Precommit, custody_dir: Path) -> OpenedLockedPayload:
    """Consume first, then open and verify one precommitted locked payload."""
    if not isinstance(payload_path, Path) or not isinstance(precommit, Precommit):
        raise ValueError("locked payload access requires a Path and Precommit")
    _verify_precommit(precommit)
    receipt = record_consumption(
        custody_dir,
        precommit.profile_sha256,
        precommit.evidence_sha256,
        precommit.protocol_sha256,
    )
    try:
        payload = payload_path.read_bytes()
    except OSError as error:
        raise ValueError(f"locked evidence payload is unreadable: {payload_path}") from error
    observed = hashlib.sha256(payload).hexdigest()
    if observed != precommit.evidence_sha256:
        raise ValueError("locked evidence payload hash differs from its precommit")
    return OpenedLockedPayload(payload, receipt)


def retire_candidate(custody_dir: Path, profile_sha256: str, evidence_sha256: str, reason: str) -> Path:
    """Append an immutable retirement record after a failed locked evaluation."""
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    evidence = _digest(evidence_sha256, "evidence")
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("calibration candidate retirement requires a non-empty reason")
    path = root / "retired" / f"{profile}.{evidence}.json"
    payload = {
        "schema_version": "calibration-retirement-v1",
        "profile_sha256": profile,
        "evidence_sha256": evidence,
        "reason": reason,
    }
    _exclusive_write(path, canonical_json_bytes(payload), "calibration candidate is already retired for this evidence")
    return path


def require_candidate_active(custody_dir: Path, profile_sha256: str, evidence_sha256: str) -> None:
    """Refuse an exact profile/evidence pair with a durable retirement record.

    Args:
        custody_dir: Local append-only safeguard root.
        profile_sha256: Exact generated profile digest.
        evidence_sha256: Exact role-specific evidence digest.

    Raises:
        ValueError: If hashes are malformed or the pair is already retired.
    """
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    evidence = _digest(evidence_sha256, "evidence")
    if (root / "retired" / f"{profile}.{evidence}.json").exists():
        raise ValueError("calibration profile/evidence pair is retired and cannot be retried")


def _verify_precommit(precommit: Precommit) -> None:
    expected = {
        "schema_version": "calibration-precommit-v1",
        "profile_sha256": _digest(precommit.profile_sha256, "profile"),
        "protocol_sha256": _digest(precommit.protocol_sha256, "protocol"),
        "evidence_sha256": _digest(precommit.evidence_sha256, "evidence"),
    }
    try:
        observed = load_strict_json_object(precommit.path.read_bytes())
    except (OSError, ValueError) as error:
        raise ValueError("calibration precommit is missing, unreadable, or invalid") from error
    if observed != expected:
        raise ValueError("calibration precommit content differs from its bound hashes")


def _exclusive_write(path: Path, payload: bytes, conflict_message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    try:
        descriptor = os.open(path, flags, 0o600)
    except FileExistsError as error:
        raise ValueError(conflict_message) from error
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        directory_descriptor = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        with suppress(OSError):
            path.unlink()
        raise


def _custody_root(value: Path) -> Path:
    if not isinstance(value, Path):
        raise ValueError("calibration custody directory must be a Path")
    value.mkdir(parents=True, exist_ok=True)
    return value


def _digest(value: str, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"calibration {label} hash must be lowercase SHA-256")
    return value
