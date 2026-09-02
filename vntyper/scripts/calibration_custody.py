"""Append-only local safeguards for one-use locked calibration evidence.

These records make accidental local reuse fail closed. They complement, but do
not replace, independent custodian controls and external evidence attestation.
"""

from __future__ import annotations

import fcntl
import hashlib
import os
from collections.abc import Callable, Iterator
from contextlib import contextmanager, suppress
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.calibration_secure_io import read_regular_path
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


@dataclass
class CandidateClaim:
    """Operation-lifetime ownership of one exact profile/evidence pair."""

    root: Path
    precommit: Precommit
    descriptor: int
    terminal_path: Path | None = None

    def retire(self, reason: str) -> Path:
        """Install the owner-authorized failure terminal and release the pair."""
        if self.terminal_path is not None:
            if self.terminal_path.parent.name == "retired":
                return self.terminal_path
            raise ValueError("completed calibration profile/evidence pair cannot be retired")
        if self.descriptor < 0:
            raise ValueError("calibration candidate claim is no longer active")
        try:
            self.terminal_path = _write_retirement(
                self.root,
                self.precommit.profile_sha256,
                self.precommit.evidence_sha256,
                reason,
            )
            return self.terminal_path
        finally:
            self._release()

    def complete(self) -> Path:
        """Install the successful terminal and release the pair."""
        if self.terminal_path is not None:
            if self.terminal_path.parent.name == "completed":
                return self.terminal_path
            raise ValueError("retired calibration profile/evidence pair cannot be completed")
        if self.descriptor < 0:
            raise ValueError("calibration candidate claim is no longer active")
        path = _completion_path(
            self.root,
            self.precommit.profile_sha256,
            self.precommit.evidence_sha256,
        )
        payload = {
            "schema_version": "calibration-completion-v1",
            "profile_sha256": self.precommit.profile_sha256,
            "protocol_sha256": self.precommit.protocol_sha256,
            "evidence_sha256": self.precommit.evidence_sha256,
        }
        _exclusive_write(path, canonical_json_bytes(payload), "calibration candidate is already completed")
        self.terminal_path = path
        self._release()
        return path

    def close(self) -> None:
        """Release an unterminated claim after a failed terminal write."""
        self._release()

    def _release(self) -> None:
        if self.descriptor >= 0:
            descriptor = self.descriptor
            self.descriptor = -1
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            except OSError:
                if self.terminal_path is None:
                    raise
            finally:
                try:
                    os.close(descriptor)
                except OSError:
                    if self.terminal_path is None:
                        raise


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


def claim_candidate(
    custody_dir: Path, profile_sha256: str, protocol_sha256: str, evidence_sha256: str
) -> CandidateClaim:
    """Own the exact pair from atomic precommit through its terminal state."""
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    protocol = _digest(protocol_sha256, "protocol")
    evidence = _digest(evidence_sha256, "evidence")
    descriptor = _acquire_pair_lock(root, profile, evidence)
    try:
        if _retirement_path(root, profile, evidence).exists():
            raise ValueError("calibration profile/evidence pair is retired and cannot be claimed")
        if _completion_path(root, profile, evidence).exists():
            raise ValueError("calibration profile/evidence pair is completed under the one-use policy")
        precommit = write_precommit(root, profile, protocol, evidence)
        return CandidateClaim(root, precommit, descriptor)
    except BaseException:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)
        raise


def open_locked_payload(
    payload_source: Path | Callable[[], bytes], precommit: Precommit, custody_dir: Path
) -> OpenedLockedPayload:
    """Consume first, then open and verify one precommitted locked payload."""
    if not isinstance(payload_source, Path) and not callable(payload_source):
        raise ValueError("locked payload access requires a Path or secure reader")
    if not isinstance(precommit, Precommit):
        raise ValueError("locked payload access requires a Precommit")
    _verify_precommit(precommit)
    receipt = record_consumption(
        custody_dir,
        precommit.profile_sha256,
        precommit.evidence_sha256,
        precommit.protocol_sha256,
    )
    try:
        payload = read_regular_path(payload_source) if isinstance(payload_source, Path) else payload_source()
    except OSError as error:
        raise ValueError("locked evidence payload is unreadable") from error
    if not isinstance(payload, bytes):
        raise ValueError("locked evidence payload reader must return bytes")
    observed = hashlib.sha256(payload).hexdigest()
    if observed != precommit.evidence_sha256:
        raise ValueError("locked evidence payload hash differs from its precommit")
    return OpenedLockedPayload(payload, receipt)


def retire_candidate(custody_dir: Path, profile_sha256: str, evidence_sha256: str, reason: str) -> Path:
    """Append an immutable retirement record after a failed locked evaluation."""
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    evidence = _digest(evidence_sha256, "evidence")
    with _pair_lock(root, profile, evidence):
        if _completion_path(root, profile, evidence).exists():
            raise ValueError("completed calibration profile/evidence pair cannot be retired")
        return _write_retirement(root, profile, evidence, reason)


def ensure_candidate_retired(custody_dir: Path, profile_sha256: str, evidence_sha256: str, reason: str) -> Path:
    """Install a retirement record, or preserve an existing exact-pair record."""
    root = _custody_root(custody_dir)
    profile = _digest(profile_sha256, "profile")
    evidence = _digest(evidence_sha256, "evidence")
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("calibration candidate retirement requires a non-empty reason")
    path = _retirement_path(root, profile, evidence)
    with _pair_lock(root, profile, evidence):
        if _completion_path(root, profile, evidence).exists():
            raise ValueError("completed calibration profile/evidence pair cannot be retired")
        if path.exists():
            observed = load_strict_json_object(read_regular_path(path))
            observed_reason = observed.get("reason")
            if (
                set(observed) != {"schema_version", "profile_sha256", "evidence_sha256", "reason"}
                or observed.get("schema_version") != "calibration-retirement-v1"
                or observed.get("profile_sha256") != profile
                or observed.get("evidence_sha256") != evidence
                or not isinstance(observed_reason, str)
                or not observed_reason.strip()
            ):
                raise ValueError("calibration candidate retirement record differs from its exact pair")
            return path
        return _write_retirement(root, profile, evidence, reason)


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
    if _retirement_path(root, profile, evidence).exists():
        raise ValueError("calibration profile/evidence pair is retired and cannot be retried")
    if _completion_path(root, profile, evidence).exists():
        raise ValueError("calibration profile/evidence pair is completed and cannot be retried")


def _verify_precommit(precommit: Precommit) -> None:
    expected = {
        "schema_version": "calibration-precommit-v1",
        "profile_sha256": _digest(precommit.profile_sha256, "profile"),
        "protocol_sha256": _digest(precommit.protocol_sha256, "protocol"),
        "evidence_sha256": _digest(precommit.evidence_sha256, "evidence"),
    }
    try:
        observed = load_strict_json_object(read_regular_path(precommit.path))
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


@contextmanager
def _pair_lock(root: Path, profile: str, evidence: str) -> Iterator[None]:
    descriptor = _acquire_pair_lock(root, profile, evidence)
    try:
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _acquire_pair_lock(root: Path, profile: str, evidence: str) -> int:
    path = root / "locks" / f"{profile}.{evidence}.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0), 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _retirement_path(root: Path, profile: str, evidence: str) -> Path:
    return root / "retired" / f"{profile}.{evidence}.json"


def _completion_path(root: Path, profile: str, evidence: str) -> Path:
    return root / "completed" / f"{profile}.{evidence}.json"


def _write_retirement(root: Path, profile: str, evidence: str, reason: str) -> Path:
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("calibration candidate retirement requires a non-empty reason")
    path = _retirement_path(root, profile, evidence)
    payload = {
        "schema_version": "calibration-retirement-v1",
        "profile_sha256": profile,
        "evidence_sha256": evidence,
        "reason": reason,
    }
    _exclusive_write(path, canonical_json_bytes(payload), "calibration candidate is already retired for this evidence")
    return path


def _digest(value: str, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"calibration {label} hash must be lowercase SHA-256")
    return value
