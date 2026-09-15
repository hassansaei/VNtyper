"""Private target-aware one-use claims, independent of the v1 custody namespace.

Immutable started/consumed records survive errors. An interrupted claim cannot be
reopened. This local safeguard complements the external exposure ledger and real
custodian separation; an operator controlling every file can restore old copies.
"""

from __future__ import annotations

import fcntl
import hashlib
import logging
import os
import re
import stat
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_exposure import (
    ExposureReceipt,
    decode_exposure_receipt,
    exposure_receipt_document,
    require_digest,
)
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader
from vntyper.scripts.calibration_target_access import (
    TargetConfirmation,
    confirmation_document,
    preflight_target_confirmation,
)
from vntyper.scripts.calibration_target_attestation import (
    TargetLockedAttestation,
    TargetValidationAttestation,
    target_locked_attestation_document,
    target_validation_attestation_document,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)
_NAMES = re.compile(
    r"[0-9a-f]{64}\.(lock|validation-started\.json|validation-result\.json|locked-started\.json|locked-consumption\.json|locked-result\.json|retired\.json)\Z"
)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _external_root(path: Path, forbidden_roots: tuple[Path, ...]) -> Path:
    if not isinstance(path, Path) or not path.is_absolute() or path.is_symlink():
        _fail("target custody requires an absolute external nonsymlink directory")
    if not isinstance(forbidden_roots, tuple) or any(not isinstance(root, Path) for root in forbidden_roots):
        _fail("target custody forbidden roots must be a Path tuple")
    resolved = path.resolve()
    if any((parent / ".git").exists() for parent in (resolved, *resolved.parents)):
        _fail("target custody must remain outside Git repositories")
    if any(resolved == root.resolve() or root.resolve() in resolved.parents for root in forbidden_roots):
        _fail("target custody must remain outside study/input/output directories")
    return resolved


def _check_root(reader: SecureDirectoryReader) -> None:
    opened = os.fstat(reader.descriptor)
    named = os.stat(reader.path, follow_symlinks=False)
    if (
        not stat.S_ISDIR(named.st_mode)
        or stat.S_IMODE(opened.st_mode) != 0o700
        or (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino)
    ):
        _fail("target custody root changed or is not private")


def _write(reader: SecureDirectoryReader, name: str, document: object) -> None:
    _check_root(reader)
    descriptor = os.open(
        name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=reader.descriptor
    )
    try:
        os.fchmod(descriptor, 0o600)
        pending = memoryview(canonical_json_bytes(document))
        while pending:
            written = os.write(descriptor, pending)
            if written <= 0:
                _fail("target custody record write did not complete")
            pending = pending[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.fsync(reader.descriptor)
    _check_root(reader)
    reader.names = reader.names | {name}


def _read(reader: SecureDirectoryReader, name: str) -> dict[str, object] | None:
    try:
        os.stat(name, dir_fd=reader.descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return None
    reader.names = reader.names | {name}
    raw = reader.read_file(name)
    document = load_strict_json_object(raw)
    if canonical_json_bytes(document) != raw:
        _fail("target custody record is corrupt or noncanonical")
    return document


def initialize_target_custody(
    path: Path,
    *,
    exposure_ledger_id: str,
    forbidden_roots: tuple[Path, ...] = (),
) -> None:
    """Exclusively initialize private external target custody for one ledger.

    Args:
        path: New external directory; existing directories are never adopted.
        exposure_ledger_id: Predeclared operator exposure-ledger identity.
        forbidden_roots: Study, input and output roots excluded from custody.

    Raises:
        ValueError: For invalid identity/location or an occupied directory.
        OSError: For durable-write failures; incomplete roots remain unusable.
    """
    path = _external_root(path, forbidden_roots)
    identity = require_digest(exposure_ledger_id, "custody ledger identity")
    try:
        path.mkdir(mode=0o700)
    except FileExistsError as error:
        raise ValueError("target custody directory already exists") from error
    path.chmod(0o700)
    with SecureDirectoryReader.open(path, set()) as reader:
        _write(
            reader, "identity.json", {"schema_version": "calibration-target-custody-v2", "exposure_ledger_id": identity}
        )
    parent = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(parent)
    finally:
        os.close(parent)


def _exposure(confirmation: TargetConfirmation, receipt: ExposureReceipt) -> None:
    exposure_receipt_document(receipt)
    expected = {
        "target": confirmation.candidate.target,
        "role": confirmation.role,
        "study_sha256": confirmation.study.sha256,
        "partition_sha256": confirmation.study.partitions.sha256,
        "evidence_sha256": confirmation.evidence_sha256,
        "exposure_ledger_id": confirmation.study.exposure_ledger_id,
    }
    if any(getattr(receipt, name) != value for name, value in expected.items()):
        _fail("target custody exposure receipt differs from confirmation")


def _start(confirmation: TargetConfirmation, receipt: ExposureReceipt) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-custody-start-v2",
        "confirmation": confirmation_document(confirmation),
        "exposure_receipt": exposure_receipt_document(receipt),
    }


def _result(
    confirmation: TargetConfirmation,
    receipt: ExposureReceipt,
    attestation: TargetValidationAttestation | TargetLockedAttestation,
) -> dict[str, object]:
    return {
        "schema_version": "calibration-target-custody-result-v2",
        "target": confirmation.candidate.target,
        "role": confirmation.role,
        "candidate_sha256": confirmation.candidate.sha256,
        "confirmation_sha256": confirmation.sha256,
        "exposure_receipt_sha256": receipt.sha256,
        "attestation_sha256": attestation.sha256,
        "status": attestation.status,
    }


def _prior_validation(reader: SecureDirectoryReader, confirmation: TargetConfirmation, prefix: str) -> None:
    validation = confirmation.validation
    if validation is None:
        _fail("locked custody requires passed validation")
    prior = preflight_target_confirmation(
        candidate=confirmation.candidate,
        study=confirmation.study,
        role="validation",
        evidence_sha256=validation.evidence_sha256,
        run_manifest_sha256=validation.run_manifest_sha256,
    )
    started = _read(reader, prefix + ".validation-started.json")
    result = _read(reader, prefix + ".validation-result.json")
    if started is None or result is None:
        _fail("locked custody requires complete local validation records")
    receipt = decode_exposure_receipt(started.get("exposure_receipt"))
    _exposure(prior, receipt)
    if (
        started != _start(prior, receipt)
        or receipt.sha256 != validation.exposure_receipt_sha256
        or result != _result(prior, receipt, validation)
        or validation.status != "passed"
    ):
        _fail("locked custody validation records differ from passed authority")


@dataclass(frozen=True)
class TargetOpenedPayload:
    """Exact locked bytes and the receipt durably written before reading them."""

    payload: bytes
    receipt: Mapping[str, object]
    receipt_sha256: str


@dataclass
class TargetCustodyClaim:
    """Process-lifetime candidate ownership; close never deletes persistent records."""

    reader: SecureDirectoryReader
    lock_descriptor: int
    confirmation: TargetConfirmation
    exposure: ExposureReceipt
    terminal: bool = False
    consumption_sha256: str | None = None
    payload_verified: bool = False

    @property
    def prefix(self) -> str:
        """Return the candidate digest used for all direct-child record names."""
        return self.confirmation.candidate.sha256

    def _active(self) -> None:
        if self.lock_descriptor < 0 or self.terminal:
            _fail("target custody claim is closed or already terminal")
        _check_root(self.reader)
        opened = os.fstat(self.lock_descriptor)
        named = os.stat(self.prefix + ".lock", dir_fd=self.reader.descriptor, follow_symlinks=False)
        if not stat.S_ISREG(named.st_mode) or (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino):
            _fail("target custody claim lock changed")

    def open_locked_payload(self, reader: Callable[[], bytes]) -> TargetOpenedPayload:
        """Consume the authority-bound payload once, then read and verify its bytes.

        Args:
            reader: Secure reader callback; never invoked before consumption fsync.

        Returns:
            Opened bytes and their aggregate-free consumption receipt.

        Raises:
            ValueError: For invalid role, repeat access or wrong payload bytes.
            OSError: From durable writes or the secure reader; consumption persists.
        """
        self._active()
        authority = self.confirmation.authority
        if self.confirmation.role != "locked-heldout" or authority is None or not callable(reader):
            _fail("target custody payload access requires locked authority and a reader")
        if self.consumption_sha256 is not None:
            _fail("target custody payload is already consumed")
        document = {
            **confirmation_document(self.confirmation),
            "schema_version": "calibration-target-consumption-v2",
            "exposure_receipt_sha256": self.exposure.sha256,
        }
        _write(self.reader, self.prefix + ".locked-consumption.json", document)
        self.consumption_sha256 = canonical_sha256(document)
        payload = reader()
        if not isinstance(payload, bytes) or hashlib.sha256(payload).hexdigest() != authority.locked_payload_sha256:
            _fail("target custody locked payload bytes differ from authority")
        self.payload_verified = True
        return TargetOpenedPayload(payload, MappingProxyType(document), self.consumption_sha256)

    def finish(self, attestation: TargetValidationAttestation | TargetLockedAttestation) -> None:
        """Record one verified completed scientific outcome; a failure stays failed.

        Args:
            attestation: Outcome binding the current candidate, evidence and receipts.

        Raises:
            ValueError: For missing consumption, wrong role or any changed binding.
            OSError: On a durable-write failure; partial records remain fail-closed.
        """
        self._active()
        if self.confirmation.role == "validation":
            if not isinstance(attestation, TargetValidationAttestation):
                _fail("validation custody requires a validation attestation")
            target_validation_attestation_document(attestation)
            phase = "validation"
        else:
            if not isinstance(attestation, TargetLockedAttestation):
                _fail("locked custody requires a locked attestation")
            target_locked_attestation_document(attestation)
            if (
                not self.payload_verified
                or self.consumption_sha256 is None
                or attestation.custody_consumption_receipt_sha256 != self.consumption_sha256
            ):
                _fail("locked custody completion requires its exact consumption receipt")
            if (
                self.confirmation.validation is None
                or self.confirmation.authority is None
                or attestation.validation_attestation_sha256 != self.confirmation.validation.sha256
                or attestation.custodian_authority_sha256 != self.confirmation.authority.sha256
            ):
                _fail("locked custody completion authority differs")
            phase = "locked"
        expected = confirmation_document(self.confirmation)
        for name in (
            "target",
            "role",
            "candidate_sha256",
            "candidate_id",
            "study_sha256",
            "protocol_sha256",
            "partition_sha256",
            "baseline_sha256",
            "evidence_sha256",
            "run_manifest_sha256",
            "exposure_ledger_id",
        ):
            if getattr(attestation, name) != expected[name]:
                _fail(f"target custody completion {name} differs")
        if attestation.exposure_receipt_sha256 != self.exposure.sha256:
            _fail("target custody completion exposure receipt differs")
        _write(
            self.reader, self.prefix + f".{phase}-result.json", _result(self.confirmation, self.exposure, attestation)
        )
        self.terminal = True

    def retire(self, reason: str) -> None:
        """Record a stable failure reason for an unfinished claimed candidate."""
        self._active()
        if reason not in {"operational-failure-after-claim", "integrity-failure-after-claim"}:
            _fail("target custody retirement reason is unsupported")
        _write(
            self.reader,
            self.prefix + ".retired.json",
            {
                "schema_version": "calibration-target-retirement-v2",
                "target": self.confirmation.candidate.target,
                "candidate_sha256": self.prefix,
                "confirmation_sha256": self.confirmation.sha256,
                "exposure_receipt_sha256": self.exposure.sha256,
                "reason": reason,
            },
        )
        self.terminal = True

    def close(self) -> None:
        """Release descriptors without changing any started or consumed state."""
        if self.lock_descriptor >= 0:
            descriptor = self.lock_descriptor
            self.lock_descriptor = -1
            try:
                os.close(descriptor)
            finally:
                self.reader.close()

    def __enter__(self) -> TargetCustodyClaim:
        self._active()
        return self

    def __exit__(self, _type, _value, _traceback) -> None:
        try:
            if not self.terminal:
                self.retire("operational-failure-after-claim")
        finally:
            self.close()


def claim_target_confirmation(
    path: Path,
    confirmation: TargetConfirmation,
    exposure: ExposureReceipt,
    *,
    forbidden_roots: tuple[Path, ...] = (),
) -> TargetCustodyClaim:
    """Claim a candidate after exposure recording and before outcome access.

    Args:
        path: Explicitly initialized external target custody directory.
        confirmation: Pure preflight lineage, revalidated before any write.
        exposure: Receipt returned by the external ledger before outcomes are read.
        forbidden_roots: All study/input/output roots excluded from custody.

    Returns:
        Context-managed exclusive claim; interrupted states cannot be retried.

    Raises:
        ValueError: For corrupt, consumed, retired or inconsistent prior state.
        OSError: On filesystem/locking failure; existing history is never removed.
    """
    confirmation_document(confirmation)
    _exposure(confirmation, exposure)
    path = _external_root(path, forbidden_roots)
    names = set(os.listdir(path))
    if "identity.json" not in names or any(name != "identity.json" and not _NAMES.fullmatch(name) for name in names):
        _fail("target custody inventory is incomplete or unknown")
    reader = SecureDirectoryReader.open(path, names)
    lock = -1
    try:
        _check_root(reader)
        identity = _read(reader, "identity.json")
        if identity != {
            "schema_version": "calibration-target-custody-v2",
            "exposure_ledger_id": confirmation.study.exposure_ledger_id,
        }:
            _fail("target custody ledger identity differs")
        prefix = confirmation.candidate.sha256
        lock = os.open(
            prefix + ".lock",
            os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC,
            0o600,
            dir_fd=reader.descriptor,
        )
        if not stat.S_ISREG(os.fstat(lock).st_mode) or os.fstat(lock).st_size != 0:
            _fail("target custody lock is not an empty regular file")
        fcntl.flock(lock, fcntl.LOCK_EX)
        claim = TargetCustodyClaim(reader, lock, confirmation, exposure)
        claim._active()
        for suffix in ("locked-started.json", "locked-consumption.json", "locked-result.json", "retired.json"):
            if _read(reader, prefix + "." + suffix) is not None:
                _fail("target custody candidate is already claimed, consumed or terminal")
        phase = "validation" if confirmation.role == "validation" else "locked"
        if phase == "validation":
            if any(
                _read(reader, prefix + "." + suffix) is not None
                for suffix in ("validation-started.json", "validation-result.json")
            ):
                _fail("target custody validation was already claimed or completed")
        else:
            _prior_validation(reader, confirmation, prefix)
        _write(reader, prefix + f".{phase}-started.json", _start(confirmation, exposure))
        return claim
    except BaseException:
        if lock >= 0:
            os.close(lock)
        reader.close()
        raise
