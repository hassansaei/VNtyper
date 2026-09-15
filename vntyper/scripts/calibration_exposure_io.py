"""Private hash-chained exposure records with a durable fixed-size head commitment."""

from __future__ import annotations

import fcntl
import logging
import os
import secrets
import stat
from pathlib import Path
from typing import BinaryIO, NoReturn

from vntyper.scripts.calibration_exposure import (
    CONFIRMATION_ROLES,
    ExposureReceipt,
    decode_exposure_identities,
    decode_exposure_receipt,
    exposure_receipt_document,
    require_digest,
)
from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

logger = logging.getLogger(__name__)
# Administrative record bound, not a sample-count or biological length cutoff.
MAX_ENTRY_BYTES = 16 * 1024 * 1024
_HEADER_FIELDS = {"schema_version", "exposure_ledger_id", "head_sequence", "head_sha256"}
_ENTRY_FIELDS = {"schema_version", "previous_sha256", "receipt", "identities", "sha256"}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _external_path(path: Path, forbidden_roots: tuple[Path, ...]) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        _fail("calibration exposure ledger requires an absolute external path")
    if path.is_symlink():
        _fail("calibration exposure ledger must not be a symlink")
    resolved = path.resolve()
    if any((parent / ".git").exists() for parent in resolved.parents):
        _fail("calibration exposure ledger must remain outside Git repositories")
    if not isinstance(forbidden_roots, tuple) or any(not isinstance(root, Path) for root in forbidden_roots):
        _fail("calibration exposure forbidden roots must be an explicit Path tuple")
    if any(resolved == root.resolve() or root.resolve() in resolved.parents for root in forbidden_roots):
        _fail("calibration exposure ledger must remain outside study and output directories")
    return resolved


def _write_all(descriptor: int, data: bytes) -> None:
    view = memoryview(data)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            _fail("calibration exposure ledger append did not complete")
        view = view[written:]
    os.fsync(descriptor)


def _anchor(ledger_id: str) -> dict[str, object]:
    return {"schema_version": "calibration-target-exposure-ledger-v2", "exposure_ledger_id": ledger_id}


def _header(ledger_id: str, sequence: int, head_sha256: str) -> bytes:
    return canonical_json_bytes({**_anchor(ledger_id), "head_sequence": f"{sequence:020d}", "head_sha256": head_sha256})


def initialize_exposure_ledger(path: Path, *, forbidden_roots: tuple[Path, ...] = ()) -> str:
    """Exclusively create an external private ledger before study declarations.

    Args:
        path: New absolute file outside repositories and all supplied output roots.
        forbidden_roots: Study/input/output roots that must not own this history.

    Returns:
        Opaque 64-hex operator-ledger identity for binding into studies/authorities.

    Raises:
        ValueError: For occupied, linked, internal or invalid ledger paths.
        OSError: On file or durable-write failure; partial files are never reused.
    """
    path = _external_path(path, forbidden_roots)
    try:
        descriptor = os.open(path, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600)
    except FileExistsError as error:
        raise ValueError("calibration exposure ledger already exists") from error
    try:
        os.fchmod(descriptor, 0o600)
        identity = secrets.token_hex(32)
        _write_all(descriptor, _header(identity, 0, canonical_sha256(_anchor(identity))))
        return identity
    finally:
        os.close(descriptor)


def _line(stream: BinaryIO) -> dict[str, object] | None:
    data = stream.readline(MAX_ENTRY_BYTES + 1)
    if not data:
        return None
    if len(data) > MAX_ENTRY_BYTES or not data.endswith(b"\n"):
        _fail("calibration exposure history is oversized or torn")
    try:
        document = load_strict_json_object(data)
        if canonical_json_bytes(document) != data:
            _fail("calibration exposure history is not canonical")
    except (ValueError, RecursionError) as error:
        raise ValueError("calibration exposure history is corrupt") from error
    return dict(document)


def _history(descriptor: int, ledger_id: str) -> tuple[int, str, set[tuple[str, str, str]]]:
    os.lseek(descriptor, 0, os.SEEK_SET)
    seen: set[tuple[str, str, str]] = set()
    sequence = 0
    with os.fdopen(os.dup(descriptor), "rb") as stream:
        header = _line(stream)
        if (
            header is None
            or set(header) != _HEADER_FIELDS
            or header["schema_version"] != "calibration-target-exposure-ledger-v2"
        ):
            _fail("calibration exposure history header is corrupt or unsupported")
        if header["exposure_ledger_id"] != ledger_id:
            _fail("calibration exposure ledger identity differs from the study")
        head_sequence = header["head_sequence"]
        if (
            not isinstance(head_sequence, str)
            or len(head_sequence) != 20
            or any(character not in "0123456789" for character in head_sequence)
            or int(head_sequence) > 2**53 - 1
        ):
            _fail("calibration exposure history head sequence is corrupt")
        head_sha256 = require_digest(header["head_sha256"], "head digest")
        previous = canonical_sha256(_anchor(ledger_id))
        while (entry := _line(stream)) is not None:
            if (
                set(entry) != _ENTRY_FIELDS
                or entry["schema_version"] != "calibration-target-exposure-entry-v2"
                or entry["previous_sha256"] != previous
            ):
                _fail("calibration exposure history chain is corrupt")
            if entry["sha256"] != canonical_sha256({key: value for key, value in entry.items() if key != "sha256"}):
                _fail("calibration exposure history entry digest differs")
            receipt = decode_exposure_receipt(entry["receipt"])
            identities = decode_exposure_identities(entry["identities"])
            sequence += 1
            if (
                receipt.sequence != sequence
                or receipt.exposure_ledger_id != ledger_id
                or receipt.membership_sha256 != canonical_sha256(entry["identities"])
            ):
                _fail("calibration exposure history receipt bindings differ")
            keys = {(receipt.target, namespace, token) for namespace, token in identities}
            if receipt.role in CONFIRMATION_ROLES and keys & seen:
                _fail("calibration exposure history contains a repeated confirmation")
            seen.update(keys)
            previous = canonical_sha256(entry)
        if sequence != int(head_sequence) or previous != head_sha256:
            _fail("calibration exposure history differs from its committed head")
    return sequence, previous, seen


def _check_open_file(descriptor: int, path: Path) -> os.stat_result:
    observed = os.fstat(descriptor)
    named = os.stat(path, follow_symlinks=False)
    if (
        not stat.S_ISREG(observed.st_mode)
        or not stat.S_ISREG(named.st_mode)
        or stat.S_IMODE(observed.st_mode) != 0o600
        or (observed.st_dev, observed.st_ino) != (named.st_dev, named.st_ino)
    ):
        _fail("calibration exposure ledger changed or is not a private regular file")
    return observed


def record_exposure(
    path: Path,
    *,
    expected_ledger_id: str,
    target: str,
    role: str,
    study_sha256: str,
    partition_sha256: str,
    evidence_sha256: str,
    identities: object,
    forbidden_roots: tuple[Path, ...] = (),
) -> ExposureReceipt:
    """Durably record audited identities before a controller reads their outcomes.

    Confirmation refuses every prior exposure for the same target and any shared
    token, regardless of study, evidence path, partition or role renaming. Other
    roles are recorded for development reuse and cannot become unseen later.
    Receipt possession never authorizes reopening labels or resets retirement.

    Args:
        path: External operator-controlled ledger initialized explicitly.
        expected_ledger_id: Identity bound into the study/authority beforehand.
        target: Callers or total-length calibration namespace.
        role: Exact authorized study role; this layer records, not grants, authority.
        study_sha256: Frozen study digest.
        partition_sha256: Frozen partition digest.
        evidence_sha256: Evidence commitment checked before access.
        identities: All audited stable identity tokens, including specimen and read/backbone.
        forbidden_roots: Every input/study/output root excluded from ledger ownership.

    Returns:
        Aggregate-free receipt after record and head fsync, while holding the file lock.

    Raises:
        ValueError: On invalid/corrupt history, identity mismatch, or prior exposure.
        OSError: On read/write failure. Torn appends remain fail-closed and are never truncated.
    """
    path = _external_path(path, forbidden_roots)
    ledger_id = require_digest(expected_ledger_id, "ledger identity")
    tokens = decode_exposure_identities(identities)
    identity_rows = [{"namespace": namespace, "sha256": token} for namespace, token in tokens]
    draft = {
        "schema_version": "calibration-target-exposure-receipt-v2",
        "target": target,
        "role": role,
        "study_sha256": study_sha256,
        "partition_sha256": partition_sha256,
        "evidence_sha256": evidence_sha256,
        "membership_sha256": canonical_sha256(identity_rows),
        "exposure_ledger_id": ledger_id,
        "sequence": 1,
    }
    decode_exposure_receipt(draft)
    descriptor = os.open(path, os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    try:
        _check_open_file(descriptor, path)
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        before = _check_open_file(descriptor, path)
        sequence, previous, seen = _history(descriptor, ledger_id)
        after = _check_open_file(descriptor, path)
        if any(getattr(before, key) != getattr(after, key) for key in ("st_size", "st_mtime_ns", "st_ctime_ns")):
            _fail("calibration exposure history changed during its read")
        if role in CONFIRMATION_ROLES and {(target, namespace, token) for namespace, token in tokens} & seen:
            _fail("calibration evidence identities were previously exposed for this target")
        receipt = decode_exposure_receipt({**draft, "sequence": sequence + 1})
        entry = {
            "schema_version": "calibration-target-exposure-entry-v2",
            "previous_sha256": previous,
            "receipt": exposure_receipt_document(receipt),
            "identities": identity_rows,
        }
        entry["sha256"] = canonical_sha256(entry)
        data = canonical_json_bytes(entry)
        if len(data) > MAX_ENTRY_BYTES:
            _fail("calibration exposure entry exceeds its administrative byte limit")
        os.lseek(descriptor, 0, os.SEEK_END)
        _write_all(descriptor, data)
        os.lseek(descriptor, 0, os.SEEK_SET)
        _write_all(descriptor, _header(ledger_id, receipt.sequence, canonical_sha256(entry)))
        final = _check_open_file(descriptor, path)
        if final.st_size != before.st_size + len(data):
            _fail("calibration exposure append size differs")
        return receipt
    finally:
        os.close(descriptor)
