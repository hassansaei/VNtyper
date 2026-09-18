"""Target-aware exposure identities and aggregate-free receipt contracts."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import NoReturn

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
ROLES = frozenset({"training", "policy-selection", "validation", "locked-heldout", "development-assessment"})
CONFIRMATION_ROLES = frozenset({"validation", "locked-heldout"})
_NAMESPACES = frozenset(
    {"specimen", "family", "named-readset", "unnamed-readset", "physical-readset", "backbone", "pair", "seed-family"}
)
_READ_IDENTITIES = frozenset({"named-readset", "unnamed-readset", "physical-readset", "backbone"})
_FIELDS = {
    "schema_version",
    "target",
    "role",
    "study_sha256",
    "partition_sha256",
    "evidence_sha256",
    "membership_sha256",
    "exposure_ledger_id",
    "sequence",
}


@dataclass(frozen=True)
class ExposureReceipt:
    """One durable pre-access declaration; never permission to reopen outcomes."""

    target: str
    role: str
    study_sha256: str
    partition_sha256: str
    evidence_sha256: str
    membership_sha256: str
    exposure_ledger_id: str
    sequence: int
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def require_digest(value: object, label: str) -> str:
    """Validate one opaque lowercase SHA256-shaped identifier."""
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        _fail(f"calibration exposure {label} must be 64 lowercase hexadecimal characters")
    return value


def decode_exposure_identities(value: object) -> tuple[tuple[str, str], ...]:
    """Validate stable identity tokens supplied by an independently audited intake.

    Args:
        value: Sorted unique namespace/SHA256 rows. Callers must supply every
            applicable specimen/family/readset/backbone token from the intake.

    Returns:
        Immutable token pairs, with a specimen and read/backbone identity required.

    Raises:
        ValueError: For missing identities, unknown namespaces or malformed tokens.
    """
    if not isinstance(value, list) or not value:
        _fail("calibration exposure requires nonempty identity rows")
    rows = []
    for row in value:
        if not isinstance(row, Mapping) or set(row) != {"namespace", "sha256"}:
            _fail("calibration exposure identity fields differ")
        namespace = row["namespace"]
        if not isinstance(namespace, str) or namespace not in _NAMESPACES:
            _fail("calibration exposure identity namespace is unknown")
        rows.append((namespace, require_digest(row["sha256"], "identity token")))
    if rows != sorted(set(rows)):
        _fail("calibration exposure identity rows must be unique and sorted")
    namespaces = {namespace for namespace, _ in rows}
    if "specimen" not in namespaces or not namespaces & _READ_IDENTITIES:
        _fail("calibration exposure requires specimen and read/backbone identities")
    return tuple(rows)


def decode_exposure_receipt(value: object) -> ExposureReceipt:
    """Decode a closed receipt without treating it as access authorization.

    Args:
        value: Canonical receipt object emitted before any outcome access.

    Returns:
        Immutable identity-bound declaration and its canonical digest.

    Raises:
        ValueError: For wrong roles, targets, identity fields or sequence numbers.
    """
    if not isinstance(value, Mapping) or set(value) != _FIELDS:
        _fail("calibration exposure receipt fields differ")
    if value["schema_version"] != "calibration-target-exposure-receipt-v2":
        _fail("calibration exposure receipt schema is unsupported")
    target, role = value["target"], value["role"]
    if not isinstance(target, str) or target not in {"length", "callers"}:
        _fail("calibration exposure target is unsupported")
    if not isinstance(role, str) or role not in ROLES:
        _fail("calibration exposure role is unsupported")
    sequence = value["sequence"]
    if type(sequence) is not int or not 1 <= sequence <= 2**53 - 1:
        _fail("calibration exposure sequence must be a positive JSON-safe integer")
    return ExposureReceipt(
        target,
        role,
        require_digest(value["study_sha256"], "study"),
        require_digest(value["partition_sha256"], "partition"),
        require_digest(value["evidence_sha256"], "evidence"),
        require_digest(value["membership_sha256"], "membership"),
        require_digest(value["exposure_ledger_id"], "ledger identity"),
        sequence,
        canonical_sha256(value),
    )


def exposure_receipt_document(receipt: ExposureReceipt) -> dict[str, object]:
    """Return the exact aggregate-free receipt after verifying its typed content.

    Args:
        receipt: Previously decoded receipt.

    Returns:
        New canonical-content mapping, excluding its own derived SHA256.

    Raises:
        ValueError: If the receipt was forged or altered.
    """
    if not isinstance(receipt, ExposureReceipt):
        _fail("calibration exposure receipt must be typed")
    value = {
        "schema_version": "calibration-target-exposure-receipt-v2",
        **{name: getattr(receipt, name) for name in _FIELDS - {"schema_version"}},
    }
    if decode_exposure_receipt(value) != receipt:
        _fail("calibration exposure receipt identity differs from content")
    return value
