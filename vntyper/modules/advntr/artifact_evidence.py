"""Verified governed evidence for recurrent adVNTR State strings."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, cast

from vntyper.scripts.canonical_json import canonical_json_bytes, canonical_sha256, load_strict_json_object

ASSERTION = "A carried-forward recurrent adVNTR State is insufficient for molecular identity."
_MODULE_DIR = Path(__file__).resolve().parent
_PACKAGED_ARTIFACT = _MODULE_DIR / "advntr_artifact_evidence.json"
_PACKAGED_DIGEST = _MODULE_DIR / "advntr_artifact_evidence.sha256"
_ROOT_FIELDS = {
    "advntr_version_upper_bound_exclusive",
    "assertion",
    "cohort_name",
    "entries",
    "exact_advntr_version",
    "schema_version",
}
_ENTRY_FIELDS = {"active", "count", "denominator", "frequency", "state", "status"}
_STATUSES = frozenset({"confirmed_artifact", "pending_renome_revalidation"})
_SHA256_RE = re.compile(r"[0-9a-f]{64}")

EvidenceStatus = Literal["confirmed_artifact", "pending_renome_revalidation"]


@dataclass(frozen=True)
class ArtifactEvidenceEntry:
    """One governed recurrent State and the facts actually retained about it."""

    state: str
    status: EvidenceStatus
    active: bool
    count: int | None
    denominator: int | None
    frequency: float | None


@dataclass(frozen=True)
class ArtifactEvidence:
    """Verified immutable artifact-evidence document and its canonical identity."""

    schema_version: int
    assertion: str
    cohort_name: str
    advntr_version_upper_bound_exclusive: str
    exact_advntr_version: str | None
    entries: tuple[ArtifactEvidenceEntry, ...]
    digest: str
    canonical_bytes: bytes

    @property
    def active_states(self) -> frozenset[str]:
        """Return the exact State strings whose governed evidence is active."""
        return frozenset(entry.state for entry in self.entries if entry.active)


def _require_exact_fields(value: dict[str, Any], expected: set[str], *, label: str) -> None:
    actual = set(value)
    if actual != expected:
        raise ValueError(f"{label} fields differ: expected {sorted(expected)}, got {sorted(actual)}")


def _optional_nonnegative_int(value: object, *, field: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"artifact evidence {field} must be a non-negative integer or null")
    return value


def _optional_frequency(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
        raise ValueError("artifact evidence frequency must be between 0 and 1 or null")
    return float(value)


def _parse_entry(value: object, *, index: int) -> ArtifactEvidenceEntry:
    if not isinstance(value, dict):
        raise ValueError(f"artifact evidence entries[{index}] must be an object")
    _require_exact_fields(value, _ENTRY_FIELDS, label=f"artifact evidence entries[{index}]")
    state = value["state"]
    status = value["status"]
    active = value["active"]
    if not isinstance(state, str) or not state:
        raise ValueError(f"artifact evidence entries[{index}].state must be a non-empty string")
    if not isinstance(status, str) or status not in _STATUSES:
        raise ValueError(f"artifact evidence entries[{index}].status is unsupported: {status}")
    if not isinstance(active, bool):
        raise ValueError(f"artifact evidence entries[{index}].active must be a Boolean")
    return ArtifactEvidenceEntry(
        state=state,
        status=cast(EvidenceStatus, status),
        active=active,
        count=_optional_nonnegative_int(value["count"], field="count"),
        denominator=_optional_nonnegative_int(value["denominator"], field="denominator"),
        frequency=_optional_frequency(value["frequency"]),
    )


def _parse_document(value: dict[str, Any], *, digest: str, canonical_bytes: bytes) -> ArtifactEvidence:
    _require_exact_fields(value, _ROOT_FIELDS, label="artifact evidence root")
    schema_version = value["schema_version"]
    assertion = value["assertion"]
    cohort_name = value["cohort_name"]
    version_upper = value["advntr_version_upper_bound_exclusive"]
    exact_version = value["exact_advntr_version"]
    raw_entries = value["entries"]
    if isinstance(schema_version, bool) or schema_version != 1:
        raise ValueError("artifact evidence schema_version must be 1")
    if assertion != ASSERTION:
        raise ValueError("artifact evidence assertion differs from the approved limited assertion")
    if cohort_name != "renome":
        raise ValueError("artifact evidence cohort_name must be renome")
    if version_upper != "2.0.4":
        raise ValueError("artifact evidence adVNTR version upper bound must be 2.0.4")
    if exact_version is not None and (not isinstance(exact_version, str) or not exact_version):
        raise ValueError("artifact evidence exact_advntr_version must be a non-empty string or null")
    if not isinstance(raw_entries, list) or not raw_entries:
        raise ValueError("artifact evidence entries must be a non-empty list")
    entries = tuple(_parse_entry(entry, index=index) for index, entry in enumerate(raw_entries))
    states = [entry.state for entry in entries]
    if len(states) != len(set(states)):
        raise ValueError("artifact evidence State strings must be unique")
    return ArtifactEvidence(
        schema_version=schema_version,
        assertion=assertion,
        cohort_name=cohort_name,
        advntr_version_upper_bound_exclusive=version_upper,
        exact_advntr_version=exact_version,
        entries=entries,
        digest=digest,
        canonical_bytes=canonical_bytes,
    )


def load_artifact_evidence(path: str | Path, *, expected_digest: str | None = None) -> ArtifactEvidence:
    """Load and verify a canonical governed artifact-evidence file.

    Args:
        path: JSON artifact path.
        expected_digest: Optional required canonical SHA-256 digest.

    Returns:
        Frozen verified artifact evidence.

    Raises:
        OSError: If the file cannot be read.
        ValueError: If bytes, digest, schema, or values violate the contract.
    """
    raw = Path(path).read_bytes()
    value = load_strict_json_object(raw)
    canonical = canonical_json_bytes(value)
    if raw != canonical:
        raise ValueError("artifact evidence must use canonical RFC 8785 bytes plus one newline")
    digest = canonical_sha256(value)
    if expected_digest is not None and digest != expected_digest:
        raise ValueError(f"artifact evidence digest mismatch: expected {expected_digest}, got {digest}")
    return _parse_document(value, digest=digest, canonical_bytes=canonical)


def load_packaged_artifact_evidence() -> ArtifactEvidence:
    """Load the packaged artifact after verifying its canonical digest sidecar.

    Returns:
        Frozen verified packaged evidence.

    Raises:
        OSError: If either packaged file cannot be read.
        ValueError: If the digest sidecar or artifact is invalid.
    """
    expected_digest = _PACKAGED_DIGEST.read_text(encoding="ascii").strip()
    if _SHA256_RE.fullmatch(expected_digest) is None:
        raise ValueError("packaged artifact evidence digest must be 64 lowercase hexadecimal characters")
    return load_artifact_evidence(_PACKAGED_ARTIFACT, expected_digest=expected_digest)
