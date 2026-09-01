"""Strict loading and canonical identity for one complete decision profile."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.decision_profile_schema import component_projection, validate_complete_inventory

ProfileKind = Literal["packaged", "explicit-custom", "generated"]
ProfileSource = Literal["package", "explicit-cli"]

_PROFILE_RESOURCE = files("vntyper").joinpath("profiles/decision_profile.json")
_DIGEST_RESOURCE = files("vntyper").joinpath("profiles/decision_profile.sha256")
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_COMPONENT_NAMES = ("advntr", "cross_match", "dominance", "kestrel", "nomenclature", "shark")


@dataclass(frozen=True)
class ResolvedDecisionProfile:
    """One validated profile with canonical snapshot bytes and component values."""

    profile_id: str
    profile_revision: str
    profile_kind: ProfileKind
    source: ProfileSource
    digest: str
    canonical_bytes: bytes
    document: Mapping[str, object]
    components: Mapping[str, object]


def _freeze(value: object) -> object:
    """Recursively freeze one JSON-compatible resolved value."""
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(child) for child in value)
    return value


def _thaw(value: object) -> object:
    """Restore immutable resolved JSON containers for schema validation."""
    if isinstance(value, Mapping):
        return {cast(str, key): _thaw(child) for key, child in value.items()}
    if isinstance(value, tuple):
        return [_thaw(child) for child in value]
    return value


def _resolved(
    document: dict[str, object],
    *,
    source: ProfileSource,
    packaged_document: Mapping[str, object] | None = None,
) -> ResolvedDecisionProfile:
    raw_kind = document.get("profile_kind")
    if source == "explicit-cli" and raw_kind == "packaged":
        raise ValueError("an explicit CLI decision profile must be explicit-custom or generated")
    mutable_packaged_document = (
        None if packaged_document is None else cast(Mapping[str, object], _thaw(packaged_document))
    )
    validate_complete_inventory(document, packaged_profile=mutable_packaged_document)
    kind = cast(ProfileKind, document["profile_kind"])
    if source == "package" and kind != "packaged":
        raise ValueError("package decision profile must have profile_kind 'packaged'")
    if source == "explicit-cli" and kind not in {"explicit-custom", "generated"}:
        raise ValueError("an explicit CLI decision profile must be explicit-custom or generated")
    canonical = canonical_json_bytes(document)
    components = MappingProxyType(
        {
            name: _freeze(component_projection(document, name, packaged_profile=mutable_packaged_document))
            for name in _COMPONENT_NAMES
        }
    )
    return ResolvedDecisionProfile(
        profile_id=cast(str, document["profile_id"]),
        profile_revision=cast(str, document["profile_revision"]),
        profile_kind=kind,
        source=source,
        digest=hashlib.sha256(canonical).hexdigest(),
        canonical_bytes=canonical,
        document=cast(Mapping[str, object], _freeze(document)),
        components=components,
    )


def parse_decision_profile(
    raw: bytes | str,
    *,
    packaged_document: Mapping[str, object],
) -> ResolvedDecisionProfile:
    """Parse and validate one complete explicitly selected decision profile.

    Args:
        raw: UTF-8 JSON bytes or text from the explicit path.
        packaged_document: Verified packaged baseline used for completeness and
            immutable-field checks.

    Returns:
        Canonical resolved explicit profile.

    Raises:
        ValueError: If JSON decoding or any closed-schema rule fails.
    """
    document = load_strict_json_object(raw)
    return _resolved(document, source="explicit-cli", packaged_document=packaged_document)


def load_packaged_decision_profile() -> ResolvedDecisionProfile:
    """Load the packaged profile after verifying canonical bytes and sidecar.

    Returns:
        Verified packaged profile.

    Raises:
        OSError: If package resources cannot be read.
        ValueError: If bytes, sidecar, schema, or values differ from contract.
    """
    raw = _PROFILE_RESOURCE.read_bytes()
    document = load_strict_json_object(raw)
    canonical = canonical_json_bytes(document)
    if raw != canonical:
        raise ValueError("packaged decision profile must use canonical RFC 8785 bytes plus one newline")
    expected_digest = _DIGEST_RESOURCE.read_text(encoding="ascii").strip()
    if _SHA256_RE.fullmatch(expected_digest) is None:
        raise ValueError("packaged decision profile digest must be 64 lowercase hexadecimal characters")
    actual_digest = hashlib.sha256(raw).hexdigest()
    if actual_digest != expected_digest:
        raise ValueError(f"packaged decision profile digest mismatch: expected {expected_digest}, got {actual_digest}")
    resolved = _resolved(document, source="package")
    if resolved.digest != expected_digest:
        raise ValueError("packaged decision profile canonical digest differs from its verified bytes")
    return resolved


def resolve_decision_profile(path: str | Path | None = None) -> ResolvedDecisionProfile:
    """Resolve exactly one packaged or explicitly selected complete profile.

    Args:
        path: Explicit complete profile path, or None for the packaged profile.

    Returns:
        Verified resolved profile. No overlays or fallback are attempted.

    Raises:
        ValueError: If an explicit path cannot be read or its contents are invalid.
    """
    packaged = load_packaged_decision_profile()
    if path is None:
        return packaged
    candidate = Path(path)
    try:
        raw = candidate.read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read explicit decision profile {candidate}: {error}") from error
    return parse_decision_profile(raw, packaged_document=packaged.document)
