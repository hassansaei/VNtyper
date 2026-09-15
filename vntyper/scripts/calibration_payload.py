"""Pure contracts binding the complete layout and bytes of a calibration payload."""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import NoReturn

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
CALLER_BUNDLE_DESCRIPTOR_PATH = "caller-bundle.json"
_CALLER_BUNDLE_FIELDS = {"schema_version", "required_callers", "components"}
_CALLER_COMPONENT_FIELDS = {"decision-profile.json", "advntr-policy.json", "background.json"}


@dataclass(frozen=True)
class PayloadFile:
    """Expected bytes of one normalized relative payload file."""

    path: str
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class PayloadManifest:
    """Complete file set and the canonical digest binding all its members."""

    files: tuple[PayloadFile, ...]
    sha256: str


@dataclass(frozen=True)
class CallerBundleDescriptor:
    """Closed caller payload composition and the digests of its component files."""

    required_callers: tuple[str, ...]
    decision_profile_sha256: str
    advntr_policy_sha256: str | None
    background_sha256: str | None
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _path(value: object) -> str:
    if not isinstance(value, str) or not value:
        _fail("payload path must be a non-empty string")
    parts = value.split("/")
    if (
        any(part in {"", ".", ".."} for part in parts)
        or ":" in value
        or "\\" in value
        or any(ord(char) < 32 or ord(char) == 127 for char in value)
        or value == "payload-manifest.json"
    ):
        _fail("payload path must be normalized, relative, and not refer to its manifest")
    return value


def _size(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        _fail("payload size_bytes must be a non-negative integer")
    return value


def _digest(value: object) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail("payload sha256 must be a lowercase SHA-256 digest")
    return value


def _optional_digest(value: object) -> str | None:
    if value is None:
        return None
    return _digest(value)


def _caller_set(value: object) -> tuple[str, ...]:
    if not isinstance(value, list) or not value or any(not isinstance(item, str) for item in value):
        _fail("required_callers must be a non-empty sorted unique list")
    callers = tuple(value)
    if callers != tuple(sorted(set(callers))) or "kestrel" not in callers or not set(callers) <= {"kestrel", "advntr"}:
        _fail("required_callers must include kestrel and optionally advntr")
    return callers


def decode_payload_manifest(value: object) -> PayloadManifest:
    """Decode a closed, sorted manifest without reading its files.

    Args:
        value: Decoded JSON array of path, size_bytes and sha256 objects.

    Returns:
        Immutable manifest whose digest includes every listed file.

    Raises:
        ValueError: If fields, layout, sizes or hashes are invalid.
    """
    if not isinstance(value, list) or not value:
        _fail("payload manifest must be a non-empty list")
    files = []
    for row in value:
        if not isinstance(row, Mapping) or set(row) != {"path", "size_bytes", "sha256"}:
            _fail("payload file fields must be path, size_bytes and sha256")
        files.append(PayloadFile(_path(row["path"]), _size(row["size_bytes"]), _digest(row["sha256"])))
    paths = tuple(item.path for item in files)
    if paths != tuple(sorted(set(paths))):
        _fail("payload layout must contain unique files in increasing path order")
    names = set(paths)
    for name in paths:
        if any(str(parent) in names for parent in PurePosixPath(name).parents):
            _fail("payload layout cannot place a file beneath another file")
    return PayloadManifest(tuple(files), canonical_sha256(value))


def decode_caller_bundle_descriptor(value: object) -> CallerBundleDescriptor:
    """Decode the closed composition descriptor for a caller calibration payload.

    Args:
        value: Decoded JSON descriptor object.

    Returns:
        Immutable caller set, component digests and descriptor digest.

    Raises:
        ValueError: If the schema, caller set or conditional component hashes are invalid.
    """
    if not isinstance(value, Mapping) or set(value) != _CALLER_BUNDLE_FIELDS:
        _fail("caller bundle descriptor fields differ from the closed contract")
    if value["schema_version"] != "caller-bundle-v2":
        _fail("caller bundle schema_version must be caller-bundle-v2")
    callers = _caller_set(value["required_callers"])
    components = value["components"]
    if not isinstance(components, Mapping) or set(components) != _CALLER_COMPONENT_FIELDS:
        _fail("caller bundle component fields differ from the closed contract")
    decision_profile = _digest(components["decision-profile.json"])
    advntr_policy = _optional_digest(components["advntr-policy.json"])
    background = _optional_digest(components["background.json"])
    if ("advntr" in callers) != (advntr_policy is not None):
        _fail("advntr-policy.json digest must be present exactly when advntr is required")
    if "advntr" not in callers and background is not None:
        _fail("background.json requires an advntr caller policy")
    return CallerBundleDescriptor(callers, decision_profile, advntr_policy, background, canonical_sha256(value))


def _manifest_rows(manifest: PayloadManifest) -> list[dict[str, object]]:
    return [{"path": item.path, "size_bytes": item.size_bytes, "sha256": item.sha256} for item in manifest.files]


def _require_payload_manifest(manifest: object) -> PayloadManifest:
    if not isinstance(manifest, PayloadManifest):
        _fail("payload manifest must be a PayloadManifest")
    if not isinstance(manifest.files, tuple) or any(not isinstance(item, PayloadFile) for item in manifest.files):
        _fail("payload manifest must use decoded immutable content")
    decoded = decode_payload_manifest(_manifest_rows(manifest))
    if decoded != manifest:
        _fail("payload manifest differs from its canonical content or digest")
    return manifest


def _caller_bundle_document(descriptor: CallerBundleDescriptor) -> dict[str, object]:
    return {
        "schema_version": "caller-bundle-v2",
        "required_callers": list(descriptor.required_callers),
        "components": {
            "decision-profile.json": descriptor.decision_profile_sha256,
            "advntr-policy.json": descriptor.advntr_policy_sha256,
            "background.json": descriptor.background_sha256,
        },
    }


def _require_caller_bundle_descriptor(descriptor: object) -> CallerBundleDescriptor:
    if not isinstance(descriptor, CallerBundleDescriptor):
        _fail("caller descriptor must be a CallerBundleDescriptor")
    if not isinstance(descriptor.required_callers, tuple):
        _fail("caller descriptor must use decoded immutable content")
    decoded = decode_caller_bundle_descriptor(_caller_bundle_document(descriptor))
    if decoded != descriptor:
        _fail("caller descriptor differs from its canonical content or digest")
    return descriptor


def payload_manifest_document(manifest: PayloadManifest) -> list[dict[str, object]]:
    """Project an immutable manifest into independent canonicalizable rows.

    Args:
        manifest: Decoded complete manifest.

    Returns:
        Fresh JSON-compatible file records.

    Raises:
        ValueError: If manifest is not a PayloadManifest.
    """
    return _manifest_rows(_require_payload_manifest(manifest))


def caller_bundle_descriptor_document(descriptor: CallerBundleDescriptor) -> dict[str, object]:
    """Project a validated caller bundle descriptor into an independent document.

    Args:
        descriptor: Decoded caller bundle descriptor.

    Returns:
        Fresh JSON-compatible descriptor content.

    Raises:
        ValueError: If the typed value is not its canonical decoded content.
    """
    return _caller_bundle_document(_require_caller_bundle_descriptor(descriptor))


def validate_payload_observations(manifest: PayloadManifest, observed: Mapping[str, tuple[int, str]]) -> None:
    """Require a reader's observed file set and bytes to match the manifest.

    The filesystem adapter must separately reject symlinks and nonregular files.
    This pure comparison cannot establish how an observation was obtained.

    Args:
        manifest: Decoded expected payload.
        observed: Relative file names mapped to observed byte count and digest.

    Raises:
        ValueError: If a file is missing, extra, malformed or changed.
    """
    manifest = _require_payload_manifest(manifest)
    if not isinstance(observed, Mapping) or set(observed) != {item.path for item in manifest.files}:
        _fail("observed payload files must match the manifest exactly")
    for item in manifest.files:
        record = observed[item.path]
        if not isinstance(record, tuple) or len(record) != 2:
            _fail("payload observation must be a (size_bytes, sha256) tuple")
        size, digest = _size(record[0]), _digest(record[1])
        if size != item.size_bytes or digest != item.sha256:
            _fail("observed payload file differs from its manifest binding")
