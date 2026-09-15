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


def payload_manifest_document(manifest: PayloadManifest) -> list[dict[str, object]]:
    """Project an immutable manifest into independent canonicalizable rows.

    Args:
        manifest: Decoded complete manifest.

    Returns:
        Fresh JSON-compatible file records.

    Raises:
        ValueError: If manifest is not a PayloadManifest.
    """
    if not isinstance(manifest, PayloadManifest):
        _fail("payload manifest must be a PayloadManifest")
    return [{"path": item.path, "size_bytes": item.size_bytes, "sha256": item.sha256} for item in manifest.files]


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
    if not isinstance(manifest, PayloadManifest):
        _fail("payload manifest must be a PayloadManifest")
    if not isinstance(observed, Mapping) or set(observed) != {item.path for item in manifest.files}:
        _fail("observed payload files must match the manifest exactly")
    for item in manifest.files:
        record = observed[item.path]
        if not isinstance(record, tuple) or len(record) != 2:
            _fail("payload observation must be a (size_bytes, sha256) tuple")
        size, digest = _size(record[0]), _digest(record[1])
        if size != item.size_bytes or digest != item.sha256:
            _fail("observed payload file differs from its manifest binding")
