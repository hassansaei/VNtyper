"""Read exact committed files after the controller authorizes outcome access.

These low-level readers grant no role access. The controller selects a role,
records its exposure, and claims confirmation custody before invoking them.
Large alignments can be verified in bounded memory without decoding outcomes.
"""

from __future__ import annotations

import hashlib
import logging
import os
import stat
from pathlib import Path
from typing import NoReturn

from vntyper.scripts.calibration_exposure import require_digest
from vntyper.scripts.calibration_target_runs import TargetRunAsset
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)
_CHUNK_BYTES = 1024 * 1024
_JSON_LIMIT = 64 * 1024 * 1024


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int]:
    return metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns


def _read(asset: TargetRunAsset, *, retain: bool, maximum_bytes: int | None) -> bytes:
    if not isinstance(asset, TargetRunAsset):
        _fail("target asset reader requires a typed commitment")
    if not isinstance(asset.path, Path) or not asset.path.is_absolute() or ".." in asset.path.parts:
        _fail("target asset path must be absolute and normalized")
    require_digest(asset.sha256, "target asset digest")
    if type(asset.size_bytes) is not int or asset.size_bytes < 0:
        _fail("target asset size must be a nonnegative integer")
    if maximum_bytes is not None and (
        type(maximum_bytes) is not int or maximum_bytes < 0 or asset.size_bytes > maximum_bytes
    ):
        _fail("target asset exceeds the declared byte limit")
    descriptor = os.open(asset.path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            _fail("target asset must be a regular file")
        if before.st_size != asset.size_bytes:
            _fail("target asset size differs from its commitment")
        digest = hashlib.sha256()
        chunks = []
        size = 0
        while True:
            chunk = os.read(descriptor, min(_CHUNK_BYTES, asset.size_bytes - size + 1))
            if not chunk:
                break
            size += len(chunk)
            if size > asset.size_bytes:
                _fail("target asset size changed during reading")
            digest.update(chunk)
            if retain:
                chunks.append(chunk)
        after = os.fstat(descriptor)
        current = asset.path.lstat()
        if _identity(before) != _identity(after) or _identity(before) != _identity(current):
            _fail("target asset changed during reading")
        if size != asset.size_bytes or digest.hexdigest() != asset.sha256:
            _fail("target asset digest differs from its commitment")
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def read_target_asset(asset: TargetRunAsset, *, maximum_bytes: int = _JSON_LIMIT) -> bytes:
    """Read one size-limited committed result from a stable file descriptor.

    Args:
        asset: Exact local path, raw SHA256, and byte size.
        maximum_bytes: Administrative memory bound, checked before opening.

    Returns:
        Verified original bytes, with no decoding or normalization.

    Raises:
        ValueError: If bytes, size, file identity, or administrative limit differ.
        OSError: If the committed regular file cannot be opened or read.
    """
    return _read(asset, retain=True, maximum_bytes=maximum_bytes)


def read_target_json(asset: TargetRunAsset, *, maximum_bytes: int = _JSON_LIMIT) -> dict[str, object]:
    """Decode verified result bytes with duplicate keys and nonfinite numbers refused.

    Args:
        asset: Exact committed JSON result.
        maximum_bytes: Administrative memory bound.

    Returns:
        Original JSON object without coercing values.

    Raises:
        ValueError: If the file differs or is not a strict JSON object.
        OSError: If the file cannot be read.
    """
    return load_strict_json_object(read_target_asset(asset, maximum_bytes=maximum_bytes))


def verify_target_asset(asset: TargetRunAsset) -> str:
    """Verify a potentially large input using bounded-memory SHA256 streaming.

    Args:
        asset: Exact committed input or reference.

    Returns:
        The verified raw file digest.

    Raises:
        ValueError: If the input differs from its commitment or changes while read.
        OSError: If the file cannot be read.
    """
    _read(asset, retain=False, maximum_bytes=None)
    return asset.sha256
