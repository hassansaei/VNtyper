"""Nonblocking, bounded reads of operator-controlled preflight files."""

from __future__ import annotations

import io
import logging
import os
import stat
from collections.abc import Mapping
from pathlib import Path
from typing import Any, NoReturn

logger = logging.getLogger(__name__)

DEFAULT_PREFLIGHT_TEXT_MAX_BYTES = 1048576


def _reject(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def configured_preflight_text_limit(config: Mapping[str, Any]) -> int:
    """Return the positive byte limit for pre-probe BED and FAI reads.

    Args:
        config: Complete or replacement pipeline configuration.

    Returns:
        The configured maximum byte count.

    Raises:
        ValueError: If the configured value is not a positive integer.
    """
    utils = config.get("utils", {})
    if not isinstance(utils, Mapping):
        _reject("utils.preflight_text_max_bytes must be a positive integer")
    value = utils.get("preflight_text_max_bytes", DEFAULT_PREFLIGHT_TEXT_MAX_BYTES)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        _reject("utils.preflight_text_max_bytes must be a positive integer")
    return value


def _open_nonblocking(path: str | Path) -> int:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NONBLOCK", 0)
    return os.open(path, flags)


def regular_file_unavailable_reason(path: str | Path, *, description: str) -> str | None:
    """Inspect one byte without allowing a FIFO or device open to block.

    Args:
        path: Candidate local input path. Symlinks to regular files remain valid.
        description: Stable operator-facing name for the input.

    Returns:
        ``None`` for a readable regular file, otherwise an actionable reason.
    """
    try:
        descriptor = _open_nonblocking(path)
    except FileNotFoundError:
        return f"{description} not found"
    except OSError as error:
        return f"{description} unreadable: {error}"
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            return f"{description} is not a regular file"
        os.read(descriptor, 1)
    except OSError as error:
        return f"{description} unreadable: {error}"
    finally:
        os.close(descriptor)
    return None


def try_read_bounded_regular_text(
    path: str | Path,
    *,
    max_bytes: int,
    description: str,
) -> tuple[str | None, str | None]:
    """Try a bounded UTF-8 read without logging an optional-input failure.

    The nonblocking open makes FIFO substitution fail closed. A descriptor-level
    type check closes the stat/open race while preserving symlinks to regular files.

    Args:
        path: Local text input path.
        max_bytes: Maximum accepted encoded size.
        description: Stable operator-facing name for the input.

    Returns:
        The decoded text and ``None``, or ``None`` and an actionable reason.
    """
    try:
        descriptor = _open_nonblocking(path)
    except FileNotFoundError:
        return None, f"{description} not found: {path}"
    except OSError as error:
        return None, f"{description} unreadable: {path}: {error}"

    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            return None, f"{description} is not a regular file: {path}"
        if metadata.st_size > max_bytes:
            return None, f"{description} exceeds configured {max_bytes}-byte preflight limit: {path}"

        chunks: list[bytes] = []
        remaining = max_bytes + 1
        while remaining:
            chunk = os.read(descriptor, min(io.DEFAULT_BUFFER_SIZE, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        encoded = b"".join(chunks)
        if len(encoded) > max_bytes:
            return None, f"{description} exceeds configured {max_bytes}-byte preflight limit: {path}"
    except OSError as error:
        return None, f"{description} unreadable: {path}: {error}"
    finally:
        os.close(descriptor)

    try:
        return encoded.decode("utf-8"), None
    except UnicodeDecodeError:
        return None, f"{description} is not valid UTF-8: {path}"


def read_bounded_regular_text(
    path: str | Path,
    *,
    max_bytes: int,
    description: str,
) -> str:
    """Read required UTF-8 text from a regular file within a byte bound.

    Args:
        path: Local text input path.
        max_bytes: Maximum accepted encoded size.
        description: Stable operator-facing name for the input.

    Returns:
        Decoded UTF-8 text.

    Raises:
        ValueError: If the path is missing, unreadable, nonregular, oversized or
            not valid UTF-8.
    """
    text, reason = try_read_bounded_regular_text(path, max_bytes=max_bytes, description=description)
    if reason is not None:
        _reject(reason)
    assert text is not None
    return text
