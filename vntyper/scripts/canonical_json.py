"""Strict JSON decoding and RFC 8785 canonical serialization."""

from __future__ import annotations

import hashlib
import json
from typing import Any

import rfc8785


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Build an object while rejecting duplicate member names."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _reject_nonfinite(constant: str) -> None:
    """Reject the non-standard numeric constants accepted by ``json.loads``."""
    raise ValueError(f"non-finite JSON constant is not allowed: {constant}")


def load_strict_json_object(raw: bytes | str) -> dict[str, Any]:
    """Decode one JSON object with duplicate names and non-finite values forbidden.

    Args:
        raw: UTF-8 JSON bytes or text.

    Returns:
        The decoded top-level object.

    Raises:
        ValueError: If JSON is invalid, contains duplicate names/non-finite values,
            or has a non-object top level.
    """
    value = json.loads(raw, object_pairs_hook=_unique_object, parse_constant=_reject_nonfinite)
    if not isinstance(value, dict):
        raise ValueError("top-level JSON value must be an object")
    return value


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize a JSON-compatible value as RFC 8785 bytes plus one newline.

    Args:
        value: JSON-compatible value to serialize.

    Returns:
        Canonical UTF-8 bytes terminated by exactly one newline.
    """
    return rfc8785.dumps(value) + b"\n"


def canonical_sha256(value: Any) -> str:
    """Return the SHA-256 hex digest of canonical newline-terminated JSON.

    Args:
        value: JSON-compatible value to digest.

    Returns:
        Lowercase SHA-256 hexadecimal digest.
    """
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()
