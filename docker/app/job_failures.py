"""Curated preflight-failure transport shared by the worker and status API."""

from __future__ import annotations

import json
import logging
import os
import stat
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

PREFLIGHT_ERROR_FILENAME = "preflight_error.json"
_PAYLOAD_KEYS = {"code", "message", "candidates"}
_MAX_ARTIFACT_BYTES = 64 * 1024


def _message_is_public(message: str) -> bool:
    return bool(message.strip()) and "/" not in message and "\\" not in message


def read_preflight_failure(output_dir: str | Path) -> dict[str, str] | None:
    """Read a valid, path-free error artifact without following a symlink.

    Args:
        output_dir: Per-job pipeline output directory.

    Returns:
        The code and public message, or ``None`` for an absent or invalid artifact.
    """
    artifact = Path(output_dir) / PREFLIGHT_ERROR_FILENAME
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(artifact, flags)
    except OSError:
        return None
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            return None
        with os.fdopen(descriptor, encoding="utf-8") as artifact_file:
            descriptor = -1
            raw = artifact_file.read(_MAX_ARTIFACT_BYTES + 1)
    except (OSError, UnicodeError):
        return None
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    if len(raw.encode("utf-8")) > _MAX_ARTIFACT_BYTES:
        return None
    try:
        payload = json.loads(raw)
    except (TypeError, ValueError):
        return None
    if not isinstance(payload, dict) or set(payload) != _PAYLOAD_KEYS or not isinstance(payload["candidates"], list):
        return None
    code = payload["code"]
    message = payload["message"]
    if not isinstance(code, str) or not code.strip() or not isinstance(message, str) or not _message_is_public(message):
        return None
    return {"code": code, "message": message}


def stored_preflight_message(redis_client: Any, job_id: str) -> str | None:
    """Return a stored curated message when both transport fields are present.

    Args:
        redis_client: Redis-compatible client for usage job hashes.
        job_id: Canonical job identifier.

    Returns:
        A path-free preflight message, or ``None`` for generic failures.
    """
    job_key = f"usage:{job_id}"
    code = redis_client.hget(job_key, "code")
    message = redis_client.hget(job_key, "message")
    if not isinstance(code, str) or not code.strip() or not isinstance(message, str):
        return None
    return message if _message_is_public(message) else None
