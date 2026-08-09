"""Run-scoped htslib reference-resolution environment policy."""

from __future__ import annotations

import logging
import os
import threading

from vntyper.scripts.reference_uri_policy import allow_ambient_reference_resolution, ref_path_remote_scheme

logger = logging.getLogger(__name__)

_REFERENCE_ENVIRONMENT_LOCK = threading.RLock()
_REFERENCE_ENVIRONMENT_STATE = threading.local()


def activate_private_reference_cache(pattern: str, *, remote_reference_path: str | None = None) -> None:
    """Restrict htslib lookup to one complete pipeline-owned M5 cache.

    Args:
        pattern: Private progressive-digest cache template.
        remote_reference_path: Optional remote-only lookup retained for header contigs
            without a local URI.

    Raises:
        RuntimeError: If no reference-resolution scope is active.
    """
    if getattr(_REFERENCE_ENVIRONMENT_STATE, "depth", 0) <= 0:
        raise RuntimeError("Private CRAM reference cache requires an active resolution scope")
    os.environ["REF_CACHE"] = pattern
    os.environ["REF_PATH"] = remote_reference_path if remote_reference_path is not None else pattern


def pin_reference_resolution(config: dict) -> str | None:
    """Pin htslib reference lookup and exclusively own the process environment.

    Args:
        config: Pipeline configuration. Missing CRAM keys use shipped defaults.

    Returns:
        The prior ``REF_PATH``, including ``None``; restore it in a ``finally``.

    Raises:
        ValueError: If the CRAM reference-resolution policy is invalid.
    """
    _REFERENCE_ENVIRONMENT_LOCK.acquire()
    try:
        previous = os.environ.get("REF_PATH")
        previous_cache = os.environ.get("REF_CACHE")
        cram_config = config.get("cram", {})
        if allow_ambient_reference_resolution(config):
            logger.warning("Ambient CRAM reference resolution is enabled and may block on a network endpoint.")
        else:
            local_ref_path = cram_config.get("local_ref_path", "%2s/%2s/%s")
            if ref_path_remote_scheme(local_ref_path) is not None:
                message = (
                    "cram.local_ref_path must not contain a remote URL when ambient reference resolution is disabled: "
                    f"{local_ref_path}"
                )
                logger.error(message)
                raise ValueError(message)
            os.environ["REF_PATH"] = local_ref_path
            os.environ.pop("REF_CACHE", None)
        depth = getattr(_REFERENCE_ENVIRONMENT_STATE, "depth", 0)
        stack = getattr(_REFERENCE_ENVIRONMENT_STATE, "cache_stack", [])
        stack.append(previous_cache)
        _REFERENCE_ENVIRONMENT_STATE.cache_stack = stack
        _REFERENCE_ENVIRONMENT_STATE.depth = depth + 1
        return previous
    except BaseException:
        _REFERENCE_ENVIRONMENT_LOCK.release()
        raise


def restore_reference_resolution(previous: str | None) -> None:
    """Restore ``REF_PATH`` and release this thread's reference environment scope.

    Args:
        previous: Original ``REF_PATH`` value, or ``None`` when it was unset.
    """
    depth = getattr(_REFERENCE_ENVIRONMENT_STATE, "depth", 0)
    try:
        if previous is None:
            os.environ.pop("REF_PATH", None)
        else:
            os.environ["REF_PATH"] = previous
        stack = getattr(_REFERENCE_ENVIRONMENT_STATE, "cache_stack", [])
        if stack:
            previous_cache = stack.pop()
            if previous_cache is None:
                os.environ.pop("REF_CACHE", None)
            else:
                os.environ["REF_CACHE"] = previous_cache
    finally:
        if depth > 0:
            _REFERENCE_ENVIRONMENT_STATE.depth = depth - 1
            _REFERENCE_ENVIRONMENT_LOCK.release()
