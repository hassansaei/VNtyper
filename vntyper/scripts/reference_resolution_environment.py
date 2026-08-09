"""Run-scoped htslib reference-resolution environment policy."""

from __future__ import annotations

import logging
import os
import threading

from vntyper.scripts.reference_uri_policy import allow_ambient_reference_resolution, ref_path_remote_scheme

logger = logging.getLogger(__name__)

_REFERENCE_ENVIRONMENT_LOCK = threading.RLock()
_REFERENCE_ENVIRONMENT_STATE = threading.local()


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
        depth = getattr(_REFERENCE_ENVIRONMENT_STATE, "depth", 0)
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
    finally:
        if depth > 0:
            _REFERENCE_ENVIRONMENT_STATE.depth = depth - 1
            _REFERENCE_ENVIRONMENT_LOCK.release()
