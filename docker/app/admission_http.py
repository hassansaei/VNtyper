"""HTTP translation for the service's shared admission registry."""

import logging
import os
from typing import NoReturn

from fastapi import HTTPException
from redis.exceptions import RedisError

from . import admission

logger = logging.getLogger(__name__)


def reserve_or_refuse(
    redis_client: admission.RedisLike,
    *,
    reservation_id: str,
    requested_bytes: int,
    paths: tuple[str, ...],
    max_jobs: int,
    minimum_free_bytes: int,
    queued_grace_seconds: int,
) -> None:
    """Reserve shared capacity or answer with a retryable refusal.

    Args:
        redis_client: Synchronous redis-py compatible client.
        reservation_id: Identifier the endpoint and worker use for release.
        requested_bytes: Conservative workload-specific headroom estimate.
        paths: Existing directories whose free space constrains admission.
        max_jobs: Maximum number of accepted unfinished submissions.
        minimum_free_bytes: Floor retained after aggregate reservations.
        queued_grace_seconds: Grace for broker-accepted tasks to start.

    Raises:
        HTTPException: 503 when capacity is exhausted or cannot be measured.
    """
    try:
        free_bytes = admission.observe_free_bytes(paths)
        outcome = admission.reserve_capacity(
            redis_client,
            reservation_id=reservation_id,
            requested_bytes=requested_bytes,
            free_bytes=free_bytes,
            max_jobs=max_jobs,
            minimum_free_bytes=minimum_free_bytes,
            queued_grace_seconds=queued_grace_seconds,
        )
    except OSError as exc:
        _raise_unavailable("storage", exc, retry_after="300")
    except (RedisError, ValueError) as exc:
        _raise_unavailable("capacity", exc, retry_after="60")

    if outcome.accepted:
        return
    if outcome.reason == "queue":
        msg = f"The job queue is full ({outcome.active_jobs} submissions reserved); retry later"
        retry_after = "60"
    else:
        msg = "Not enough free storage to accept a new submission; retry later"
        retry_after = "300"
    logger.error(msg)
    raise HTTPException(status_code=503, detail=msg, headers={"Retry-After": retry_after})


def total_file_bytes_or_refuse(paths: tuple[str, ...]) -> int:
    """Return aggregate member bytes or a retryable storage refusal.

    Args:
        paths: Existing cohort archives to size for the reservation.

    Returns:
        int: Aggregate size in bytes.

    Raises:
        HTTPException: 503 if a member disappears or cannot be inspected.
    """
    try:
        return sum(os.path.getsize(path) for path in paths)
    except OSError as exc:
        _raise_unavailable("storage", exc, retry_after="300")


def _raise_unavailable(kind: str, exc: Exception, *, retry_after: str) -> NoReturn:
    """Raise the common path-safe 503 response for dependency failures."""
    msg = f"Submission {kind} is temporarily unavailable; retry later"
    logger.error(f"{msg}: {exc}")
    raise HTTPException(status_code=503, detail=msg, headers={"Retry-After": retry_after}) from exc
