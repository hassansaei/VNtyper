"""Atomic queue and aggregate-storage admission reservations."""

import logging
import shutil
import time
from collections.abc import Callable
from dataclasses import dataclass
from threading import Event, Thread
from types import TracebackType
from typing import Any, Protocol

from redis.exceptions import RedisError, WatchError

logger = logging.getLogger(__name__)

ADMISSION_RESERVATIONS = "vntyper:admission:reservations"
ADMISSION_DEADLINES = "vntyper:admission:deadlines"
DEFAULT_QUEUED_GRACE_SECONDS = 7 * 24 * 60 * 60


class RedisLike(Protocol):
    """The small redis-py surface used by admission reservations."""

    def pipeline(self) -> Any:
        """Return a transactional pipeline."""

    def lrem(self, name: str, count: int, value: str) -> Any:
        """Remove a value from a Redis list."""


@dataclass(frozen=True)
class AdmissionOutcome:
    """Result of one atomic capacity-reservation attempt."""

    accepted: bool
    reason: str | None = None
    active_jobs: int = 0
    limiting_free_bytes: int | None = None


class CapacityLease:
    """Rollback owner for a reservation until broker dispatch succeeds."""

    def __init__(self, redis_client: RedisLike, reservation_id: str) -> None:
        """Bind a lease to an existing reservation.

        Args:
            redis_client: Synchronous redis-py compatible client.
            reservation_id: Identifier already stored by ``reserve_capacity``.
        """
        self._redis_client = redis_client
        self._reservation_id = reservation_id
        self._committed = False

    def __enter__(self) -> "CapacityLease":
        """Return this lease to the endpoint's dispatch block."""
        return self

    def commit(self) -> None:
        """Transfer reservation ownership to the dispatched worker."""
        self._committed = True

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Release an uncommitted reservation without masking the refusal."""
        if self._committed:
            return
        try:
            release_capacity(self._redis_client, self._reservation_id)
        except RedisError as release_error:
            logger.error(f"Failed to roll back capacity reservation {self._reservation_id}: {release_error}")


class ReservationHeartbeat:
    """Renew a running worker's capacity lease until task cleanup begins."""

    def __init__(
        self,
        redis_client: RedisLike,
        reservation_id: str,
        *,
        active_lease_seconds: int,
        heartbeat_seconds: int,
        clock: Callable[[], float] = time.time,
    ) -> None:
        """Configure a worker lease heartbeat.

        Args:
            redis_client: Synchronous redis-py compatible client.
            reservation_id: Identifier already reserved by the API.
            active_lease_seconds: Expiry horizon renewed by every heartbeat.
            heartbeat_seconds: Delay between background renewals.
            clock: Time source, injectable for deterministic tests.

        Raises:
            ValueError: If the timing configuration cannot safely renew.
        """
        _validate_lease_seconds(active_lease_seconds, heartbeat_seconds)
        self._redis_client = redis_client
        self._reservation_id = reservation_id
        self._active_lease_seconds = active_lease_seconds
        self._heartbeat_seconds = heartbeat_seconds
        self._clock = clock
        self._stop = Event()
        self._thread: Thread | None = None

    def start(self) -> None:
        """Renew immediately and start periodic renewal while the reservation is active."""
        if not renew_capacity(
            self._redis_client,
            self._reservation_id,
            active_lease_seconds=self._active_lease_seconds,
            clock=self._clock,
        ):
            msg = f"Capacity reservation {self._reservation_id} was absent when its worker started"
            logger.error(msg)
            raise ValueError(msg)
        thread = Thread(target=self._run, name=f"admission-heartbeat-{self._reservation_id}", daemon=True)
        thread.start()
        self._thread = thread

    def stop(self) -> None:
        """Stop periodic renewal without delaying task cleanup."""
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5)

    def _run(self) -> None:
        """Renew until stopped, released, or pruned after a prolonged outage."""
        while not self._stop.wait(self._heartbeat_seconds):
            try:
                renewed = renew_capacity(
                    self._redis_client,
                    self._reservation_id,
                    active_lease_seconds=self._active_lease_seconds,
                    clock=self._clock,
                )
            except RedisError as exc:
                logger.error(f"Failed to renew capacity reservation {self._reservation_id}: {exc}")
                continue
            if not renewed:
                return


def observe_free_bytes(paths: tuple[str, ...]) -> tuple[int, ...]:
    """Return a free-byte observation for every configured path.

    Args:
        paths: Existing directories whose backing volumes must have capacity.

    Returns:
        tuple[int, ...]: One free-byte observation per path.

    Raises:
        OSError: If a path or its filesystem cannot be inspected.
        ValueError: If no paths were supplied.
    """
    if not paths:
        msg = "At least one capacity path is required"
        logger.error(msg)
        raise ValueError(msg)

    return tuple(shutil.disk_usage(path).free for path in paths)


def reserve_capacity(
    redis_client: RedisLike,
    *,
    reservation_id: str,
    requested_bytes: int,
    free_bytes: tuple[int, ...],
    max_jobs: int,
    minimum_free_bytes: int,
    queued_grace_seconds: int = DEFAULT_QUEUED_GRACE_SECONDS,
    clock: Callable[[], float] = time.time,
) -> AdmissionOutcome:
    """Atomically prune expired leases and reserve a slot plus heuristic bytes.

    Redis ``WATCH`` makes the read/decision/write cycle retry when another API
    process changes the shared reservation state. Each stored value is an
    operational headroom estimate, not a bound on what the task can write.
    Observed filesystem free space is reduced by existing estimates and this
    request before the configured floor is applied.

    Args:
        redis_client: Synchronous redis-py compatible client.
        reservation_id: Stable identifier released by the worker.
        requested_bytes: Conservative workload-specific headroom estimate.
        free_bytes: Current free bytes observed at each configured volume path.
        max_jobs: Maximum number of accepted unfinished submissions.
        minimum_free_bytes: Floor retained after all reservations.
        queued_grace_seconds: Time an accepted, not-yet-started task may wait.
            The default is deliberately conservative, but operators must set it
            above their intended maximum broker backlog.
        clock: Time source, injectable for deterministic tests.

    Returns:
        AdmissionOutcome: Accepted, queue-full, or disk-full decision.

    Raises:
        ValueError: If an internal limit or request is invalid.
        redis.RedisError: If the reservation registry is unavailable.
    """
    _validate_request(
        reservation_id,
        requested_bytes,
        free_bytes,
        max_jobs,
        minimum_free_bytes,
        queued_grace_seconds,
    )
    while True:
        now = clock()
        deadline = now + queued_grace_seconds
        with redis_client.pipeline() as pipeline:
            try:
                pipeline.watch(ADMISSION_RESERVATIONS, ADMISSION_DEADLINES)
                raw_reservations = pipeline.hgetall(ADMISSION_RESERVATIONS)
                expired_ids = set(pipeline.zrangebyscore(ADMISSION_DEADLINES, "-inf", now))
                active_reservations = {key: value for key, value in raw_reservations.items() if key not in expired_ids}
                if reservation_id in active_reservations:
                    pipeline.unwatch()
                    return AdmissionOutcome(accepted=True, active_jobs=len(active_reservations))

                active_jobs = len(active_reservations)
                if active_jobs >= max_jobs:
                    _commit_prune(pipeline, expired_ids)
                    return AdmissionOutcome(accepted=False, reason="queue", active_jobs=active_jobs)

                reserved_bytes = sum(int(value) for value in active_reservations.values())
                limiting_free = min(free_bytes)
                if limiting_free - reserved_bytes - requested_bytes < minimum_free_bytes:
                    _commit_prune(pipeline, expired_ids)
                    return AdmissionOutcome(
                        accepted=False,
                        reason="disk",
                        active_jobs=active_jobs,
                        limiting_free_bytes=limiting_free,
                    )

                pipeline.multi()
                if expired_ids:
                    pipeline.hdel(ADMISSION_RESERVATIONS, *expired_ids)
                    pipeline.zrem(ADMISSION_DEADLINES, *expired_ids)
                pipeline.hset(ADMISSION_RESERVATIONS, reservation_id, requested_bytes)
                pipeline.zadd(ADMISSION_DEADLINES, {reservation_id: deadline})
                pipeline.execute()
                return AdmissionOutcome(accepted=True, active_jobs=active_jobs + 1)
            except WatchError:
                continue


def release_capacity(redis_client: RedisLike, reservation_id: str) -> None:
    """Idempotently release a worker's capacity reservation.

    Args:
        redis_client: Synchronous redis-py compatible client.
        reservation_id: Identifier originally passed to ``reserve_capacity``.
    """
    with redis_client.pipeline() as pipeline:
        pipeline.multi()
        pipeline.hdel(ADMISSION_RESERVATIONS, reservation_id)
        pipeline.zrem(ADMISSION_DEADLINES, reservation_id)
        pipeline.execute()


def renew_capacity(
    redis_client: RedisLike,
    reservation_id: str,
    *,
    active_lease_seconds: int,
    clock: Callable[[], float] = time.time,
) -> bool:
    """Atomically extend an unexpired worker lease without resurrecting it.

    The reservation hash and deadline index are one authoritative record. A
    missing half or a deadline at or before the current time is pruned in the
    same watched transaction and reported as absent.

    Args:
        redis_client: Synchronous redis-py compatible client.
        reservation_id: Identifier originally reserved by the API.
        active_lease_seconds: New expiry horizon from the current time.
        clock: Time source, injectable for deterministic tests.

    Returns:
        bool: True when an existing, unexpired reservation was renewed.

    Raises:
        ValueError: If the identifier or lease duration is invalid.
        redis.RedisError: If the reservation registry is unavailable.
    """
    if not reservation_id:
        msg = "Capacity reservation identifier must not be empty"
        logger.error(msg)
        raise ValueError(msg)
    if active_lease_seconds <= 0:
        msg = "Active capacity lease seconds must be positive"
        logger.error(msg)
        raise ValueError(msg)

    while True:
        now = clock()
        deadline = now + active_lease_seconds
        with redis_client.pipeline() as pipeline:
            try:
                pipeline.watch(ADMISSION_RESERVATIONS, ADMISSION_DEADLINES)
                exists = pipeline.hexists(ADMISSION_RESERVATIONS, reservation_id)
                current_deadline = pipeline.zscore(ADMISSION_DEADLINES, reservation_id)
                if not exists or current_deadline is None or current_deadline <= now:
                    if not exists and current_deadline is None:
                        pipeline.unwatch()
                        return False
                    pipeline.multi()
                    pipeline.hdel(ADMISSION_RESERVATIONS, reservation_id)
                    pipeline.zrem(ADMISSION_DEADLINES, reservation_id)
                    pipeline.execute()
                    return False
                pipeline.multi()
                pipeline.zadd(ADMISSION_DEADLINES, {reservation_id: deadline})
                pipeline.execute()
                return True
            except WatchError:
                continue


def cohort_reservation_bytes(member_bytes: int, base_bytes: int, member_factor: int) -> int:
    """Calculate the cohort work-area heuristic from explicit settings.

    The member factor reserves headroom for the owned member snapshot and the
    result archive. The fixed base covers reports and temporary output. This is
    admission planning only; it does not enforce a runtime write ceiling.

    Args:
        member_bytes: Aggregate bytes in cohort member archives.
        base_bytes: Fixed report and temporary-output headroom.
        member_factor: Number of member-sized working copies to reserve.

    Returns:
        int: Workload-specific reservation estimate.

    Raises:
        ValueError: If any input is negative.
    """
    if member_bytes < 0 or base_bytes < 0 or member_factor < 0:
        msg = "Cohort reservation budget inputs must not be negative"
        logger.error(msg)
        raise ValueError(msg)
    return base_bytes + member_factor * member_bytes


def release_worker_bookkeeping(
    redis_client: RedisLike,
    task_id: str,
    reservation_id: str,
    task_description: str,
    log: logging.Logger,
) -> None:
    """Best-effort release of display and admission bookkeeping.

    The two operations are isolated so failure of the legacy queue-position
    display cannot prevent release of authoritative admission capacity, and
    neither cleanup error can replace the task's actual result.

    Args:
        redis_client: Synchronous redis-py compatible client.
        task_id: Celery identifier stored in the display list.
        reservation_id: Endpoint identifier stored in the admission hash.
        task_description: Human-readable task identifier for log messages.
        log: Owning task module's logger.
    """
    try:
        redis_client.lrem("vntyper_job_queue", 0, task_id)
        log.info(f"Removed {task_description} from vntyper_job_queue")
    except RedisError as exc:
        log.error(f"Error removing {task_description} from vntyper_job_queue: {exc}")

    try:
        release_capacity(redis_client, reservation_id)
        log.info(f"Released capacity reservation for {reservation_id}")
    except RedisError as exc:
        log.error(f"Error releasing capacity reservation for {reservation_id}: {exc}")


def _validate_request(
    reservation_id: str,
    requested_bytes: int,
    free_bytes: tuple[int, ...],
    max_jobs: int,
    minimum_free_bytes: int,
    queued_grace_seconds: int,
) -> None:
    """Reject invalid internal admission inputs before touching Redis."""
    if not reservation_id:
        msg = "Capacity reservation identifier must not be empty"
    elif requested_bytes < 0:
        msg = "Capacity reservation bytes must not be negative"
    elif not free_bytes or any(value < 0 for value in free_bytes):
        msg = "Filesystem free-byte observations must be present and non-negative"
    elif max_jobs <= 0:
        msg = "Maximum queued jobs must be positive"
    elif minimum_free_bytes < 0:
        msg = "Minimum free disk bytes must not be negative"
    elif queued_grace_seconds <= 0:
        msg = "Queued capacity grace seconds must be positive"
    else:
        return
    logger.error(msg)
    raise ValueError(msg)


def _commit_prune(pipeline: Any, expired_ids: set[str]) -> None:
    """Commit expiration cleanup while preserving the watched decision."""
    if not expired_ids:
        pipeline.unwatch()
        return
    pipeline.multi()
    pipeline.hdel(ADMISSION_RESERVATIONS, *expired_ids)
    pipeline.zrem(ADMISSION_DEADLINES, *expired_ids)
    pipeline.execute()


def _validate_lease_seconds(active_lease_seconds: int, heartbeat_seconds: int) -> None:
    """Require renewal to occur before an active lease can expire."""
    if active_lease_seconds <= 0:
        msg = "Active capacity lease seconds must be positive"
    elif heartbeat_seconds <= 0 or heartbeat_seconds >= active_lease_seconds:
        msg = "Capacity heartbeat seconds must be positive and below the active lease"
    else:
        return
    logger.error(msg)
    raise ValueError(msg)
