"""Redis-backed aggregation for the usage-statistics endpoint."""

import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from typing import Protocol, TypedDict

logger = logging.getLogger(__name__)


RedisScalar = str | bytes
RedisCursor = int | str | bytes
_SCAN_COUNT = 1_000
_FETCH_BATCH_SIZE = 100


class CumulativeStatistics(TypedDict):
    """Cumulative usage statistics since counter inception."""

    total_jobs: int
    unique_users: int
    since: str | None
    job_statuses: dict[str, int]


class UsageStatistics(TypedDict):
    """Aggregate values returned by the usage-statistics endpoint."""

    total_jobs: int
    unique_users: int
    job_statuses: dict[str, int]
    cumulative: CumulativeStatistics


class UsagePipeline(Protocol):
    """Redis pipeline operations used for batched usage-hash reads."""

    def hgetall(self, key: RedisScalar) -> "UsagePipeline":
        """Queue one hash read."""

    def execute(self) -> list[Mapping[RedisScalar, RedisScalar]]:
        """Execute queued reads and return their results in command order."""


class UsageStore(Protocol):
    """Redis operations required to aggregate stored usage hashes."""

    def scan(
        self,
        cursor: int,
        *,
        match: str,
        count: int,
    ) -> tuple[RedisCursor, list[RedisScalar]]:
        """Return one cursor page of usage keys."""

    def pipeline(self, *, transaction: bool) -> UsagePipeline:
        """Return a command pipeline."""


def _key_identity(key: RedisScalar) -> str:
    """Normalise text and byte responses for duplicate detection."""
    return key.decode("utf-8") if isinstance(key, bytes) else key


def _field(data: Mapping[RedisScalar, RedisScalar], name: str) -> str | None:
    """Read and normalise one hash field from decoded or raw Redis data."""
    value = data.get(name)
    if value is None:
        value = data.get(name.encode("utf-8"))
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def _fetch_hashes(store: UsageStore, keys: Sequence[RedisScalar]) -> list[Mapping[RedisScalar, RedisScalar]]:
    """Fetch a bounded key batch in one Redis round trip."""
    pipeline = store.pipeline(transaction=False)
    for key in keys:
        pipeline.hgetall(key)
    return pipeline.execute()


def aggregate_usage_statistics(store: UsageStore) -> UsageStatistics:
    """Aggregate usage hashes without blocking Redis key discovery.

    Args:
        store: Redis-compatible usage database.

    Returns:
        Counts of jobs, distinct users and job statuses.
    """
    unique_users: set[str | None] = set()
    job_statuses: Counter[str] = Counter()
    seen_keys: set[str] = set()
    seen_cursors: set[int] = set()
    pending_keys: list[RedisScalar] = []
    cursor = 0

    def consume_pending() -> None:
        """Add one pipelined batch to the aggregate."""
        for data in _fetch_hashes(store, pending_keys):
            unique_users.add(_field(data, "user_hash"))
            status = _field(data, "status")
            if status is None:
                status = "unknown"
            job_statuses[status] += 1
        pending_keys.clear()

    while True:
        next_cursor, keys = store.scan(cursor, match="usage:*", count=_SCAN_COUNT)
        for key in keys:
            identity = _key_identity(key)
            if identity.startswith(("usage:cumulative:", "usage:counted:")):
                continue
            if identity in seen_keys:
                continue
            seen_keys.add(identity)
            pending_keys.append(key)
            if len(pending_keys) == _FETCH_BATCH_SIZE:
                consume_pending()

        cursor = int(next_cursor)
        if cursor == 0:
            break
        if cursor in seen_cursors:
            msg = f"Redis usage SCAN repeated cursor {cursor}; refusing to loop indefinitely."
            logger.error(msg)
            raise RuntimeError(msg)
        seen_cursors.add(cursor)

    if pending_keys:
        consume_pending()

    # Read cumulative counters
    cumulative_jobs = 0
    cumulative_users = 0
    cumulative_since: str | None = None
    cumulative_completed = 0
    cumulative_failed = 0

    getter = getattr(store, "get", None)
    if getter is not None:
        try:
            val = getter("usage:cumulative:jobs")
            cumulative_jobs = int(val) if val is not None else 0
        except (ValueError, TypeError):
            cumulative_jobs = 0

        try:
            val = getter("usage:cumulative:since")
            cumulative_since = (
                val.decode("utf-8") if isinstance(val, bytes) else (str(val) if val is not None else None)
            )
        except (UnicodeDecodeError, AttributeError, ValueError, TypeError):
            cumulative_since = None

        try:
            val = getter("usage:cumulative:completed")
            cumulative_completed = int(val) if val is not None else 0
        except (ValueError, TypeError):
            cumulative_completed = 0

        try:
            val = getter("usage:cumulative:failed")
            cumulative_failed = int(val) if val is not None else 0
        except (ValueError, TypeError):
            cumulative_failed = 0

    pfcount_fn = getattr(store, "pfcount", None)
    if pfcount_fn is not None:
        try:
            cumulative_users = int(pfcount_fn("usage:cumulative:users") or 0)
        except (ValueError, TypeError, AttributeError):
            cumulative_users = 0

    # Seed lower bound from active rolling window if unseeded and cumulative is 0
    setter = getattr(store, "set", None)
    if cumulative_jobs == 0 and len(seen_keys) > 0 and setter is not None:
        try:
            was_unseeded = store.set("usage:cumulative:seeded", "1", nx=True)
            if was_unseeded:
                store.set("usage:cumulative:jobs", str(len(seen_keys)))
                cumulative_jobs = len(seen_keys)
                if job_statuses.get("completed"):
                    store.set("usage:cumulative:completed", str(job_statuses["completed"]))
                    cumulative_completed = job_statuses["completed"]
                if job_statuses.get("failed"):
                    store.set("usage:cumulative:failed", str(job_statuses["failed"]))
                    cumulative_failed = job_statuses["failed"]
                pfadd_fn = getattr(store, "pfadd", None)
                if pfadd_fn is not None:
                    for u in unique_users:
                        if u:
                            pfadd_fn("usage:cumulative:users", u)
                    if pfcount_fn is not None:
                        cumulative_users = int(pfcount_fn("usage:cumulative:users") or 0)
        except (ValueError, TypeError, AttributeError, RuntimeError, OSError) as seed_err:
            logger.warning("Could not seed cumulative usage statistics: %s", seed_err)

    if cumulative_since is None and setter is not None:
        now_iso = datetime.now(timezone.utc).isoformat()
        try:
            setnx_fn = getattr(store, "setnx", None)
            if setnx_fn is not None:
                setnx_fn("usage:cumulative:since", now_iso)
            else:
                store.set("usage:cumulative:since", now_iso, nx=True)
            cumulative_since = now_iso
        except (ValueError, TypeError, AttributeError, RuntimeError, OSError):
            pass

    return {
        "total_jobs": len(seen_keys),
        "unique_users": len(unique_users),
        "job_statuses": dict(job_statuses),
        "cumulative": {
            "total_jobs": cumulative_jobs,
            "unique_users": cumulative_users,
            "since": cumulative_since,
            "job_statuses": {
                "completed": cumulative_completed,
                "failed": cumulative_failed,
            },
        },
    }
