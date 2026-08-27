"""Redis-backed aggregation for the usage-statistics endpoint."""

import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Protocol, TypedDict

logger = logging.getLogger(__name__)


RedisScalar = str | bytes
RedisCursor = int | str | bytes
_SCAN_COUNT = 1_000
_FETCH_BATCH_SIZE = 100


class UsageStatistics(TypedDict):
    """Aggregate values returned by the usage-statistics endpoint."""

    total_jobs: int
    unique_users: int
    job_statuses: dict[str, int]


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

    return {
        "total_jobs": len(seen_keys),
        "unique_users": len(unique_users),
        "job_statuses": dict(job_statuses),
    }
