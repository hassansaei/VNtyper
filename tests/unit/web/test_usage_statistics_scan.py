"""Usage statistics must scan and batch-read the shared Redis database."""

from collections.abc import Mapping

import pytest
from app.usage_statistics import RedisScalar, aggregate_usage_statistics

pytestmark = pytest.mark.unit


class _PagingPipeline:
    """Small Redis-pipeline stand-in that executes queued hash reads together."""

    def __init__(self, store: "_PagingUsageStore") -> None:
        self._store = store
        self._keys: list[RedisScalar] = []

    def hgetall(self, key: RedisScalar) -> "_PagingPipeline":
        """Queue one hash read and preserve Redis pipeline chaining."""
        self._keys.append(key)
        return self

    def execute(self) -> list[Mapping[RedisScalar, RedisScalar]]:
        """Return queued hashes in command order."""
        self._store.pipeline_executions += 1
        self._store.pipeline_batch_sizes.append(len(self._keys))
        return [self._store.record(key) for key in self._keys]


class _PagingUsageStore:
    """Deterministic, cursor-paged Redis stand-in with mixed response types."""

    def __init__(self) -> None:
        self.scan_calls: list[int] = []
        self.pipeline_executions = 0
        self.pipeline_batch_sizes: list[int] = []
        self._pages: dict[int, tuple[int | bytes, list[RedisScalar]]] = {
            0: (b"5", [b"usage:job-1", "usage:job-2"]),
            5: (7, [b"usage:job-2", b"usage:job-3"]),
            7: (b"0", []),
        }
        self._records: dict[str, Mapping[RedisScalar, RedisScalar]] = {
            "usage:job-1": {b"user_hash": b"u1", b"status": b"completed"},
            "usage:job-2": {"user_hash": "u2", "status": "failed"},
            "usage:job-3": {b"user_hash": b"u1"},
        }

    def scan(self, cursor: int, *, match: str, count: int) -> tuple[int | bytes, list[RedisScalar]]:
        """Return one page and verify the bounded usage-prefix query."""
        assert match == "usage:*"
        assert 0 < count <= 1_000
        self.scan_calls.append(cursor)
        return self._pages[cursor]

    def pipeline(self, *, transaction: bool) -> _PagingPipeline:
        """Return a non-transactional batch reader."""
        assert transaction is False
        return _PagingPipeline(self)

    def hgetall(self, _key: RedisScalar) -> Mapping[RedisScalar, RedisScalar]:
        """Reject the N+1 direct-read pattern."""
        raise AssertionError("usage hashes must be fetched through a pipeline")

    def record(self, key: RedisScalar) -> Mapping[RedisScalar, RedisScalar]:
        """Look up a record after normalising Redis key response types."""
        normalized = key.decode("utf-8") if isinstance(key, bytes) else key
        return self._records[normalized]


class _RepeatedCursorStore(_PagingUsageStore):
    """Broken Redis stand-in that never advances its non-zero cursor."""

    def scan(self, cursor: int, *, match: str, count: int) -> tuple[int | bytes, list[RedisScalar]]:
        """Repeat cursor five and fail loudly if the caller keeps looping."""
        assert match == "usage:*"
        assert 0 < count <= 1_000
        self.scan_calls.append(cursor)
        if len(self.scan_calls) > 3:
            raise AssertionError("usage aggregation did not stop on a repeated cursor")
        return b"5", []


class _StatusEdgeStore(_PagingUsageStore):
    """Usage hashes that distinguish an absent status from an empty one."""

    def __init__(self) -> None:
        super().__init__()
        self._pages = {0: (0, ["usage:missing", "usage:empty"])}
        self._records = {
            "usage:missing": {"user_hash": "u1"},
            "usage:empty": {"user_hash": "u2", "status": ""},
        }


class _LargeUsageStore(_PagingUsageStore):
    """One SCAN page large enough to require multiple fetch pipelines."""

    def __init__(self) -> None:
        super().__init__()
        text_keys = [f"usage:job-{index}" for index in range(205)]
        keys: list[RedisScalar] = list(text_keys)
        pages: dict[int, tuple[int | bytes, list[RedisScalar]]] = {0: (0, keys)}
        self._pages = pages
        self._records = {
            key: {"user_hash": f"user-{index}", "status": "completed"} for index, key in enumerate(text_keys)
        }


class _ScanFailureStore(_PagingUsageStore):
    """Redis stand-in whose incremental key discovery fails."""

    def scan(self, cursor: int, *, match: str, count: int) -> tuple[int | bytes, list[RedisScalar]]:
        """Raise the connection failure that the endpoint must not hide."""
        raise ConnectionError("usage Redis unavailable")


def test_usage_statistics_are_computed_without_calling_keys(client, fake_redis, monkeypatch) -> None:
    """The route preserves its aggregate without a blocking ``KEYS`` request.

    Args:
        client: In-process FastAPI test client.
        fake_redis: In-process Redis stand-in used by the application.
        monkeypatch: Standard pytest fixture.
    """
    fake_redis.hset("usage:job-1", mapping={"user_hash": "u1", "status": "completed"})
    fake_redis.hset("usage:job-2", mapping={"user_hash": "u2", "status": "failed"})
    fake_redis.hset("usage:job-3", mapping={"user_hash": "u1", "status": "completed"})

    def _refuse_keys(*_args, **_kwargs):
        raise AssertionError("KEYS blocks the shared Redis instance; use SCAN")

    monkeypatch.setattr(fake_redis, "keys", _refuse_keys)

    response = client.get("/usage-statistics/")

    assert response.status_code == 200
    assert response.json() == {
        "total_jobs": 3,
        "unique_users": 2,
        "job_statuses": {"completed": 2, "failed": 1},
    }


def test_usage_statistics_with_no_recorded_usage_are_all_zeroes(client, fake_redis, monkeypatch) -> None:
    """An empty incremental scan returns the unchanged zero-valued response.

    Args:
        client: In-process FastAPI test client.
        fake_redis: In-process Redis stand-in used by the application.
        monkeypatch: Standard pytest fixture.
    """

    def _refuse_keys(*_args, **_kwargs):
        raise AssertionError("KEYS blocks the shared Redis instance; use SCAN")

    monkeypatch.setattr(fake_redis, "keys", _refuse_keys)

    response = client.get("/usage-statistics/")

    assert response.status_code == 200
    assert response.json() == {"total_jobs": 0, "unique_users": 0, "job_statuses": {}}


def test_usage_statistics_scan_pages_are_deduplicated_and_batch_fetched() -> None:
    """Mixed cursor/key types preserve counts without direct per-key reads."""
    store = _PagingUsageStore()

    result = aggregate_usage_statistics(store)

    assert result == {
        "total_jobs": 3,
        "unique_users": 2,
        "job_statuses": {"completed": 1, "failed": 1, "unknown": 1},
    }
    assert store.scan_calls == [0, 5, 7]
    assert store.pipeline_executions == 1
    assert store.pipeline_batch_sizes == [3]
    assert list(result["job_statuses"]) == ["completed", "failed", "unknown"]


def test_usage_statistics_reject_a_repeated_nonzero_scan_cursor() -> None:
    """A broken SCAN response cannot trap a request in an infinite loop."""
    store = _RepeatedCursorStore()

    with pytest.raises(RuntimeError, match="Redis usage SCAN repeated cursor 5"):
        aggregate_usage_statistics(store)

    assert store.scan_calls == [0, 5]


def test_usage_statistics_preserve_an_empty_status_value() -> None:
    """Only a missing status maps to ``unknown``; an empty value stays empty."""
    result = aggregate_usage_statistics(_StatusEdgeStore())

    assert result == {
        "total_jobs": 2,
        "unique_users": 2,
        "job_statuses": {"unknown": 1, "": 1},
    }


def test_usage_statistics_bound_each_pipeline_fetch_batch() -> None:
    """A large SCAN page is split into bounded Redis pipeline executions."""
    store = _LargeUsageStore()

    result = aggregate_usage_statistics(store)

    assert result == {
        "total_jobs": 205,
        "unique_users": 205,
        "job_statuses": {"completed": 205},
    }
    assert sum(store.pipeline_batch_sizes) == 205
    assert len(store.pipeline_batch_sizes) > 1
    assert max(store.pipeline_batch_sizes) <= 100


def test_usage_statistics_propagate_scan_failures() -> None:
    """A Redis discovery failure is not converted into misleading zeroes."""
    store = _ScanFailureStore()

    with pytest.raises(ConnectionError, match="usage Redis unavailable"):
        aggregate_usage_statistics(store)

    assert store.pipeline_executions == 0
