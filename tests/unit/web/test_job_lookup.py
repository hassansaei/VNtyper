"""The job-to-task lookup must hand routes a text task identifier.

The three job clients are built with `decode_responses=True`, so at runtime
`GET` returns text. redis-py 8.1 cannot express that in its types: `Redis.get`
is declared `bytes | str | None` for every client, so a route that formats the
value into a log line trips `str-bytes-safe`, and one that hands it to Celery
would pass bytes through unchanged if a raw client were ever wired in.
`job_lookup.stored_task_id` narrows the value once, the way
`usage_statistics._field` does for hash fields.
"""

import uuid
from datetime import datetime, timezone
from types import SimpleNamespace

import fakeredis
import pytest
from app.job_lookup import job_timing, stored_task_id, submitted_at_key

pytestmark = pytest.mark.unit

JOB_ID = str(uuid.uuid4())


def test_returns_the_stored_task_id_as_text() -> None:
    """A decoding client's answer is returned unchanged."""
    store = fakeredis.FakeRedis(decode_responses=True)
    store.set(JOB_ID, "task-1")

    assert stored_task_id(store, JOB_ID) == "task-1"


def test_decodes_a_raw_byte_response() -> None:
    """A non-decoding client's bytes are decoded rather than passed on."""
    store = fakeredis.FakeRedis(decode_responses=False)
    store.set(JOB_ID, b"task-1")

    assert stored_task_id(store, JOB_ID) == "task-1"


def test_returns_none_when_no_task_is_stored() -> None:
    """An unknown job has no task, and the caller sees `None`, not a falsy string."""
    store = fakeredis.FakeRedis(decode_responses=True)

    assert stored_task_id(store, JOB_ID) is None


def test_job_status_hands_celery_a_text_task_id_even_from_a_raw_client(
    monkeypatch: pytest.MonkeyPatch, client, web_app
) -> None:
    """The status route resolves the task by text identifier whatever the client decodes.

    Args:
        monkeypatch: Standard pytest fixture.
        client: TestClient over the patched app.
        web_app: The patched `app.main` module.
    """
    raw_store = fakeredis.FakeRedis(decode_responses=False)
    raw_store.set(JOB_ID, b"task-1")
    monkeypatch.setattr(web_app, "redis_client", raw_store)
    seen: list[object] = []

    def _async_result(task_id):
        seen.append(task_id)
        return SimpleNamespace(status="PENDING", info=None)

    monkeypatch.setattr(web_app, "AsyncResult", _async_result)

    response = client.get(f"/job-status/{JOB_ID}/")

    assert response.status_code == 200
    assert response.json() == {
        "job_id": JOB_ID,
        "status": "pending",
        "error": None,
        "submitted_at": None,
        "started_at": None,
        "completed_at": None,
    }
    assert seen == ["task-1"]


def test_job_timing_reports_the_three_recorded_instants_in_utc() -> None:
    """Submission comes from the job db, start from the usage record, finish from Celery."""
    store = fakeredis.FakeRedis(decode_responses=True)
    store.set(submitted_at_key(JOB_ID), "2026-11-10T09:00:00+00:00")
    store.hset(f"usage:{JOB_ID}", mapping={"timestamp": "2026-11-10T10:00:05+01:00"})
    done = SimpleNamespace(status="SUCCESS", date_done=datetime(2026, 11, 10, 9, 3, 0))  # noqa: DTZ001 - naive is read as UTC

    assert job_timing(store, store, JOB_ID, done) == {
        "submitted_at": "2026-11-10T09:00:00+00:00",
        "started_at": "2026-11-10T09:00:05+00:00",
        "completed_at": "2026-11-10T09:03:00+00:00",
    }


def test_job_timing_leaves_unrecorded_and_unfinished_instants_empty() -> None:
    """A running task has no completion time even if Celery carries a stale date_done."""
    store = fakeredis.FakeRedis(decode_responses=True)
    running = SimpleNamespace(status="STARTED", date_done=datetime(2026, 11, 10, 9, 3, 0, tzinfo=timezone.utc))

    assert job_timing(store, store, JOB_ID, running) == {
        "submitted_at": None,
        "started_at": None,
        "completed_at": None,
    }


def test_job_status_exposes_the_timing_fields(monkeypatch: pytest.MonkeyPatch, client, web_app, fake_redis) -> None:
    """The status route adds the instants without changing the existing fields."""
    fake_redis.set(JOB_ID, "task-1")
    fake_redis.set(submitted_at_key(JOB_ID), "2026-11-10T09:00:00+00:00")
    monkeypatch.setattr(
        web_app,
        "AsyncResult",
        lambda task_id: SimpleNamespace(status="SUCCESS", info=None, date_done="2026-11-10T09:03:00+00:00"),
    )

    body = client.get(f"/job-status/{JOB_ID}/").json()

    assert body["status"] == "completed"
    assert body["submitted_at"] == "2026-11-10T09:00:00+00:00"
    assert body["completed_at"] == "2026-11-10T09:03:00+00:00"
