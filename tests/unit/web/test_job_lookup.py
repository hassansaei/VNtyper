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
from types import SimpleNamespace

import fakeredis
import pytest
from app.job_lookup import stored_task_id

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
    assert response.json() == {"job_id": JOB_ID, "status": "pending", "error": None}
    assert seen == ["task-1"]
