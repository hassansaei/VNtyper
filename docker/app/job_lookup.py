"""Reading the job-to-task mapping back from the job database.

The API stores each job identifier as a key whose value is the Celery task id,
and three routes read it back. The job clients are built with
``decode_responses=True``, so at runtime the value is text -- but redis-py 8.1
declares ``Redis.get`` as ``bytes | str | None`` for every client, since the
decoding flag is not part of the type. `stored_task_id` narrows the answer once,
the way `usage_statistics._field` does for hash fields: a byte answer is
decoded rather than passed on, so a route never formats ``b'...'`` into a log
line or hands bytes to Celery, and an absent key is ``None`` rather than a
falsy value of some other type.
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from typing import Protocol

logger = logging.getLogger(__name__)


class JobStore(Protocol):
    """The one Redis operation the lookup needs."""

    def get(self, name: str) -> object:
        """Return whatever is stored under ``name``, or None when there is none.

        Declared as ``object`` on purpose: redis-py's own return annotation for
        ``GET`` differs between the versions this project has pinned, and the
        lookup narrows the answer itself rather than trusting either.
        """


class HashStore(Protocol):
    """The one Redis hash operation the timing lookup needs."""

    def hget(self, name: str, key: str) -> object:
        """Return one hash field, or None when it is absent."""


def stored_task_id(store: JobStore, job_id: str) -> str | None:
    """Return the Celery task id recorded for a job, as text.

    Args:
        store: The job database client.
        job_id: The canonical job identifier, already validated by the route.

    Returns:
        The task identifier, or None when no task is stored for the job.
    """
    value = store.get(job_id)
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, str):
        return value
    return None


def submitted_at_key(job_id: str) -> str:
    """Return the job-database key holding a job's submission instant."""
    return f"{job_id}:submitted_at"


def _iso_utc(value: object) -> str | None:
    """Render a stored instant as ISO 8601 UTC; naive datetimes are taken as UTC."""
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        try:
            value = datetime.fromisoformat(value)
        except ValueError:
            return None
    if not isinstance(value, datetime):
        return None
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc).isoformat()


def job_timing(job_store: JobStore, usage_store: HashStore, job_id: str, task_result: object) -> dict[str, str | None]:
    """Return the submitted/started/completed instants already recorded for a job.

    Args:
        job_store: The job database (holds ``<job_id>:submitted_at``).
        usage_store: The usage database (``usage:<job_id>`` ``timestamp`` is written
            when the worker starts the task).
        job_id: The canonical job identifier.
        task_result: The Celery ``AsyncResult``; its ``date_done`` is the finish time
            of a succeeded or failed task.

    Returns:
        ``submitted_at``, ``started_at`` and ``completed_at`` as ISO 8601 UTC strings,
        each None when it was not recorded (or has expired).
    """
    finished = getattr(task_result, "status", None) in ("SUCCESS", "FAILURE")
    return {
        "submitted_at": _iso_utc(job_store.get(submitted_at_key(job_id))),
        "started_at": _iso_utc(usage_store.hget(f"usage:{job_id}", "timestamp")),
        "completed_at": _iso_utc(getattr(task_result, "date_done", None)) if finished else None,
    }
