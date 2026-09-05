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
