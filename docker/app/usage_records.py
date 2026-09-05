"""The usage hash a worker writes when a job starts.

Both Celery tasks record the same four fields when they begin -- a hash of the
caller, a timestamp, the job identifier and ``"started"`` -- and the cohort task
adds its analysis type and cohort. `usage_statistics.py` reads these hashes back
to count jobs and unique users, so the field names and the hash input are part
of a stored contract, not a formatting choice: change them here, deliberately,
and nowhere else.

The record is declared with the key and value types redis-py publishes rather
than as ``dict[str, str]``. redis-py 8.1 types ``hset(mapping=...)`` as
``Mapping[FieldT, EncodableT]``, and ``Mapping`` is invariant in its key type,
so a mapping whose keys are inferred as plain ``str`` is rejected even though
``str`` is one arm of ``FieldT``. Naming the aliases from ``redis.typing`` keeps
this module tracking whatever the installed client accepts instead of restating
the union here.
"""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime, timezone

from redis.typing import EncodableT, FieldT

logger = logging.getLogger(__name__)

UsageRecord = dict[FieldT, EncodableT]


def client_hash(client_ip: str | None, user_agent: str | None) -> str:
    """Return the caller hash usage statistics count unique users by.

    The input is the address and agent joined by a hyphen, formatted with
    ``str()`` so a missing value hashes as the word ``None`` -- the cohort task
    has always handed ``None`` through and the stored hashes depend on it.

    Args:
        client_ip: The caller's address, or None when the task was given none.
        user_agent: The caller's user agent, or None when the task was given none.

    Returns:
        The hex SHA-256 digest of ``"<ip>-<agent>"``.
    """
    return hashlib.sha256(f"{client_ip}-{user_agent}".encode()).hexdigest()


def started_usage_record(
    job_id: str,
    *,
    client_ip: str | None,
    user_agent: str | None,
    analysis_type: str | None = None,
    cohort_id: str | None = None,
    now: datetime | None = None,
) -> UsageRecord:
    """Build the hash a task writes under ``usage:<job_id>`` when it starts.

    Args:
        job_id: The job or analysis identifier the hash is keyed by.
        client_ip: The caller's address, hashed rather than stored.
        user_agent: The caller's user agent, hashed rather than stored.
        analysis_type: The cohort task's analysis label; omitted when None.
        cohort_id: The cohort the analysis belongs to; omitted when None.
        now: The instant to stamp; defaults to the current aware UTC time.

    Returns:
        The mapping to hand to ``hset(mapping=...)``.
    """
    stamped = now if now is not None else datetime.now(timezone.utc)
    record: UsageRecord = {
        "user_hash": client_hash(client_ip, user_agent),
        "timestamp": stamped.isoformat(),
        "job_id": job_id,
        "status": "started",
    }
    if analysis_type is not None:
        record["analysis_type"] = analysis_type
    if cohort_id is not None:
        record["cohort_id"] = cohort_id
    return record
