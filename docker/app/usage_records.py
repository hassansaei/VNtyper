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


def record_cumulative_job_started(
    store: any,
    job_id: str,
    user_hash: str | None,
    now: datetime | None = None,
    ttl_seconds: int = 30 * 86400,
) -> bool:
    """Record job start in cumulative counters, guarded against task retries.

    Returns True if this was the first start attempt counted, False if already counted.
    """
    guard_key = f"usage:counted:started:{job_id}"
    setter = getattr(store, "set", None)
    if setter is not None:
        try:
            was_set = store.set(guard_key, "1", nx=True, ex=ttl_seconds)
            if not was_set:
                return False
        except TypeError:
            setnx_fn = getattr(store, "setnx", None)
            if setnx_fn is not None and not setnx_fn(guard_key, "1"):
                return False
            expire_fn = getattr(store, "expire", None)
            if expire_fn is not None:
                expire_fn(guard_key, ttl_seconds)

    stamped = now if now is not None else datetime.now(timezone.utc)
    since_val = stamped.isoformat()

    # SETNX usage:cumulative:since
    setnx = getattr(store, "setnx", None)
    if setnx is not None:
        setnx("usage:cumulative:since", since_val)
    elif setter is not None:
        store.set("usage:cumulative:since", since_val, nx=True)

    # INCR usage:cumulative:jobs
    incr = getattr(store, "incr", None)
    if incr is not None:
        incr("usage:cumulative:jobs")

    # PFADD usage:cumulative:users
    pfadd = getattr(store, "pfadd", None)
    if pfadd is not None and user_hash:
        pfadd("usage:cumulative:users", user_hash)

    return True


def record_cumulative_job_completed(
    store: any,
    job_id: str,
    ttl_seconds: int = 30 * 86400,
) -> bool:
    """Record job completion in cumulative counters, guarded against retries."""
    guard_key = f"usage:counted:completed:{job_id}"
    setter = getattr(store, "set", None)
    if setter is not None:
        try:
            was_set = store.set(guard_key, "1", nx=True, ex=ttl_seconds)
            if not was_set:
                return False
        except TypeError:
            setnx_fn = getattr(store, "setnx", None)
            if setnx_fn is not None and not setnx_fn(guard_key, "1"):
                return False
            expire_fn = getattr(store, "expire", None)
            if expire_fn is not None:
                expire_fn(guard_key, ttl_seconds)

    incr = getattr(store, "incr", None)
    if incr is not None:
        incr("usage:cumulative:completed")
    return True


def record_cumulative_job_failed(
    store: any,
    job_id: str,
    ttl_seconds: int = 30 * 86400,
) -> bool:
    """Record job failure in cumulative counters, guarded against duplicate failure reporting."""
    guard_key = f"usage:counted:failed:{job_id}"
    setter = getattr(store, "set", None)
    if setter is not None:
        try:
            was_set = store.set(guard_key, "1", nx=True, ex=ttl_seconds)
            if not was_set:
                return False
        except TypeError:
            setnx_fn = getattr(store, "setnx", None)
            if setnx_fn is not None and not setnx_fn(guard_key, "1"):
                return False
            expire_fn = getattr(store, "expire", None)
            if expire_fn is not None:
                expire_fn(guard_key, ttl_seconds)

    incr = getattr(store, "incr", None)
    if incr is not None:
        incr("usage:cumulative:failed")
    return True
