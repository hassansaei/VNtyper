"""The started-job usage record the worker writes must be one Redis-typed mapping.

redis-py 8.1 declares `hset(mapping=...)` as `Mapping[FieldT, EncodableT]`, and
`Mapping` is invariant in its key type, so an inferred `dict[str, str]` is
rejected by mypy even though every key is a `str`. `usage_records` is the one
place that names the type Redis accepts, so the two worker tasks share one
record shape instead of each carrying an ad-hoc dict.
"""

import hashlib
from datetime import datetime, timezone

import fakeredis
import pytest
from app.usage_records import client_hash, started_usage_record

pytestmark = pytest.mark.unit

FIXED_NOW = datetime(2026, 9, 5, 12, 30, 45, tzinfo=timezone.utc)


def test_client_hash_is_the_sha256_of_ip_and_agent_joined_by_a_hyphen() -> None:
    """The hash is stable across releases: it is what usage statistics count."""
    assert client_hash("127.0.0.1", "pytest") == hashlib.sha256(b"127.0.0.1-pytest").hexdigest()


def test_client_hash_keeps_hashing_missing_values_as_the_word_none() -> None:
    """The cohort task hands `None` through; the historical hash input is `"None-None"`."""
    assert client_hash(None, None) == hashlib.sha256(b"None-None").hexdigest()


def test_started_record_carries_exactly_hash_timestamp_job_and_status() -> None:
    """A pipeline job's record has four fields and nothing else."""
    record = started_usage_record("job-1", client_ip="127.0.0.1", user_agent="pytest", now=FIXED_NOW)

    assert record == {
        "user_hash": hashlib.sha256(b"127.0.0.1-pytest").hexdigest(),
        "timestamp": "2026-09-05T12:30:45+00:00",
        "job_id": "job-1",
        "status": "started",
    }


def test_started_record_adds_cohort_fields_only_when_they_are_given() -> None:
    """The cohort task names its analysis type and cohort; the pipeline task does not."""
    cohort_record = started_usage_record(
        "analysis",
        client_ip=None,
        user_agent=None,
        now=FIXED_NOW,
        analysis_type="cohort_analysis",
        cohort_id="cohort-1",
    )

    assert cohort_record["analysis_type"] == "cohort_analysis"
    assert cohort_record["cohort_id"] == "cohort-1"
    assert set(cohort_record) == {"user_hash", "timestamp", "job_id", "status", "analysis_type", "cohort_id"}


def test_started_record_timestamp_defaults_to_an_aware_utc_now() -> None:
    """Without `now`, the timestamp is an aware UTC instant taken at call time."""
    before = datetime.now(timezone.utc)
    record = started_usage_record("job-1", client_ip="127.0.0.1", user_agent="pytest")
    after = datetime.now(timezone.utc)

    stamped = datetime.fromisoformat(str(record["timestamp"]))
    assert stamped.tzinfo is not None
    assert before <= stamped <= after


def test_started_record_round_trips_through_a_real_redis_hash() -> None:
    """The record is what `hset(mapping=...)` accepts, at runtime and for mypy."""
    store = fakeredis.FakeRedis(decode_responses=True)
    record = started_usage_record("job-1", client_ip="127.0.0.1", user_agent="pytest", now=FIXED_NOW)

    store.hset("usage:job-1", mapping=record)

    assert store.hgetall("usage:job-1") == record
