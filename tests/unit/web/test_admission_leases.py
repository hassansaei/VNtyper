"""Expiry, renewal and heuristic byte budgets for admission reservations."""

import logging
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Event, Lock

import pytest

pytestmark = pytest.mark.unit


def _reserve(admission, redis_client, reservation_id: str, *, now: float, requested_bytes: int = 100):
    """Reserve with small deterministic limits for lease tests."""
    return admission.reserve_capacity(
        redis_client,
        reservation_id=reservation_id,
        requested_bytes=requested_bytes,
        free_bytes=(300,),
        max_jobs=1,
        minimum_free_bytes=100,
        queued_grace_seconds=60,
        clock=lambda: now,
    )


def test_expired_queued_reservation_is_pruned_during_the_next_reserve(fake_redis) -> None:
    """A broker-accepted task that never starts eventually stops consuming capacity."""
    from app import admission

    assert _reserve(admission, fake_redis, "abandoned", now=100).accepted

    replacement = _reserve(admission, fake_redis, "replacement", now=160)

    assert replacement.accepted
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {"replacement": "100"}
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "abandoned") is None


def test_reservation_resamples_the_clock_after_a_watch_retry_crosses_queued_grace(fake_redis) -> None:
    """A contended acceptance starts its grace period from the successful retry."""
    from app import admission
    from redis.exceptions import WatchError

    real_pipeline = fake_redis.pipeline
    first_execute = True

    class RetryOncePipeline:
        def __init__(self):
            self._pipeline = real_pipeline()

        def __enter__(self):
            self._pipeline.__enter__()
            return self

        def __exit__(self, *args):
            return self._pipeline.__exit__(*args)

        def __getattr__(self, name):
            return getattr(self._pipeline, name)

        def execute(self):
            nonlocal first_execute
            if first_execute:
                first_execute = False
                raise WatchError
            return self._pipeline.execute()

    fake_redis.pipeline = RetryOncePipeline
    clock_values = iter((0.0, 2.0))

    outcome = admission.reserve_capacity(
        fake_redis,
        reservation_id="contended",
        requested_bytes=100,
        free_bytes=(300,),
        max_jobs=1,
        minimum_free_bytes=100,
        queued_grace_seconds=1,
        clock=lambda: next(clock_values),
    )

    assert outcome.accepted
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "contended") == 3


def test_active_worker_renewal_prevents_expiry(fake_redis) -> None:
    """Heartbeat renewal keeps valid long-running work authoritative."""
    from app import admission

    assert _reserve(admission, fake_redis, "running", now=100).accepted
    assert admission.renew_capacity(fake_redis, "running", active_lease_seconds=30, clock=lambda: 150)

    blocked = _reserve(admission, fake_redis, "other", now=179)

    assert not blocked.accepted
    assert blocked.reason == "queue"
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "running") == 180


def test_renewal_prunes_and_rejects_a_lease_at_its_deadline_without_new_admission(fake_redis) -> None:
    """An overdue queued lease cannot revive itself when no new request triggered pruning."""
    from app import admission

    assert _reserve(admission, fake_redis, "overdue", now=100).accepted

    assert not admission.renew_capacity(fake_redis, "overdue", active_lease_seconds=30, clock=lambda: 160)
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "overdue") is None


def test_expired_renewal_retries_when_an_active_renewal_changes_the_watched_deadline(fake_redis) -> None:
    """A stale prune decision cannot delete a lease renewed after that decision's read."""
    from app import admission

    assert _reserve(admission, fake_redis, "racing", now=100).accepted
    real_pipeline = fake_redis.pipeline
    expired_read = Event()
    active_written = Event()
    read_lock = Lock()
    reads = 0
    watched_keys: list[tuple[str, ...]] = []

    class SynchronizedPipeline:
        def __init__(self):
            self._pipeline = real_pipeline()

        def __enter__(self):
            self._pipeline.__enter__()
            return self

        def __exit__(self, *args):
            return self._pipeline.__exit__(*args)

        def __getattr__(self, name):
            return getattr(self._pipeline, name)

        def watch(self, *keys):
            with read_lock:
                watched_keys.append(keys)
            return self._pipeline.watch(*keys)

        def zscore(self, key, member):
            nonlocal reads
            result = self._pipeline.zscore(key, member)
            with read_lock:
                reads += 1
                first_read = reads == 1
            if first_read:
                expired_read.set()
                assert active_written.wait(timeout=5)
            return result

    fake_redis.pipeline = SynchronizedPipeline

    with ThreadPoolExecutor(max_workers=1) as pool:
        expired = pool.submit(
            admission.renew_capacity,
            fake_redis,
            "racing",
            active_lease_seconds=30,
            clock=lambda: 161,
        )
        assert expired_read.wait(timeout=5)
        assert admission.renew_capacity(fake_redis, "racing", active_lease_seconds=30, clock=lambda: 150)
        active_written.set()
        assert expired.result(timeout=5)

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {"racing": "100"}
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "racing") == 191
    assert watched_keys
    assert all(admission.ADMISSION_DEADLINES in keys for keys in watched_keys)


def test_renewal_resamples_the_clock_after_a_watch_retry_crosses_the_deadline(fake_redis) -> None:
    """Contention cannot make renewal rely on a pre-expiry timestamp after expiry."""
    from app import admission
    from redis.exceptions import WatchError

    assert _reserve(admission, fake_redis, "contended", now=100).accepted
    real_pipeline = fake_redis.pipeline
    first_execute = True

    class RetryOncePipeline:
        def __init__(self):
            self._pipeline = real_pipeline()

        def __enter__(self):
            self._pipeline.__enter__()
            return self

        def __exit__(self, *args):
            return self._pipeline.__exit__(*args)

        def __getattr__(self, name):
            return getattr(self._pipeline, name)

        def execute(self):
            nonlocal first_execute
            if first_execute:
                first_execute = False
                raise WatchError
            return self._pipeline.execute()

    fake_redis.pipeline = RetryOncePipeline
    clock_values = iter((159.0, 161.0))

    assert not admission.renew_capacity(
        fake_redis,
        "contended",
        active_lease_seconds=30,
        clock=lambda: next(clock_values),
    )
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "contended") is None


def test_worker_heartbeat_renews_periodically_with_a_deterministic_clock(
    fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The background owner performs renewal beyond its immediate startup claim."""
    from app import admission

    assert _reserve(admission, fake_redis, "running", now=90).accepted
    clock_values = iter((100.0, 110.0))

    class WaitOnce:
        def __init__(self):
            self.calls = 0

        def wait(self, _seconds):
            self.calls += 1
            return self.calls > 1

        def set(self):
            return None

    class ImmediateThread:
        def __init__(self, *, target, **_kwargs):
            self.target = target

        def start(self):
            self.target()

        def join(self, timeout=None):
            return None

    monkeypatch.setattr(admission, "Event", WaitOnce)
    monkeypatch.setattr(admission, "Thread", ImmediateThread)
    heartbeat = admission.ReservationHeartbeat(
        fake_redis,
        "running",
        active_lease_seconds=30,
        heartbeat_seconds=10,
        clock=lambda: next(clock_values),
    )

    heartbeat.start()

    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, "running") == 140


def test_worker_refuses_to_run_after_its_reservation_was_pruned(fake_redis) -> None:
    """A late broker delivery cannot execute outside the authoritative cap."""
    from app import admission

    heartbeat = admission.ReservationHeartbeat(
        fake_redis,
        "already-pruned",
        active_lease_seconds=30,
        heartbeat_seconds=10,
        clock=lambda: 100,
    )

    with pytest.raises(ValueError, match="was absent when its worker started"):
        heartbeat.start()


def test_heartbeat_stop_is_safe_after_thread_start_fails(fake_redis, monkeypatch: pytest.MonkeyPatch) -> None:
    """A thread bootstrap error remains primary when task cleanup stops the heartbeat."""
    from app import admission

    assert _reserve(admission, fake_redis, "running", now=100).accepted

    class FailingThread:
        def __init__(self, **_kwargs):
            pass

        def start(self):
            raise RuntimeError("heartbeat thread failed to start")

        def join(self, timeout=None):
            pytest.fail("an unstarted heartbeat thread must not be joined")

    monkeypatch.setattr(admission, "Thread", FailingThread)
    heartbeat = admission.ReservationHeartbeat(
        fake_redis,
        "running",
        active_lease_seconds=30,
        heartbeat_seconds=10,
        clock=lambda: 110,
    )

    with pytest.raises(RuntimeError, match="heartbeat thread failed to start"):
        heartbeat.start()

    heartbeat.stop()


def test_background_heartbeat_logs_redis_errors_and_reconciles_absence(
    fake_redis, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A transient renewal error does not kill the heartbeat loop or resurrect a lease."""
    from app import admission
    from redis.exceptions import ConnectionError as RedisConnectionError

    outcomes = iter((True, RedisConnectionError("redis unavailable"), False))
    calls = 0

    def fake_renew(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        outcome = next(outcomes)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    class NeverStopped:
        def wait(self, _seconds):
            return False

        def set(self):
            return None

    class ImmediateThread:
        def __init__(self, *, target, **_kwargs):
            self.target = target

        def start(self):
            self.target()

        def join(self, timeout=None):
            return None

    monkeypatch.setattr(admission, "renew_capacity", fake_renew)
    monkeypatch.setattr(admission, "Event", NeverStopped)
    monkeypatch.setattr(admission, "Thread", ImmediateThread)
    caplog.set_level(logging.ERROR, logger="app.admission")
    heartbeat = admission.ReservationHeartbeat(
        fake_redis,
        "running",
        active_lease_seconds=30,
        heartbeat_seconds=10,
    )

    heartbeat.start()

    assert calls == 3
    assert "Failed to renew capacity reservation running: redis unavailable" in caplog.text


def test_killed_worker_reservation_expires_after_its_active_lease(fake_redis) -> None:
    """A worker that stops heartbeating cannot retain capacity forever."""
    from app import admission

    assert _reserve(admission, fake_redis, "killed", now=100).accepted
    assert admission.renew_capacity(fake_redis, "killed", active_lease_seconds=30, clock=lambda: 150)

    replacement = _reserve(admission, fake_redis, "replacement", now=180)

    assert replacement.accepted
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {"replacement": "100"}


def test_renewal_does_not_resurrect_an_expired_or_released_reservation(fake_redis) -> None:
    """A late heartbeat cannot recreate state removed by prune or cleanup."""
    from app import admission

    assert _reserve(admission, fake_redis, "job-a", now=100).accepted
    admission.release_capacity(fake_redis, "job-a")

    assert not admission.renew_capacity(fake_redis, "job-a", active_lease_seconds=30, clock=lambda: 110)
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}
    assert fake_redis.zcard(admission.ADMISSION_DEADLINES) == 0


@pytest.mark.parametrize(
    ("reservation_id", "lease_seconds", "message"),
    [("", 30, "identifier must not be empty"), ("job-a", 0, "lease seconds must be positive")],
)
def test_renewal_rejects_invalid_inputs_before_redis(
    fake_redis, reservation_id: str, lease_seconds: int, message: str
) -> None:
    """Invalid worker lease inputs fail closed without creating state."""
    from app import admission

    with pytest.raises(ValueError, match=message):
        admission.renew_capacity(fake_redis, reservation_id, active_lease_seconds=lease_seconds)

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


@pytest.mark.parametrize("active_seconds,heartbeat_seconds", [(0, 1), (30, 0), (30, 30)])
def test_heartbeat_rejects_timing_that_cannot_renew_before_expiry(
    fake_redis, active_seconds: int, heartbeat_seconds: int
) -> None:
    """A deployment cannot configure a non-renewable active lease."""
    from app import admission

    with pytest.raises(ValueError):
        admission.ReservationHeartbeat(
            fake_redis,
            "job-a",
            active_lease_seconds=active_seconds,
            heartbeat_seconds=heartbeat_seconds,
        )


def test_prune_and_reserve_remain_atomic_under_concurrency(fake_redis) -> None:
    """Two contenders cannot both replace the same expired queue slot."""
    from app import admission

    assert _reserve(admission, fake_redis, "stale", now=100).accepted
    real_pipeline = fake_redis.pipeline
    read_barrier = Barrier(2)
    read_lock = Lock()
    reads = 0

    class SynchronizedPipeline:
        def __init__(self):
            self._pipeline = real_pipeline()

        def __enter__(self):
            self._pipeline.__enter__()
            return self

        def __exit__(self, *args):
            return self._pipeline.__exit__(*args)

        def __getattr__(self, name):
            return getattr(self._pipeline, name)

        def hgetall(self, key):
            nonlocal reads
            result = self._pipeline.hgetall(key)
            with read_lock:
                reads += 1
                this_read = reads
            if this_read <= 2:
                read_barrier.wait()
            return result

    fake_redis.pipeline = SynchronizedPipeline

    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = list(pool.map(lambda job_id: _reserve(admission, fake_redis, job_id, now=161), ("a", "b")))

    assert sorted(outcome.accepted for outcome in outcomes) == [False, True]
    assert fake_redis.hlen(admission.ADMISSION_RESERVATIONS) == 1
    assert fake_redis.zcard(admission.ADMISSION_DEADLINES) == 1


def test_concurrent_disk_reservations_cannot_both_spend_the_same_headroom(fake_redis) -> None:
    """The byte decision and write share the same WATCH transaction."""
    from app import admission

    real_pipeline = fake_redis.pipeline
    read_barrier = Barrier(2)
    read_lock = Lock()
    reads = 0

    class SynchronizedPipeline:
        def __init__(self):
            self._pipeline = real_pipeline()

        def __enter__(self):
            self._pipeline.__enter__()
            return self

        def __exit__(self, *args):
            return self._pipeline.__exit__(*args)

        def __getattr__(self, name):
            return getattr(self._pipeline, name)

        def hgetall(self, key):
            nonlocal reads
            result = self._pipeline.hgetall(key)
            with read_lock:
                reads += 1
                this_read = reads
            if this_read <= 2:
                read_barrier.wait()
            return result

    fake_redis.pipeline = SynchronizedPipeline

    def reserve(job_id: str):
        return admission.reserve_capacity(
            fake_redis,
            reservation_id=job_id,
            requested_bytes=100,
            free_bytes=(250,),
            max_jobs=10,
            minimum_free_bytes=100,
            queued_grace_seconds=60,
            clock=lambda: 0,
        )

    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = list(pool.map(reserve, ("a", "b")))

    assert sorted((outcome.accepted, outcome.reason) for outcome in outcomes) == [(False, "disk"), (True, None)]
    assert fake_redis.hlen(admission.ADMISSION_RESERVATIONS) == 1


@pytest.mark.parametrize(
    ("member_bytes", "base_bytes", "member_factor", "expected"),
    [(0, 1_000, 2, 1_000), (250, 1_000, 2, 1_500), (250, 1_000, 3, 1_750)],
)
def test_cohort_reservation_budget_arithmetic(
    member_bytes: int, base_bytes: int, member_factor: int, expected: int
) -> None:
    """Cohort reservations include fixed work area plus snapshot/archive headroom."""
    from app.admission import cohort_reservation_bytes

    assert cohort_reservation_bytes(member_bytes, base_bytes, member_factor) == expected


@pytest.mark.parametrize("values", [(-1, 1, 1), (1, -1, 1), (1, 1, -1)])
def test_cohort_reservation_budget_rejects_negative_inputs(values: tuple[int, int, int]) -> None:
    """Misconfigured heuristic coefficients fail closed."""
    from app.admission import cohort_reservation_bytes

    with pytest.raises(ValueError, match="must not be negative"):
        cohort_reservation_bytes(*values)


def test_pipeline_budget_default_reserves_input_and_two_working_copies(monkeypatch: pytest.MonkeyPatch) -> None:
    """The shipped heuristic is three upload ceilings, not an asserted runtime bound."""
    from app.config import pipeline_reservation_bytes

    monkeypatch.delenv("PIPELINE_RESERVATION_BYTES", raising=False)
    assert pipeline_reservation_bytes(1_000) == 3_000
    monkeypatch.setenv("PIPELINE_RESERVATION_BYTES", "4500")
    assert pipeline_reservation_bytes(1_000) == 4_500


def test_admission_lease_defaults_are_conservative_and_internally_valid() -> None:
    """Queued grace is long, while heartbeat cadence stays below the active lease."""
    from app.config import Settings

    assert Settings.ADMISSION_QUEUED_GRACE_SECONDS == 7 * 24 * 60 * 60
    assert 0 < Settings.ADMISSION_HEARTBEAT_SECONDS < Settings.ADMISSION_ACTIVE_LEASE_SECONDS
