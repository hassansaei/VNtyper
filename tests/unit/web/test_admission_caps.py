"""Atomic admission limits for pipeline and cohort submissions."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier, Lock
from types import SimpleNamespace

import pytest
import yaml
from redis.exceptions import ConnectionError as RedisConnectionError

pytestmark = pytest.mark.unit


def _submit(client):
    """Submit one minimal pipeline request."""
    return client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"x", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )


def _capacity_module():
    """Import the production helper while keeping the initial RED a test failure."""
    try:
        from app import admission
    except ImportError:
        return None
    return admission


def test_a_full_reservation_queue_refuses_before_writing(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The admission registry, not display bookkeeping, owns the queue cap."""
    monkeypatch_values = {
        "MAX_QUEUED_JOBS": 1,
        "MIN_FREE_DISK_BYTES": 0,
        "MAX_UPLOAD_BYTES": 1,
    }
    for name, value in monkeypatch_values.items():
        monkeypatch.setattr(web_app.settings, name, value)
    fake_redis.hset("vntyper:admission:reservations", "running-job", 1)

    response = _submit(client)

    assert response.status_code == 503
    assert response.headers["Retry-After"] == "60"
    assert list((tmp_path / "input").iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()


def test_stale_display_entries_do_not_consume_authoritative_capacity(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The legacy queue-position list is deliberately not the admission count."""
    monkeypatch.setattr(web_app.settings, "MAX_QUEUED_JOBS", 1)
    fake_redis.rpush("vntyper_job_queue", "already-finished-task")

    response = _submit(client)

    assert response.status_code == 200
    web_app.run_vntyper_job.delay.assert_called_once()


def test_concurrent_reservations_cannot_both_take_one_remaining_slot(fake_redis) -> None:
    """WATCH/MULTI must turn two simultaneous observations into one admission."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    real_pipeline = fake_redis.pipeline
    read_barrier = Barrier(2)
    read_lock = Lock()
    read_count = 0

    class SynchronizedPipeline:
        """Force both transactions to read the watched empty hash once."""

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
            nonlocal read_count
            result = self._pipeline.hgetall(key)
            with read_lock:
                read_count += 1
                this_read = read_count
            if this_read <= 2:
                read_barrier.wait()
            return result

    fake_redis.pipeline = SynchronizedPipeline

    def reserve(identifier: str):
        return admission.reserve_capacity(
            fake_redis,
            reservation_id=identifier,
            requested_bytes=100,
            free_bytes=(300,),
            max_jobs=1,
            minimum_free_bytes=100,
        )

    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = list(pool.map(reserve, ("job-a", "job-b")))

    assert sorted(outcome.accepted for outcome in outcomes) == [False, True]
    assert fake_redis.hlen(admission.ADMISSION_RESERVATIONS) == 1


@pytest.mark.parametrize(
    ("free_bytes", "accepted"),
    [
        ((201,), True),
        ((200,), True),
        ((199,), False),
        ((500, 199), False),
    ],
)
def test_disk_decision_retains_the_floor_after_all_reservations(
    fake_redis, free_bytes: tuple[int, ...], accepted: bool
) -> None:
    """The boundary is floor + aggregate reservations, on every volume."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    outcome = admission.reserve_capacity(
        fake_redis,
        reservation_id="job-a",
        requested_bytes=100,
        free_bytes=free_bytes,
        max_jobs=10,
        minimum_free_bytes=100,
    )

    assert outcome.accepted is accepted


def test_release_is_idempotent_and_frees_the_reserved_bytes(fake_redis) -> None:
    """Worker retries and endpoint compensation may release the same id safely."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    first = admission.reserve_capacity(
        fake_redis,
        reservation_id="job-a",
        requested_bytes=100,
        free_bytes=(300,),
        max_jobs=1,
        minimum_free_bytes=100,
    )
    assert first.accepted

    admission.release_capacity(fake_redis, "job-a")
    admission.release_capacity(fake_redis, "job-a")

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}
    assert fake_redis.zcard(admission.ADMISSION_DEADLINES) == 0
    second = admission.reserve_capacity(
        fake_redis,
        reservation_id="job-b",
        requested_bytes=100,
        free_bytes=(300,),
        max_jobs=1,
        minimum_free_bytes=100,
    )
    assert second.accepted


def test_reserving_the_same_identifier_twice_does_not_double_count(fake_redis) -> None:
    """A retried API operation recognizes its existing reservation."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    kwargs = {
        "reservation_id": "job-a",
        "requested_bytes": 100,
        "free_bytes": (300,),
        "max_jobs": 1,
        "minimum_free_bytes": 100,
    }

    first = admission.reserve_capacity(fake_redis, **kwargs)
    second = admission.reserve_capacity(fake_redis, **kwargs)

    assert first.accepted and second.accepted
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {"job-a": "100"}


def test_lease_logs_a_redis_failure_while_preserving_the_dispatch_failure(
    fake_redis, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """Rollback bookkeeping cannot replace the broker exception being raised."""
    import logging

    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    def fail_release(*_args, **_kwargs):
        raise RedisConnectionError("redis unavailable")

    monkeypatch.setattr(admission, "release_capacity", fail_release)
    caplog.set_level(logging.ERROR, logger="app.admission")

    with pytest.raises(RuntimeError, match="broker unavailable"), admission.CapacityLease(fake_redis, "job-a"):
        raise RuntimeError("broker unavailable")

    assert "Failed to roll back capacity reservation job-a: redis unavailable" in caplog.text


def test_volume_observation_inspects_both_configured_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    """Path-specific quotas cannot hide behind matching device identifiers."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    calls: list[str] = []

    def disk_usage(path: str):
        calls.append(path)
        return SimpleNamespace(free={"/input": 1234, "/output": 567}[path])

    monkeypatch.setattr(admission.shutil, "disk_usage", disk_usage)

    assert admission.observe_free_bytes(("/input", "/output")) == (1234, 567)
    assert calls == ["/input", "/output"]


def test_filesystem_observation_failure_refuses_before_writing(
    client, web_app, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """An unreadable volume is unavailable capacity, not permission to enqueue."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    def fail_stat(_paths):
        raise OSError("stat failed")

    monkeypatch.setattr(admission, "observe_free_bytes", fail_stat)

    response = _submit(client)

    assert response.status_code == 503
    assert response.headers["Retry-After"] == "300"
    assert list((tmp_path / "input").iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()


def test_low_disk_refuses_before_writing(client, web_app, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """One byte below floor plus reservation is a retryable storage refusal."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    available = web_app.settings.MIN_FREE_DISK_BYTES + web_app.settings.PIPELINE_RESERVATION_BYTES - 1
    monkeypatch.setattr(admission, "observe_free_bytes", lambda _paths: (available,))

    response = _submit(client)

    assert response.status_code == 503
    assert response.headers["Retry-After"] == "300"
    assert list((tmp_path / "input").iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()


def test_pipeline_disk_boundary_accepts_the_configured_heuristic(
    client, web_app, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Admission accepts exactly floor plus the pipeline headroom setting."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    available = web_app.settings.MIN_FREE_DISK_BYTES + web_app.settings.PIPELINE_RESERVATION_BYTES
    monkeypatch.setattr(admission, "observe_free_bytes", lambda _paths: (available,))

    response = _submit(client)

    assert response.status_code == 200


def test_redis_failure_refuses_before_writing(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Admission fails closed when its shared atomic registry is unreachable."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    def fail_pipeline(*_args, **_kwargs):
        raise RedisConnectionError("redis unavailable")

    monkeypatch.setattr(fake_redis, "pipeline", fail_pipeline)

    response = _submit(client)

    assert response.status_code == 503
    assert response.headers["Retry-After"] == "60"
    assert list((tmp_path / "input").iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()


def test_enqueue_failure_rolls_back_the_pipeline_reservation(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A broker refusal cannot consume an admission slot permanently."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    def fail_enqueue(**_kwargs):
        raise RuntimeError("broker unavailable")

    monkeypatch.setattr(web_app.run_vntyper_job, "delay", fail_enqueue)

    with pytest.raises(RuntimeError, match="broker unavailable"):
        _submit(client)

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


def test_accepted_pipeline_submission_transfers_the_reservation_to_the_worker(client, web_app, fake_redis) -> None:
    """A successful broker dispatch keeps its slot until worker cleanup."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    response = _submit(client)

    assert response.status_code == 200
    job_id = response.json()["job_id"]
    assert fake_redis.hget(admission.ADMISSION_RESERVATIONS, job_id) == str(3 * 1024 * 1024 * 1024)
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, job_id) is not None
    assert web_app.run_vntyper_job.delay.call_args.kwargs["capacity_reserved"] is True


@pytest.mark.parametrize("failure", [ValueError("upload failed"), OSError("volume write failed")])
def test_pipeline_upload_or_write_failure_rolls_back_capacity(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, failure: Exception
) -> None:
    """A failure after admission but before dispatch cannot leak the lease."""
    from app import admission

    def fail_write(*_args, **_kwargs):
        raise failure

    monkeypatch.setattr(web_app, "save_upload_bounded", fail_write)

    if isinstance(failure, ValueError):
        assert _submit(client).status_code == 413
    else:
        with pytest.raises(OSError, match="volume write failed"):
            _submit(client)

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}
    assert fake_redis.zcard(admission.ADMISSION_DEADLINES) == 0
    web_app.run_vntyper_job.delay.assert_not_called()


def test_cohort_submission_uses_the_shared_cap_and_rolls_back_dispatch_failure(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Cohort analyses cannot bypass or leak the shared reservation registry."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    created = client.post("/create-cohort/", data={"passphrase": "correct-horse-battery-staple"})
    cohort_id = created.json()["cohort_id"]
    fake_redis.sadd(f"cohort:{cohort_id}:jobs", "member-a")
    (tmp_path / "output" / "member-a.zip").write_bytes(b"member-result")

    def fail_enqueue(**_kwargs):
        raise RuntimeError("broker unavailable")

    monkeypatch.setattr(web_app.run_cohort_analysis_job, "delay", fail_enqueue)

    with pytest.raises(RuntimeError, match="broker unavailable"):
        client.post(
            "/cohort-analysis/",
            data={"cohort_id": cohort_id, "passphrase": "correct-horse-battery-staple"},
        )

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


def test_cohort_submission_cannot_bypass_a_full_shared_queue(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The capacity cap applies equally to joint cohort analysis tasks."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    created = client.post("/create-cohort/", data={"passphrase": "correct-horse-battery-staple"})
    cohort_id = created.json()["cohort_id"]
    fake_redis.sadd(f"cohort:{cohort_id}:jobs", "member-a")
    (tmp_path / "output" / "member-a.zip").write_bytes(b"member-result")
    monkeypatch.setattr(web_app.settings, "MAX_QUEUED_JOBS", 1)
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, "running-job", 1)

    response = client.post(
        "/cohort-analysis/",
        data={"cohort_id": cohort_id, "passphrase": "correct-horse-battery-staple"},
    )

    assert response.status_code == 503
    assert response.headers["Retry-After"] == "60"
    web_app.run_cohort_analysis_job.delay.assert_not_called()


def test_accepted_cohort_reserves_snapshot_archive_and_fixed_headroom(
    client, web_app, fake_redis, tmp_path: Path
) -> None:
    """The endpoint applies the reviewed cohort heuristic to actual member bytes."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    created = client.post("/create-cohort/", data={"passphrase": "correct-horse-battery-staple"})
    cohort_id = created.json()["cohort_id"]
    fake_redis.sadd(f"cohort:{cohort_id}:jobs", "member-a")
    member_bytes = b"member-result"
    (tmp_path / "output" / "member-a.zip").write_bytes(member_bytes)

    response = client.post(
        "/cohort-analysis/",
        data={"cohort_id": cohort_id, "passphrase": "correct-horse-battery-staple"},
    )

    assert response.status_code == 200
    analysis_id = response.json()["analysis_job_id"]
    expected = web_app.settings.COHORT_RESERVATION_BASE_BYTES + 2 * len(member_bytes)
    assert fake_redis.hget(admission.ADMISSION_RESERVATIONS, analysis_id) == str(expected)
    assert web_app.run_cohort_analysis_job.delay.call_args.kwargs["capacity_reserved"] is True


def test_cohort_member_stat_failure_refuses_without_dispatch(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A member disappearing during admission fails closed before reservation."""
    created = client.post("/create-cohort/", data={"passphrase": "correct-horse-battery-staple"})
    cohort_id = created.json()["cohort_id"]
    fake_redis.sadd(f"cohort:{cohort_id}:jobs", "member-a")
    (tmp_path / "output" / "member-a.zip").write_bytes(b"member-result")

    def fail_size(_path):
        raise OSError("member disappeared")

    monkeypatch.setattr(web_app.os.path, "getsize", fail_size)

    response = client.post(
        "/cohort-analysis/",
        data={"cohort_id": cohort_id, "passphrase": "correct-horse-battery-staple"},
    )

    assert response.status_code == 503
    assert response.headers["Retry-After"] == "300"
    web_app.run_cohort_analysis_job.delay.assert_not_called()


def test_the_caps_have_documented_defaults() -> None:
    """The public configuration exposes the reviewed queue and disk values."""
    from app.config import Settings

    assert Settings.MAX_QUEUED_JOBS == 100
    assert Settings.MIN_FREE_DISK_BYTES == 2 * Settings.MAX_UPLOAD_BYTES


def test_default_disk_floor_tracks_an_overridden_upload_ceiling(monkeypatch: pytest.MonkeyPatch) -> None:
    """Compose's empty optional value must preserve the derived 2x default."""
    from app.config import minimum_free_disk_bytes

    monkeypatch.delenv("MIN_FREE_DISK_BYTES", raising=False)
    assert minimum_free_disk_bytes(200) == 400
    monkeypatch.setenv("MIN_FREE_DISK_BYTES", "")
    assert minimum_free_disk_bytes(200) == 400
    monkeypatch.setenv("MIN_FREE_DISK_BYTES", "123")
    assert minimum_free_disk_bytes(200) == 123

    monkeypatch.setenv("MIN_FREE_DISK_BYTES", "not-an-integer")
    with pytest.raises(ValueError):
        minimum_free_disk_bytes(200)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"reservation_id": ""}, "identifier must not be empty"),
        ({"requested_bytes": -1}, "bytes must not be negative"),
        ({"free_bytes": ()}, "observations must be present"),
        ({"free_bytes": (-1,)}, "observations must be present"),
        ({"max_jobs": 0}, "must be positive"),
        ({"minimum_free_bytes": -1}, "must not be negative"),
        ({"queued_grace_seconds": 0}, "grace seconds must be positive"),
    ],
)
def test_invalid_internal_capacity_inputs_fail_before_redis(
    fake_redis, overrides: dict[str, object], message: str
) -> None:
    """Misconfigured limits fail closed instead of corrupting the registry."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"
    kwargs = {
        "reservation_id": "job-a",
        "requested_bytes": 1,
        "free_bytes": (100,),
        "max_jobs": 2,
        "minimum_free_bytes": 1,
        "queued_grace_seconds": 60,
    }
    kwargs.update(overrides)

    with pytest.raises(ValueError, match=message):
        admission.reserve_capacity(fake_redis, **kwargs)

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


def test_volume_observation_requires_at_least_one_path() -> None:
    """An empty volume set is an internal error, never unlimited capacity."""
    admission = _capacity_module()
    assert admission is not None, "app.admission has not been implemented"

    with pytest.raises(ValueError, match="At least one capacity path is required"):
        admission.observe_free_bytes(())


def test_compose_passes_the_admission_knobs_to_api_and_worker() -> None:
    """Values set in docker/.env must reach both service processes."""
    compose_path = Path(__file__).resolve().parents[3] / "docker" / "docker-compose.yml"
    compose = yaml.safe_load(compose_path.read_text())

    for service_name in ("api", "worker"):
        environment = compose["services"][service_name]["environment"]
        assert "MAX_UPLOAD_BYTES=${MAX_UPLOAD_BYTES:-1073741824}" in environment
        assert "MAX_QUEUED_JOBS=${MAX_QUEUED_JOBS:-100}" in environment
        assert "MIN_FREE_DISK_BYTES=${MIN_FREE_DISK_BYTES:-}" in environment
        assert "PIPELINE_RESERVATION_BYTES=${PIPELINE_RESERVATION_BYTES:-}" in environment
        assert "COHORT_RESERVATION_BASE_BYTES=${COHORT_RESERVATION_BASE_BYTES:-}" in environment
        assert "COHORT_MEMBER_RESERVATION_FACTOR=${COHORT_MEMBER_RESERVATION_FACTOR:-2}" in environment
        assert "ADMISSION_QUEUED_GRACE_SECONDS=${ADMISSION_QUEUED_GRACE_SECONDS:-604800}" in environment
        assert "ADMISSION_ACTIVE_LEASE_SECONDS=${ADMISSION_ACTIVE_LEASE_SECONDS:-86400}" in environment
        assert "ADMISSION_HEARTBEAT_SECONDS=${ADMISSION_HEARTBEAT_SECONDS:-300}" in environment
