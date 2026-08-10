"""How `delete_old_results()` decides that a result archive has aged out.

The task compares a per-file timestamp against a cutoff derived from
`MAX_RESULT_AGE_DAYS`. Both sides of that comparison have to be built the same
way: Python refuses to order a timezone-aware datetime against a naive one and
raises `TypeError`, so a change that converts one side and not the other turns
scheduled cleanup into a task that dies on its first archive and never frees
disk again. These tests pin both ends of the decision - past the limit is
deleted, inside the limit is kept - so the comparison stays evaluable and keeps
its meaning.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.tasks` is importable here.
"""

import logging
import time
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_SECONDS_PER_DAY = 86400


@pytest.fixture
def expiry_task(monkeypatch: pytest.MonkeyPatch, fake_redis, tmp_path: Path):
    """Run `delete_old_results` against a scratch output tree and fake Redis.

    `os.path.getctime` is stubbed rather than the files backdated: creation time
    is an inode field on Linux and `os.utime` cannot set it, so a real file
    cannot be made to look old.

    Args:
        monkeypatch: Standard pytest fixture; restores every patch at teardown.
        fake_redis: In-process Redis stand-in from conftest.
        tmp_path: Scratch directory used as the output tree.

    Returns:
        Callable[[dict[str, float]], None]: Called with a mapping of file name
            to age in seconds, it creates those files and runs the task.
    """
    from app import tasks

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, attr, fake_redis)
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(tmp_path))

    def _invoke(ages_in_seconds: dict[str, float]) -> None:
        """Create the named archives with the given ages and run the task.

        Args:
            ages_in_seconds: File name to how long ago it should appear to have
                been created.
        """
        now = time.time()
        creation_times = {}
        for name, age in ages_in_seconds.items():
            (tmp_path / name).write_bytes(b"result data")
            creation_times[str(tmp_path / name)] = now - age

        monkeypatch.setattr(tasks.os.path, "getctime", lambda path: creation_times[path])
        tasks.delete_old_results.run()

    return _invoke


def test_an_archive_past_the_age_limit_is_deleted(expiry_task, tmp_path: Path) -> None:
    """An archive older than `MAX_RESULT_AGE_DAYS` is removed.

    Args:
        expiry_task: The task fixture above.
        tmp_path: The scratch output tree the fixture configured.
    """
    from app import tasks

    too_old = (tasks.settings.MAX_RESULT_AGE_DAYS + 1) * _SECONDS_PER_DAY
    expiry_task({"stale.zip": too_old})

    assert not (tmp_path / "stale.zip").exists()


def test_an_archive_inside_the_age_limit_is_kept(expiry_task, tmp_path: Path) -> None:
    """An archive younger than `MAX_RESULT_AGE_DAYS` survives the sweep.

    Args:
        expiry_task: The task fixture above.
        tmp_path: The scratch output tree the fixture configured.
    """
    from app import tasks

    still_fresh = (tasks.settings.MAX_RESULT_AGE_DAYS - 1) * _SECONDS_PER_DAY
    expiry_task({"recent.zip": still_fresh})

    assert (tmp_path / "recent.zip").exists()


def test_the_sweep_leaves_files_that_are_not_result_archives_alone(expiry_task, tmp_path: Path) -> None:
    """Only `.zip` archives are candidates, however old the other files are.

    Args:
        expiry_task: The task fixture above.
        tmp_path: The scratch output tree the fixture configured.
    """
    from app import tasks

    too_old = (tasks.settings.MAX_RESULT_AGE_DAYS + 1) * _SECONDS_PER_DAY
    expiry_task({"stale.zip": too_old, "stale.log": too_old})

    assert not (tmp_path / "stale.zip").exists()
    assert (tmp_path / "stale.log").exists()


def test_cleanup_continues_after_one_delete_error(
    monkeypatch: pytest.MonkeyPatch, fake_redis, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """One failed expired-cohort deletion does not stop its remaining cleanup."""
    from app import tasks

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, attr, fake_redis)
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(tasks.os, "listdir", lambda _path: [])
    redis_exists = fake_redis.exists
    monkeypatch.setattr(fake_redis, "scan_iter", lambda _pattern: iter(["cohort:expired"]))
    monkeypatch.setattr(fake_redis, "exists", lambda key: False if key == "cohort:expired" else redis_exists(key))
    monkeypatch.setattr(fake_redis, "smembers", lambda _key: ["blocked-job", "later-job"])

    blocked_archive = tmp_path / "blocked-job.zip"
    later_archive = tmp_path / "later-job.zip"
    blocked_archive.write_bytes(b"blocked")
    later_archive.write_bytes(b"later")
    cohort_jobs_key = "cohort:expired:jobs"
    fake_redis.sadd(cohort_jobs_key, "blocked-job", "later-job")
    for job_id in ("blocked-job", "later-job"):
        fake_redis.set(job_id, "result")
        fake_redis.set(f"celery-task-meta-{job_id}", "metadata")

    attempted: list[str] = []

    def remove(path: str) -> None:
        attempted.append(path)
        if path == str(blocked_archive):
            raise OSError("blocked")
        Path(path).unlink()

    monkeypatch.setattr(tasks.os, "remove", remove)

    with caplog.at_level(logging.ERROR, logger=tasks.logger.name):
        assert tasks.delete_old_results.run() is None

    assert attempted == [str(blocked_archive), str(later_archive)]
    assert blocked_archive.exists()
    assert not later_archive.exists()
    assert not redis_exists(cohort_jobs_key)
    for job_id in ("blocked-job", "later-job"):
        assert fake_redis.get(job_id) is None
        assert fake_redis.get(f"celery-task-meta-{job_id}") is None
    errors = [
        record for record in caplog.records if record.name == tasks.logger.name and record.levelno >= logging.ERROR
    ]
    assert [record.levelno for record in errors] == [logging.ERROR]
    assert errors[0].getMessage() == f"Error deleting file {blocked_archive}: blocked"
