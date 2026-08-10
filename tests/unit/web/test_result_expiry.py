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
    """One failed result deletion does not stop the scheduled cleanup sweep."""
    from app import tasks

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, attr, fake_redis)
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(tmp_path))
    monkeypatch.setattr(tasks.os, "listdir", lambda _path: ["blocked.zip", "later.zip"])
    monkeypatch.setattr(tasks.os.path, "isfile", lambda _path: True)
    too_old = (tasks.settings.MAX_RESULT_AGE_DAYS + 1) * _SECONDS_PER_DAY
    monkeypatch.setattr(tasks.os.path, "getctime", lambda _path: time.time() - too_old)

    attempted: list[str] = []

    def remove(path: str) -> None:
        attempted.append(path)
        if path.endswith("blocked.zip"):
            raise OSError("blocked")

    monkeypatch.setattr(tasks.os, "remove", remove)

    with caplog.at_level(logging.ERROR, logger=tasks.logger.name):
        assert tasks.delete_old_results.run() is None

    assert attempted == [str(tmp_path / "blocked.zip"), str(tmp_path / "later.zip")]
    errors = [
        record for record in caplog.records if record.name == tasks.logger.name and record.levelno >= logging.ERROR
    ]
    assert [record.levelno for record in errors] == [logging.ERROR]
    assert "blocked" in caplog.text
