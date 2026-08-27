"""What the joint cohort analysis task is allowed to delete when it finishes.

`run_cohort_analysis_job` reads the result archive of every member of a cohort.
Its cleanup block therefore runs with a list of paths that belong to other jobs,
and the distinction between "the scratch file this task wrote" and "the inputs
this task was given" is the whole subject of this module: cleanup may remove the
former and must leave the latter alone. Member archives are what
`/download/{job_id}/` serves and what a later analysis reads again; they are
aged out centrally by `delete_old_results()`, not by whoever happened to read
them last.

The task is replaced by a MagicMock in the `web_app` fixture, so it is exercised
here directly instead. Authorization for the endpoint that enqueues it lives in
`test_cohort_authz.py`.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.tasks` is importable here.
"""

import os
import zipfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


@pytest.fixture
def cohort_analysis_task(monkeypatch: pytest.MonkeyPatch, fake_redis):
    """Run `run_cohort_analysis_job`'s body with Redis and subprocess neutralised.

    It is invoked through `.run()` inside a pushed request rather than through
    `.apply()`: `.apply()` still routes the result through Celery's configured
    result backend, and this app's broker URL is built in `celery_app.py`, which
    belongs to another task in this wave.

    Args:
        monkeypatch: Standard pytest fixture; restores every patch at teardown.
        fake_redis: In-process Redis stand-in from conftest.

    Returns:
        Callable[..., None]: Invokes the task body with the given keyword
            arguments, with a request context in place.
    """
    from app import tasks

    for attr in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, attr, fake_redis)

    def _fake_run(command, *args, **kwargs):
        """Stand in for the `vntyper cohort` subprocess without running it."""
        return None

    monkeypatch.setattr(tasks.subprocess, "run", _fake_run)

    def _invoke(**kwargs) -> None:
        """Call the task body as a worker would.

        Args:
            **kwargs: Arguments forwarded to `run_cohort_analysis_job`.
        """
        tasks.run_cohort_analysis_job.push_request(id="analysis-task-1")
        try:
            tasks.run_cohort_analysis_job.run(**kwargs)
        finally:
            tasks.run_cohort_analysis_job.pop_request()

    return _invoke


def test_cohort_analysis_preserves_every_member_result(cohort_analysis_task, tmp_path: Path) -> None:
    """Running an analysis leaves the members' own result archives in place.

    Each member's `.zip` is what `/download/{job_id}/` serves and what a second
    analysis would read, so consuming them turns one analysis request into the
    permanent loss of every member's results.

    Args:
        cohort_analysis_task: The task fixture above.
        tmp_path: Scratch directory standing in for the output tree.
    """
    zip_paths = []
    for job_id in ("job-a", "job-b"):
        zip_path = tmp_path / f"{job_id}.zip"
        zip_path.write_bytes(b"result data for " + job_id.encode())
        zip_paths.append(str(zip_path))

    cohort_analysis_task(
        cohort_id="cohort-1",
        zip_paths=zip_paths,
        output_dir=str(tmp_path / "analysis"),
    )

    for job_id, member_result in zip(("job-a", "job-b"), zip_paths, strict=True):
        assert os.path.exists(member_result), f"{job_id} lost its results to the cohort analysis"
        assert Path(member_result).read_bytes() == b"result data for " + job_id.encode()


def test_cohort_analysis_still_removes_its_own_scratch_file(cohort_analysis_task, tmp_path: Path) -> None:
    """The listing the task writes for itself is still cleaned up.

    Args:
        cohort_analysis_task: The task fixture above.
        tmp_path: Scratch directory standing in for the output tree.
    """
    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    output_dir = tmp_path / "analysis"

    cohort_analysis_task(
        cohort_id="cohort-1",
        zip_paths=[str(zip_path)],
        output_dir=str(output_dir),
    )

    assert not (output_dir / "cohort_input.txt").exists()
    with zipfile.ZipFile(f"{output_dir}.zip") as archive:
        assert "cohort_input.txt" not in archive.namelist()
        assert str(zip_path).encode() not in Path(f"{output_dir}.zip").read_bytes()


def test_cohort_analysis_surfaces_a_setup_failure(cohort_analysis_task, monkeypatch, tmp_path: Path) -> None:
    """A failure before the scratch file exists reports itself, not a NameError.

    The cleanup block runs whether or not the task got as far as naming its
    scratch file, so it must not assume the name is bound.

    Args:
        cohort_analysis_task: The task fixture above.
        monkeypatch: Standard pytest fixture.
        tmp_path: Scratch directory standing in for the output tree.
    """
    from app import tasks

    def _refuse(*args, **kwargs):
        """Fail at the first filesystem call the task makes."""
        raise OSError("no space left on device")

    monkeypatch.setattr(tasks, "cohort_workspace", _refuse)

    with pytest.raises(OSError, match="no space left on device"):
        cohort_analysis_task(
            cohort_id="cohort-1",
            zip_paths=[str(tmp_path / "job-a.zip")],
            output_dir=str(tmp_path / "analysis"),
        )
