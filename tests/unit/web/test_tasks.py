"""Coverage of the Celery worker's failure, retry and cleanup paths.

`docker/app/tasks.py` runs the pipeline out-of-process (via `subprocess.run`)
and talks to three separate Redis databases. Every path that runs when
something goes wrong -- a bad exit code, an exception before the job even
starts, a filesystem error while tidying up -- was previously exercised only
by accident, if at all: the existing sibling modules
(`test_index_handoff.py`, `test_cohort_analysis_task.py`,
`test_result_expiry.py`) each cover one narrow slice of this file for their
own purposes and leave the rest alone. This module targets what they leave
uncovered: the `subprocess.CalledProcessError` branches, the failure-email
content, the `try`/`except`/`finally` cleanup block's own error handling, and
the parts of `delete_old_results` and `run_cohort_analysis_job` nothing else
reaches.

Every Celery task here is invoked through `.run()` (or `.push_request()` /
`.run()` / `.pop_request()` for bound tasks) rather than `.delay()` or
`.apply()`: both of the latter route through the broker URL built in
`celery_app.py`, which this suite must not touch. The three Redis clients are
replaced with `unittest.mock.MagicMock` rather than `fakeredis`, because most
assertions here are about which call was made with which arguments (a status
transition, a specific key), which a mock's call history states directly
rather than requiring the state to be read back out of a store.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.tasks` is importable here.
"""

import inspect
import logging
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import ANY, MagicMock, patch

import pytest

from vntyper.scripts.utils import validate_bam_file

pytestmark = pytest.mark.unit

BAM_BYTES = b"alignment-bytes"
INDEX_BYTES = b"index-bytes"


# ---------------------------------------------------------------------------
# Shared fixtures and helpers
# ---------------------------------------------------------------------------


@pytest.fixture
def redis_mocks(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Replace the three module-level Redis clients with independent mocks.

    Args:
        monkeypatch: Standard pytest fixture; restores every patch at teardown.

    Returns:
        SimpleNamespace: `.queue`, `.cohort` and `.usage`, matching
            `tasks.redis_client`, `tasks.redis_cohort_client` and
            `tasks.redis_usage_client` respectively.
    """
    from app import tasks

    mocks = SimpleNamespace(
        queue=MagicMock(name="redis_client"),
        cohort=MagicMock(name="redis_cohort_client"),
        usage=MagicMock(name="redis_usage_client"),
    )
    monkeypatch.setattr(tasks, "redis_client", mocks.queue)
    monkeypatch.setattr(tasks, "redis_cohort_client", mocks.cohort)
    monkeypatch.setattr(tasks, "redis_usage_client", mocks.usage)
    return mocks


@pytest.fixture
def no_email_task(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace `tasks.send_email_task` with a mock, standing in for Celery.

    `run_vntyper_job` calls `send_email_task.delay(...)` directly on the
    module-level task object. Left real, that call would try to reach the
    broker `celery_app.py` builds from the environment.

    Args:
        monkeypatch: Standard pytest fixture; restores the patch at teardown.

    Returns:
        MagicMock: The stand-in task object; assert on `.delay.call_args`.
    """
    from app import tasks

    mock_task = MagicMock(name="send_email_task")
    monkeypatch.setattr(tasks, "send_email_task", mock_task)
    return mock_task


def _job_kwargs(tmp_path: Path, bam_path: Path, job_id: str = "job-1", **overrides) -> dict:
    """Build a full keyword-argument set for `run_vntyper_job`.

    Args:
        tmp_path: Scratch directory standing in for the job tree.
        bam_path: The alignment path to submit.
        job_id: The job identifier; becomes the basename of `output_dir`.
        **overrides: Any field to replace in the default set.

    Returns:
        dict: Keyword arguments ready to hand to `run_vntyper_job.run(...)`.
    """
    kwargs = {
        "bam_path": str(bam_path),
        "output_dir": str(tmp_path / "output" / job_id),
        "thread": 1,
        "reference_assembly": "hg38",
        "fast_mode": False,
        "keep_intermediates": False,
        "archive_results": False,
        "email": None,
        "cohort_key": None,
        "client_ip": "127.0.0.1",
        "user_agent": "pytest",
        "advntr_mode": False,
        "index_path": None,
    }
    kwargs.update(overrides)
    return kwargs


def _make_job_input(tmp_path: Path, job_id: str = "job-1", index_name: str | None = None) -> tuple[Path, Path | None]:
    """Lay out a job input directory the way the endpoint would have.

    Args:
        tmp_path: Scratch directory standing in for the input tree.
        job_id: The job identifier; becomes the input subdirectory name.
        index_name: The stored index's filename, or None to store no index.

    Returns:
        tuple: The alignment path, and the index path or None.
    """
    job_input_dir = tmp_path / "input" / job_id
    job_input_dir.mkdir(parents=True)
    bam_path = job_input_dir / "sample.bam"
    bam_path.write_bytes(BAM_BYTES)
    if index_name is None:
        return bam_path, None
    index_path = job_input_dir / index_name
    index_path.write_bytes(INDEX_BYTES)
    return bam_path, index_path


def _invoke_vntyper_job(tasks, **kwargs) -> None:
    """Call `run_vntyper_job`'s body as a worker would, with a request context.

    Args:
        tasks: The imported `app.tasks` module.
        **kwargs: Arguments forwarded to `run_vntyper_job`.
    """
    tasks.run_vntyper_job.push_request(id="task-1")
    try:
        tasks.run_vntyper_job.run(**kwargs)
    finally:
        tasks.run_vntyper_job.pop_request()


def _invoke_cohort_job(tasks, **kwargs) -> None:
    """Call `run_cohort_analysis_job`'s body as a worker would.

    Args:
        tasks: The imported `app.tasks` module.
        **kwargs: Arguments forwarded to `run_cohort_analysis_job`.
    """
    tasks.run_cohort_analysis_job.push_request(id="analysis-1")
    try:
        tasks.run_cohort_analysis_job.run(**kwargs)
    finally:
        tasks.run_cohort_analysis_job.pop_request()


def _subprocess_stub(
    monkeypatch: pytest.MonkeyPatch,
    tasks,
    *,
    index_error: Exception | None = None,
    pipeline_error: Exception | None = None,
) -> list:
    """Replace `tasks.subprocess.run` with a recorder that can be told to fail.

    Standing in for `samtools index`, it also writes the `.bai` file that
    command would have produced, so "was an index built" is an assertion
    about the commands issued rather than a file that appeared by magic.

    Args:
        monkeypatch: Standard pytest fixture; restores the patch at teardown.
        tasks: The imported `app.tasks` module.
        index_error: If given, raised instead of running `samtools index`.
        pipeline_error: If given, raised instead of running `vntyper pipeline`.

    Returns:
        list: The argument vector of every command issued, in order.
    """
    commands: list[list[str]] = []

    def _run(command, *args, **kwargs):
        commands.append(list(command))
        if command[:2] == ["samtools", "index"]:
            if index_error is not None:
                raise index_error
            Path(f"{command[2]}.bai").write_bytes(b"generated-index")
            return None
        if pipeline_error is not None:
            raise pipeline_error
        return None

    monkeypatch.setattr(tasks.subprocess, "run", _run)
    return commands


def _quickcheck_the_way_the_pipeline_does(command: list[str]) -> None:
    """Reproduce the pipeline's own startup write, by running the code that makes it.

    A stand-in that only records commands makes the input directory look tidier
    than any real run ever leaves it, so it cannot tell a fixed worker from an
    unfixed one. This one calls the production `validate_bam_file` -- the
    function that decides where the `samtools quickcheck` log goes -- with
    `run_command` replaced so nothing is executed and the log is simply created
    where that decision points.

    The `inspect.signature` check is what makes the test discriminate: before
    #201 the function had no `log_dir` parameter at all, so the shim omits it,
    the log lands beside the input alignment exactly as it did in production, and
    the cleanup assertions below fail. It does *not* prove that `pipeline.py`
    passes `log_dir` -- that is a separate call site, pinned by
    `tests/unit/test_pipeline_cwd.py`; what it proves is what the worker's
    cleanup sees when it does.

    Args:
        command: The argument vector the task asked for.
    """
    if "pipeline" not in command:
        return
    flag = "--cram" if "--cram" in command else "--bam"
    alignment = command[command.index(flag) + 1]
    output_dir = command[command.index("-o") + 1]
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    kwargs = {"log_dir": output_dir} if "log_dir" in inspect.signature(validate_bam_file).parameters else {}

    def _create_the_log(_command, log_file, critical=False, cwd=None):
        Path(log_file).write_text("quickcheck output\n")
        return True

    with patch("vntyper.scripts.utils.run_command", _create_the_log):
        validate_bam_file(alignment, **kwargs)


def _subprocess_stub_that_runs_the_pipelines_own_startup(monkeypatch: pytest.MonkeyPatch, tasks) -> list:
    """`_subprocess_stub`, plus the one production write the pipeline makes at startup.

    Args:
        monkeypatch: Standard pytest fixture; restores the patch at teardown.
        tasks: The imported `app.tasks` module.

    Returns:
        list: The argument vector of every command issued, in order.
    """
    commands: list[list[str]] = []

    def _run(command, *args, **kwargs):
        commands.append(list(command))
        if command[:2] == ["samtools", "index"]:
            Path(f"{command[2]}.bai").write_bytes(b"generated-index")
            return None
        _quickcheck_the_way_the_pipeline_does(list(command))
        return None

    monkeypatch.setattr(tasks.subprocess, "run", _run)
    return commands


# ---------------------------------------------------------------------------
# send_email_task
# ---------------------------------------------------------------------------


def test_send_email_task_sends_and_logs_on_success(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A successful send calls `send_email` with the given fields and logs it.

    Args:
        monkeypatch: Standard pytest fixture.
        caplog: Captures the task's log record.
    """
    from app import tasks

    sent_with: list[dict] = []
    monkeypatch.setattr(tasks, "send_email", lambda **kwargs: sent_with.append(kwargs))
    caplog.set_level(logging.INFO, logger="app.tasks")

    tasks.send_email_task.run(to_email="user@example.com", subject="Hello", content="<p>hi</p>")

    assert sent_with == [{"to_email": "user@example.com", "subject": "Hello", "content": "<p>hi</p>"}]
    assert "Email sent to user@example.com with subject 'Hello'" in caplog.text


def test_send_email_task_reraises_the_smtp_failure_outside_a_worker(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A `send_email` failure is logged, then surfaced through `self.retry`.

    `self.retry()` is Celery's own machinery: called outside of a worker
    (no `push_request()` here, matching how this task is actually invoked --
    directly, never `.delay()`d onto itself) it re-raises the original
    exception rather than scheduling a retry. This test pins that observed
    behaviour; it does not assert Celery's retry scheduling itself.

    Args:
        monkeypatch: Standard pytest fixture.
        caplog: Captures the task's log record.
    """
    from app import tasks

    failure = RuntimeError("smtp connection refused")
    monkeypatch.setattr(tasks, "send_email", lambda **kwargs: (_ for _ in ()).throw(failure))
    caplog.set_level(logging.ERROR, logger="app.tasks")

    with pytest.raises(RuntimeError, match="smtp connection refused"):
        tasks.send_email_task.run(to_email="user@example.com", subject="Hello", content="<p>hi</p>")

    assert "Failed to send email to user@example.com: smtp connection refused" in caplog.text


# ---------------------------------------------------------------------------
# run_vntyper_job: a failure with nothing to do with subprocess
# ---------------------------------------------------------------------------


def test_a_failure_before_the_pipeline_even_starts_is_marked_failed_and_cleaned_up(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """A job that raises before touching the filesystem still cleans up.

    `index_path` is resolved to its fallback before the `try` block
    specifically so cleanup has a value to use even when the job fails at
    its very first Redis call. This test exercises exactly that: the usage
    bookkeeping call itself raises, before any subprocess is invoked.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    subprocess_mock = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", subprocess_mock)
    redis_mocks.usage.hset.side_effect = [RuntimeError("redis unavailable"), None]

    with pytest.raises(RuntimeError, match="redis unavailable"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    subprocess_mock.assert_not_called()
    assert redis_mocks.usage.hset.call_args_list == [
        (("usage:job-1",), {"mapping": ANY}),
        (("usage:job-1", "status", "failed"), {}),
    ]
    redis_mocks.queue.lrem.assert_called_once_with("vntyper_job_queue", 0, "task-1")
    assert not bam_path.exists(), "cleanup must still remove the alignment after an early failure"


# ---------------------------------------------------------------------------
# run_vntyper_job: subprocess exits non-zero
# ---------------------------------------------------------------------------


def test_index_generation_failure_stops_the_job_before_the_pipeline_runs(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """`samtools index` exiting non-zero fails the job without running the pipeline.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    index_error = subprocess.CalledProcessError(1, ["samtools", "index", str(bam_path)])
    commands = _subprocess_stub(monkeypatch, tasks, index_error=index_error)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, email="user@example.com"))

    assert commands == [["samtools", "index", str(bam_path)]], "the pipeline must not run after the index step fails"
    assert ("usage:job-1", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]
    assert no_email_task.delay.call_count == 0, "only a pipeline failure emails; an index failure is not a pipeline run"


def test_pipeline_failure_marks_the_job_failed_and_sends_a_failure_email(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """`vntyper pipeline` exiting non-zero fails the job and emails the failure.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])
    commands = _subprocess_stub(monkeypatch, tasks, pipeline_error=pipeline_error)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, email="user@example.com", cohort_key=None))

    assert commands[0] == ["samtools", "index", str(bam_path)], "the index is still built before the pipeline runs"
    assert commands[1][:6] == ["conda", "run", "-n", "vntyper", "vntyper", "pipeline"]
    assert ("usage:job-1", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]
    no_email_task.delay.assert_called_once()
    email_kwargs = no_email_task.delay.call_args.kwargs
    assert email_kwargs["to_email"] == "user@example.com"
    assert email_kwargs["subject"] == "VNtyper Job Failed"
    assert "Job ID <strong>job-1</strong>" in email_kwargs["content"]
    assert "Cohort ID" not in email_kwargs["content"]
    assert not bam_path.exists(), "cleanup must still run after a pipeline failure"


def test_pipeline_failure_email_names_the_cohort_when_one_is_set(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """The failure email names the cohort when the job belongs to one.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])
    _subprocess_stub(monkeypatch, tasks, pipeline_error=pipeline_error)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(
            tasks, **_job_kwargs(tmp_path, bam_path, email="user@example.com", cohort_key="cohort:my-cohort")
        )

    email_kwargs = no_email_task.delay.call_args.kwargs
    assert "Cohort ID: <strong>my-cohort</strong>" in email_kwargs["content"]


def test_pipeline_failure_without_an_email_sends_no_notification(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """No `email` argument means no failure notification is enqueued.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])
    _subprocess_stub(monkeypatch, tasks, pipeline_error=pipeline_error)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, email=None))

    no_email_task.delay.assert_not_called()


def test_archive_failure_is_marked_failed_and_the_directory_is_left_in_place(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A `shutil.make_archive` failure fails the job without removing `output_dir`.

    `shutil.rmtree(output_dir)` runs immediately after `make_archive` in the
    same `try`; if the archive step itself fails, the results directory must
    still be there afterwards -- otherwise a failed archive attempt loses
    the very results it was trying to package.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    monkeypatch.setattr(tasks.shutil, "make_archive", MagicMock(side_effect=OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert output_dir.exists(), "a failed archive attempt must not remove the results it could not package"
    assert ("usage:job-1", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]


# ---------------------------------------------------------------------------
# run_vntyper_job: success paths
# ---------------------------------------------------------------------------


def test_a_successful_job_is_marked_completed_and_cleans_up_its_inputs(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A clean run marks the job completed and removes its own input files.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    redis_mocks.usage.hset.assert_any_call("usage:job-1", "status", "completed")
    assert not any(c.args == ("usage:job-1", "status", "failed") for c in redis_mocks.usage.hset.call_args_list)
    redis_mocks.cohort.sadd.assert_not_called()
    no_email_task.delay.assert_not_called()
    assert not bam_path.exists()
    assert not (tmp_path / "input" / "job-1").exists()
    redis_mocks.queue.lrem.assert_called_once_with("vntyper_job_queue", 0, "task-1")


def test_a_clean_job_removes_its_input_directory_and_fires_no_leftover_warning(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """#201: the `os.rmdir` at the end of cleanup is reachable again.

    The pipeline used to write `<alignment>.quickcheck.log` into the job's input
    directory, so cleanup always found that directory non-empty: the `os.rmdir`
    never ran in production, one directory and inode leaked per job, and the
    "still holds files" warning fired on 100% of jobs -- which is exactly when a
    warning stops being able to report a genuinely unexpected leftover.

    The subprocess stand-in here runs the pipeline's own startup validation
    rather than merely recording the command, so this test fails against the
    unfixed code instead of passing vacuously.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub_that_runs_the_pipelines_own_startup(monkeypatch, tasks)
    caplog.set_level(logging.WARNING, logger="app.tasks")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    job_input_dir = tmp_path / "input" / "job-1"
    assert not job_input_dir.exists(), (
        f"the input directory was left behind, holding {sorted(p.name for p in job_input_dir.iterdir())}"
    )
    assert "still holds files" not in caplog.text
    redis_mocks.usage.hset.assert_any_call("usage:job-1", "status", "completed")
    # Non-vacuity: the stand-in really did run the pipeline's validation, and the
    # log it produced is still there to diagnose a corrupt upload with -- it is
    # kept rather than deleted, deliberately.
    assert (tmp_path / "output" / "job-1" / "sample.bam.quickcheck.log").exists()


def test_the_leftover_warning_still_fires_for_something_the_job_did_not_create(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The guard must keep working -- it must only stop firing on the normal path.

    Same realistic stand-in as the test above, plus one file the job has no claim
    to. Without this, "no warning on a clean job" would also be satisfied by
    deleting the guard.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    bystander = bam_path.parent / "someone_elses.txt"
    bystander.write_bytes(b"not this job's to delete")
    _subprocess_stub_that_runs_the_pipelines_own_startup(monkeypatch, tasks)
    caplog.set_level(logging.WARNING, logger="app.tasks")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert bystander.read_bytes() == b"not this job's to delete"
    assert f"Input directory {bam_path.parent} still holds files and was left in place" in caplog.text


def test_a_successful_job_with_email_sends_a_completion_notice(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """Providing `email` sends a completion notice naming the job and download link.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, email="user@example.com"))

    no_email_task.delay.assert_called_once()
    email_kwargs = no_email_task.delay.call_args.kwargs
    assert email_kwargs["to_email"] == "user@example.com"
    assert email_kwargs["subject"] == "VNtyper Job Completed Successfully"
    assert "Job ID: <strong>job-1</strong>" in email_kwargs["content"]
    assert f"{tasks.settings.API_BASE_URL}/api/download/job-1/" in email_kwargs["content"]


def test_a_successful_job_with_a_cohort_joins_it_and_extends_its_retention(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A completed job in a cohort is added to the cohort's job set, and the
    cohort's retention is extended -- both only on success.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)
    retention_spy = MagicMock(name="extend_cohort_retention")
    monkeypatch.setattr(tasks, "extend_cohort_retention", retention_spy)

    _invoke_vntyper_job(
        tasks, **_job_kwargs(tmp_path, bam_path, email="user@example.com", cohort_key="cohort:my-cohort")
    )

    redis_mocks.cohort.sadd.assert_called_once_with("cohort:my-cohort:jobs", "job-1")
    retention_spy.assert_called_once_with(
        redis_mocks.cohort, "cohort:my-cohort", tasks.settings.COHORT_RETENTION_DAYS * 86400
    )
    email_kwargs = no_email_task.delay.call_args.kwargs
    assert "Cohort ID: <strong>my-cohort</strong>" in email_kwargs["content"]


def test_optional_flags_are_all_appended_and_a_successful_archive_replaces_the_directory(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """Every optional flag reaches the command, and a successful archive zips
    the results and removes the original directory.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    commands = _subprocess_stub(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")

    _invoke_vntyper_job(
        tasks,
        **_job_kwargs(
            tmp_path,
            bam_path,
            fast_mode=True,
            keep_intermediates=True,
            archive_results=True,
            advntr_mode=True,
        ),
    )

    pipeline_command = commands[1]
    assert pipeline_command == [
        "conda",
        "run",
        "-n",
        "vntyper",
        "vntyper",
        "pipeline",
        "--bam",
        str(bam_path),
        "-o",
        str(output_dir),
        "--threads",
        "1",
        "--reference-assembly",
        "hg38",
        "--fast-mode",
        "--keep-intermediates",
        "--archive-results",
        "--extra-modules",
        "advntr",
        "--advntr-max-coverage",
        "300",
    ]
    assert not output_dir.exists(), "a successful archive removes the original results directory"
    assert Path(f"{output_dir}.zip").exists(), "a successful archive leaves the zip behind"


# ---------------------------------------------------------------------------
# run_vntyper_job: index_path=None compatibility default
# ---------------------------------------------------------------------------


def test_index_path_none_falls_back_to_the_conventional_bai_name_and_skips_rebuilding(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """With `index_path=None`, an index already at `f"{bam_path}.bai"` is used as-is.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    conventional_index = Path(f"{bam_path}.bai")
    conventional_index.write_bytes(INDEX_BYTES)
    commands = _subprocess_stub(monkeypatch, tasks)

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, index_path=None))

    assert not any(c[:2] == ["samtools", "index"] for c in commands), (
        "an index already present at the fallback path must not be rebuilt"
    )
    assert not conventional_index.exists(), "cleanup must still remove the index it fell back to"


def test_index_path_none_still_builds_an_index_when_none_exists(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """With `index_path=None` and no index on disk, one is generated at the
    conventional name.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    commands = _subprocess_stub(monkeypatch, tasks)

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, index_path=None))

    assert commands[0] == ["samtools", "index", str(bam_path)]
    assert not Path(f"{bam_path}.bai").exists(), "the generated index is cleaned up like any other job input"


# ---------------------------------------------------------------------------
# run_vntyper_job: cleanup edge cases
# ---------------------------------------------------------------------------


def test_cleanup_tolerates_input_paths_and_a_directory_that_were_never_created(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """Cleanup does not raise when the alignment, its index, and the input
    directory itself never existed on disk.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path = tmp_path / "input" / "job-1" / "sample.bam"  # never created; no job_input_dir either

    def _run(command, *args, **kwargs):
        return None

    monkeypatch.setattr(tasks.subprocess, "run", _run)

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    redis_mocks.usage.hset.assert_any_call("usage:job-1", "status", "completed")


def test_cleanup_leaves_a_nonempty_input_directory_and_warns_instead_of_deleting_it(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """An input directory holding something other than this job's own files
    is left in place, with a warning, rather than removed.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    unexpected_file = bam_path.parent / "unexpected.txt"
    unexpected_file.write_bytes(b"not this job's to delete")
    _subprocess_stub(monkeypatch, tasks)
    caplog.set_level(logging.WARNING, logger="app.tasks")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert bam_path.parent.exists(), "a directory holding an unexpected file must not be removed"
    assert unexpected_file.exists()
    assert f"Input directory {bam_path.parent} still holds files and was left in place" in caplog.text


def test_cleanup_logs_but_does_not_raise_when_removing_the_alignment_fails(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A filesystem error deleting the alignment is logged, not raised -- it
    must not mask a job that otherwise completed (or failed) cleanly.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)
    real_remove = tasks.os.remove

    def _flaky_remove(path, *args, **kwargs):
        if path == str(bam_path):
            raise OSError("permission denied")
        return real_remove(path, *args, **kwargs)

    monkeypatch.setattr(tasks.os, "remove", _flaky_remove)
    caplog.set_level(logging.ERROR, logger="app.tasks")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))  # must not raise

    assert bam_path.exists(), "the removal genuinely failed; the file is still there"
    assert f"Error deleting input file {bam_path}: permission denied" in caplog.text


def test_cleanup_logs_but_does_not_raise_when_removing_the_input_directory_fails(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A filesystem error removing the now-empty input directory is logged,
    not raised.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    job_input_dir = bam_path.parent
    _subprocess_stub(monkeypatch, tasks)
    real_rmdir = tasks.os.rmdir

    def _flaky_rmdir(path, *args, **kwargs):
        if path == str(job_input_dir):
            raise OSError("directory busy")
        return real_rmdir(path, *args, **kwargs)

    monkeypatch.setattr(tasks.os, "rmdir", _flaky_rmdir)
    caplog.set_level(logging.ERROR, logger="app.tasks")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))  # must not raise

    assert job_input_dir.exists(), "the removal genuinely failed; the directory is still there"
    assert f"Error deleting input directory {job_input_dir}: directory busy" in caplog.text


# ---------------------------------------------------------------------------
# delete_old_results: paths test_result_expiry.py does not reach
# ---------------------------------------------------------------------------


def test_delete_old_results_logs_but_does_not_raise_when_removal_of_an_aged_archive_fails(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A filesystem error deleting an aged-out archive is logged, and the
    sweep still completes instead of dying on its first failure.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory used as the output tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(tmp_path))
    redis_mocks.cohort.scan_iter.return_value = iter([])
    stale = tmp_path / "stale.zip"
    stale.write_bytes(b"result data")
    too_old = (tasks.settings.MAX_RESULT_AGE_DAYS + 1) * 86400
    import time

    monkeypatch.setattr(tasks.os.path, "getctime", lambda path: time.time() - too_old)
    monkeypatch.setattr(tasks.os, "remove", MagicMock(side_effect=OSError("permission denied")))
    caplog.set_level(logging.ERROR, logger="app.tasks")

    tasks.delete_old_results.run()  # must not raise

    assert stale.exists(), "the removal genuinely failed; the archive is still there"
    assert f"Error deleting file {stale}: permission denied" in caplog.text


def test_delete_old_results_removes_data_for_a_cohort_that_has_expired(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """A cohort key that no longer `exists()` (its Redis TTL has passed) has
    its members' results and Redis records removed, one archive is kept
    because it was too recent for the age sweep, and a missing member
    archive does not stop the other Redis records from being cleared.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory used as the output tree.
    """
    from app import tasks

    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(tmp_path))
    fresh_zip = tmp_path / "job-a.zip"
    fresh_zip.write_bytes(b"result data")  # freshly created: the age sweep must leave it alone
    redis_mocks.cohort.scan_iter.return_value = iter(["cohort:expired"])
    redis_mocks.cohort.exists.return_value = False
    redis_mocks.cohort.smembers.return_value = {"job-a", "job-b"}

    tasks.delete_old_results.run()

    assert not fresh_zip.exists(), "an expired cohort's member archive is removed even though it was not aged out"
    redis_mocks.client_calls = redis_mocks.queue.delete.call_args_list
    deleted = {c.args[0] for c in redis_mocks.queue.delete.call_args_list}
    assert deleted == {"job-a", "celery-task-meta-job-a", "job-b", "celery-task-meta-job-b"}
    redis_mocks.cohort.delete.assert_called_once_with("cohort:expired:jobs")


# ---------------------------------------------------------------------------
# run_cohort_analysis_job: paths test_cohort_analysis_task.py does not reach
# ---------------------------------------------------------------------------


def test_cohort_analysis_archive_failure_is_marked_failed_and_reraised(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """A `shutil.make_archive` failure after the `vntyper cohort` subprocess
    succeeds still fails the job and updates the usage record.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the output tree.
    """
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    monkeypatch.setattr(tasks.subprocess, "run", lambda *a, **k: None)
    monkeypatch.setattr(tasks.shutil, "make_archive", MagicMock(side_effect=OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        _invoke_cohort_job(
            tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(tmp_path / "analysis")
        )

    assert ("usage:analysis", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]


def test_cohort_analysis_skips_retention_extension_when_no_cohort_id_is_given(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """An empty `cohort_id` means there is no cohort record to extend.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the output tree.
    """
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    monkeypatch.setattr(tasks.subprocess, "run", lambda *a, **k: None)
    retention_spy = MagicMock(name="extend_cohort_retention")
    monkeypatch.setattr(tasks, "extend_cohort_retention", retention_spy)

    _invoke_cohort_job(tasks, cohort_id="", zip_paths=[str(zip_path)], output_dir=str(tmp_path / "analysis"))

    retention_spy.assert_not_called()


def test_cohort_analysis_logs_but_does_not_raise_when_deleting_its_scratch_file_fails(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A filesystem error removing the task's own `cohort_input.txt` scratch
    file is logged, not raised.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the output tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    monkeypatch.setattr(tasks.subprocess, "run", lambda *a, **k: None)
    monkeypatch.setattr(tasks.os, "remove", MagicMock(side_effect=OSError("permission denied")))
    caplog.set_level(logging.ERROR, logger="app.tasks")

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(tmp_path / "analysis"))

    input_file = tmp_path / "analysis" / "cohort_input.txt"
    assert input_file.exists(), "the removal genuinely failed; the scratch file is still there"
    assert f"Error deleting cohort input file {input_file}: permission denied" in caplog.text


def test_cohort_analysis_logs_but_does_not_raise_when_removing_the_empty_output_dir_fails(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A filesystem error removing the now-empty output directory is logged,
    not raised.

    The output directory becomes empty once the task deletes its own scratch
    file, because the members' `.zip` archives it read are outside `output_dir`.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the output tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    monkeypatch.setattr(tasks.subprocess, "run", lambda *a, **k: None)
    output_dir = tmp_path / "analysis"
    real_rmdir = tasks.os.rmdir

    def _flaky_rmdir(path, *args, **kwargs):
        if path == str(output_dir):
            raise OSError("directory busy")
        return real_rmdir(path, *args, **kwargs)

    monkeypatch.setattr(tasks.os, "rmdir", _flaky_rmdir)
    caplog.set_level(logging.ERROR, logger="app.tasks")

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert output_dir.exists(), "the removal genuinely failed; the directory is still there"
    assert f"Error deleting directory {output_dir}: directory busy" in caplog.text
