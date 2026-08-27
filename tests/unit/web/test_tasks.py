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

import hashlib
import inspect
import json
import logging
import os
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import ANY, MagicMock, patch

import pytest
from redis.exceptions import ConnectionError as RedisConnectionError

from vntyper.scripts.utils import validate_bam_file

pytestmark = pytest.mark.unit

BAM_BYTES = b"alignment-bytes"
INDEX_BYTES = b"index-bytes"
PATIENT_BYTES = b"patient-alignment-bytes-that-must-never-enter-an-archive"


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
    reservation_pipeline = mocks.queue.pipeline.return_value.__enter__.return_value
    reservation_pipeline.hexists.return_value = True
    reservation_pipeline.zscore.return_value = float("inf")
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
    output_dir = tmp_path / "output" / job_id
    output_dir.mkdir(parents=True, exist_ok=True)
    supplied_index = overrides.get("index_path")
    identity_index = Path(supplied_index) if supplied_index else None
    input_metadata = bam_path.parent.stat()
    output_metadata = output_dir.stat()
    alignment_metadata = bam_path.stat()
    index_metadata = identity_index.stat() if identity_index is not None and identity_index.exists() else None
    kwargs = {
        "bam_path": str(bam_path),
        "output_dir": str(output_dir),
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
        "capacity_reserved": True,
        "workspace_identity": {
            "input_dir": [input_metadata.st_dev, input_metadata.st_ino],
            "output_dir": [output_metadata.st_dev, output_metadata.st_ino],
            "alignment": [alignment_metadata.st_dev, alignment_metadata.st_ino],
            "alignment_sha256": hashlib.sha256(bam_path.read_bytes()).hexdigest(),
            "index": None if index_metadata is None else [index_metadata.st_dev, index_metadata.st_ino],
            "index_sha256": (
                None
                if index_metadata is None or identity_index is None
                else hashlib.sha256(identity_index.read_bytes()).hexdigest()
            ),
        },
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
    pipeline_error: Exception | None = None,
) -> list:
    """Replace `tasks.subprocess.run` with a recorder that can be told to fail.

    Args:
        monkeypatch: Standard pytest fixture; restores the patch at teardown.
        tasks: The imported `app.tasks` module.
        pipeline_error: If given, raised instead of running `vntyper pipeline`.

    Returns:
        list: The argument vector of every command issued, in order.
    """
    commands: list[list[str]] = []

    def _run(command, *args, **kwargs):
        commands.append(list(command))
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
# run_vntyper_job: pipeline subprocess exits non-zero
# ---------------------------------------------------------------------------


def test_missing_index_is_deferred_to_pipeline_preflight_without_writing_beside_upload(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """The worker must let pipeline preflight build its run-local index.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked email task.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    commands = _subprocess_stub(monkeypatch, tasks)

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert len(commands) == 1
    assert commands[0][:7] == ["conda", "run", "--no-capture-output", "-n", "vntyper", "vntyper", "pipeline"]
    assert not any(command[:2] == ["samtools", "index"] for command in commands)


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

    assert commands[0][:7] == ["conda", "run", "--no-capture-output", "-n", "vntyper", "vntyper", "pipeline"]
    assert ("usage:job-1", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]
    no_email_task.delay.assert_called_once()
    email_kwargs = no_email_task.delay.call_args.kwargs
    assert email_kwargs["to_email"] == "user@example.com"
    assert email_kwargs["subject"] == "VNtyper Job Failed"
    assert "Job ID <strong>job-1</strong>" in email_kwargs["content"]
    assert "Cohort ID" not in email_kwargs["content"]
    assert not bam_path.exists(), "cleanup must still run after a pipeline failure"


def test_pipeline_failure_stores_the_curated_preflight_code_and_message(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A nonzero pipeline exit transports its reviewed artifact into the job hash.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked email task.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    artifact = {
        "code": "reference_unresolved",
        "message": "Unable to resolve CRAM reference: contig=chr1, M5=digest.",
        "candidates": [["cli", "full.fa", "probe exited non-zero"]],
    }
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])

    def _run_pipeline_and_emit_current_artifact(command, *args, **kwargs):
        (output_dir / "preflight_error.json").write_text(json.dumps(artifact), encoding="utf-8")
        raise pipeline_error

    monkeypatch.setattr(tasks.subprocess, "run", _run_pipeline_and_emit_current_artifact)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    redis_mocks.usage.hset.assert_any_call(
        "usage:job-1",
        mapping={"code": "reference_unresolved", "message": artifact["message"]},
    )


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


def test_pipeline_failure_email_escapes_values_without_exposing_argv_or_worker_paths(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The HTML notification carries only escaped, recipient-actionable facts."""
    from app import tasks

    job_id = 'job-<b>&"1'
    cohort_id = 'cohort-<img src=x onerror="alert(1)">&'
    bam_path, _ = _make_job_input(tmp_path, job_id=job_id)
    private_command = [
        "conda-secret",
        "--token",
        "/opt/vntyper/input/private/sample.bam",
        "-o",
        "/opt/vntyper/output/private",
    ]

    class HostileReturnCode(int):
        def __str__(self) -> str:
            return '<status&"9>'

        def __format__(self, format_spec: str) -> str:
            del format_spec
            return str(self)

    pipeline_error = subprocess.CalledProcessError(HostileReturnCode(9), private_command)
    _subprocess_stub(monkeypatch, tasks, pipeline_error=pipeline_error)

    with caplog.at_level(logging.ERROR), pytest.raises(subprocess.CalledProcessError) as raised:
        _invoke_vntyper_job(
            tasks,
            **_job_kwargs(
                tmp_path,
                bam_path,
                job_id=job_id,
                email="user@example.com",
                cohort_key=f"cohort:{cohort_id}",
            ),
        )

    assert raised.value is pipeline_error
    content = no_email_task.delay.call_args.kwargs["content"]
    assert "job-&lt;b&gt;&amp;&quot;1" in content
    assert "cohort-&lt;img src=x onerror=&quot;alert(1)&quot;&gt;&amp;" in content
    assert "exited with status &lt;status&amp;&quot;9&gt;" in content
    assert "<pre>" not in content
    assert "conda-secret" not in content
    assert "--token" not in content
    assert str(bam_path) not in content
    assert str(tmp_path) not in content
    assert "/opt/vntyper/input" not in content
    assert "/opt/vntyper/output" not in content
    assert "conda-secret" in caplog.text
    assert "/opt/vntyper/input/private/sample.bam" in caplog.text


def test_pipeline_timeout_email_uses_the_actual_deadline_without_adding_a_timeout(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A runner timeout is reported without this task inventing a pipeline deadline."""
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    timeout_error = subprocess.TimeoutExpired(
        cmd=["vntyper", "pipeline", "--bam", "/opt/vntyper/input/private/sample.bam"],
        timeout=3723.5,
        output="partial output from /opt/vntyper/output/private",
        stderr="failed while reading /opt/vntyper/input/private/sample.bam",
    )
    runner = MagicMock(name="subprocess.run", side_effect=timeout_error)
    monkeypatch.setattr(tasks.subprocess, "run", runner)

    with pytest.raises(subprocess.TimeoutExpired) as raised:
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, email="user@example.com"))

    assert raised.value is timeout_error
    assert runner.call_args.kwargs == {"check": True}
    assert ("usage:job-1", "status", "failed") in [call.args for call in redis_mocks.usage.hset.call_args_list]
    assert not bam_path.exists()
    content = no_email_task.delay.call_args.kwargs["content"]
    assert "did not finish within 3723.5 seconds" in content
    assert "['vntyper', 'pipeline'" not in content
    assert "/opt/vntyper/input" not in content
    assert "/opt/vntyper/output" not in content
    assert "partial output" not in content
    assert str(bam_path) not in content
    assert "Cohort ID" not in content


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
    """An archive creation failure fails the job without removing `output_dir`.

    `shutil.rmtree(output_dir)` runs immediately after the archive helper in the
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
    monkeypatch.setattr(tasks, "create_safe_archive", MagicMock(side_effect=OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert output_dir.exists(), "a failed archive attempt must not remove the results it could not package"
    assert ("usage:job-1", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]


def test_result_directory_cleanup_failure_removes_public_alias_without_following_target(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A failed post-archive cleanup cannot leave a current public download."""
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    patient = tmp_path / "external-patient.bam"
    patient.write_bytes(PATIENT_BYTES)

    def replace_archive_and_fail(_workspace) -> None:
        archive = Path(f"{output_dir}.zip")
        archive.unlink()
        archive.symlink_to(patient)
        raise OSError("result directory busy")

    monkeypatch.setattr(tasks.PipelineJobWorkspace, "remove_output", replace_archive_and_fail)

    with pytest.raises(OSError, match="result directory busy"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    public_archive = Path(f"{output_dir}.zip")
    assert not public_archive.exists() and not public_archive.is_symlink()
    assert patient.read_bytes() == PATIENT_BYTES
    assert ("usage:job-1", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]


def test_partial_result_cleanup_quarantines_one_complete_archive(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
) -> None:
    """Partial source deletion cannot destroy both complete representations."""
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"complete result")

    def partially_remove_then_fail(_workspace) -> None:
        (output_dir / "result.txt").unlink()
        raise OSError("partial result cleanup")

    monkeypatch.setattr(tasks.PipelineJobWorkspace, "remove_output", partially_remove_then_fail)

    with pytest.raises(OSError, match="partial result cleanup"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert not Path(f"{output_dir}.zip").exists()
    quarantines = list(output_dir.parent.glob(".job-1.zip.failed-*"))
    assert len(quarantines) == 1
    with zipfile.ZipFile(quarantines[0]) as archive:
        assert archive.read("result.txt") == b"complete result"


def test_worker_archive_refuses_a_real_symlink_without_reading_or_deleting_its_target(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """The worker's second archive pass must fail closed on an alignment view.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        no_email_task: The mocked `send_email_task`.
        tmp_path: Scratch directory standing in for the job tree.
    """
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    _subprocess_stub(monkeypatch, tasks)
    patient_alignment = tmp_path / "patient-source.bam"
    patient_alignment.write_bytes(PATIENT_BYTES)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    alignment_view = output_dir / "alignment_view.bam"
    alignment_view.symlink_to(patient_alignment)

    with pytest.raises(ValueError, match="symbolic link.*alignment_view\\.bam"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert output_dir.exists(), "a rejected archive must leave the result tree for diagnosis"
    assert alignment_view.is_symlink()
    assert not Path(f"{output_dir}.zip").exists(), "an unsafe archive must never be installed"
    assert patient_alignment.read_bytes() == PATIENT_BYTES
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
    caplog.set_level(logging.WARNING, logger="app.pipeline_job_workspace")

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
    caplog.set_level(logging.WARNING, logger="app.pipeline_job_workspace")

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
    # `cohort_retention_days()`, not the raw `COHORT_RETENTION_DAYS`: a cohort must not
    # outlive the archives it lists, so the configured value is clamped to
    # `MAX_RESULT_AGE_DAYS`. Asserting the raw constant here would pass while the task
    # extended a cohort past its members' expiry, which is the defect being closed.
    retention_spy.assert_called_once_with(
        redis_mocks.cohort, "cohort:my-cohort", tasks.settings.cohort_retention_days() * 86400
    )
    assert tasks.settings.cohort_retention_days() <= tasks.settings.MAX_RESULT_AGE_DAYS
    email_kwargs = no_email_task.delay.call_args.kwargs
    assert "Cohort ID: <strong>my-cohort</strong>" in email_kwargs["content"]
    # The completion email stated no deadline at all, which became misleading once the window
    # shortened. It says "about N days" rather than naming a hard cutoff, because
    # MAX_RESULT_AGE_DAYS is a cleanup-eligibility threshold and not an enforced deadline:
    # the sweep runs daily and /download/ performs no age check, so a job finishing just
    # after the sweep is retrievable for nearly N+1 days.
    # Whitespace-normalised: the body is a wrapped triple-quoted f-string, so phrases span
    # newlines and a raw substring test would pass or fail on the wrap position.
    content = " ".join(email_kwargs["content"].split())
    assert f"about {tasks.settings.MAX_RESULT_AGE_DAYS} days" in content
    assert "stops working after" not in content, "the service enforces no such deadline"
    assert "anyone with this link can retrieve them" in content, "the link is a capability, and must say so"


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

    pipeline_command = commands[0]
    assert pipeline_command == [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "vntyper",
        "vntyper",
        "pipeline",
        "--bam",
        pipeline_command[8],
        "-o",
        pipeline_command[10],
        "--threads",
        "1",
        "--reference-assembly",
        "hg38",
        "--fast-mode",
        "--keep-intermediates",
        "--extra-modules",
        "advntr",
        "--advntr-max-coverage",
        "300",
    ]
    assert Path(pipeline_command[8]).parent.parent == Path("/tmp")
    assert not Path(pipeline_command[8]).is_symlink()
    assert Path(pipeline_command[8]).name == bam_path.name
    assert pipeline_command[10].startswith(f"/proc/{os.getpid()}/fd/")
    assert not output_dir.exists(), "a successful archive removes the original results directory"
    assert Path(f"{output_dir}.zip").exists(), "a successful archive leaves the zip behind"


def test_failed_retry_removes_stale_public_archive_without_following_its_target(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A stale symlink disappears before the subprocess while its patient target is untouched."""
    from app import tasks

    bam_path, _ = _make_job_input(tmp_path)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    patient = tmp_path / "external-patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    stale = Path(f"{output_dir}.zip")
    stale.symlink_to(patient)
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])
    _subprocess_stub(monkeypatch, tasks, pipeline_error=pipeline_error)

    with pytest.raises(subprocess.CalledProcessError):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert not stale.exists() and not stale.is_symlink()
    assert patient.read_bytes() == PATIENT_BYTES


# ---------------------------------------------------------------------------
# run_vntyper_job: index_path=None compatibility default
# ---------------------------------------------------------------------------


def test_index_path_none_does_not_adopt_an_unreported_conventional_bai(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """With `index_path=None`, a conventional neighbour is outside worker authority.

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
    assert conventional_index.read_bytes() == INDEX_BYTES


def test_index_path_none_defers_a_missing_index_to_pipeline_preflight(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """With no uploaded index, the worker still writes nothing beside the input.

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

    assert not any(command[:2] == ["samtools", "index"] for command in commands)
    assert not Path(f"{bam_path}.bai").exists()


# ---------------------------------------------------------------------------
# run_vntyper_job: cleanup edge cases
# ---------------------------------------------------------------------------


def test_worker_rejects_a_task_without_the_api_workspace_identity(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, no_email_task: MagicMock, tmp_path: Path
) -> None:
    """A path-only task cannot bypass the API-to-worker identity boundary.

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

    with pytest.raises(RuntimeError, match="missing workspace identity token"):
        _invoke_vntyper_job(
            tasks,
            bam_path=str(bam_path),
            output_dir=str(tmp_path / "output" / "job-1"),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
        )

    redis_mocks.usage.hset.assert_any_call("usage:job-1", "status", "failed")


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
    caplog.set_level(logging.WARNING, logger="app.pipeline_job_workspace")

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
    from app import pipeline_job_workspace

    real_unlink = pipeline_job_workspace.os.unlink
    input_identity = (bam_path.parent.stat().st_dev, bam_path.parent.stat().st_ino)

    def _flaky_unlink(path, *args, **kwargs):
        dir_fd = kwargs.get("dir_fd")
        if path == bam_path.name and dir_fd is not None:
            metadata = os.fstat(dir_fd)
            if (metadata.st_dev, metadata.st_ino) == input_identity:
                raise OSError("permission denied")
        return real_unlink(path, *args, **kwargs)

    monkeypatch.setattr(pipeline_job_workspace.os, "unlink", _flaky_unlink)
    caplog.set_level(logging.ERROR, logger="app.pipeline_job_workspace")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))  # must not raise

    assert bam_path.exists(), "the removal genuinely failed; the file is still there"
    assert f"Error deleting input file {bam_path}: permission denied" in caplog.text


def test_input_cleanup_logs_one_error_and_attempts_every_owned_path(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """One blocked owned input does not prevent cleanup of its other indexes."""
    from app import tasks

    attempted: list[str] = []
    monkeypatch.setattr(tasks.os.path, "exists", lambda _path: True)

    def remove(path: str) -> None:
        attempted.append(path)
        if path.endswith("sample.bam"):
            raise OSError("blocked")

    monkeypatch.setattr(tasks.os, "remove", remove)

    with caplog.at_level(logging.ERROR, logger=tasks.logger.name):
        tasks.remove_job_input_files("/jobs/1/sample.bam", "/jobs/1/sample.bam.bai")

    assert attempted == list(
        dict.fromkeys(
            (
                "/jobs/1/sample.bam",
                "/jobs/1/sample.bam.bai",
                *tasks.derived_index_paths("/jobs/1/sample.bam"),
            )
        )
    )
    errors = [
        record for record in caplog.records if record.name == tasks.logger.name and record.levelno >= logging.ERROR
    ]
    assert [record.levelno for record in errors] == [logging.ERROR]
    assert "blocked" in caplog.text


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
    """An archive creation failure after the `vntyper cohort` subprocess
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
    monkeypatch.setattr(tasks, "create_safe_archive", MagicMock(side_effect=OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        _invoke_cohort_job(
            tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(tmp_path / "analysis")
        )

    assert ("usage:analysis", "status", "failed") in [c.args for c in redis_mocks.usage.hset.call_args_list]


def test_cohort_analysis_refuses_a_pre_existing_output_symlink(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """A planted cohort destination cannot redirect child output into another directory."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"member archive")
    victim = tmp_path / "victim"
    victim.mkdir()
    (victim / "keep.txt").write_bytes(PATIENT_BYTES)
    output_dir = tmp_path / "analysis"
    output_dir.symlink_to(victim, target_is_directory=True)

    def write_child_result(command, *args, **kwargs):
        child_output = Path(command[command.index("-o") + 1])
        (child_output / "cohort_result.tsv").write_bytes(b"redirected child output")

    monkeypatch.setattr(tasks.subprocess, "run", write_child_result)

    with pytest.raises(RuntimeError, match="freshly created"):
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(output_dir))

    assert (victim / "keep.txt").read_bytes() == PATIENT_BYTES
    assert not (victim / "cohort_result.tsv").exists()
    assert output_dir.is_symlink()
    assert not Path(f"{output_dir}.zip").exists()


def test_cohort_analysis_refuses_a_planted_input_symlink_in_an_existing_output_directory(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """An occupied output directory is rejected before a planted child link is opened."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"member archive")
    victim = tmp_path / "victim.txt"
    victim.write_bytes(PATIENT_BYTES)
    output_dir = tmp_path / "analysis"
    output_dir.mkdir()
    planted_input = output_dir / "cohort_input.txt"
    planted_input.symlink_to(victim)
    child = MagicMock(name="cohort subprocess")
    monkeypatch.setattr(tasks.subprocess, "run", child)

    with pytest.raises(RuntimeError, match="freshly created"):
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(output_dir))

    assert victim.read_bytes() == PATIENT_BYTES
    assert planted_input.is_symlink()
    assert planted_input.resolve() == victim
    child.assert_not_called()
    assert not Path(f"{output_dir}.zip").exists()


def test_cohort_workspace_never_adopts_a_replacement_in_the_creation_window(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """Replacing a just-created public entry cannot redirect cohort task writes."""
    from app import cohort_workspace as cohort_workspace_module
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"member archive")
    output_dir = tmp_path / "analysis"
    displaced_entry = tmp_path / "analysis-created-by-task"
    replacement = tmp_path / "replacement"
    replacement.mkdir()
    (replacement / "keep.txt").write_bytes(PATIENT_BYTES)
    real_open = cohort_workspace_module.os.open
    swapped = False

    def swap_public_entry_during_first_open(path, flags, *args, **kwargs):
        nonlocal swapped
        is_public_entry = path == output_dir.name and kwargs.get("dir_fd") is not None
        is_atomic_reservation = bool(flags & cohort_workspace_module.os.O_CREAT)
        if is_public_entry and is_atomic_reservation:
            descriptor = real_open(path, flags, *args, **kwargs)
        elif is_public_entry:
            descriptor = None
        else:
            return real_open(path, flags, *args, **kwargs)
        if not swapped:
            swapped = True
            parent_descriptor = kwargs["dir_fd"]
            cohort_workspace_module.os.rename(
                output_dir.name,
                displaced_entry.name,
                src_dir_fd=parent_descriptor,
                dst_dir_fd=parent_descriptor,
            )
            cohort_workspace_module.os.rename(
                replacement.name,
                output_dir.name,
                src_dir_fd=parent_descriptor,
                dst_dir_fd=parent_descriptor,
            )
        if descriptor is not None:
            return descriptor
        return real_open(path, flags, *args, **kwargs)

    def write_child_result(command, *args, **kwargs):
        child_output = Path(command[command.index("-o") + 1])
        (child_output / "cohort_result.tsv").write_bytes(b"task-owned result")

    monkeypatch.setattr(cohort_workspace_module.os, "open", swap_public_entry_during_first_open)
    monkeypatch.setattr(tasks.subprocess, "run", write_child_result)

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(output_dir))

    assert displaced_entry.is_file(), "the atomically opened reservation remains bound after its name is moved"
    assert output_dir.is_dir(), "the replacement stays under the name chosen by the concurrent actor"
    assert (output_dir / "keep.txt").read_bytes() == PATIENT_BYTES
    assert not (output_dir / "cohort_result.tsv").exists()
    with zipfile.ZipFile(f"{output_dir}.zip") as archive:
        assert archive.read("cohort_result.tsv") == b"task-owned result"


def test_fresh_workspace_refusal_preserves_an_existing_public_archive(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """A redelivered task cannot erase its prior archive before ownership is established."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"member archive")
    output_dir = tmp_path / "analysis"
    output_dir.mkdir()
    (output_dir / "prior-result.tsv").write_bytes(b"prior result directory")
    archive = Path(f"{output_dir}.zip")
    prior_archive = b"complete prior public archive"
    archive.write_bytes(prior_archive)
    child = MagicMock(name="cohort subprocess")
    monkeypatch.setattr(tasks.subprocess, "run", child)

    with pytest.raises(RuntimeError, match="freshly created"):
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(output_dir))

    assert archive.read_bytes() == prior_archive
    assert (output_dir / "prior-result.tsv").read_bytes() == b"prior result directory"
    child.assert_not_called()


def test_cohort_analysis_keeps_child_writes_bound_when_output_directory_is_swapped(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """Replacing the public directory name cannot redirect writes or archive reads."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"member archive")
    output_dir = tmp_path / "analysis"
    displaced_output = tmp_path / "analysis-displaced"
    victim = tmp_path / "victim"
    victim.mkdir()
    (victim / "keep.txt").write_bytes(PATIENT_BYTES)

    def swap_output_before_child_writes(command, *args, **kwargs):
        output_dir.rename(displaced_output)
        output_dir.symlink_to(victim, target_is_directory=True)
        child_output = Path(command[command.index("-o") + 1])
        (child_output / "cohort_result.tsv").write_bytes(b"bound child output")

    monkeypatch.setattr(tasks.subprocess, "run", swap_output_before_child_writes)

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(output_dir))

    assert (victim / "keep.txt").read_bytes() == PATIENT_BYTES
    assert not (victim / "cohort_result.tsv").exists()
    assert displaced_output.is_file(), "the moved public reservation is never used as task output"
    with zipfile.ZipFile(f"{output_dir}.zip") as archive:
        assert archive.read("cohort_result.tsv") == b"bound child output"


def test_cohort_workspace_cleanup_uses_the_bound_staging_directory_after_its_name_is_swapped(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """Cleanup empties task staging through its fd and never traverses a replacement."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"member archive")
    output_dir = tmp_path / "analysis"
    replacement = tmp_path / "replacement"
    replacement.mkdir()
    (replacement / "keep.txt").write_bytes(PATIENT_BYTES)
    displaced_paths: list[Path] = []
    replacement_paths: list[Path] = []

    def swap_staging_before_child_writes(command, *args, **kwargs):
        child_output = Path(command[command.index("-o") + 1])
        staging_path = Path(child_output.readlink())
        displaced = staging_path.with_name(f"{staging_path.name}-displaced")
        staging_path.rename(displaced)
        replacement.rename(staging_path)
        displaced_paths.append(displaced)
        replacement_paths.append(staging_path)
        (child_output / "cohort_result.tsv").write_bytes(b"bound child output")

    monkeypatch.setattr(tasks.subprocess, "run", swap_staging_before_child_writes)

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(output_dir))

    assert (replacement_paths[0] / "keep.txt").read_bytes() == PATIENT_BYTES
    assert not (replacement_paths[0] / "cohort_result.tsv").exists()
    assert displaced_paths[0].is_dir()
    assert list(displaced_paths[0].iterdir()) == []
    with zipfile.ZipFile(f"{output_dir}.zip") as archive:
        assert archive.read("cohort_result.tsv") == b"bound child output"


def test_cohort_analysis_consumes_a_bound_snapshot_when_member_path_is_replaced(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """The cohort subprocess consumes task-owned bytes, never a later path replacement."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"original member archive")
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    consumed: list[bytes] = []
    consumed_names: list[str] = []

    def replace_member_before_consumption(command, *args, **kwargs):
        member.unlink()
        member.symlink_to(patient)
        input_file = Path(command[command.index("--input-file") + 1])
        listed = Path(input_file.read_text().strip())
        consumed_names.append(listed.name)
        consumed.append(listed.read_bytes())

    monkeypatch.setattr(tasks.subprocess, "run", replace_member_before_consumption)

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(member)], output_dir=str(tmp_path / "analysis"))

    assert consumed == [b"original member archive"]
    assert consumed_names == ["job-a.zip"]
    assert member.is_symlink()
    assert patient.read_bytes() == PATIENT_BYTES


def test_cohort_analysis_consumes_the_copied_bytes_when_snapshot_directory_name_is_replaced(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """The cohort child reopens snapshots through the directory inode, not its old name."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"original member archive")
    output_dir = tmp_path / "analysis"
    snapshot_dirs: list[Path] = []
    displaced_snapshot_dirs: list[Path] = []
    consumed = tmp_path / "consumed.bin"
    real_subprocess_run = subprocess.run

    def replace_snapshot_directory_before_child_opens_input(command, *args, **kwargs):
        input_file = Path(command[command.index("--input-file") + 1])
        listed = Path(input_file.read_text().strip())
        snapshot_dir = Path(listed.parent.readlink())
        displaced_snapshot_dir = snapshot_dir.with_name(f"{snapshot_dir.name}-displaced")
        snapshot_dirs.append(snapshot_dir)
        displaced_snapshot_dirs.append(displaced_snapshot_dir)
        snapshot_dir.rename(displaced_snapshot_dir)
        snapshot_dir.mkdir()
        (snapshot_dir / member.name).write_bytes(PATIENT_BYTES)
        reader = (
            "from pathlib import Path; import sys; "
            "listed = Path(sys.argv[1]).read_text().strip(); "
            "Path(sys.argv[2]).write_bytes(Path(listed).read_bytes())"
        )
        real_subprocess_run([sys.executable, "-c", reader, str(input_file), str(consumed)], check=True)

    monkeypatch.setattr(tasks.subprocess, "run", replace_snapshot_directory_before_child_opens_input)

    with pytest.raises(ValueError, match="snapshot directory changed"):
        _invoke_cohort_job(
            tasks,
            cohort_id="cohort-1",
            zip_paths=[str(member)],
            output_dir=str(output_dir),
        )

    assert consumed.read_bytes() == b"original member archive"
    assert not (displaced_snapshot_dirs[0] / member.name).exists()
    assert not snapshot_dirs[0].exists()
    assert not Path(f"{output_dir}.zip").exists()


def test_cohort_snapshot_cleanup_never_masks_the_subprocess_failure_after_directory_replacement(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A replaced snapshot name is secondary when the cohort child already failed."""
    from app import tasks

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"original member archive")
    output_dir = tmp_path / "analysis"
    snapshot_dirs: list[Path] = []
    displaced_snapshot_dirs: list[Path] = []
    primary_failure = subprocess.CalledProcessError(9, ["vntyper", "cohort"])

    def replace_snapshot_directory_then_fail(command, *args, **kwargs):
        input_file = Path(command[command.index("--input-file") + 1])
        listed = Path(input_file.read_text().strip())
        snapshot_dir = Path(listed.parent.readlink())
        displaced_snapshot_dir = snapshot_dir.with_name(f"{snapshot_dir.name}-displaced")
        snapshot_dirs.append(snapshot_dir)
        displaced_snapshot_dirs.append(displaced_snapshot_dir)
        snapshot_dir.rename(displaced_snapshot_dir)
        snapshot_dir.mkdir()
        (snapshot_dir / member.name).write_bytes(PATIENT_BYTES)
        raise primary_failure

    monkeypatch.setattr(tasks.subprocess, "run", replace_snapshot_directory_then_fail)
    caplog.set_level(logging.ERROR, logger="app.archive_delivery")

    with pytest.raises(subprocess.CalledProcessError) as raised:
        _invoke_cohort_job(
            tasks,
            cohort_id="cohort-1",
            zip_paths=[str(member)],
            output_dir=str(output_dir),
        )

    assert raised.value is primary_failure
    assert not (displaced_snapshot_dirs[0] / member.name).exists()
    assert not snapshot_dirs[0].exists()
    assert "snapshot directory changed" in caplog.text


def test_cohort_snapshot_copy_failure_removes_the_partial_task_owned_copy(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failed copy leaves neither its partial bytes nor its private directory behind."""
    from app import archive_delivery

    member = tmp_path / "job-a.zip"
    member.write_bytes(b"original member archive")
    snapshot_dir = tmp_path / ".cohort-members"
    primary_failure = OSError("snapshot copy failed")

    def copy_part_then_fail(source, target, *args, **kwargs):
        target.write(source.read(5))
        raise primary_failure

    monkeypatch.setattr(archive_delivery.shutil, "copyfileobj", copy_part_then_fail)

    with pytest.raises(OSError) as raised:
        archive_delivery.snapshot_owned_archives([str(member)], snapshot_dir)

    assert raised.value is primary_failure
    assert not snapshot_dir.exists()


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


def test_cohort_analysis_fails_before_archiving_when_deleting_its_scratch_file_fails(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A path-bearing scratch file must never be packaged when cleanup fails.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the output tree.
        caplog: Captures the task's log record.
    """
    from app import tasks
    from app.cohort_workspace import CohortWorkspace

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    monkeypatch.setattr(tasks.subprocess, "run", lambda *a, **k: None)
    monkeypatch.setattr(CohortWorkspace, "unlink", MagicMock(side_effect=OSError("permission denied")))
    caplog.set_level(logging.ERROR, logger="app.tasks")

    with pytest.raises(OSError, match="permission denied"):
        _invoke_cohort_job(
            tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(tmp_path / "analysis")
        )

    assert not (tmp_path / "analysis").exists()
    assert not (tmp_path / "analysis.zip").exists()
    assert "permission denied" in caplog.text


@pytest.mark.parametrize("cleanup_target", ["staging", "reservation"])
def test_cohort_analysis_cleanup_failure_marks_failed_preserves_archive_and_allows_retry(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    cleanup_target: str,
) -> None:
    """A transient workspace-cleanup error fails safely and leaves the next attempt runnable.

    Args:
        monkeypatch: Standard pytest fixture.
        redis_mocks: The three mocked Redis clients.
        tmp_path: Scratch directory standing in for the output tree.
        caplog: Captures the task's log record.
    """
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")

    def write_cohort_result(command, *args, **kwargs):
        child_output = Path(command[command.index("-o") + 1])
        (child_output / "cohort_result.tsv").write_bytes(b"complete cohort result")

    monkeypatch.setattr(tasks.subprocess, "run", write_cohort_result)
    output_dir = tmp_path / "analysis"
    from app import cohort_workspace as cohort_workspace_module

    cleanup_failure = OSError(f"transient {cleanup_target} cleanup failure")
    if cleanup_target == "staging":
        real_clear_directory = cohort_workspace_module._clear_directory
        cleanup_calls = 0

        def _flaky_clear_directory(descriptor):
            nonlocal cleanup_calls
            cleanup_calls += 1
            if cleanup_calls == 1:
                raise cleanup_failure
            return real_clear_directory(descriptor)

        monkeypatch.setattr(cohort_workspace_module, "_clear_directory", _flaky_clear_directory)
    else:
        real_unlink = cohort_workspace_module.os.unlink
        cleanup_calls = 0

        def _flaky_unlink(path, *args, **kwargs):
            nonlocal cleanup_calls
            if path == output_dir.name and kwargs.get("dir_fd") is not None:
                cleanup_calls += 1
                if cleanup_calls == 1:
                    raise cleanup_failure
            return real_unlink(path, *args, **kwargs)

        monkeypatch.setattr(cohort_workspace_module.os, "unlink", _flaky_unlink)
    caplog.set_level(logging.ERROR, logger="app.cohort_workspace")

    with pytest.raises(OSError) as raised:
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert raised.value is cleanup_failure
    assert ("usage:analysis", "status", "failed") in [call.args for call in redis_mocks.usage.hset.call_args_list]
    assert not output_dir.exists(), "the transient cleanup must converge before descriptors are released"
    assert not Path(f"{output_dir}.zip").exists(), "a failed task must not leave a public downloadable archive"
    quarantines = list(tmp_path.glob(".analysis.zip.failed-*"))
    assert len(quarantines) == 1
    with zipfile.ZipFile(quarantines[0]) as archive:
        assert archive.read("cohort_result.tsv") == b"complete cohort result"
    assert "transient" in caplog.text

    _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert Path(f"{output_dir}.zip").is_file()
    assert not output_dir.exists()


def test_cohort_workspace_cleanup_never_masks_an_active_subprocess_failure(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The task error remains authoritative when descriptor cleanup also fails transiently."""
    from app import cohort_workspace as cohort_workspace_module
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    output_dir = tmp_path / "analysis"
    primary_failure = subprocess.CalledProcessError(17, ["vntyper", "cohort"])
    cleanup_failure = OSError("transient staging cleanup failure")
    real_clear_directory = cohort_workspace_module._clear_directory
    cleanup_calls = 0

    def fail_subprocess(*args, **kwargs):
        raise primary_failure

    def _flaky_clear_directory(descriptor):
        nonlocal cleanup_calls
        cleanup_calls += 1
        if cleanup_calls == 1:
            raise cleanup_failure
        return real_clear_directory(descriptor)

    monkeypatch.setattr(tasks.subprocess, "run", fail_subprocess)
    monkeypatch.setattr(cohort_workspace_module, "_clear_directory", _flaky_clear_directory)
    caplog.set_level(logging.ERROR)

    with pytest.raises(subprocess.CalledProcessError) as raised:
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert raised.value is primary_failure
    assert ("usage:analysis", "status", "failed") in [call.args for call in redis_mocks.usage.hset.call_args_list]
    assert not output_dir.exists()
    assert not list(tmp_path.glob(".analysis.cohort-*"))
    assert "workspace cleanup also failed" in caplog.text


def test_pipeline_worker_releases_its_capacity_reservation(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, fake_redis, tmp_path: Path
) -> None:
    """Every completed pipeline task returns its shared queue and byte capacity."""
    from app import admission, tasks

    monkeypatch.setattr(tasks, "redis_client", fake_redis)
    _subprocess_stub(monkeypatch, tasks)
    bam_path, _ = _make_job_input(tmp_path)
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, "job-1", 100)
    fake_redis.zadd(admission.ADMISSION_DEADLINES, {"job-1": time.time() + 60})

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


def test_cohort_worker_releases_its_capacity_reservation(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, fake_redis, tmp_path: Path
) -> None:
    """Cohort tasks release through the same idempotent worker-cleanup path."""
    from app import admission, tasks

    monkeypatch.setattr(tasks, "redis_client", fake_redis)
    monkeypatch.setattr(tasks.subprocess, "run", lambda *args, **kwargs: None)
    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, "analysis", 100)

    _invoke_cohort_job(
        tasks,
        cohort_id="cohort-1",
        zip_paths=[str(zip_path)],
        output_dir=str(tmp_path / "analysis"),
    )

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


@pytest.mark.parametrize("failing_call", ["hset", "expire"])
def test_pipeline_worker_startup_failure_releases_capacity_without_masking_primary(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    fake_redis,
    tmp_path: Path,
    failing_call: str,
) -> None:
    """Either initial usage write remains inside pipeline cleanup ownership."""
    from app import admission, tasks

    monkeypatch.setattr(tasks, "redis_client", fake_redis)
    bam_path, _ = _make_job_input(tmp_path)
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, "job-1", 100)
    fake_redis.zadd(admission.ADMISSION_DEADLINES, {"job-1": time.time() + 60})
    primary = RuntimeError(f"initial {failing_call} failed")
    if failing_call == "hset":
        redis_mocks.usage.hset.side_effect = [primary, None]
    else:
        redis_mocks.usage.expire.side_effect = primary

    with pytest.raises(RuntimeError, match=f"initial {failing_call} failed"):
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


@pytest.mark.parametrize("failing_call", ["hset", "expire"])
def test_cohort_worker_startup_failure_releases_capacity_without_masking_primary(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    fake_redis,
    tmp_path: Path,
    failing_call: str,
) -> None:
    """Cohort initial usage writes are covered by the task's finally block."""
    from app import admission, tasks

    monkeypatch.setattr(tasks, "redis_client", fake_redis)
    archive = tmp_path / "member.zip"
    archive.write_bytes(b"member")
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, "analysis", 100)
    fake_redis.zadd(admission.ADMISSION_DEADLINES, {"analysis": time.time() + 60})
    primary = RuntimeError(f"initial {failing_call} failed")
    if failing_call == "hset":
        redis_mocks.usage.hset.side_effect = [primary, None]
    else:
        redis_mocks.usage.expire.side_effect = primary

    with pytest.raises(RuntimeError, match=f"initial {failing_call} failed"):
        _invoke_cohort_job(
            tasks,
            cohort_id="cohort-1",
            zip_paths=[str(archive)],
            output_dir=str(tmp_path / "analysis"),
            capacity_reserved=True,
        )

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


@pytest.mark.parametrize("job_type", ["pipeline", "cohort"])
def test_worker_execution_failure_stops_heartbeat_and_releases_capacity(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    fake_redis,
    tmp_path: Path,
    job_type: str,
) -> None:
    """Both task bodies own heartbeat and admission release on execution errors."""
    from app import admission, tasks

    events: list[str] = []

    class FakeHeartbeat:
        def __init__(self, *_args, **_kwargs):
            events.append("created")

        def start(self):
            events.append("started")

        def stop(self):
            events.append("stopped")

    monkeypatch.setattr(admission, "ReservationHeartbeat", FakeHeartbeat)
    monkeypatch.setattr(tasks, "redis_client", fake_redis)
    failure = subprocess.CalledProcessError(1, ["vntyper", job_type])
    monkeypatch.setattr(tasks.subprocess, "run", lambda *_args, **_kwargs: (_ for _ in ()).throw(failure))
    reservation_id = "job-1" if job_type == "pipeline" else "analysis"
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, reservation_id, 100)

    with pytest.raises(subprocess.CalledProcessError):
        if job_type == "pipeline":
            bam_path, _ = _make_job_input(tmp_path)
            _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))
        else:
            archive = tmp_path / "member.zip"
            archive.write_bytes(b"member")
            _invoke_cohort_job(
                tasks,
                cohort_id="cohort-1",
                zip_paths=[str(archive)],
                output_dir=str(tmp_path / "analysis"),
                capacity_reserved=True,
            )

    assert events == ["created", "started", "stopped"]
    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}


@pytest.mark.parametrize("job_type", ["pipeline", "cohort"])
def test_heartbeat_thread_start_failure_preserves_error_and_releases_capacity(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    fake_redis,
    tmp_path: Path,
    job_type: str,
) -> None:
    """A failed heartbeat bootstrap cannot mask itself or skip either task's release."""
    from app import admission, tasks

    class FailingThread:
        def __init__(self, **_kwargs):
            pass

        def start(self):
            raise RuntimeError("heartbeat thread failed to start")

        def join(self, timeout=None):
            pytest.fail("an unstarted heartbeat thread must not be joined")

    monkeypatch.setattr(admission, "Thread", FailingThread)
    monkeypatch.setattr(tasks, "redis_client", fake_redis)
    reservation_id = "job-1" if job_type == "pipeline" else "analysis"
    now = time.time()
    fake_redis.hset(admission.ADMISSION_RESERVATIONS, reservation_id, 100)
    fake_redis.zadd(admission.ADMISSION_DEADLINES, {reservation_id: now + 60})

    with pytest.raises(RuntimeError, match="heartbeat thread failed to start"):
        if job_type == "pipeline":
            bam_path, _ = _make_job_input(tmp_path)
            _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))
        else:
            archive = tmp_path / "member.zip"
            archive.write_bytes(b"member")
            _invoke_cohort_job(
                tasks,
                cohort_id="cohort-1",
                zip_paths=[str(archive)],
                output_dir=str(tmp_path / "analysis"),
                capacity_reserved=True,
            )

    assert fake_redis.hgetall(admission.ADMISSION_RESERVATIONS) == {}
    assert fake_redis.zscore(admission.ADMISSION_DEADLINES, reservation_id) is None


def test_pipeline_worker_logs_a_capacity_release_failure_without_failing_the_result(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A Redis outage during release cannot turn completed work into a failure."""
    from app import admission, tasks

    _subprocess_stub(monkeypatch, tasks)
    bam_path, _ = _make_job_input(tmp_path)

    def fail_release(*_args, **_kwargs):
        raise RedisConnectionError("redis unavailable")

    monkeypatch.setattr(admission, "release_capacity", fail_release)
    caplog.set_level(logging.ERROR, logger="app.tasks")

    _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path))

    assert "Error releasing capacity reservation for job-1: redis unavailable" in caplog.text


def test_cohort_worker_logs_a_capacity_release_failure_without_failing_the_result(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The cohort cleanup path treats admission release as best-effort too."""
    from app import admission, tasks

    monkeypatch.setattr(tasks.subprocess, "run", lambda *args, **kwargs: None)
    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")

    def fail_release(*_args, **_kwargs):
        raise RedisConnectionError("redis unavailable")

    monkeypatch.setattr(admission, "release_capacity", fail_release)
    caplog.set_level(logging.ERROR, logger="app.tasks")

    _invoke_cohort_job(
        tasks,
        cohort_id="cohort-1",
        zip_paths=[str(zip_path)],
        output_dir=str(tmp_path / "analysis"),
    )

    assert "Error releasing capacity reservation for analysis: redis unavailable" in caplog.text
