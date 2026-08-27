"""Secondary failure bookkeeping must never replace a pipeline exception."""

import hashlib
import logging
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

pytestmark = pytest.mark.unit


@pytest.fixture
def worker(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Provide isolated Redis, subprocess, and email boundaries for the task."""
    from app import failure_reporting, tasks

    queue = MagicMock(name="redis_client")
    cohort = MagicMock(name="redis_cohort_client")
    usage = MagicMock(name="redis_usage_client")
    email = MagicMock(name="send_email_task")
    monkeypatch.setattr(tasks, "redis_client", queue)
    monkeypatch.setattr(tasks, "redis_cohort_client", cohort)
    monkeypatch.setattr(tasks, "redis_usage_client", usage)
    monkeypatch.setattr(tasks, "send_email_task", email)
    return SimpleNamespace(
        tasks=tasks,
        reporting=failure_reporting,
        queue=queue,
        cohort=cohort,
        usage=usage,
        email=email,
    )


def _input_and_kwargs(tmp_path: Path, **overrides: object) -> tuple[Path, dict[str, object]]:
    """Create one worker-owned BAM and its complete task argument mapping."""
    input_dir = tmp_path / "input" / "job-1"
    input_dir.mkdir(parents=True)
    bam_path = input_dir / "sample.bam"
    bam_path.write_bytes(b"alignment")
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    input_metadata = input_dir.stat()
    output_metadata = output_dir.stat()
    alignment_metadata = bam_path.stat()
    kwargs: dict[str, object] = {
        "bam_path": str(bam_path),
        "output_dir": str(output_dir),
        "thread": 1,
        "reference_assembly": "hg38",
        "fast_mode": False,
        "keep_intermediates": False,
        "archive_results": False,
        "email": "user@example.com",
        "cohort_key": None,
        "client_ip": "127.0.0.1",
        "user_agent": "pytest",
        "advntr_mode": False,
        "index_path": None,
        "capacity_reserved": False,
        "workspace_identity": {
            "input_dir": [input_metadata.st_dev, input_metadata.st_ino],
            "output_dir": [output_metadata.st_dev, output_metadata.st_ino],
            "alignment": [alignment_metadata.st_dev, alignment_metadata.st_ino],
            "alignment_sha256": hashlib.sha256(bam_path.read_bytes()).hexdigest(),
            "index": None,
            "index_sha256": None,
        },
    }
    kwargs.update(overrides)
    return bam_path, kwargs


def _invoke(worker: SimpleNamespace, kwargs: dict[str, object]) -> None:
    """Run the bound Celery task with a deterministic request identifier."""
    worker.tasks.run_vntyper_job.push_request(id="task-1")
    try:
        worker.tasks.run_vntyper_job.run(**kwargs)
    finally:
        worker.tasks.run_vntyper_job.pop_request()


def _raise_primary(monkeypatch: pytest.MonkeyPatch, worker: SimpleNamespace, primary: BaseException) -> None:
    """Make only the external pipeline execution raise ``primary``."""
    monkeypatch.setattr(
        worker.tasks.subprocess,
        "run",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(primary),
    )


@pytest.mark.parametrize(
    "primary",
    [
        pytest.param(subprocess.CalledProcessError(9, ["vntyper", "pipeline"]), id="nonzero"),
        pytest.param(subprocess.TimeoutExpired(["vntyper", "pipeline"], 123.5), id="timeout"),
    ],
)
def test_failed_status_write_preserves_each_pipeline_exception(
    monkeypatch: pytest.MonkeyPatch,
    worker: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    primary: BaseException,
) -> None:
    """A Redis outage on the inner failed-status write cannot become the task error."""
    bam_path, kwargs = _input_and_kwargs(tmp_path)
    secondary = RuntimeError("failed-status Redis unavailable")
    worker.usage.hset.side_effect = [None, secondary, None]
    _raise_primary(monkeypatch, worker, primary)

    with caplog.at_level(logging.ERROR, logger=worker.tasks.logger.name), pytest.raises(type(primary)) as raised:
        _invoke(worker, kwargs)

    assert raised.value is primary
    assert "Error recording failed status for job-1: failed-status Redis unavailable" in caplog.text
    assert not bam_path.exists()


def test_preflight_metadata_read_failure_preserves_primary_and_does_not_skip_email(
    monkeypatch: pytest.MonkeyPatch,
    worker: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Unreadable optional metadata is logged before notification continues."""
    bam_path, kwargs = _input_and_kwargs(tmp_path, cohort_key="cohort:cohort-1")
    primary = subprocess.CalledProcessError(17, ["vntyper", "pipeline"])
    monkeypatch.setattr(
        worker.reporting,
        "read_preflight_failure",
        MagicMock(side_effect=OSError("metadata unreadable")),
    )
    _raise_primary(monkeypatch, worker, primary)

    with caplog.at_level(logging.ERROR, logger=worker.tasks.logger.name), pytest.raises(type(primary)) as raised:
        _invoke(worker, kwargs)

    assert raised.value is primary
    worker.email.delay.assert_called_once()
    assert "Cohort ID: <strong>cohort-1</strong>" in worker.email.delay.call_args.kwargs["content"]
    assert "Error recording preflight failure for job-1: metadata unreadable" in caplog.text
    assert not bam_path.exists()


def test_preflight_metadata_write_failure_preserves_primary_and_does_not_skip_email(
    monkeypatch: pytest.MonkeyPatch,
    worker: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A rejected metadata hash update is independent of failure notification."""
    bam_path, kwargs = _input_and_kwargs(tmp_path)
    primary = subprocess.TimeoutExpired(["vntyper", "pipeline"], 456.5)
    metadata_error = RuntimeError("metadata Redis unavailable")
    worker.usage.hset.side_effect = [None, None, metadata_error, None]
    monkeypatch.setattr(
        worker.reporting,
        "read_preflight_failure",
        MagicMock(return_value={"code": "reference_unresolved", "message": "Unable to resolve reference."}),
    )
    _raise_primary(monkeypatch, worker, primary)

    with caplog.at_level(logging.ERROR, logger=worker.tasks.logger.name), pytest.raises(type(primary)) as raised:
        _invoke(worker, kwargs)

    assert raised.value is primary
    worker.email.delay.assert_called_once()
    assert "Error recording preflight failure for job-1: metadata Redis unavailable" in caplog.text
    assert not bam_path.exists()


def test_failure_email_rendering_failure_preserves_timeout(
    monkeypatch: pytest.MonkeyPatch,
    worker: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A renderer bug cannot replace the exact timeout raised by the runner."""
    bam_path, kwargs = _input_and_kwargs(tmp_path, cohort_key="cohort:cohort-1")
    primary = subprocess.TimeoutExpired(["vntyper", "pipeline"], 789.5)
    monkeypatch.setattr(
        worker.reporting,
        "build_failure_email_content",
        MagicMock(side_effect=ValueError("renderer failed")),
    )
    _raise_primary(monkeypatch, worker, primary)

    with caplog.at_level(logging.ERROR, logger=worker.tasks.logger.name), pytest.raises(type(primary)) as raised:
        _invoke(worker, kwargs)

    assert raised.value is primary
    worker.email.delay.assert_not_called()
    assert "Error queueing failure email for job-1: renderer failed" in caplog.text
    assert not bam_path.exists()


@pytest.mark.parametrize(
    ("primary", "cohort_key"),
    [
        pytest.param(subprocess.CalledProcessError(9, ["vntyper", "pipeline"]), None, id="nonzero"),
        pytest.param(
            subprocess.TimeoutExpired(["vntyper", "pipeline"], 321.5),
            "cohort:cohort-1",
            id="timeout-cohort",
        ),
    ],
)
def test_failure_email_enqueue_failure_preserves_each_pipeline_exception(
    monkeypatch: pytest.MonkeyPatch,
    worker: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    primary: BaseException,
    cohort_key: str | None,
) -> None:
    """A broker failure is logged while the runner's exact exception survives."""
    bam_path, kwargs = _input_and_kwargs(tmp_path, cohort_key=cohort_key)
    worker.email.delay.side_effect = RuntimeError("broker unavailable")
    _raise_primary(monkeypatch, worker, primary)

    with caplog.at_level(logging.ERROR, logger=worker.tasks.logger.name), pytest.raises(type(primary)) as raised:
        _invoke(worker, kwargs)

    assert raised.value is primary
    assert "Error queueing failure email for job-1: broker unavailable" in caplog.text
    assert ("usage:job-1", "status", "failed") in [call.args for call in worker.usage.hset.call_args_list]
    assert not bam_path.exists()
