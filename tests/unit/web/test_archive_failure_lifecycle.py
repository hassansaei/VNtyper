"""Failed archive publication and quarantine lifecycle coverage."""

import hashlib
import logging
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

pytestmark = pytest.mark.unit

BAM_BYTES = b"alignment-bytes"
PATIENT_BYTES = b"patient-alignment-bytes-that-must-never-enter-an-archive"


@pytest.fixture
def redis_mocks(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Replace the worker's three Redis clients with independent mocks."""
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
    """Replace the email Celery task so tests never reach its broker."""
    from app import tasks

    mock_task = MagicMock(name="send_email_task")
    monkeypatch.setattr(tasks, "send_email_task", mock_task)
    return mock_task


def _job_kwargs(tmp_path: Path, bam_path: Path, **overrides) -> dict:
    """Build the complete argument set for an individual worker run."""
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_metadata = bam_path.parent.stat()
    output_metadata = output_dir.stat()
    alignment_metadata = bam_path.stat()
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
    return kwargs


def _make_job_input(tmp_path: Path) -> Path:
    """Create the individual worker's alignment input."""
    job_input_dir = tmp_path / "input" / "job-1"
    job_input_dir.mkdir(parents=True)
    bam_path = job_input_dir / "sample.bam"
    bam_path.write_bytes(BAM_BYTES)
    return bam_path


def _invoke_vntyper_job(tasks, **kwargs) -> None:
    """Invoke the bound individual task with a worker request context."""
    tasks.run_vntyper_job.push_request(id="task-1")
    try:
        tasks.run_vntyper_job.run(**kwargs)
    finally:
        tasks.run_vntyper_job.pop_request()


def _invoke_cohort_job(tasks, **kwargs) -> None:
    """Invoke the bound cohort task with a worker request context."""
    tasks.run_cohort_analysis_job.push_request(id="analysis-1")
    try:
        tasks.run_cohort_analysis_job.run(**kwargs)
    finally:
        tasks.run_cohort_analysis_job.pop_request()


def _stub_pipeline(monkeypatch: pytest.MonkeyPatch, tasks) -> None:
    """Replace the individual pipeline subprocess with a successful no-op."""
    monkeypatch.setattr(tasks.subprocess, "run", lambda *args, **kwargs: None)


def _write_cohort_result(monkeypatch: pytest.MonkeyPatch, tasks, output_dir: Path) -> None:
    """Replace the cohort subprocess with one complete scientific result."""

    def write_cohort_result(command, *args, **kwargs):
        child_output = Path(command[command.index("-o") + 1])
        (child_output / "cohort_result.tsv").write_bytes(b"complete cohort result")

    monkeypatch.setattr(tasks.subprocess, "run", write_cohort_result)


def test_first_quarantine_failure_retries_without_leaving_a_regular_public_archive(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A transient quarantine failure cannot leave a failed job downloadable."""
    from app import tasks

    bam_path = _make_job_input(tmp_path)
    _stub_pipeline(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    cleanup_failure = OSError("result directory busy")

    def partially_remove_then_fail(_workspace) -> None:
        (output_dir / "result.txt").unlink()
        raise cleanup_failure

    real_quarantine = tasks.quarantine_archive
    quarantine_calls = 0

    def fail_first_quarantine(*args, **kwargs):
        nonlocal quarantine_calls
        quarantine_calls += 1
        if quarantine_calls == 1:
            raise OSError("first quarantine denied")
        return real_quarantine(*args, **kwargs)

    monkeypatch.setattr(tasks.PipelineJobWorkspace, "remove_output", partially_remove_then_fail)
    monkeypatch.setattr(tasks, "quarantine_archive", fail_first_quarantine)
    caplog.set_level(logging.ERROR, logger="app.tasks")

    with pytest.raises(OSError) as raised:
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert raised.value is cleanup_failure
    assert "first quarantine denied" in caplog.text
    public_archive = Path(f"{output_dir}.zip")
    assert not public_archive.exists()
    quarantines = list(output_dir.parent.glob(".job-1.zip.failed-*"))
    assert len(quarantines) == 1
    with zipfile.ZipFile(quarantines[0]) as archive:
        assert archive.read("result.txt") == b"result data"


@pytest.mark.parametrize("fatal_step", ["completed", "cohort", "email"])
def test_post_publish_failure_quarantines_the_individual_complete_archive(
    fatal_step: str,
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
) -> None:
    """Bookkeeping failures hide but preserve the complete scientific output."""
    from app import tasks

    bam_path = _make_job_input(tmp_path)
    _stub_pipeline(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    failure = RuntimeError(f"{fatal_step} failed")

    if fatal_step == "completed":

        def fail_completed(*args, **kwargs):
            if args == ("usage:job-1", "status", "completed"):
                raise failure
            return None

        redis_mocks.usage.hset.side_effect = fail_completed
    elif fatal_step == "cohort":
        redis_mocks.cohort.sadd.side_effect = failure
    else:
        no_email_task.delay.side_effect = failure

    with pytest.raises(RuntimeError) as raised:
        _invoke_vntyper_job(
            tasks,
            **_job_kwargs(
                tmp_path,
                bam_path,
                archive_results=True,
                email="user@example.com" if fatal_step == "email" else None,
                cohort_key="cohort:my-cohort" if fatal_step == "cohort" else None,
            ),
        )

    assert raised.value is failure
    assert not Path(f"{output_dir}.zip").exists()
    quarantines = list(output_dir.parent.glob(".job-1.zip.failed-*"))
    assert len(quarantines) == 1
    with zipfile.ZipFile(quarantines[0]) as archive:
        assert archive.read("result.txt") == b"result data"
    assert ("usage:job-1", "status", "failed") in [call.args for call in redis_mocks.usage.hset.call_args_list]


def test_post_publish_alias_cleanup_does_not_follow_target_or_mask_the_individual_primary_error(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """An archive alias is unlinked without touching its target or masking failure."""
    from app import tasks

    bam_path = _make_job_input(tmp_path)
    _stub_pipeline(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    patient = tmp_path / "external-patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    primary = RuntimeError("completed state unavailable")

    def fail_completed(*args, **kwargs):
        if args == ("usage:job-1", "status", "completed"):
            archive = Path(f"{output_dir}.zip")
            archive.unlink()
            archive.symlink_to(patient)
            raise primary
        return None

    redis_mocks.usage.hset.side_effect = fail_completed
    caplog.set_level(logging.ERROR, logger="app.tasks")

    with pytest.raises(RuntimeError) as raised:
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert raised.value is primary
    assert "symbolic link" in caplog.text
    public_archive = Path(f"{output_dir}.zip")
    assert not public_archive.exists() and not public_archive.is_symlink()
    assert patient.read_bytes() == PATIENT_BYTES


def test_post_publish_hardlink_cleanup_removes_only_the_individual_public_name(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    no_email_task: MagicMock,
    tmp_path: Path,
) -> None:
    """A hard-linked public alias is revoked without changing its external name."""
    from app import tasks

    bam_path = _make_job_input(tmp_path)
    _stub_pipeline(monkeypatch, tasks)
    output_dir = tmp_path / "output" / "job-1"
    output_dir.mkdir(parents=True)
    (output_dir / "result.txt").write_bytes(b"result data")
    external = tmp_path / "external-patient.bam"
    external.write_bytes(PATIENT_BYTES)
    external_inode = external.stat().st_ino
    primary = RuntimeError("completed state unavailable")

    def fail_completed(*args, **kwargs):
        if args == ("usage:job-1", "status", "completed"):
            archive = Path(f"{output_dir}.zip")
            archive.unlink()
            archive.hardlink_to(external)
            raise primary
        return None

    redis_mocks.usage.hset.side_effect = fail_completed

    with pytest.raises(RuntimeError) as raised:
        _invoke_vntyper_job(tasks, **_job_kwargs(tmp_path, bam_path, archive_results=True))

    assert raised.value is primary
    assert not Path(f"{output_dir}.zip").exists()
    assert external.read_bytes() == PATIENT_BYTES
    assert external.stat().st_ino == external_inode


def test_cohort_completed_state_failure_quarantines_the_complete_archive(
    monkeypatch: pytest.MonkeyPatch, redis_mocks: SimpleNamespace, tmp_path: Path
) -> None:
    """A failed cohort transition hides but preserves its complete output."""
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    output_dir = tmp_path / "analysis"
    _write_cohort_result(monkeypatch, tasks, output_dir)
    failure = RuntimeError("completed state unavailable")

    def fail_completed(*args, **kwargs):
        if args == ("usage:analysis", "status", "completed"):
            raise failure
        return None

    redis_mocks.usage.hset.side_effect = fail_completed

    with pytest.raises(RuntimeError) as raised:
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert raised.value is failure
    assert not (tmp_path / "analysis.zip").exists()
    quarantines = list(tmp_path.glob(".analysis.zip.failed-*"))
    assert len(quarantines) == 1
    with zipfile.ZipFile(quarantines[0]) as archive:
        assert archive.read("cohort_result.tsv") == b"complete cohort result"
    assert zip_path.read_bytes() == b"result data"
    assert ("usage:analysis", "status", "failed") in [call.args for call in redis_mocks.usage.hset.call_args_list]


def test_cohort_post_publish_alias_cleanup_does_not_follow_target_or_mask_the_primary_error(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A cohort archive alias is removed without following it or masking failure."""
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    output_dir = tmp_path / "analysis"
    patient = tmp_path / "external-patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    _write_cohort_result(monkeypatch, tasks, output_dir)
    primary = RuntimeError("completed state unavailable")

    def fail_completed(*args, **kwargs):
        if args == ("usage:analysis", "status", "completed"):
            archive = Path(f"{output_dir}.zip")
            archive.unlink()
            archive.symlink_to(patient)
            raise primary
        return None

    redis_mocks.usage.hset.side_effect = fail_completed
    caplog.set_level(logging.ERROR, logger="app.tasks")

    with pytest.raises(RuntimeError) as raised:
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert raised.value is primary
    assert "symbolic link" in caplog.text
    public_archive = tmp_path / "analysis.zip"
    assert not public_archive.exists() and not public_archive.is_symlink()
    assert patient.read_bytes() == PATIENT_BYTES


def test_cohort_post_publish_hardlink_cleanup_removes_only_the_public_name(
    monkeypatch: pytest.MonkeyPatch,
    redis_mocks: SimpleNamespace,
    tmp_path: Path,
) -> None:
    """A cohort hard-link alias is revoked without changing its external name."""
    from app import tasks

    zip_path = tmp_path / "job-a.zip"
    zip_path.write_bytes(b"result data")
    output_dir = tmp_path / "analysis"
    external = tmp_path / "external-patient.bam"
    external.write_bytes(PATIENT_BYTES)
    external_inode = external.stat().st_ino
    _write_cohort_result(monkeypatch, tasks, output_dir)
    primary = RuntimeError("completed state unavailable")

    def fail_completed(*args, **kwargs):
        if args == ("usage:analysis", "status", "completed"):
            archive = Path(f"{output_dir}.zip")
            archive.unlink()
            archive.hardlink_to(external)
            raise primary
        return None

    redis_mocks.usage.hset.side_effect = fail_completed

    with pytest.raises(RuntimeError) as raised:
        _invoke_cohort_job(tasks, cohort_id="cohort-1", zip_paths=[str(zip_path)], output_dir=str(output_dir))

    assert raised.value is primary
    assert not Path(f"{output_dir}.zip").exists()
    assert external.read_bytes() == PATIENT_BYTES
    assert external.stat().st_ino == external_inode
