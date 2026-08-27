"""Output binding and real pipeline-contract handoff tests."""

import json
import os
import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tests.unit.web.pipeline_job_inode_handoff_support import _digest, _identity, _invoke_job

pytestmark = pytest.mark.unit


def test_descriptor_bound_paths_pass_the_real_pipeline_ownership_boundary(tmp_path: Path) -> None:
    """The worker handoff stays compatible with VNtyper's pre-write guard."""
    from app.pipeline_job_workspace import open_pipeline_job_workspace

    from tests.support.pipeline_harness import MINIMAL_CONFIG
    from vntyper.scripts.alignment_target_io import protect_alignment_inputs

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"accepted-alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(b"accepted-alignment"),
        "index": None,
        "index_sha256": None,
    }
    workspace = open_pipeline_job_workspace(str(alignment), None, str(output_dir), identity)
    try:
        protect_alignment_inputs(
            workspace.bound_output_path,
            workspace.bound_alignment_path,
            "bam",
            None,
            None,
            MINIMAL_CONFIG,
            "hg38",
        )
    finally:
        workspace.remove_views()
        workspace.close()


def test_child_failure_still_clears_an_output_inode_displaced_during_execution(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Detached partial results are cleared while the child error remains primary."""
    from app import tasks

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"accepted-alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(b"accepted-alignment"),
        "index": None,
        "index_sha256": None,
    }
    displaced_output = output_dir.with_name("job-1-displaced")
    child_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks, "report_pipeline_failure", MagicMock())

    def displace_then_fail(command, *, check: bool) -> None:
        assert check is True
        output_dir.rename(displaced_output)
        output_dir.mkdir()
        (output_dir / "sentinel.txt").write_text("replacement", encoding="utf-8")
        child_output = Path(command[command.index("-o") + 1])
        partial_dir = child_output / "partial-stage"
        partial_dir.mkdir()
        (partial_dir / "partial-result.txt").write_text("partial", encoding="utf-8")
        raise child_error

    monkeypatch.setattr(tasks.subprocess, "run", displace_then_fail)

    with pytest.raises(subprocess.CalledProcessError) as caught:
        _invoke_job(
            tasks,
            bam_path=str(alignment),
            output_dir=str(output_dir),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
            workspace_identity=identity,
        )

    assert caught.value is child_error
    assert (output_dir / "sentinel.txt").read_text(encoding="utf-8") == "replacement"
    assert list(displaced_output.iterdir()) == []


def test_bound_output_persists_a_real_curated_preflight_failure(tmp_path: Path) -> None:
    """The procfd output spelling remains writable by the real preflight transport."""
    from app.pipeline_job_workspace import open_pipeline_job_workspace

    from vntyper.scripts.preflight_error_io import PreflightErrorContext, persist_preflight_failure

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(b"alignment"),
        "index": None,
        "index_sha256": None,
    }
    payload = {
        "code": "alignment_preflight_failed",
        "message": "Alignment preflight failed before processing; inspect the server logs for the job.",
        "candidates": [],
    }
    workspace = open_pipeline_job_workspace(str(alignment), None, str(output_dir), identity)
    original = RuntimeError("preflight failed")
    try:
        with (
            pytest.raises(RuntimeError) as caught,
            persist_preflight_failure(PreflightErrorContext(workspace.bound_output_path, payload=payload)),
        ):
            raise original
        assert caught.value is original
        assert json.loads((output_dir / "preflight_error.json").read_text(encoding="utf-8")) == payload
    finally:
        workspace.remove_views()
        workspace.close()


@pytest.mark.parametrize("replacement_kind", ["regular", "symlink"])
def test_cleanup_preserves_unbound_conventional_index_replacements(
    replacement_kind: str,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """No index token means cleanup owns no conventional index pathname."""
    from app import tasks

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "index": None,
    }
    replacement = input_dir / "sample.bam.bai"
    external = tmp_path / "external-index.bai"
    external.write_bytes(b"external-index")

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))

    def install_replacement(_command, *, check: bool) -> None:
        assert check is True
        if replacement_kind == "symlink":
            replacement.symlink_to(external)
        else:
            replacement.write_bytes(b"replacement-index")

    monkeypatch.setattr(tasks.subprocess, "run", install_replacement)

    _invoke_job(
        tasks,
        bam_path=str(alignment),
        output_dir=str(output_dir),
        thread=1,
        reference_assembly="hg38",
        fast_mode=False,
        keep_intermediates=False,
        archive_results=False,
        workspace_identity=identity,
    )

    assert os.path.lexists(replacement)
    if replacement_kind == "symlink":
        assert replacement.is_symlink()
        assert external.read_bytes() == b"external-index"
    else:
        assert replacement.read_bytes() == b"replacement-index"


@pytest.mark.parametrize("identity_token", [[], "not-a-mapping"])
def test_non_mapping_workspace_tokens_fail_through_the_task_contract(
    identity_token,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Malformed broker payloads use the logged refusal and normal failed-job path."""
    from app import tasks

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"alignment")
    queue = MagicMock(name="redis_client")
    usage = MagicMock(name="redis_usage_client")
    monkeypatch.setattr(tasks, "redis_client", queue)
    monkeypatch.setattr(tasks, "redis_cohort_client", MagicMock(name="redis_cohort_client"))
    monkeypatch.setattr(tasks, "redis_usage_client", usage)
    pipeline = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", pipeline)
    caplog.set_level("ERROR")

    with pytest.raises(RuntimeError, match="invalid workspace identity token"):
        _invoke_job(
            tasks,
            bam_path=str(alignment),
            output_dir=str(output_dir),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
            workspace_identity=identity_token,
        )

    pipeline.assert_not_called()
    usage.hset.assert_any_call("usage:job-1", "status", "failed")
    queue.lrem.assert_called_once_with("vntyper_job_queue", 0, "task-inode")
    assert "Refusing pipeline job: invalid workspace identity token" in caplog.text


def test_output_displacement_during_completion_bookkeeping_fails_the_task(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A swap after the child check cannot turn detached results into success."""
    from app import tasks

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "index": None,
    }
    displaced_output = output_dir.with_name("job-1-displaced")
    usage = MagicMock(name="redis_usage_client")

    for name in ("redis_client", "redis_cohort_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks, "redis_usage_client", usage)

    def write_result(command, *, check: bool) -> None:
        assert check is True
        child_output = Path(command[command.index("-o") + 1])
        (child_output / "result.txt").write_text("owned", encoding="utf-8")

    def displace_on_completed(_key, *args, **_kwargs) -> None:
        if args == ("status", "completed"):
            output_dir.rename(displaced_output)
            output_dir.mkdir()
            (output_dir / "sentinel.txt").write_text("replacement", encoding="utf-8")

    monkeypatch.setattr(tasks.subprocess, "run", write_result)
    usage.hset.side_effect = displace_on_completed

    with pytest.raises(RuntimeError, match="output directory identity changed during completion"):
        _invoke_job(
            tasks,
            bam_path=str(alignment),
            output_dir=str(output_dir),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
            workspace_identity=identity,
        )

    assert (output_dir / "sentinel.txt").read_text(encoding="utf-8") == "replacement"
    assert list(displaced_output.iterdir()) == []
    usage.hset.assert_any_call("usage:job-1", "status", "failed")


@pytest.mark.parametrize(
    ("alignment_name", "index_name", "file_format"),
    [
        ("sample.bam", "sample.bam.bai", "bam"),
        ("sample.bam", "sample.bai", "bam"),
        ("sample.cram", "sample.cram.crai", "cram"),
        ("sample.cram", "sample.crai", "cram"),
    ],
)
def test_real_ownership_guard_and_index_resolver_accept_every_uploaded_index_name(
    alignment_name: str,
    index_name: str,
    file_format: str,
    tmp_path: Path,
) -> None:
    """The bound view satisfies both production ownership and resolver contracts."""
    from app.pipeline_job_workspace import open_pipeline_job_workspace

    from tests.support.pipeline_harness import MINIMAL_CONFIG
    from vntyper.scripts.alignment_index import resolve_any_index, resolve_bam_index
    from vntyper.scripts.alignment_target_io import protect_alignment_inputs

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / alignment_name
    index = input_dir / index_name
    alignment.write_bytes(b"alignment")
    index.write_bytes(b"index")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(b"alignment"),
        "index": _identity(index),
        "index_sha256": _digest(b"index"),
    }
    workspace = open_pipeline_job_workspace(str(alignment), str(index), str(output_dir), identity)
    try:
        protect_alignment_inputs(
            workspace.bound_output_path,
            workspace.bound_alignment_path,
            file_format,
            None,
            None,
            MINIMAL_CONFIG,
            "hg38",
        )
        resolved = (
            resolve_bam_index(workspace.bound_alignment_path)
            if file_format == "bam"
            else resolve_any_index(workspace.bound_alignment_path, file_format)
        )
        assert resolved is not None
        assert workspace.bound_index_path is not None
        assert os.path.samefile(resolved, workspace.bound_index_path)
    finally:
        workspace.remove_views()
        workspace.close()


def test_worker_holds_input_and_output_descriptors_through_execution_and_cleanup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Mid-run path swaps neither substitute input nor redirect output/cleanup."""
    from app import tasks

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"accepted-alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "index": None,
    }
    displaced_input = input_dir.with_name("job-1-displaced")
    displaced_output = output_dir.with_name("job-1-displaced")

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))

    def swap_during_pipeline(command, *, check: bool) -> None:
        assert check is True
        input_dir.rename(displaced_input)
        input_dir.mkdir()
        (input_dir / "sample.bam").write_bytes(b"substituted-alignment")
        output_dir.rename(displaced_output)
        output_dir.mkdir()
        (output_dir / "sentinel.txt").write_text("replacement", encoding="utf-8")

        flag = "--cram" if "--cram" in command else "--bam"
        child_alignment = Path(command[command.index(flag) + 1])
        child_output = Path(command[command.index("-o") + 1])
        assert child_alignment.read_bytes() == b"accepted-alignment"
        assert child_alignment.parent.parent == Path("/tmp")
        assert not child_alignment.is_symlink()
        assert str(child_output).startswith(f"/proc/{os.getpid()}/fd/")
        (child_output / "pipeline-result.txt").write_text("owned", encoding="utf-8")

    monkeypatch.setattr(tasks.subprocess, "run", swap_during_pipeline)

    with pytest.raises(RuntimeError, match="output directory identity changed"):
        _invoke_job(
            tasks,
            bam_path=str(alignment),
            output_dir=str(output_dir),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
            workspace_identity=identity,
        )

    assert (input_dir / "sample.bam").read_bytes() == b"substituted-alignment"
    assert (output_dir / "sentinel.txt").read_text(encoding="utf-8") == "replacement"
    assert not (displaced_input / "sample.bam").exists()
    assert list(displaced_output.iterdir()) == []


def test_descriptor_cleanup_does_not_mask_the_pipeline_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A secondary bound-workspace cleanup error preserves the child error."""
    from app import tasks

    input_dir = tmp_path / "input" / "job-1"
    output_dir = tmp_path / "output" / "job-1"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    alignment = input_dir / "sample.bam"
    alignment.write_bytes(b"accepted-alignment")
    identity = {
        "input_dir": _identity(input_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "index": None,
    }
    pipeline_error = subprocess.CalledProcessError(1, ["vntyper", "pipeline"])

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks, "report_pipeline_failure", MagicMock())
    monkeypatch.setattr(tasks.subprocess, "run", MagicMock(side_effect=pipeline_error))
    real_remove_views = tasks.PipelineJobWorkspace.remove_views

    def remove_views_then_report_failure(workspace) -> None:
        real_remove_views(workspace)
        raise OSError("cleanup failed")

    monkeypatch.setattr(tasks.PipelineJobWorkspace, "remove_views", remove_views_then_report_failure)

    with pytest.raises(subprocess.CalledProcessError) as caught:
        _invoke_job(
            tasks,
            bam_path=str(alignment),
            output_dir=str(output_dir),
            thread=1,
            reference_assembly="hg38",
            fast_mode=False,
            keep_intermediates=False,
            archive_results=False,
            workspace_identity=identity,
        )

    assert caught.value is pipeline_error
