"""Identity-token validation and descriptor cleanup handoff tests."""

import os
from contextlib import suppress
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tests.unit.web.pipeline_job_inode_handoff_support import _digest, _identity

pytestmark = pytest.mark.unit


def test_worker_rejects_a_malformed_identity_token_before_opening_paths(tmp_path: Path) -> None:
    """Broker payloads cannot replace the expected two-integer identity shape."""
    from app.pipeline_job_workspace import open_pipeline_job_workspace

    with pytest.raises(RuntimeError, match="invalid input directory identity token"):
        open_pipeline_job_workspace(
            str(tmp_path / "input" / "sample.bam"),
            None,
            str(tmp_path / "output"),
            {"input_dir": "not-an-identity", "output_dir": [1, 2], "alignment": [1, 3], "index": None},
        )


@pytest.mark.parametrize("identity_token", [[], "not-a-mapping", {"input_dir": "malformed"}])
def test_early_invalid_workspace_tokens_reclaim_only_the_matching_uuid_spool_job(
    identity_token,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """An early token refusal cannot strand accepted bytes in the protected spool."""
    from app import tasks

    job_id = "12345678-1234-4234-8234-123456789abc"
    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / job_id
    output_root = tmp_path / "output"
    output_dir = output_root / job_id
    shared_input = tmp_path / "input"
    spool_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    shared_input.mkdir()
    alignment = spool_dir / "sample.bam"
    index = spool_dir / "sample.bam.bai"
    alignment.write_bytes(b"accepted-alignment")
    index.write_bytes(b"accepted-index")

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks.settings, "DEFAULT_HANDOFF_SPOOL_DIR", str(spool_root), raising=False)
    monkeypatch.setattr(tasks.settings, "DEFAULT_INPUT_DIR", str(shared_input))
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(output_root))
    pipeline = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", pipeline)

    tasks.run_vntyper_job.push_request(id="task-invalid-token")
    try:
        with pytest.raises(RuntimeError, match="invalid .*identity token"):
            tasks.run_vntyper_job.run(
                bam_path=str(alignment),
                index_path=str(index),
                output_dir=str(output_dir),
                thread=1,
                reference_assembly="hg38",
                fast_mode=False,
                keep_intermediates=False,
                archive_results=False,
                workspace_identity=identity_token,
            )
    finally:
        tasks.run_vntyper_job.pop_request()

    pipeline.assert_not_called()
    assert not alignment.exists()
    assert not index.exists()
    assert not spool_dir.exists()


def test_missing_index_tokens_refuse_and_reclaim_the_supplied_uuid_spool_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A supplied index cannot be omitted from snapshot and cleanup ownership."""
    from app import tasks

    job_id = "12345678-1234-4234-8234-123456789abc"
    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / job_id
    output_root = tmp_path / "output"
    output_dir = output_root / job_id
    shared_input = tmp_path / "input"
    spool_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    shared_input.mkdir()
    alignment = spool_dir / "sample.bam"
    index = spool_dir / "sample.bam.bai"
    alignment.write_bytes(b"accepted-alignment")
    index.write_bytes(b"accepted-index")
    identity = {
        "input_dir": _identity(spool_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(b"accepted-alignment"),
        "index": None,
        "index_sha256": None,
    }

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks.settings, "DEFAULT_HANDOFF_SPOOL_DIR", str(spool_root), raising=False)
    monkeypatch.setattr(tasks.settings, "DEFAULT_INPUT_DIR", str(shared_input))
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(output_root))
    pipeline = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", pipeline)

    tasks.run_vntyper_job.push_request(id="task-index-token")
    try:
        with pytest.raises(RuntimeError, match="index path and identity token disagree"):
            tasks.run_vntyper_job.run(
                bam_path=str(alignment),
                index_path=str(index),
                output_dir=str(output_dir),
                thread=1,
                reference_assembly="hg38",
                fast_mode=False,
                keep_intermediates=False,
                archive_results=False,
                workspace_identity=identity,
            )
    finally:
        tasks.run_vntyper_job.pop_request()

    pipeline.assert_not_called()
    assert not alignment.exists()
    assert not index.exists()
    assert not spool_dir.exists()


def test_invalid_token_cleanup_refuses_paths_outside_or_for_a_different_spool_job(
    tmp_path: Path,
) -> None:
    """The fallback capability cannot delete another job or an arbitrary path."""
    from app.pipeline_job_workspace import reclaim_unopened_spool_inputs

    spool_root = tmp_path / "handoff"
    first_job = spool_root / "12345678-1234-4234-8234-123456789abc"
    second_job = spool_root / "87654321-4321-4321-8321-cba987654321"
    outside = tmp_path / "outside"
    for directory in (first_job, second_job, outside):
        directory.mkdir(parents=True)
    protected = second_job / "sample.bam"
    arbitrary = outside / "sample.bam"
    protected.write_bytes(b"other-job")
    arbitrary.write_bytes(b"outside")

    with pytest.raises(RuntimeError, match="does not match the output job"):
        reclaim_unopened_spool_inputs(
            str(protected),
            None,
            str(tmp_path / "output" / first_job.name),
            str(spool_root),
        )
    with pytest.raises(RuntimeError, match="outside the protected handoff spool"):
        reclaim_unopened_spool_inputs(
            str(arbitrary),
            None,
            str(tmp_path / "output" / first_job.name),
            str(spool_root),
        )

    assert protected.read_bytes() == b"other-job"
    assert arbitrary.read_bytes() == b"outside"


def test_fallback_cleanup_closes_all_descriptors_without_masking_an_active_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A cleanup failure remains primary while both retained handles are closed."""
    from app import pipeline_job_workspace

    job_id = "12345678-1234-4234-8234-123456789abc"
    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / job_id
    spool_dir.mkdir(parents=True)
    alignment = spool_dir / "sample.bam"
    alignment.write_bytes(b"accepted")
    real_open = pipeline_job_workspace.os.open
    real_close = pipeline_job_workspace.os.close
    opened: list[int] = []
    closed: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_stat(*_args, **_kwargs):
        raise OSError("primary cleanup failure")

    def failing_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)
        raise OSError(f"close failure {descriptor}")

    monkeypatch.setattr(pipeline_job_workspace.os, "open", recording_open)
    monkeypatch.setattr(pipeline_job_workspace.os, "stat", failing_stat)
    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)
    try:
        with pytest.raises(OSError, match="primary cleanup failure"):
            pipeline_job_workspace.reclaim_unopened_spool_inputs(
                str(alignment),
                None,
                str(tmp_path / "output" / job_id),
                str(spool_root),
            )
    finally:
        for descriptor in opened:
            with suppress(OSError):
                real_close(descriptor)

    assert closed == list(reversed(opened))


def test_fallback_cleanup_attempts_all_closes_and_surfaces_the_first_close_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Successful cleanup reports its first close error only after closing both handles."""
    from app import pipeline_job_workspace

    job_id = "12345678-1234-4234-8234-123456789abc"
    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / job_id
    spool_dir.mkdir(parents=True)
    alignment = spool_dir / "sample.bam"
    alignment.write_bytes(b"accepted")
    real_open = pipeline_job_workspace.os.open
    real_close = pipeline_job_workspace.os.close
    opened: list[int] = []
    closed: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)
        label = "job" if descriptor == opened[-1] else "root"
        raise OSError(f"{label} close failure")

    monkeypatch.setattr(pipeline_job_workspace.os, "open", recording_open)
    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)

    with pytest.raises(OSError, match="job close failure"):
        pipeline_job_workspace.reclaim_unopened_spool_inputs(
            str(alignment),
            None,
            str(tmp_path / "output" / job_id),
            str(spool_root),
        )

    assert closed == list(reversed(opened))


def test_worker_rejects_a_replaced_alignment_file_before_creating_private_views(tmp_path: Path) -> None:
    """A file-level replacement is refused even when both directory names remain stable."""
    from app.pipeline_job_workspace import open_pipeline_job_workspace

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
    original = input_dir / "original.bam"
    alignment.rename(original)
    alignment.write_bytes(b"replacement-alignment")

    with pytest.raises(RuntimeError, match="alignment identity changed before worker start"):
        open_pipeline_job_workspace(str(alignment), None, str(output_dir), identity)

    assert original.read_bytes() == b"accepted-alignment"
    assert alignment.read_bytes() == b"replacement-alignment"


def test_open_directory_closes_after_fstat_failure_without_masking_it(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A directory fstat failure owns the error while its descriptor is still closed."""
    from app import pipeline_job_workspace

    directory = tmp_path / "directory"
    directory.mkdir()
    directory_metadata = directory.stat()
    directory_identity = (directory_metadata.st_dev, directory_metadata.st_ino)
    real_open = pipeline_job_workspace.os.open
    real_close = pipeline_job_workspace.os.close
    opened: list[int] = []
    closes: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_fstat(_descriptor: int):
        raise OSError("directory fstat failed")

    def failing_close(descriptor: int) -> None:
        closes.append(descriptor)
        real_close(descriptor)
        raise OSError("directory close failed")

    monkeypatch.setattr(pipeline_job_workspace.os, "open", recording_open)
    monkeypatch.setattr(pipeline_job_workspace.os, "fstat", failing_fstat)
    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)
    try:
        with pytest.raises(OSError, match="directory fstat failed"):
            pipeline_job_workspace._open_directory(str(directory), directory_identity, "input")
    finally:
        for descriptor in opened:
            with suppress(OSError):
                real_close(descriptor)

    assert closes == opened


def test_open_file_closes_after_fstat_failure_without_masking_it(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A file fstat failure owns the error while its descriptor is still closed."""
    from app import pipeline_job_workspace

    directory = tmp_path / "directory"
    directory.mkdir()
    source = directory / "sample.bam"
    source.write_bytes(b"alignment")
    source_metadata = source.stat()
    source_identity = (source_metadata.st_dev, source_metadata.st_ino)
    parent_descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    real_open = pipeline_job_workspace.os.open
    real_close = pipeline_job_workspace.os.close
    opened: list[int] = []
    closes: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_fstat(_descriptor: int):
        raise OSError("file fstat failed")

    def failing_close(descriptor: int) -> None:
        closes.append(descriptor)
        real_close(descriptor)
        raise OSError("file close failed")

    monkeypatch.setattr(pipeline_job_workspace.os, "open", recording_open)
    monkeypatch.setattr(pipeline_job_workspace.os, "fstat", failing_fstat)
    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)
    try:
        with pytest.raises(OSError, match="file fstat failed"):
            pipeline_job_workspace._open_file(parent_descriptor, str(source), source_identity, "alignment")
    finally:
        real_close(parent_descriptor)
        for descriptor in opened:
            with suppress(OSError):
                real_close(descriptor)

    assert closes == opened


def test_open_file_identity_refusal_survives_a_secondary_close_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Closing a substituted file cannot replace the stable identity refusal."""
    from app import pipeline_job_workspace

    directory = tmp_path / "directory"
    directory.mkdir()
    source = directory / "sample.bam"
    source.write_bytes(b"alignment")
    parent_descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    real_close = pipeline_job_workspace.os.close

    def failing_close(descriptor: int) -> None:
        real_close(descriptor)
        raise OSError("file close failed")

    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)
    try:
        with pytest.raises(RuntimeError, match="alignment identity changed before worker start"):
            pipeline_job_workspace._open_file(parent_descriptor, str(source), (0, 0), "alignment")
    finally:
        real_close(parent_descriptor)


def test_snapshot_digest_refusal_survives_a_secondary_close_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A destination close error cannot mask the byte-integrity refusal."""
    from app import pipeline_job_workspace

    source = tmp_path / "source.bam"
    destination = tmp_path / "snapshot"
    source.write_bytes(b"alignment")
    destination.mkdir()
    source_descriptor = os.open(source, os.O_RDONLY)
    destination_descriptor = os.open(destination, os.O_RDONLY | os.O_DIRECTORY)
    real_open = pipeline_job_workspace.os.open
    real_close = pipeline_job_workspace.os.close
    snapshot_descriptors: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        snapshot_descriptors.append(descriptor)
        return descriptor

    def failing_close(descriptor: int) -> None:
        real_close(descriptor)
        if descriptor in snapshot_descriptors:
            raise OSError("snapshot close failed")

    monkeypatch.setattr(pipeline_job_workspace.os, "open", recording_open)
    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)
    try:
        with pytest.raises(RuntimeError, match="alignment SHA-256 changed before worker start"):
            pipeline_job_workspace._copy_snapshot(
                source_descriptor,
                destination_descriptor,
                "sample.bam",
                _digest(b"different"),
                "alignment",
            )
    finally:
        real_close(source_descriptor)
        real_close(destination_descriptor)


def test_pipeline_recursive_cleanup_preserves_primary_while_closing_ancestors(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Nested output cleanup closes cannot replace the recursive primary error."""
    from app import pipeline_job_workspace

    root = tmp_path / "root"
    trigger = root / "nested" / "deeper" / "trigger"
    trigger.parent.mkdir(parents=True)
    trigger.write_bytes(b"data")
    real_open = pipeline_job_workspace.os.open
    real_stat = pipeline_job_workspace.os.stat
    real_close = pipeline_job_workspace.os.close
    root_descriptor = real_open(root, pipeline_job_workspace.os.O_RDONLY | pipeline_job_workspace.os.O_DIRECTORY)
    opened: list[int] = []
    closed: list[int] = []

    def recording_open(*args, **kwargs) -> int:
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_stat(path, *args, **kwargs):
        if path == "trigger":
            raise OSError("primary recursive cleanup failure")
        return real_stat(path, *args, **kwargs)

    def failing_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)
        raise OSError(f"secondary close failure {descriptor}")

    monkeypatch.setattr(pipeline_job_workspace.os, "open", recording_open)
    monkeypatch.setattr(pipeline_job_workspace.os, "stat", failing_stat)
    monkeypatch.setattr(pipeline_job_workspace.os, "close", failing_close)
    try:
        with pytest.raises(OSError, match="primary recursive cleanup failure"):
            pipeline_job_workspace._clear_directory(root_descriptor)
    finally:
        real_close(root_descriptor)
        for descriptor in opened:
            with suppress(OSError):
                real_close(descriptor)

    assert closed == list(reversed(opened))


def test_workspace_close_attempts_every_descriptor_after_one_close_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    """One descriptor-close error cannot leak the remaining retained handles."""
    from app import pipeline_job_workspace

    workspace = pipeline_job_workspace.PipelineJobWorkspace(
        input_descriptor=11,
        output_descriptor=12,
        view_descriptor=15,
        alignment_descriptor=13,
        index_descriptor=14,
        alignment_path="/input/sample.bam",
        index_path="/input/sample.bam.bai",
        output_dir="/output/job-1",
        view_dir="/private/view",
        input_identity=(1, 1),
        output_identity=(1, 2),
        view_identity=(1, 5),
        alignment_identity=(1, 3),
        index_identity=(1, 4),
        alignment_view_name="sample.bam",
        index_view_name="sample.bam.bai",
    )
    closed: list[int] = []

    def fail_first_close(descriptor: int) -> None:
        closed.append(descriptor)
        if descriptor == 14:
            raise OSError("first close failed")

    monkeypatch.setattr(pipeline_job_workspace.os, "close", fail_first_close)

    with pytest.raises(OSError, match="first close failed"):
        workspace.close()

    assert closed == [14, 13, 15, 12, 11]
