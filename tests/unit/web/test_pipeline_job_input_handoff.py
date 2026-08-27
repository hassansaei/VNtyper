"""Input identity, digest, and worker-snapshot handoff tests."""

import stat
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tests.unit.web.pipeline_job_inode_handoff_support import _digest, _identity, _invoke_job

pytestmark = pytest.mark.unit


def test_api_queues_digest_bound_spool_paths_and_never_writes_shared_input(client, web_app) -> None:
    """Accepted research data exists only in the protected API/worker spool."""
    shared_input = Path(web_app.DEFAULT_INPUT_DIR)
    displaced_input = shared_input.with_name("input-displaced")
    shared_input.rename(displaced_input)
    shared_input.mkdir(exist_ok=True)
    payload = b"accepted-alignment"

    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", payload, "application/octet-stream")},
    )

    assert response.status_code == 200
    queued = web_app.run_vntyper_job.delay.call_args.kwargs
    queued_alignment = Path(queued["bam_path"])
    assert queued_alignment.parent.parent == Path(web_app.DEFAULT_HANDOFF_SPOOL_DIR)
    assert queued_alignment.read_bytes() == payload
    assert queued["workspace_identity"]["alignment_sha256"] == _digest(payload)
    assert list(shared_input.iterdir()) == []
    assert list(displaced_input.iterdir()) == []


def test_upload_receipt_carries_lowercase_sha256(tmp_path: Path) -> None:
    """The byte digest is captured by the same copy that creates the handoff file."""
    import io

    from app.uploads import UploadReceipt, save_upload_bounded

    payload = b"exact-upload-bytes"
    destination = tmp_path / "sample.bam"

    receipt = save_upload_bounded(io.BytesIO(payload), str(destination), 100, capture_identity=True)

    assert isinstance(receipt, UploadReceipt)
    assert receipt.sha256 == _digest(payload)
    assert receipt.sha256.islower()


def test_same_inode_spool_mutation_is_rejected_and_reclaimed(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A queued digest detects changed alignment bytes even when inode identity survives."""
    from app import tasks

    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / "job-1"
    output_root = tmp_path / "output"
    output_dir = output_root / "job-1"
    shared_input = tmp_path / "input"
    spool_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    shared_input.mkdir(exist_ok=True)
    accepted = b"accepted-alignment"
    alignment = spool_dir / "sample.bam"
    alignment.write_bytes(accepted)
    identity = {
        "input_dir": _identity(spool_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(accepted),
        "index": None,
        "index_sha256": None,
    }
    original_inode = alignment.stat().st_ino
    alignment.write_bytes(b"mutated-alignment")
    assert alignment.stat().st_ino == original_inode

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks.settings, "DEFAULT_HANDOFF_SPOOL_DIR", str(spool_root), raising=False)
    monkeypatch.setattr(tasks.settings, "DEFAULT_INPUT_DIR", str(shared_input))
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(output_root))
    pipeline = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", pipeline)

    with pytest.raises(RuntimeError, match="alignment SHA-256 changed before worker start"):
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

    pipeline.assert_not_called()
    assert not alignment.exists()
    assert not spool_dir.exists() or list(spool_dir.iterdir()) == []


def test_index_digest_mutation_is_rejected_before_the_child(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The uploaded alignment index has an independent byte-integrity token."""
    from app import tasks

    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / "job-1"
    output_root = tmp_path / "output"
    output_dir = output_root / "job-1"
    shared_input = tmp_path / "input"
    spool_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    shared_input.mkdir(exist_ok=True)
    alignment_bytes = b"accepted-alignment"
    index_bytes = b"accepted-index"
    alignment = spool_dir / "sample.bam"
    index = spool_dir / "sample.bam.bai"
    alignment.write_bytes(alignment_bytes)
    index.write_bytes(index_bytes)
    identity = {
        "input_dir": _identity(spool_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(alignment_bytes),
        "index": _identity(index),
        "index_sha256": _digest(index_bytes),
    }
    index.write_bytes(b"mutated-index")

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks.settings, "DEFAULT_HANDOFF_SPOOL_DIR", str(spool_root), raising=False)
    monkeypatch.setattr(tasks.settings, "DEFAULT_INPUT_DIR", str(shared_input))
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(output_root))
    pipeline = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", pipeline)

    with pytest.raises(RuntimeError, match="index SHA-256 changed before worker start"):
        _invoke_job(
            tasks,
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

    pipeline.assert_not_called()
    assert not alignment.exists()
    assert not index.exists()


def test_worker_uses_real_single_link_tmp_snapshots_despite_hostile_tmpdir(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The child receives immutable-by-scope real files outside every shared mount."""
    from app import tasks

    spool_root = tmp_path / "handoff"
    spool_dir = spool_root / "job-1"
    output_root = tmp_path / "output"
    output_dir = output_root / "job-1"
    shared_input = tmp_path / "input"
    spool_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    shared_input.mkdir()
    alignment_bytes = b"accepted-alignment"
    index_bytes = b"accepted-index"
    alignment = spool_dir / "sample.bam"
    index = spool_dir / "sample.bam.bai"
    alignment.write_bytes(alignment_bytes)
    index.write_bytes(index_bytes)
    identity = {
        "input_dir": _identity(spool_dir),
        "output_dir": _identity(output_dir),
        "alignment": _identity(alignment),
        "alignment_sha256": _digest(alignment_bytes),
        "index": _identity(index),
        "index_sha256": _digest(index_bytes),
    }

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    monkeypatch.setattr(tasks.settings, "DEFAULT_HANDOFF_SPOOL_DIR", str(spool_root), raising=False)
    monkeypatch.setattr(tasks.settings, "DEFAULT_INPUT_DIR", str(shared_input))
    monkeypatch.setattr(tasks.settings, "DEFAULT_OUTPUT_DIR", str(output_root))
    monkeypatch.setenv("TMPDIR", str(shared_input))

    def inspect_snapshot(command, *, check: bool) -> None:
        assert check is True
        child_alignment = Path(command[command.index("--bam") + 1])
        snapshot_dir = child_alignment.parent
        assert snapshot_dir.parent == Path("/tmp")
        assert child_alignment.read_bytes() == alignment_bytes
        assert child_alignment.stat().st_nlink == 1
        assert stat.S_ISREG(child_alignment.stat().st_mode)
        entries = list(snapshot_dir.iterdir())
        assert entries
        assert all(not entry.is_symlink() for entry in entries)
        resolved_index = child_alignment.with_suffix(child_alignment.suffix + ".bai")
        assert resolved_index.read_bytes() == index_bytes
        assert resolved_index.stat().st_nlink == 1

    monkeypatch.setattr(tasks.subprocess, "run", inspect_snapshot)

    _invoke_job(
        tasks,
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


def test_upload_stays_in_the_created_inode_but_a_replaced_name_refuses_handoff(
    client,
    web_app,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A replacement cannot capture bytes or be accepted as the worker handoff."""
    real_save = web_app.save_upload_bounded
    created_input: Path | None = None
    replacement_input: Path | None = None

    def swap_then_save(source, destination: str, max_bytes: int, *args, **kwargs):
        nonlocal created_input, replacement_input
        public_input = Path(destination).parent
        created_input = public_input.with_name(f"{public_input.name}-created")
        public_input.rename(created_input)
        public_input.mkdir()
        replacement_input = public_input
        return real_save(source, destination, max_bytes, *args, **kwargs)

    monkeypatch.setattr(web_app, "save_upload_bounded", swap_then_save)

    with pytest.raises(RuntimeError, match="input directory identity changed before handoff"):
        client.post(
            "/run-job/",
            files={"bam_file": ("sample.bam", b"created-inode-bytes", "application/octet-stream")},
        )

    assert created_input is not None
    assert replacement_input is not None
    assert list(created_input.iterdir()) == []
    assert list(replacement_input.iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()
    web_app.run_vntyper_job.apply_async.assert_not_called()


def test_worker_rejects_a_replaced_input_directory_before_starting_the_pipeline(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The task cannot adopt a different real directory under the queued path."""
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

    displaced = input_dir.with_name("job-1-displaced")
    input_dir.rename(displaced)
    input_dir.mkdir()
    replacement = input_dir / "sample.bam"
    replacement.write_bytes(b"substituted-alignment")

    for name in ("redis_client", "redis_cohort_client", "redis_usage_client"):
        monkeypatch.setattr(tasks, name, MagicMock(name=name))
    pipeline = MagicMock(name="subprocess.run")
    monkeypatch.setattr(tasks.subprocess, "run", pipeline)

    with pytest.raises(RuntimeError, match="input directory identity changed"):
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

    pipeline.assert_not_called()
    assert replacement.read_bytes() == b"substituted-alignment"
    assert (displaced / "sample.bam").read_bytes() == b"accepted-alignment"


def test_worker_reports_a_symlinked_input_directory_as_a_closed_refusal(tmp_path: Path) -> None:
    """An input symlink is rejected with the worker's stable refusal contract."""
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
    displaced = input_dir.with_name("job-1-displaced")
    input_dir.rename(displaced)
    input_dir.symlink_to(displaced, target_is_directory=True)

    with pytest.raises(RuntimeError, match="workspace could not be safely bound"):
        open_pipeline_job_workspace(str(alignment), None, str(output_dir), identity)
