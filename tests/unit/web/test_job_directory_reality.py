"""A job's directories must be real directories this service created.

The containment check in ``app.uploads.safe_upload_path`` compares two
realpaths. When the per-job directory is a symlink, both sides resolve to the
link's target, so the workspace must refuse planted entries and reclaim a
post-creation swap without following it.

The reproduced escape requires write access to the shared input volume plus a
fresh-UUID guess, or a race between directory creation and the upload open.
Nothing in the HTTP surface creates symlinks; these checks defend the service's
filesystem-containment promise rather than describing a remote exploit.
"""

import io
import os
import stat
from pathlib import Path
from uuid import UUID

import pytest

pytestmark = pytest.mark.unit

from app import job_workspace as job_workspace_module  # noqa: E402
from app import uploads as uploads_module  # noqa: E402
from app.job_workspace import job_workspace  # noqa: E402
from app.uploads import save_upload_bounded  # noqa: E402


def test_workspace_refuses_a_symlinked_job_input_directory(tmp_path: Path) -> None:
    """A pre-planted symlink under the input root is refused, not adopted."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    input_root = tmp_path / "input"
    input_root.mkdir()
    (input_root / "job-1").symlink_to(elsewhere)

    with pytest.raises(RuntimeError, match="job-1"):  # noqa: SIM117
        with job_workspace(str(input_root), str(tmp_path / "output"), "job-1"):
            raise AssertionError("the workspace must never yield a symlinked directory")

    assert list(elsewhere.iterdir()) == []
    assert not (tmp_path / "output" / "job-1").exists()


def test_workspace_refuses_a_pre_existing_job_directory(tmp_path: Path) -> None:
    """A fresh per-request id must never find its input name already present."""
    input_root = tmp_path / "input"
    input_root.mkdir()
    (input_root / "job-1").mkdir()

    with pytest.raises(RuntimeError, match="job-1"):  # noqa: SIM117
        with job_workspace(str(input_root), str(tmp_path / "output"), "job-1"):
            raise AssertionError("an already-existing entry is evidence of interference")


def test_workspace_reclaims_input_when_output_name_is_already_taken(tmp_path: Path) -> None:
    """An output collision refuses the job and removes its newly made input."""
    output_root = tmp_path / "output"
    output_root.mkdir()
    (output_root / "job-1").mkdir()

    with pytest.raises(RuntimeError, match="job-1"):  # noqa: SIM117
        with job_workspace(str(tmp_path / "input"), str(output_root), "job-1"):
            raise AssertionError("an occupied output name must prevent admission")

    assert not (tmp_path / "input" / "job-1").exists()
    assert (output_root / "job-1").is_dir()


def test_reclaim_leaves_a_swapped_symlink_and_does_not_touch_its_target(tmp_path: Path) -> None:
    """A replacement link is neither removed nor followed during bound cleanup."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (elsewhere / "keep.txt").write_text("keep", encoding="utf-8")
    input_root = tmp_path / "input"

    with pytest.raises(RuntimeError, match="refused"):  # noqa: SIM117
        with job_workspace(str(input_root), str(tmp_path / "output"), "job-1") as (job_input, _):
            os.rmdir(job_input)
            Path(job_input).symlink_to(elsewhere)
            msg = "the submission was refused"
            raise RuntimeError(msg)

    entry = input_root / "job-1"
    assert entry.is_symlink(), "a replacement public entry is not owned by the refused request"
    assert (elsewhere / "keep.txt").read_text(encoding="utf-8") == "keep"


def test_reclaim_failure_never_masks_the_submission_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A failed symlink unlink neither follows it nor masks the request error."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    input_root = tmp_path / "input"

    with pytest.raises(RuntimeError, match="original refusal"):  # noqa: SIM117
        with job_workspace(str(input_root), str(tmp_path / "output"), "job-1") as (job_input, _):
            os.rmdir(job_input)
            Path(job_input).symlink_to(elsewhere)

            def _fail_unlink(_path: str) -> None:
                raise PermissionError("cannot unlink")

            real_rmtree = job_workspace_module.shutil.rmtree

            def _reject_link_rmtree(path: str, **kwargs) -> None:
                if path == job_input:
                    raise AssertionError("a symlink must never be handed to rmtree")
                real_rmtree(path, **kwargs)

            monkeypatch.setattr(job_workspace_module.os, "unlink", _fail_unlink)
            monkeypatch.setattr(job_workspace_module.shutil, "rmtree", _reject_link_rmtree)
            msg = "original refusal"
            raise RuntimeError(msg)

    assert (input_root / "job-1").is_symlink()
    assert list(elsewhere.iterdir()) == []


def test_endpoint_never_stores_an_upload_through_a_symlinked_job_directory(
    web_app, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A pinned request is refused before it can write through a planted link."""
    from fastapi.testclient import TestClient

    pinned = UUID("b7c41d90-2e58-4a63-9f07-15c8de3b6a24")
    monkeypatch.setattr(web_app, "uuid4", lambda: pinned)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (tmp_path / "handoff" / str(pinned)).symlink_to(elsewhere)

    client = TestClient(web_app.app, raise_server_exceptions=False)
    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", b"bamdata", "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 500
    assert list(elsewhere.iterdir()) == []
    web_app.run_vntyper_job.delay.assert_not_called()
    web_app.run_vntyper_job.apply_async.assert_not_called()


# ---------------------------------------------------------------------------
# save_upload_bounded: the write itself refuses an unreal parent or entry.
# ---------------------------------------------------------------------------


def test_save_refuses_a_symlinked_parent_directory(tmp_path: Path) -> None:
    """Bytes are never written through a job directory that is really a link."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    input_root = tmp_path / "input"
    input_root.mkdir()
    job_dir = input_root / "job-1"
    job_dir.symlink_to(elsewhere)

    with pytest.raises(RuntimeError, match="job directory"):
        save_upload_bounded(io.BytesIO(b"x"), str(job_dir / "sample.bam"), 10)

    assert list(elsewhere.iterdir()) == []


def test_save_refuses_a_symlinked_destination_entry(tmp_path: Path) -> None:
    """An existing destination link is refused without changing its target."""
    victim = tmp_path / "victim.bam"
    victim.write_bytes(b"original")
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    (job_dir / "sample.bam").symlink_to(victim)

    with pytest.raises(RuntimeError, match="new file"):
        save_upload_bounded(io.BytesIO(b"x"), str(job_dir / "sample.bam"), 10)

    assert victim.read_bytes() == b"original"
    assert (job_dir / "sample.bam").is_symlink()


def test_save_refuses_an_existing_regular_destination(tmp_path: Path) -> None:
    """An existing regular file is never silently overwritten."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    destination.write_bytes(b"original")

    with pytest.raises(RuntimeError, match="new file"):
        save_upload_bounded(io.BytesIO(b"replacement"), str(destination), 20)

    assert destination.read_bytes() == b"original"


def test_save_still_writes_exact_bytes_with_a_readable_non_executable_mode(tmp_path: Path) -> None:
    """A fresh destination preserves the ordinary upload bytes and file mode."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"

    written = save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert written == 7
    assert destination.read_bytes() == b"payload"
    assert stat.S_IMODE(destination.stat().st_mode) == 0o644


def test_save_accepts_a_symlink_above_the_final_job_directory(tmp_path: Path) -> None:
    """A mounted or symlinked root remains valid when the job entry is real."""
    real_input_root = tmp_path / "real-input"
    real_input_root.mkdir()
    mounted_input_root = tmp_path / "mounted-input"
    mounted_input_root.symlink_to(real_input_root)
    job_dir = mounted_input_root / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"

    written = save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert written == 7
    assert (real_input_root / "job-1" / "sample.bam").read_bytes() == b"payload"


def test_a_refused_oversize_copy_reclaims_through_the_directory_descriptor(tmp_path: Path) -> None:
    """The size refusal leaves no partial destination in the real parent."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()

    with pytest.raises(ValueError, match="maximum accepted size"):
        save_upload_bounded(io.BytesIO(b"too many bytes"), str(job_dir / "sample.bam"), 4)

    assert list(job_dir.iterdir()) == []


def test_a_write_failure_reclaims_the_partial_destination(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A disk write error closes and removes the file created for the upload."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_fdopen = os.fdopen
    writer_closed = False

    class PartlyFailingWriter:
        """Real descriptor owner that fails after one byte reaches the file."""

        def __init__(self, descriptor: int) -> None:
            self._handle = real_fdopen(descriptor, "wb")

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            self.close()

        def close(self) -> None:
            nonlocal writer_closed
            self._handle.close()
            writer_closed = True

        def write(self, chunk: bytes) -> None:
            self._handle.write(chunk[:1])
            raise OSError("simulated full volume")

    monkeypatch.setattr(uploads_module.os, "fdopen", lambda descriptor, _mode: PartlyFailingWriter(descriptor))

    with pytest.raises(OSError, match="simulated full volume"):
        save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert writer_closed
    assert not destination.exists()


def test_fdopen_failure_closes_both_descriptors_and_reclaims_the_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Failure to wrap the created file cannot leak either raw descriptor."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_open = os.open
    opened_descriptors: list[int] = []

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        opened_descriptors.append(descriptor)
        return descriptor

    def refuse_fdopen(_descriptor: int, _mode: str):
        raise OSError("cannot construct file object")

    monkeypatch.setattr(uploads_module.os, "open", recording_open)
    monkeypatch.setattr(uploads_module.os, "fdopen", refuse_fdopen)

    with pytest.raises(OSError, match="cannot construct file object"):
        save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert len(opened_descriptors) == 2
    for descriptor in opened_descriptors:
        with pytest.raises(OSError):
            os.fstat(descriptor)
    assert not destination.exists()


def test_file_close_failure_does_not_prevent_parent_descriptor_cleanup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failing file close still reclaims the file and closes its parent."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_open = os.open
    real_fdopen = os.fdopen
    opened_descriptors: list[int] = []

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        opened_descriptors.append(descriptor)
        return descriptor

    class CloseFailingWriter:
        """File wrapper that closes its descriptor, then reports a close error."""

        def __init__(self, descriptor: int) -> None:
            self._handle = real_fdopen(descriptor, "wb")

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            self.close()

        def close(self) -> None:
            self._handle.close()
            raise OSError("simulated close failure")

        def write(self, chunk: bytes) -> None:
            self._handle.write(chunk)

    monkeypatch.setattr(uploads_module.os, "open", recording_open)
    monkeypatch.setattr(uploads_module.os, "fdopen", lambda descriptor, _mode: CloseFailingWriter(descriptor))

    with pytest.raises(OSError, match="simulated close failure"):
        save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert len(opened_descriptors) == 2
    for descriptor in opened_descriptors:
        with pytest.raises(OSError):
            os.fstat(descriptor)
    assert not destination.exists()


def test_file_close_failure_does_not_replace_an_active_read_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed close cannot hide the read failure that caused cleanup."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_fdopen = os.fdopen

    class ReadFailingSource:
        def read(self, _size: int) -> bytes:
            raise OSError("primary read failure")

    class CloseFailingWriter:
        def __init__(self, descriptor: int) -> None:
            self._handle = real_fdopen(descriptor, "wb")

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            self.close()

        def close(self) -> None:
            self._handle.close()
            raise OSError("secondary file close failure")

        def write(self, chunk: bytes) -> None:
            self._handle.write(chunk)

    monkeypatch.setattr(uploads_module.os, "fdopen", lambda descriptor, _mode: CloseFailingWriter(descriptor))

    with pytest.raises(OSError, match="primary read failure"):
        save_upload_bounded(ReadFailingSource(), str(destination), 10)

    assert not destination.exists()


def test_file_close_failure_does_not_replace_an_active_write_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed close cannot hide the write failure that caused cleanup."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_fdopen = os.fdopen

    class WriteAndCloseFailingWriter:
        def __init__(self, descriptor: int) -> None:
            self._handle = real_fdopen(descriptor, "wb")

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            self.close()

        def close(self) -> None:
            self._handle.close()
            raise OSError("secondary file close failure")

        def write(self, chunk: bytes) -> None:
            self._handle.write(chunk[:1])
            raise OSError("primary write failure")

    monkeypatch.setattr(uploads_module.os, "fdopen", lambda descriptor, _mode: WriteAndCloseFailingWriter(descriptor))

    with pytest.raises(OSError, match="primary write failure"):
        save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert not destination.exists()


def test_file_close_failure_does_not_replace_an_active_size_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed close cannot hide the bounded-copy size refusal."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_fdopen = os.fdopen

    class CloseFailingWriter:
        def __init__(self, descriptor: int) -> None:
            self._handle = real_fdopen(descriptor, "wb")

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            self.close()

        def close(self) -> None:
            self._handle.close()
            raise OSError("secondary file close failure")

        def write(self, chunk: bytes) -> None:
            self._handle.write(chunk)

    monkeypatch.setattr(uploads_module.os, "fdopen", lambda descriptor, _mode: CloseFailingWriter(descriptor))

    with pytest.raises(ValueError, match="maximum accepted size"):
        save_upload_bounded(io.BytesIO(b"oversize"), str(destination), 2)

    assert not destination.exists()


def test_parent_close_failure_does_not_replace_an_active_destination_open_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed parent close cannot hide refusal of an occupied destination."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    destination.write_bytes(b"original")
    real_open = os.open
    real_close = os.close
    parent_descriptor: int | None = None

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal parent_descriptor
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if dir_fd is None:
            parent_descriptor = descriptor
        return descriptor

    def close_parent_with_error(descriptor: int) -> None:
        real_close(descriptor)
        if descriptor == parent_descriptor:
            raise OSError("secondary parent close failure")

    monkeypatch.setattr(uploads_module.os, "open", recording_open)
    monkeypatch.setattr(uploads_module.os, "close", close_parent_with_error)

    with pytest.raises(RuntimeError, match="destination could not be created"):
        save_upload_bounded(io.BytesIO(b"replacement"), str(destination), 20)

    assert destination.read_bytes() == b"original"


def test_parent_close_failure_surfaces_when_the_copy_itself_succeeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A parent close error is reported when no earlier failure is active."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    destination = job_dir / "sample.bam"
    real_open = os.open
    real_close = os.close
    parent_descriptor: int | None = None

    def recording_open(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal parent_descriptor
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if dir_fd is None:
            parent_descriptor = descriptor
        return descriptor

    def close_parent_with_error(descriptor: int) -> None:
        real_close(descriptor)
        if descriptor == parent_descriptor:
            raise OSError("parent close failure")

    monkeypatch.setattr(uploads_module.os, "open", recording_open)
    monkeypatch.setattr(uploads_module.os, "close", close_parent_with_error)

    with pytest.raises(OSError, match="parent close failure"):
        save_upload_bounded(io.BytesIO(b"payload"), str(destination), 10)

    assert destination.read_bytes() == b"payload"


def test_post_open_parent_swap_cannot_redirect_creation_or_partial_cleanup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Creation and cleanup stay on the opened directory after a path swap."""
    job_dir = tmp_path / "job-1"
    job_dir.mkdir()
    original_directory = tmp_path / "original-job-1"
    swapped_target = tmp_path / "swapped-target"
    swapped_target.mkdir()
    real_open = os.open
    parent_was_swapped = False

    def swap_after_parent_open(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal parent_was_swapped
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if path == str(job_dir) and dir_fd is None:
            job_dir.rename(original_directory)
            job_dir.symlink_to(swapped_target)
            parent_was_swapped = True
        return descriptor

    monkeypatch.setattr(uploads_module.os, "open", swap_after_parent_open)

    with pytest.raises(ValueError, match="maximum accepted size"):
        save_upload_bounded(io.BytesIO(b"oversize"), str(job_dir / "sample.bam"), 2)

    assert parent_was_swapped
    assert list(original_directory.iterdir()) == []
    assert list(swapped_target.iterdir()) == []
