"""Identity-bound spool input, worker-local snapshots, and bound output.

The shipped boundary covers an actor that can mutate the legacy shared input
mount. Outputs use a service-private named volume mounted only into the API and
worker; an operator who replaces it with a shared bind weakens this boundary.
Arbitrary same-UID code inside the worker namespace remains out of scope: it
can access the worker's descriptors and local ``/tmp`` snapshots directly.
"""

from __future__ import annotations

import hashlib
import logging
import os
import stat
import sys
import tempfile
from collections.abc import Mapping, Sequence
from contextlib import suppress
from uuid import UUID

logger = logging.getLogger(__name__)

_CLOSE_ON_EXEC = getattr(os, "O_CLOEXEC", 0)
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | _CLOSE_ON_EXEC
_FILE_FLAGS = os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | _CLOSE_ON_EXEC
_SNAPSHOT_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | _CLOSE_ON_EXEC
_COPY_CHUNK_BYTES = 1024 * 1024
_SNAPSHOT_ROOT = "/tmp"


def _expected_identity(value: object, label: str) -> tuple[int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
        or any(not isinstance(part, int) or isinstance(part, bool) or part < 0 for part in value)
    ):
        msg = f"Refusing pipeline job: invalid {label} identity token"
        logger.error(msg)
        raise RuntimeError(msg)
    return value[0], value[1]


def _expected_digest(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        msg = f"Refusing pipeline job: invalid {label} SHA-256 token"
        logger.error(msg)
        raise RuntimeError(msg)
    return value


def _same_identity(metadata: os.stat_result, expected: tuple[int, int]) -> bool:
    return (metadata.st_dev, metadata.st_ino) == expected


def _close_after_failure(descriptor: int, label: str) -> None:
    """Close a partially-opened resource without replacing its primary error."""
    try:
        os.close(descriptor)
    except OSError as error:
        logger.warning(f"Closing refused pipeline {label} also failed: {error}")


def _open_directory(path: str, expected: tuple[int, int], label: str) -> int:
    descriptor = os.open(path, _DIRECTORY_FLAGS)
    try:
        metadata = os.fstat(descriptor)
    except BaseException:
        _close_after_failure(descriptor, f"{label} directory")
        raise
    if not stat.S_ISDIR(metadata.st_mode) or not _same_identity(metadata, expected):
        _close_after_failure(descriptor, f"{label} directory")
        msg = f"Refusing pipeline job: {label} directory identity changed before worker start"
        logger.error(msg)
        raise RuntimeError(msg)
    return descriptor


def _open_file(directory_descriptor: int, path: str, expected: tuple[int, int], label: str) -> int:
    descriptor = os.open(os.path.basename(path), _FILE_FLAGS, dir_fd=directory_descriptor)
    try:
        metadata = os.fstat(descriptor)
    except BaseException:
        _close_after_failure(descriptor, label)
        raise
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1 or not _same_identity(metadata, expected):
        _close_after_failure(descriptor, label)
        msg = f"Refusing pipeline job: {label} identity changed before worker start"
        logger.error(msg)
        raise RuntimeError(msg)
    return descriptor


def _clear_directory(descriptor: int) -> None:
    for name in os.listdir(descriptor):
        metadata = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        if not stat.S_ISDIR(metadata.st_mode):
            os.unlink(name, dir_fd=descriptor)
            continue
        child_descriptor = os.open(name, _DIRECTORY_FLAGS, dir_fd=descriptor)
        try:
            opened = os.fstat(child_descriptor)
            if not _same_identity(opened, (metadata.st_dev, metadata.st_ino)):
                raise OSError(f"pipeline output child {name} changed while being opened")
            _clear_directory(child_descriptor)
            current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            if not _same_identity(current, (opened.st_dev, opened.st_ino)):
                raise OSError(f"pipeline output child {name} changed during cleanup")
            os.rmdir(name, dir_fd=descriptor)
        except BaseException:
            _close_after_failure(child_descriptor, "output child")
            raise
        else:
            os.close(child_descriptor)


def _validate_spool_path(alignment_path: str, spool_root: str | None) -> None:
    if spool_root is None:
        return
    root = os.path.realpath(spool_root)
    parent = os.path.realpath(os.path.dirname(alignment_path))
    if os.path.commonpath((root, parent)) != root:
        msg = "Refusing pipeline job: alignment path is outside the protected handoff spool"
        logger.error(msg)
        raise RuntimeError(msg)


def _validate_snapshot_root(shared_roots: Sequence[str]) -> None:
    snapshot_root = os.path.realpath(_SNAPSHOT_ROOT)
    snapshot_metadata = os.stat(_SNAPSHOT_ROOT, follow_symlinks=False)
    if not stat.S_ISDIR(snapshot_metadata.st_mode):
        msg = "Refusing pipeline job: worker snapshot root is not a real directory"
        logger.error(msg)
        raise RuntimeError(msg)
    for configured_root in shared_roots:
        root = os.path.realpath(configured_root)
        try:
            root_metadata = os.stat(configured_root)
        except OSError as error:
            msg = f"Refusing pipeline job: shared root cannot be verified: {error}"
            logger.error(msg)
            raise RuntimeError(msg) from error
        aliases_snapshot_root = _same_identity(
            root_metadata,
            (snapshot_metadata.st_dev, snapshot_metadata.st_ino),
        )
        if aliases_snapshot_root or os.path.commonpath((snapshot_root, root)) == root:
            msg = "Refusing pipeline job: worker snapshot root aliases a shared service root"
            logger.error(msg)
            raise RuntimeError(msg)


def _copy_snapshot(
    source_descriptor: int,
    destination_directory: int,
    name: str,
    expected_sha256: str,
    label: str,
) -> None:
    destination = os.open(name, _SNAPSHOT_FLAGS, 0o400, dir_fd=destination_directory)
    digest = hashlib.sha256()
    try:
        os.lseek(source_descriptor, 0, os.SEEK_SET)
        while True:
            chunk = os.read(source_descriptor, _COPY_CHUNK_BYTES)
            if not chunk:
                break
            digest.update(chunk)
            view = memoryview(chunk)
            while view:
                written = os.write(destination, view)
                view = view[written:]
        os.fsync(destination)
        metadata = os.fstat(destination)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise OSError(f"worker-local {label} snapshot is not a single-link regular file")
        if digest.hexdigest() != expected_sha256:
            msg = f"Refusing pipeline job: {label} SHA-256 changed before worker start"
            logger.error(msg)
            raise RuntimeError(msg)
    except BaseException:
        _close_after_failure(destination, f"{label} snapshot")
        raise
    else:
        os.close(destination)


def _remove_owned_entries(
    descriptor: int,
    entries: Sequence[tuple[str | None, tuple[int, int] | None]],
) -> None:
    for path, expected in entries:
        if path is None or expected is None:
            continue
        name = os.path.basename(path)
        try:
            current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            if _same_identity(current, expected):
                os.unlink(name, dir_fd=descriptor)
        except FileNotFoundError:
            pass


def _canonical_job_id(value: str) -> bool:
    try:
        parsed = UUID(value)
    except ValueError:
        return False
    return parsed.version == 4 and str(parsed) == value


def reclaim_unopened_spool_inputs(
    alignment_path: str,
    index_path: str | None,
    output_dir: str,
    spool_root: str,
) -> None:
    """Reclaim a refused task's exact UUID-scoped protected-spool files.

    This fallback is used only when a malformed broker token prevents the
    identity-bound workspace from opening. The protected spool is not mounted
    into the scoped shared-volume actor, so its exact UUID directory is a
    narrow cleanup capability; arbitrary and cross-job paths are refused.
    """
    job_id = os.path.basename(os.path.normpath(output_dir))
    alignment_directory = os.path.dirname(os.path.abspath(alignment_path))
    expected_directory = os.path.join(os.path.abspath(spool_root), job_id)
    if not _canonical_job_id(job_id):
        msg = "Refusing fallback cleanup: output job name is not a canonical UUID"
        logger.error(msg)
        raise RuntimeError(msg)
    if alignment_directory != expected_directory:
        root_prefix = os.path.abspath(spool_root) + os.sep
        if not alignment_directory.startswith(root_prefix):
            msg = "Refusing fallback cleanup: alignment path is outside the protected handoff spool"
        else:
            msg = "Refusing fallback cleanup: alignment path does not match the output job"
        logger.error(msg)
        raise RuntimeError(msg)
    paths = (alignment_path,) if index_path is None else (alignment_path, index_path)
    if any(
        os.path.dirname(os.path.abspath(path)) != expected_directory or os.path.basename(path) in {"", ".", ".."}
        for path in paths
    ):
        msg = "Refusing fallback cleanup: input paths do not share the exact protected spool job directory"
        logger.error(msg)
        raise RuntimeError(msg)

    root_descriptor = os.open(os.path.abspath(spool_root), _DIRECTORY_FLAGS)
    job_descriptor = -1
    try:
        try:
            job_descriptor = os.open(job_id, _DIRECTORY_FLAGS, dir_fd=root_descriptor)
        except FileNotFoundError:
            return
        for path in paths:
            name = os.path.basename(path)
            try:
                metadata = os.stat(name, dir_fd=job_descriptor, follow_symlinks=False)
            except FileNotFoundError:
                continue
            if stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1:
                os.unlink(name, dir_fd=job_descriptor)
        if not os.listdir(job_descriptor):
            os.rmdir(job_id, dir_fd=root_descriptor)
    finally:
        active_error = sys.exc_info()[1]
        first_close_error: OSError | None = None
        for descriptor in (job_descriptor, root_descriptor):
            if descriptor < 0:
                continue
            try:
                os.close(descriptor)
            except OSError as error:
                if first_close_error is None:
                    first_close_error = error
        if first_close_error is not None:
            if active_error is None:
                raise first_close_error
            logger.warning(f"Closing fallback cleanup descriptors also failed: {first_close_error}")


class PipelineJobWorkspace:
    """Retained spool/output handles and worker-local input snapshots."""

    def __init__(
        self,
        input_descriptor: int,
        output_descriptor: int,
        view_descriptor: int,
        alignment_descriptor: int,
        index_descriptor: int | None,
        alignment_path: str,
        index_path: str | None,
        output_dir: str,
        view_dir: str,
        input_identity: tuple[int, int],
        output_identity: tuple[int, int],
        view_identity: tuple[int, int],
        alignment_identity: tuple[int, int],
        index_identity: tuple[int, int] | None,
        alignment_view_name: str,
        index_view_name: str | None,
    ) -> None:
        self.input_descriptor = input_descriptor
        self.output_descriptor = output_descriptor
        self.view_descriptor = view_descriptor
        self.alignment_descriptor = alignment_descriptor
        self.index_descriptor = index_descriptor
        self.bound_alignment_path = os.path.join(view_dir, alignment_view_name)
        self.bound_output_path = f"/proc/{os.getpid()}/fd/{output_descriptor}/."
        self.bound_index_path = None if index_view_name is None else os.path.join(view_dir, index_view_name)
        self._alignment_path = alignment_path
        self._index_path = index_path
        self._output_dir = output_dir
        self._view_dir = view_dir
        self._input_identity = input_identity
        self._output_identity = output_identity
        self._view_identity = view_identity
        self._alignment_identity = alignment_identity
        self._index_identity = index_identity
        self._alignment_view_name = alignment_view_name
        self._index_view_name = index_view_name
        self._views_removed = False
        self.output_displaced = False
        self._closed = False

    def require_current_output(self, phase: str = "worker execution") -> None:
        """Refuse success when the public output name lost its captured identity."""
        try:
            current = os.stat(self._output_dir, follow_symlinks=False)
        except OSError as error:
            self.output_displaced = True
            msg = f"Refusing pipeline job: output directory identity changed during {phase}: {error}"
            logger.error(msg)
            raise RuntimeError(msg) from error
        if not stat.S_ISDIR(current.st_mode) or not _same_identity(current, self._output_identity):
            self.output_displaced = True
            msg = f"Refusing pipeline job: output directory identity changed during {phase}"
            logger.error(msg)
            raise RuntimeError(msg)

    def remove_inputs(self) -> None:
        """Remove only spool entries named by API identity tokens."""
        for path, expected in (
            (self._alignment_path, self._alignment_identity),
            (self._index_path, self._index_identity),
        ):
            if path is None or expected is None:
                continue
            try:
                current = os.stat(os.path.basename(path), dir_fd=self.input_descriptor, follow_symlinks=False)
                if _same_identity(current, expected):
                    os.unlink(os.path.basename(path), dir_fd=self.input_descriptor)
            except FileNotFoundError:
                pass
            except OSError as error:
                logger.error(f"Error deleting input file {path}: {error}")
        public_input = os.path.dirname(self._alignment_path)
        try:
            current = os.stat(public_input, follow_symlinks=False)
            if _same_identity(current, self._input_identity):
                if os.listdir(self.input_descriptor):
                    logger.warning(f"Input directory {public_input} still holds files and was left in place")
                else:
                    os.rmdir(public_input)
                    logger.info(f"Deleted empty handoff directory: {public_input}")
        except OSError as error:
            logger.error(f"Error deleting input directory {public_input}: {error}")

    def remove_views(self) -> None:
        """Remove worker-local real-file snapshots without following public paths."""
        if self._views_removed:
            return
        for name in (self._index_view_name, self._alignment_view_name):
            if name is not None:
                with suppress(FileNotFoundError):
                    os.unlink(name, dir_fd=self.view_descriptor)
        try:
            current = os.stat(self._view_dir, follow_symlinks=False)
        except OSError:
            pass
        else:
            if stat.S_ISDIR(current.st_mode) and _same_identity(current, self._view_identity):
                os.rmdir(self._view_dir)
        self._views_removed = True

    def discard_detached_output(self) -> None:
        """Clear results from a captured output inode that lost its public name."""
        try:
            current = os.stat(self._output_dir, follow_symlinks=False)
        except OSError:
            self.output_displaced = True
        else:
            if not stat.S_ISDIR(current.st_mode) or not _same_identity(current, self._output_identity):
                self.output_displaced = True
        if self.output_displaced:
            _clear_directory(self.output_descriptor)

    def remove_output(self) -> None:
        """Clear and remove the still-current public output directory."""
        self.require_current_output()
        _clear_directory(self.output_descriptor)
        self.require_current_output()
        os.rmdir(self._output_dir)

    def close(self) -> None:
        """Close every retained descriptor and report the first close error."""
        if self._closed:
            return
        self._closed = True
        first_error: OSError | None = None
        for descriptor in (
            self.index_descriptor,
            self.alignment_descriptor,
            self.view_descriptor,
            self.output_descriptor,
            self.input_descriptor,
        ):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError as error:
                    if first_error is None:
                        first_error = error
        if first_error is not None:
            raise first_error


def open_pipeline_job_workspace(
    alignment_path: str,
    index_path: str | None,
    output_dir: str,
    identity: Mapping[str, object] | None,
    *,
    spool_root: str | None = None,
    shared_roots: Sequence[str] = (),
) -> PipelineJobWorkspace:
    """Bind protected spool inputs, verify digests, and create local snapshots."""
    if identity is None:
        msg = "Refusing pipeline job: missing workspace identity token"
        logger.error(msg)
        raise RuntimeError(msg)
    if not isinstance(identity, Mapping):
        msg = "Refusing pipeline job: invalid workspace identity token"
        logger.error(msg)
        raise RuntimeError(msg)

    input_descriptor = output_descriptor = view_descriptor = alignment_descriptor = -1
    index_descriptor: int | None = None
    view_dir: str | None = None
    view_identity: tuple[int, int] | None = None
    alignment_view_name: str | None = None
    index_view_name: str | None = None
    input_identity: tuple[int, int] | None = None
    alignment_identity: tuple[int, int] | None = None
    expected_index_identity: tuple[int, int] | None = None
    try:
        _validate_spool_path(alignment_path, spool_root)
        _validate_snapshot_root(shared_roots)
        input_identity = _expected_identity(identity.get("input_dir"), "input directory")
        output_identity = _expected_identity(identity.get("output_dir"), "output directory")
        alignment_identity = _expected_identity(identity.get("alignment"), "alignment")
        alignment_sha256 = _expected_digest(identity.get("alignment_sha256"), "alignment")
        input_descriptor = _open_directory(os.path.dirname(alignment_path), input_identity, "input")
        output_descriptor = _open_directory(output_dir, output_identity, "output")
        alignment_descriptor = _open_file(input_descriptor, alignment_path, alignment_identity, "alignment")

        index_identity = identity.get("index")
        if (index_path is None) != (index_identity is None and identity.get("index_sha256") is None):
            msg = "Refusing pipeline job: index path and identity token disagree"
            logger.error(msg)
            raise RuntimeError(msg)
        expected_index_sha256: str | None = None
        if index_identity is not None:
            if index_path is None:
                msg = "Refusing pipeline job: index identity has no index path"
                logger.error(msg)
                raise RuntimeError(msg)
            expected_index_identity = _expected_identity(index_identity, "index")
            expected_index_sha256 = _expected_digest(identity.get("index_sha256"), "index")
            index_descriptor = _open_file(input_descriptor, index_path, expected_index_identity, "index")
        elif identity.get("index_sha256") is not None:
            msg = "Refusing pipeline job: index SHA-256 has no index identity"
            logger.error(msg)
            raise RuntimeError(msg)

        view_dir = tempfile.mkdtemp(prefix="vntyper-pipeline-snapshot-", dir=_SNAPSHOT_ROOT)
        created_view = os.stat(view_dir, follow_symlinks=False)
        view_identity = (created_view.st_dev, created_view.st_ino)
        view_descriptor = _open_directory(view_dir, view_identity, "worker snapshot")
        alignment_view_name = os.path.basename(alignment_path)
        _copy_snapshot(
            alignment_descriptor,
            view_descriptor,
            alignment_view_name,
            alignment_sha256,
            "alignment",
        )
        if index_descriptor is not None and expected_index_sha256 is not None and index_path is not None:
            index_view_name = os.path.basename(index_path)
            _copy_snapshot(index_descriptor, view_descriptor, index_view_name, expected_index_sha256, "index")

        return PipelineJobWorkspace(
            input_descriptor,
            output_descriptor,
            view_descriptor,
            alignment_descriptor,
            index_descriptor,
            alignment_path,
            index_path,
            output_dir,
            view_dir,
            input_identity,
            output_identity,
            view_identity,
            alignment_identity,
            expected_index_identity,
            alignment_view_name,
            index_view_name,
        )
    except (OSError, RuntimeError) as error:
        if input_descriptor >= 0 and input_identity is not None and alignment_identity is not None:
            with suppress(OSError):
                _remove_owned_entries(
                    input_descriptor,
                    ((alignment_path, alignment_identity), (index_path, expected_index_identity)),
                )
            public_input = os.path.dirname(alignment_path)
            with suppress(OSError):
                current = os.stat(public_input, follow_symlinks=False)
                if _same_identity(current, input_identity) and not os.listdir(input_descriptor):
                    os.rmdir(public_input)
        if view_descriptor >= 0:
            with suppress(OSError):
                _clear_directory(view_descriptor)
        if view_dir is not None and view_identity is not None:
            with suppress(OSError):
                current_view = os.stat(view_dir, follow_symlinks=False)
                if stat.S_ISDIR(current_view.st_mode) and _same_identity(current_view, view_identity):
                    os.rmdir(view_dir)
        for descriptor in (
            index_descriptor,
            alignment_descriptor,
            view_descriptor,
            output_descriptor,
            input_descriptor,
        ):
            if descriptor is not None and descriptor >= 0:
                with suppress(OSError):
                    os.close(descriptor)
        if isinstance(error, RuntimeError):
            raise
        msg = f"Refusing pipeline job: workspace could not be safely bound: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error
