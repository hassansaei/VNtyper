"""Descriptor-bound reads of exclusively owned result archives."""

from __future__ import annotations

import logging
import os
import shutil
import stat
import tempfile
import zipfile
from pathlib import Path
from typing import Literal

from starlette.responses import FileResponse

logger = logging.getLogger(__name__)

_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
_FILE_FLAGS = os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)


def _same_identity(left: os.stat_result, right: os.stat_result) -> bool:
    return (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)


def open_owned_regular(path: str | Path) -> tuple[int, os.stat_result]:
    """Open one exclusively owned regular file without following its final path.

    Args:
        path: File to bind.

    Returns:
        Its open descriptor and metadata.

    Raises:
        OSError: If the path cannot be opened without following links.
        ValueError: If it is not a singly linked regular file.
    """
    archive_path = Path(path)
    parent_descriptor = os.open(archive_path.parent, _DIRECTORY_FLAGS)
    descriptor = -1
    try:
        descriptor = os.open(archive_path.name, _FILE_FLAGS, dir_fd=parent_descriptor)
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise ValueError(f"Unsafe result archive is not an exclusively owned regular file: {archive_path}")
        return descriptor, metadata
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        raise
    finally:
        os.close(parent_descriptor)


def _unlink_if_same(path: Path, expected: os.stat_result) -> None:
    parent_descriptor = os.open(path.parent, _DIRECTORY_FLAGS)
    try:
        try:
            current = os.stat(path.name, dir_fd=parent_descriptor, follow_symlinks=False)
        except FileNotFoundError:
            return
        if _same_identity(expected, current):
            os.unlink(path.name, dir_fd=parent_descriptor)
    finally:
        os.close(parent_descriptor)


def _close_snapshot_directory(
    directory: Path,
    parent_descriptor: int,
    directory_descriptor: int,
    directory_metadata: os.stat_result,
    names: tuple[str, ...],
) -> None:
    """Remove bound snapshot members and close their directory descriptors."""
    cleanup_failure: BaseException | None = None

    def record_failure(error: BaseException) -> None:
        nonlocal cleanup_failure
        if cleanup_failure is None:
            cleanup_failure = error
        else:
            logger.error(f"Multiple cohort snapshot cleanup operations failed: {error}")

    for name in names:
        try:
            os.unlink(name, dir_fd=directory_descriptor)
        except FileNotFoundError:
            pass
        except BaseException as error:
            record_failure(error)

    try:
        current = os.stat(directory.name, dir_fd=parent_descriptor, follow_symlinks=False)
        if not _same_identity(directory_metadata, current):
            raise ValueError(f"Cohort member snapshot directory changed during analysis: {directory}")
        os.rmdir(directory.name, dir_fd=parent_descriptor)
    except BaseException as error:
        record_failure(error)

    for descriptor in (directory_descriptor, parent_descriptor):
        try:
            os.close(descriptor)
        except BaseException as error:
            record_failure(error)

    if cleanup_failure is not None:
        raise cleanup_failure


class OwnedArchiveSnapshots:
    """Task-owned archive copies kept reachable through an open directory."""

    def __init__(
        self,
        directory: Path,
        parent_descriptor: int,
        directory_descriptor: int,
        directory_metadata: os.stat_result,
        names: tuple[str, ...],
    ) -> None:
        self.paths = tuple(f"/proc/{os.getpid()}/fd/{directory_descriptor}/{name}" for name in names)
        self._directory = directory
        self._parent_descriptor = parent_descriptor
        self._directory_descriptor = directory_descriptor
        self._directory_metadata = directory_metadata
        self._names = names
        self._closed = False

    def __enter__(self) -> OwnedArchiveSnapshots:
        """Keep the snapshots alive for the surrounding operation."""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> Literal[False]:
        """Release snapshots without replacing an error from their consumer."""
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc_value is None:
                raise
            logger.error(f"Cohort subprocess failed and snapshot cleanup also failed: {cleanup_error}")
        return False

    def close(self) -> None:
        """Delete the bound copies and close their anchoring descriptors."""
        if self._closed:
            return
        self._closed = True
        _close_snapshot_directory(
            self._directory,
            self._parent_descriptor,
            self._directory_descriptor,
            self._directory_metadata,
            self._names,
        )


class OwnedFileResponse(FileResponse):
    """A ``FileResponse`` bound to an already validated descriptor."""

    def __init__(
        self,
        path: str | Path,
        *,
        media_type: str,
        filename: str,
        remove_after_send: bool = False,
    ) -> None:
        self.owned_path = Path(path)
        self.descriptor, metadata = open_owned_regular(self.owned_path)
        self._owned_metadata = metadata
        self._remove_after_send = remove_after_send
        super().__init__(
            f"/proc/self/fd/{self.descriptor}",
            media_type=media_type,
            filename=filename,
            stat_result=metadata,
        )

    async def __call__(self, scope, receive, send) -> None:
        """Stream the bound inode, then close it and optionally remove its exact path."""
        safe_scope = dict(scope)
        safe_scope["extensions"] = {
            name: value for name, value in scope.get("extensions", {}).items() if name != "http.response.pathsend"
        }
        response_failure: BaseException | None = None
        try:
            await super().__call__(safe_scope, receive, send)
        except BaseException as error:
            response_failure = error
            raise
        finally:
            cleanup_failure: BaseException | None = None
            try:
                os.close(self.descriptor)
            except BaseException as error:
                cleanup_failure = error
            if self._remove_after_send:
                try:
                    _unlink_if_same(self.owned_path, self._owned_metadata)
                except BaseException as error:
                    if cleanup_failure is None:
                        cleanup_failure = error
                    else:
                        logger.error(f"Response descriptor and scratch cleanup both failed: {error}")
            if cleanup_failure is not None:
                if response_failure is None:
                    raise cleanup_failure
                logger.error(f"Response failed and archive delivery cleanup also failed: {cleanup_failure}")


def write_owned_zip_member(archive: zipfile.ZipFile, path: str | Path, arcname: str) -> None:
    """Copy one safely bound result archive into another ZIP.

    Args:
        archive: Destination ZIP.
        path: Existing result archive.
        arcname: Member name in the destination.

    Raises:
        OSError: If the source cannot be safely opened or changes ownership.
        ValueError: If the source is not exclusively owned.
    """
    try:
        descriptor, _ = open_owned_regular(path)
    except (OSError, ValueError) as error:
        raise ValueError(f"Missing or unsafe cohort member archive: {path}") from error
    try:
        with tempfile.SpooledTemporaryFile(max_size=8 * 1024 * 1024) as staged:
            with os.fdopen(os.dup(descriptor), "rb") as source:
                shutil.copyfileobj(source, staged)
            if os.fstat(descriptor).st_nlink != 1:
                raise ValueError(f"Unsafe result archive gained another hard link while being read: {path}")
            staged.seek(0)
            with archive.open(arcname, "w", force_zip64=True) as target:
                shutil.copyfileobj(staged, target)
    finally:
        os.close(descriptor)


def snapshot_owned_archives(paths: list[str], snapshot_dir: str | Path) -> OwnedArchiveSnapshots:
    """Copy member archives and return the owner of their anchored lifetime.

    Args:
        paths: Member archive paths supplied to the cohort task.
        snapshot_dir: Newly created private directory for the copies.

    Returns:
        Owner exposing descriptor-bound paths in input order.
    """
    directory = Path(snapshot_dir)
    directory.mkdir(mode=0o700)
    parent_descriptor = os.open(directory.parent, _DIRECTORY_FLAGS)
    directory_descriptor = -1
    names: list[str] = []
    try:
        directory_descriptor = os.open(directory.name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
        directory_metadata = os.fstat(directory_descriptor)
        for path in paths:
            source_descriptor, _ = open_owned_regular(path)
            name = Path(path).name
            target_descriptor = -1
            try:
                target_descriptor = os.open(
                    name,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
                    0o600,
                    dir_fd=directory_descriptor,
                )
                names.append(name)
                with (
                    os.fdopen(os.dup(source_descriptor), "rb") as source,
                    os.fdopen(os.dup(target_descriptor), "wb") as target,
                ):
                    shutil.copyfileobj(source, target)
                if os.fstat(source_descriptor).st_nlink != 1:
                    raise ValueError(f"Unsafe result archive gained another hard link while being read: {path}")
            finally:
                if target_descriptor >= 0:
                    os.close(target_descriptor)
                os.close(source_descriptor)
    except BaseException:
        if directory_descriptor >= 0:
            try:
                _close_snapshot_directory(
                    directory,
                    parent_descriptor,
                    directory_descriptor,
                    os.fstat(directory_descriptor),
                    tuple(names),
                )
            except BaseException as cleanup_error:
                logger.error(f"Snapshot creation failed and cleanup also failed: {cleanup_error}")
        else:
            os.close(parent_descriptor)
        raise
    return OwnedArchiveSnapshots(
        directory,
        parent_descriptor,
        directory_descriptor,
        directory_metadata,
        tuple(names),
    )
