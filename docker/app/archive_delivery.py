"""Descriptor-bound reads of exclusively owned result archives."""

from __future__ import annotations

import logging
import os
import shutil
import stat
import zipfile
from pathlib import Path

from starlette.responses import FileResponse

logger = logging.getLogger(__name__)

_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
_FILE_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)


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
        with os.fdopen(os.dup(descriptor), "rb") as source, archive.open(arcname, "w", force_zip64=True) as target:
            shutil.copyfileobj(source, target)
        if os.fstat(descriptor).st_nlink != 1:
            raise ValueError(f"Unsafe result archive gained another hard link while being read: {path}")
    finally:
        os.close(descriptor)


def snapshot_owned_archives(paths: list[str], snapshot_dir: str | Path) -> list[str]:
    """Copy member archives from bound descriptors into task-owned files.

    Args:
        paths: Member archive paths supplied to the cohort task.
        snapshot_dir: Newly created private directory for the copies.

    Returns:
        Paths of the task-owned snapshots, in input order.
    """
    directory = Path(snapshot_dir)
    directory.mkdir(mode=0o700)
    directory_descriptor = os.open(directory, _DIRECTORY_FLAGS)
    snapshots: list[str] = []
    try:
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
                with (
                    os.fdopen(os.dup(source_descriptor), "rb") as source,
                    os.fdopen(os.dup(target_descriptor), "wb") as target,
                ):
                    shutil.copyfileobj(source, target)
                if os.fstat(source_descriptor).st_nlink != 1:
                    raise ValueError(f"Unsafe result archive gained another hard link while being read: {path}")
                snapshots.append(str(directory / name))
            finally:
                if target_descriptor >= 0:
                    os.close(target_descriptor)
                os.close(source_descriptor)
    except BaseException:
        shutil.rmtree(directory)
        raise
    finally:
        os.close(directory_descriptor)
    return snapshots
