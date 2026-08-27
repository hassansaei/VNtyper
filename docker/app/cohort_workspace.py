"""Fresh, descriptor-bound workspace for one cohort analysis."""

from __future__ import annotations

import logging
import os
import stat
import tempfile
from contextlib import suppress
from pathlib import Path
from types import TracebackType
from typing import Literal

logger = logging.getLogger(__name__)

_CLOSE_ON_EXEC = getattr(os, "O_CLOEXEC", 0)
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | _CLOSE_ON_EXEC
_RESERVATION_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | _CLOSE_ON_EXEC


def _same_identity(left: os.stat_result, right: os.stat_result) -> bool:
    return (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)


def _clear_directory(descriptor: int) -> None:
    """Remove directory contents through an already-bound descriptor."""
    for name in os.listdir(descriptor):
        entry = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        if not stat.S_ISDIR(entry.st_mode):
            os.unlink(name, dir_fd=descriptor)
            continue

        child_descriptor = os.open(name, _DIRECTORY_FLAGS, dir_fd=descriptor)
        try:
            child = os.fstat(child_descriptor)
            if not _same_identity(entry, child):
                raise OSError(f"cohort staging child {name} changed while being opened")
            _clear_directory(child_descriptor)
            current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            if not _same_identity(child, current):
                raise OSError(f"cohort staging child {name} changed during cleanup")
            os.rmdir(name, dir_fd=descriptor)
        finally:
            os.close(child_descriptor)


class CohortWorkspace:
    """An exclusively reserved public name and private task output directory."""

    def __init__(
        self,
        path: Path,
        parent_descriptor: int,
        reservation_descriptor: int,
        reservation_metadata: os.stat_result,
        staging_path: Path,
        staging_name: str,
        descriptor: int,
        metadata: os.stat_result,
    ) -> None:
        self.path = path
        self.parent_descriptor = parent_descriptor
        self.reservation_descriptor = reservation_descriptor
        self.reservation_metadata = reservation_metadata
        self.staging_path = staging_path
        self.staging_name = staging_name
        self.descriptor = descriptor
        self.metadata = metadata
        self.bound_path = f"/proc/{os.getpid()}/fd/{descriptor}"
        self._closed = False

    def __enter__(self) -> CohortWorkspace:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> Literal[False]:
        try:
            self.close()
        except BaseException as cleanup_error:
            if exc_value is None:
                raise
            logger.error(f"Cohort analysis failed and workspace cleanup also failed: {cleanup_error}")
        return False

    def child_path(self, name: str) -> str:
        """Return a child path anchored through the open workspace descriptor."""
        return f"{self.bound_path}/{name}"

    def unlink(self, name: str) -> None:
        """Unlink one workspace entry without resolving the staging directory name."""
        os.unlink(name, dir_fd=self.descriptor)

    def close(self) -> None:
        """Remove private output and release the public name only while still owned.

        Raises:
            BaseException: The first cleanup failure, after retrying it once and
                attempting every independent cleanup operation.
        """
        if self._closed:
            return
        self._closed = True
        cleanup_failure: BaseException | None = None

        def record_cleanup_failure(error: BaseException, message: str) -> None:
            nonlocal cleanup_failure
            logger.error(f"{message}: {error}")
            if cleanup_failure is None:
                cleanup_failure = error
            else:
                logger.error(f"Multiple cohort workspace cleanup operations failed: {error}")

        def remove_staging() -> None:
            _clear_directory(self.descriptor)
            current_staging = os.stat(
                self.staging_name,
                dir_fd=self.parent_descriptor,
                follow_symlinks=False,
            )
            if stat.S_ISDIR(current_staging.st_mode) and _same_identity(self.metadata, current_staging):
                os.rmdir(self.staging_name, dir_fd=self.parent_descriptor)
            else:
                logger.error(f"Error deleting cohort staging directory {self.staging_path}: directory changed")

        def remove_reservation() -> None:
            current_reservation = os.stat(
                self.path.name,
                dir_fd=self.parent_descriptor,
                follow_symlinks=False,
            )
            if _same_identity(self.reservation_metadata, current_reservation):
                os.unlink(self.path.name, dir_fd=self.parent_descriptor)
            else:
                logger.error(f"Error deleting cohort output reservation {self.path}: entry changed")

        for attempt in range(2):
            try:
                remove_staging()
                break
            except BaseException as error:
                suffix = "" if attempt == 0 else " after retry"
                record_cleanup_failure(error, f"Error deleting cohort staging directory {self.staging_path}{suffix}")

        for attempt in range(2):
            try:
                remove_reservation()
                break
            except BaseException as error:
                suffix = "" if attempt == 0 else " after retry"
                record_cleanup_failure(error, f"Error deleting cohort output reservation {self.path}{suffix}")

        for descriptor in (self.descriptor, self.reservation_descriptor, self.parent_descriptor):
            try:
                os.close(descriptor)
            except BaseException as error:
                if cleanup_failure is None:
                    cleanup_failure = error
                else:
                    logger.error(f"Multiple cohort workspace cleanup operations failed: {error}")
        if cleanup_failure is not None:
            raise cleanup_failure


def cohort_workspace(output_dir: str | Path) -> CohortWorkspace:
    """Reserve a fresh public name and bind a private cohort output directory.

    Args:
        output_dir: Public task output name to reserve exclusively.

    Returns:
        CohortWorkspace: Owner of the reservation and private staging directory.

    Raises:
        RuntimeError: If the public name cannot be reserved or private staging
            cannot be freshly created and bound without following links.
    """
    path = Path(output_dir)
    parent_descriptor = -1
    reservation_descriptor = -1
    descriptor = -1
    reservation_metadata: os.stat_result | None = None
    expected_staging_metadata: os.stat_result | None = None
    metadata: os.stat_result | None = None
    staging_path: Path | None = None
    staging_name: str | None = None
    try:
        parent_descriptor = os.open(path.parent, _DIRECTORY_FLAGS)
        reservation_descriptor = os.open(
            path.name,
            _RESERVATION_FLAGS,
            0o600,
            dir_fd=parent_descriptor,
        )
        reservation_metadata = os.fstat(reservation_descriptor)
        if not stat.S_ISREG(reservation_metadata.st_mode):
            raise OSError("created cohort output reservation is not a regular file")

        staging_path = Path(tempfile.mkdtemp(prefix=f".{path.name}.cohort-", dir=path.parent))
        staging_name = staging_path.name
        expected_staging_metadata = os.stat(staging_name, dir_fd=parent_descriptor, follow_symlinks=False)
        descriptor = os.open(staging_name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
        metadata = os.fstat(descriptor)
        if not stat.S_ISDIR(expected_staging_metadata.st_mode) or not _same_identity(
            expected_staging_metadata, metadata
        ):
            raise OSError("private cohort staging directory changed while being opened")
    except OSError as error:
        if descriptor >= 0:
            if metadata is not None:
                with suppress(OSError):
                    _clear_directory(descriptor)
            os.close(descriptor)
        if staging_name is not None and parent_descriptor >= 0:
            with suppress(OSError):
                current_staging = os.stat(staging_name, dir_fd=parent_descriptor, follow_symlinks=False)
                owned_staging = metadata or expected_staging_metadata
                if owned_staging is not None and _same_identity(owned_staging, current_staging):
                    os.rmdir(staging_name, dir_fd=parent_descriptor)
        if reservation_descriptor >= 0:
            reservation_metadata = reservation_metadata or os.fstat(reservation_descriptor)
            if parent_descriptor >= 0:
                with suppress(OSError):
                    current = os.stat(path.name, dir_fd=parent_descriptor, follow_symlinks=False)
                    if _same_identity(reservation_metadata, current):
                        os.unlink(path.name, dir_fd=parent_descriptor)
            os.close(reservation_descriptor)
        if parent_descriptor >= 0:
            os.close(parent_descriptor)
        msg = f"Refusing cohort analysis {path.name}: its output directory could not be freshly created: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error
    assert reservation_metadata is not None
    assert staging_path is not None
    assert staging_name is not None
    assert metadata is not None
    return CohortWorkspace(
        path,
        parent_descriptor,
        reservation_descriptor,
        reservation_metadata,
        staging_path,
        staging_name,
        descriptor,
        metadata,
    )
