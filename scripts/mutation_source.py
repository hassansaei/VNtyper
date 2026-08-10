"""Descriptor-anchored source reads and atomic mutation installation."""

from __future__ import annotations

import contextlib
import os
import signal
import stat
import uuid
from collections.abc import Iterator
from pathlib import Path

from mutation_workspace_fs import RootCapability, parent_directory_fd

WorkspaceRoot = Path | RootCapability


@contextlib.contextmanager
def _defer_termination() -> Iterator[None]:
    """Defer SIGINT/SIGTERM across descriptor ownership transitions."""
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
    try:
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous)


def read_source_bytes(root: WorkspaceRoot, relative: str) -> bytes:
    """Read one regular source file without following its final component.

    Args:
        root: Open workspace capability or workspace directory to borrow.
        relative: Safe workspace-relative source name.

    Returns:
        The complete source bytes.

    Raises:
        RuntimeError: If the selected entry is not a regular file.
        OSError: If anchored traversal or reading fails.
    """
    with parent_directory_fd(root, relative, create=False) as (parent_fd, name):
        try:
            descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent_fd)
        except OSError as exc:
            raise RuntimeError(f"mutation target is not a regular file: {relative}") from exc
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise RuntimeError(f"mutation target is not a regular file: {relative}")
            chunks: list[bytes] = []
            while chunk := os.read(descriptor, 1024 * 1024):
                chunks.append(chunk)
            return b"".join(chunks)
        finally:
            os.close(descriptor)


def atomic_replace_source(root: WorkspaceRoot, relative: str, content: bytes) -> None:
    """Atomically install complete source bytes beneath a pinned workspace root.

    The existing entry is never opened for writing. A same-parent temporary regular
    file is fully written and fsynced, then replaces only the selected directory name.
    This breaks a hardlink instead of modifying the shared inode and replaces a final
    symlink entry instead of following it.

    Args:
        root: Open workspace capability or workspace directory to borrow.
        relative: Safe workspace-relative source name.
        content: Complete bytes to install.

    Raises:
        RuntimeError: If the selected entry is not a regular file.
        OSError: If creation, writing, syncing, replacement, or cleanup fails.
    """
    with parent_directory_fd(root, relative, create=False) as (parent_fd, name):
        target_status = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if not stat.S_ISREG(target_status.st_mode):
            raise RuntimeError(f"mutation target is not a regular file: {relative}")

        temporary_name = ""
        temporary_fd: int | None = None
        try:
            while temporary_fd is None:
                candidate = f".{name}.mutation-{uuid.uuid4().hex}"
                with _defer_termination():
                    try:
                        temporary_fd = os.open(
                            candidate,
                            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
                            stat.S_IMODE(target_status.st_mode),
                            dir_fd=parent_fd,
                        )
                    except FileExistsError:
                        continue
                    temporary_name = candidate

            offset = 0
            while offset < len(content):
                written = os.write(temporary_fd, content[offset:])
                if written <= 0:
                    raise OSError("source write made no progress")
                offset += written
            os.fchmod(temporary_fd, stat.S_IMODE(target_status.st_mode))
            os.fsync(temporary_fd)
            with _defer_termination():
                os.close(temporary_fd)
                temporary_fd = None
            os.replace(
                temporary_name,
                name,
                src_dir_fd=parent_fd,
                dst_dir_fd=parent_fd,
            )
            temporary_name = ""
            os.fsync(parent_fd)
        finally:
            try:
                if temporary_fd is not None:
                    with _defer_termination():
                        os.close(temporary_fd)
                        temporary_fd = None
            finally:
                if temporary_name:
                    with _defer_termination(), contextlib.suppress(FileNotFoundError):
                        os.unlink(temporary_name, dir_fd=parent_fd)
                        temporary_name = ""
