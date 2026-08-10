"""Atomic output installation and path separation for the mutation harness."""

from __future__ import annotations

import contextlib
import os
import signal
import stat
import uuid
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TextIO

_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)


@contextlib.contextmanager
def _defer_termination() -> Iterator[None]:
    """Defer SIGINT/SIGTERM across descriptor ownership transitions."""
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
    try:
        yield
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous)


def _canonical_identity(path: Path) -> tuple[Path, tuple[int, int] | None]:
    """Return a canonical name and existing inode identity for one path."""
    canonical = path.resolve(strict=False)
    try:
        status = path.stat()
    except FileNotFoundError:
        return canonical, None
    return canonical, (status.st_dev, status.st_ino)


def _overlap(left: tuple[Path, tuple[int, int] | None], right: tuple[Path, tuple[int, int] | None]) -> bool:
    """Return whether two paths name or contain the same storage."""
    left_path, left_identity = left
    right_path, right_identity = right
    return (
        left_path == right_path
        or left_path.is_relative_to(right_path)
        or right_path.is_relative_to(left_path)
        or (left_identity is not None and left_identity == right_identity)
    )


def validate_disjoint_paths(
    artifacts: Sequence[tuple[str, Path | None]],
    protected_sources: Sequence[tuple[str, Path]] = (),
) -> None:
    """Require artifact paths to be distinct and disjoint from selected sources.

    Args:
        artifacts: Named inputs and output destinations involved in one operation.
        protected_sources: Named mutation sources that no artifact may alias or contain.

    Raises:
        ValueError: If two artifacts overlap by name, ancestry, symlink resolution, or
            existing inode, or if an artifact overlaps a selected source.
    """
    present = [(label, path, _canonical_identity(path)) for label, path in artifacts if path is not None]
    for index, (left_label, _left_path, left_state) in enumerate(present):
        for right_label, _right_path, right_state in present[index + 1 :]:
            if _overlap(left_state, right_state):
                raise ValueError(f"mutation artifact paths must be distinct: {left_label} and {right_label}")

    protected = [(label, path, _canonical_identity(path)) for label, path in protected_sources]
    for artifact_label, _artifact_path, artifact_state in present:
        for source_label, _source_path, source_state in protected:
            if _overlap(artifact_state, source_state):
                raise ValueError(f"mutation artifact overlaps selected source: {artifact_label} and {source_label}")


def atomic_write_text(destination: Path, content: str) -> None:
    """Install complete UTF-8 text at ``destination`` atomically and durably.

    Args:
        destination: Final path to replace.
        content: Complete report content to install.

    Raises:
        OSError: If writing, flushing, syncing, replacement, or cleanup fails.
    """
    parent_fd = os.open(destination.parent, _DIRECTORY_FLAGS)
    temporary_name = ""
    raw_fd: int | None = None
    stream_fd: int | None = None
    stream: TextIO | None = None
    try:
        try:
            previous_status = os.stat(destination.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            previous_mode = 0o600
        else:
            previous_mode = stat.S_IMODE(previous_status.st_mode) if stat.S_ISREG(previous_status.st_mode) else 0o600

        while raw_fd is None:
            candidate = f".{destination.name}.mutation-{uuid.uuid4().hex}"
            with _defer_termination():
                try:
                    raw_fd = os.open(
                        candidate,
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
                        previous_mode,
                        dir_fd=parent_fd,
                    )
                except FileExistsError:
                    continue
                temporary_name = candidate

        os.fchmod(raw_fd, previous_mode)
        with _defer_termination():
            stream_fd = os.dup(raw_fd)
        try:
            with _defer_termination():
                stream = os.fdopen(stream_fd, "w", encoding="utf-8")
                stream_fd = None
        except BaseException:
            try:
                if stream is not None:
                    stream.close()
                    stream = None
            finally:
                if stream_fd is not None:
                    with _defer_termination():
                        os.close(stream_fd)
                        stream_fd = None
            raise
        assert stream is not None
        with stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        stream = None

        with _defer_termination():
            os.close(raw_fd)
            raw_fd = None
        os.replace(
            temporary_name,
            destination.name,
            src_dir_fd=parent_fd,
            dst_dir_fd=parent_fd,
        )
        temporary_name = ""
        os.fsync(parent_fd)
    finally:
        try:
            if stream_fd is not None:
                with _defer_termination():
                    os.close(stream_fd)
                    stream_fd = None
        finally:
            try:
                if raw_fd is not None:
                    with _defer_termination():
                        os.close(raw_fd)
                        raw_fd = None
            finally:
                try:
                    if temporary_name:
                        with _defer_termination(), contextlib.suppress(FileNotFoundError):
                            os.unlink(temporary_name, dir_fd=parent_fd)
                            temporary_name = ""
                finally:
                    with _defer_termination():
                        os.close(parent_fd)
