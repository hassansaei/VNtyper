"""Descriptor-relative no-follow reads for calibration custody boundaries."""

from __future__ import annotations

import logging
import os
import stat
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_DIRECTORY = getattr(os, "O_DIRECTORY", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)


@dataclass
class SecureDirectoryReader:
    """One pinned directory descriptor with closed regular-file inventory."""

    path: Path
    descriptor: int
    names: frozenset[str]

    @classmethod
    def open(cls, path: Path, expected_names: set[str]) -> SecureDirectoryReader:
        """Open and pin a nonsymlink directory with an exact regular-file inventory."""
        if not isinstance(path, Path) or not _NOFOLLOW:
            raise ValueError("secure calibration directory reads require Path and O_NOFOLLOW support")
        try:
            descriptor = os.open(path, os.O_RDONLY | _DIRECTORY | _CLOEXEC | _NOFOLLOW)
        except OSError as error:
            raise ValueError(f"secure calibration directory is missing, unreadable, or a symlink: {path}") from error
        try:
            if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise ValueError("secure calibration import root must be a directory")
            names = frozenset(os.listdir(descriptor))
            if names != frozenset(expected_names):
                raise ValueError("secure calibration import inventory differs")
            for name in names:
                metadata = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
                if not stat.S_ISREG(metadata.st_mode):
                    raise ValueError("secure calibration import entries must be regular non-symlink files")
        except BaseException:
            os.close(descriptor)
            raise
        return cls(path, descriptor, names)

    def read_files(self, names: tuple[str, ...]) -> dict[str, bytes]:
        """Open all named regular files first, then read each from that same descriptor."""
        descriptors: dict[str, int] = {}
        try:
            for name in names:
                if name not in self.names:
                    raise ValueError(f"secure calibration import file is undeclared: {name}")
                descriptor = os.open(name, os.O_RDONLY | _CLOEXEC | _NOFOLLOW, dir_fd=self.descriptor)
                if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                    os.close(descriptor)
                    raise ValueError("secure calibration import entry changed from a regular file")
                descriptors[name] = descriptor
            return {name: _read_descriptor(descriptor) for name, descriptor in descriptors.items()}
        except OSError as error:
            raise ValueError("secure calibration import file changed, is unreadable, or is a symlink") from error
        finally:
            for descriptor in descriptors.values():
                os.close(descriptor)

    def read_file(self, name: str) -> bytes:
        """Read one regular file relative to the pinned directory without following links."""
        return self.read_files((name,))[name]

    def close(self) -> None:
        """Close the pinned directory descriptor exactly once."""
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1

    def __enter__(self) -> SecureDirectoryReader:
        return self

    def __exit__(self, _type, _value, _traceback) -> None:
        self.close()


def read_regular_path(path: Path) -> bytes:
    """Read one exact regular path from the same no-follow file descriptor."""
    try:
        descriptor = os.open(path, os.O_RDONLY | _CLOEXEC | _NOFOLLOW)
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError("secure calibration payload must be a regular file")
        return _read_descriptor(descriptor)
    except OSError as error:
        raise ValueError("secure calibration payload is unreadable or is a symlink") from error
    finally:
        if "descriptor" in locals():
            os.close(descriptor)


def _read_descriptor(descriptor: int) -> bytes:
    chunks: list[bytes] = []
    while True:
        chunk = os.read(descriptor, 1024 * 1024)
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)
