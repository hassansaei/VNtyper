"""Pinned regular-file sources for calibration read fingerprinting."""

from __future__ import annotations

import errno
import hashlib
import logging
import os
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, NoReturn

logger = logging.getLogger(__name__)

_BUFFER_SIZE = 1024 * 1024
_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)
_NONBLOCK = getattr(os, "O_NONBLOCK", 0)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int]:
    return metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns


@dataclass
class PinnedReadSource:
    """One no-follow regular file held open for hashing and parsing."""

    path: Path
    label: str
    descriptor: int
    identity: tuple[int, int, int, int, int]

    @classmethod
    def open(cls, path: Path, label: str) -> PinnedReadSource:
        """Open one regular path without following its final symlink."""
        if not isinstance(path, Path) or not isinstance(label, str) or not label or not _NOFOLLOW or not _NONBLOCK:
            _fail("calibration pinned reads require a Path, label, O_NOFOLLOW, and O_NONBLOCK support")
        try:
            descriptor = os.open(path, os.O_RDONLY | _NONBLOCK | _CLOEXEC | _NOFOLLOW)
        except OSError as error:
            if error.errno == errno.ELOOP:
                _fail(f"{label} is unreadable or a symlink")
            raise RuntimeError(f"{label} could not be read") from None
        try:
            metadata = os.fstat(descriptor)
            if not stat.S_ISREG(metadata.st_mode):
                _fail(f"{label} must be a regular file and must not be a symlink")
            return cls(path, label, descriptor, _identity(metadata))
        except BaseException:
            os.close(descriptor)
            raise

    def rewind(self) -> None:
        """Move the pinned regular source to its first byte."""
        try:
            os.lseek(self.descriptor, 0, os.SEEK_SET)
        except OSError:
            raise RuntimeError(f"{self.label} could not be read") from None

    def duplicate_descriptor(self) -> int:
        """Return an owned duplicate positioned at the first byte."""
        self.rewind()
        try:
            return os.dup(self.descriptor)
        except OSError:
            raise RuntimeError(f"{self.label} could not be read") from None

    def digest(self) -> str:
        """Hash complete pinned bytes and reject mutation or short reads."""
        self.rewind()
        digest = hashlib.sha256()
        observed = 0
        try:
            while True:
                chunk = os.read(self.descriptor, _BUFFER_SIZE)
                if not chunk:
                    break
                digest.update(chunk)
                observed += len(chunk)
        except OSError:
            raise RuntimeError(f"{self.label} could not be read") from None
        self.verify(observed_size=observed)
        return digest.hexdigest()

    def copy_and_digest(self, target: BinaryIO) -> str:
        """Copy pinned bytes to a private snapshot while hashing that same stream."""
        self.rewind()
        digest = hashlib.sha256()
        observed = 0
        try:
            while True:
                chunk = os.read(self.descriptor, _BUFFER_SIZE)
                if not chunk:
                    break
                target.write(chunk)
                digest.update(chunk)
                observed += len(chunk)
            target.flush()
        except OSError:
            raise RuntimeError(f"{self.label} could not be read") from None
        self.verify(observed_size=observed)
        return digest.hexdigest()

    def verify(self, *, observed_size: int | None = None) -> None:
        """Require descriptor and no-follow path identities to remain unchanged."""
        try:
            descriptor_metadata = os.fstat(self.descriptor)
            path_metadata = os.stat(self.path, follow_symlinks=False)
        except OSError:
            raise RuntimeError(f"{self.label} changed during read fingerprinting") from None
        if (
            not stat.S_ISREG(path_metadata.st_mode)
            or _identity(descriptor_metadata) != self.identity
            or _identity(path_metadata) != self.identity
            or (observed_size is not None and observed_size != self.identity[2])
        ):
            raise RuntimeError(f"{self.label} changed during read fingerprinting")

    def close(self) -> None:
        """Close the pinned descriptor exactly once."""
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1

    def __enter__(self) -> PinnedReadSource:
        return self

    def __exit__(self, _type, _value, _traceback) -> None:
        self.close()
