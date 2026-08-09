"""Lifetime-owned descriptor binding for one generated alignment index."""

from __future__ import annotations

import logging
import os
import stat
from contextlib import suppress
from pathlib import Path

logger = logging.getLogger(__name__)


class AlignmentIndexBinding:
    """Keep the exact generated index inode available to child processes."""

    def __init__(self, index_path: str | Path, fallback_path: str | Path) -> None:
        """Open an index and select a descriptor or verified-hardlink view.

        Args:
            index_path: Fresh temporary index to bind before publication.
            fallback_path: Reserved same-filesystem path used without procfs.

        Raises:
            RuntimeError: If the index is not regular or neither view can be retained.
        """
        self._descriptor: int | None = None
        self._identity: tuple[int, int] | None = None
        self._fallback_path: Path | None = None
        self._fallback_identity: tuple[int, int] | None = None
        self.consumer_path: str | None = None
        try:
            descriptor = os.open(index_path, os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0))
            self._descriptor = descriptor
            descriptor_stat = os.fstat(descriptor)
            if not stat.S_ISREG(descriptor_stat.st_mode):
                raise OSError("generated index is not a regular file")
            self._identity = (descriptor_stat.st_dev, descriptor_stat.st_ino)
            proc_path = f"/proc/{os.getpid()}/fd/{descriptor}"
            try:
                proc_stat = os.stat(proc_path)
            except OSError:
                self._install_hardlink(index_path, Path(fallback_path))
            else:
                if (proc_stat.st_dev, proc_stat.st_ino) != self._identity:
                    self._install_hardlink(index_path, Path(fallback_path))
                else:
                    self.consumer_path = proc_path
        except OSError as error:
            self._close_descriptor()
            message = f"Unable to retain the fresh generated index {index_path}: {error}"
            logger.error(message)
            raise RuntimeError(message) from error

    @property
    def is_open(self) -> bool:
        """Return whether the generated-index descriptor is retained."""
        return self._descriptor is not None

    def _install_hardlink(self, index_path: str | Path, fallback_path: Path) -> None:
        if os.path.lexists(fallback_path):
            raise OSError(f"reserved generated-index binding already exists: {fallback_path}")
        created_identity: tuple[int, int] | None = None
        try:
            os.link(index_path, fallback_path, follow_symlinks=False)
            linked_stat = os.lstat(fallback_path)
            linked_identity = (linked_stat.st_dev, linked_stat.st_ino)
            created_identity = linked_identity
            if not stat.S_ISREG(linked_stat.st_mode) or linked_identity != self._identity:
                raise OSError("hardlink does not identify the already-open generated index")
        except OSError:
            if created_identity is not None:
                with suppress(OSError):
                    current = os.lstat(fallback_path)
                    if (current.st_dev, current.st_ino) == created_identity:
                        os.unlink(fallback_path)
            raise
        self._fallback_path = fallback_path
        self._fallback_identity = linked_identity
        self.consumer_path = str(fallback_path)

    def _remove_fallback(self) -> None:
        if self._fallback_path is None or self._fallback_identity is None:
            return
        try:
            current_stat = os.stat(self._fallback_path, follow_symlinks=False)
        except FileNotFoundError:
            self._fallback_path = None
            self._fallback_identity = None
            return
        except OSError as error:
            message = f"Unable to inspect generated-index binding {self._fallback_path} before release: {error}"
            logger.error(message)
            raise RuntimeError(message) from error
        if (
            not stat.S_ISREG(current_stat.st_mode)
            or (current_stat.st_dev, current_stat.st_ino) != self._fallback_identity
        ):
            message = f"Refusing to remove generated-index binding because it was replaced: {self._fallback_path}"
            logger.error(message)
            raise RuntimeError(message)
        try:
            os.unlink(self._fallback_path)
        except OSError as error:
            message = f"Unable to remove generated-index binding {self._fallback_path} before release: {error}"
            logger.error(message)
            raise RuntimeError(message) from error
        self._fallback_path = None
        self._fallback_identity = None

    def _close_descriptor(self) -> None:
        descriptor = self._descriptor
        if descriptor is None:
            return
        self._descriptor = None
        with suppress(OSError):
            os.close(descriptor)

    def close(self) -> None:
        """Remove an exact fallback entry, then release the index descriptor."""
        self._remove_fallback()
        self._close_descriptor()

    def __del__(self) -> None:
        """Preserve the descriptor if safe fallback cleanup is refused."""
        try:
            self.close()
        except Exception as error:
            with suppress(Exception):
                logger.error(f"Preserving generated-index descriptor because cleanup was refused: {error}")
