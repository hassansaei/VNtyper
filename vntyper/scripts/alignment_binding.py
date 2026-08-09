"""Lifetime-owned descriptor bindings for run-local alignment views."""

from __future__ import annotations

import logging
import os
import stat
from contextlib import suppress

logger = logging.getLogger(__name__)


class AlignmentBinding:
    """Hold one alignment inode open for the lifetime of its pipeline consumers."""

    def __init__(self, input_path: str) -> None:
        """Open and validate the descriptor exposed through this process's procfs.

        Args:
            input_path: Alignment pathname to bind. A symlink is followed once here.

        Raises:
            RuntimeError: If the file or its stable procfs descriptor path is unavailable.
        """
        self._descriptor: int | None = None
        self.input_path = input_path
        try:
            descriptor = os.open(input_path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
            self._descriptor = descriptor
            descriptor_stat = os.fstat(descriptor)
            if not stat.S_ISREG(descriptor_stat.st_mode):
                raise OSError("alignment is not a regular file")
            self.view_target = f"/proc/{os.getpid()}/fd/{descriptor}"
            target_stat = os.stat(self.view_target)
            if (target_stat.st_dev, target_stat.st_ino) != (descriptor_stat.st_dev, descriptor_stat.st_ino):
                raise OSError("procfs descriptor path does not resolve to the opened alignment")
        except OSError as error:
            self.close()
            message = f"Unable to create a stable run-local alignment binding for {input_path}: {error}"
            logger.error(message)
            raise RuntimeError(message) from error

    @property
    def is_open(self) -> bool:
        """Return whether the bound descriptor is still available to consumers."""
        return self._descriptor is not None

    def close(self) -> None:
        """Release the descriptor exactly once."""
        descriptor = self._descriptor
        if descriptor is None:
            return
        self._descriptor = None
        with suppress(OSError):
            os.close(descriptor)

    def __del__(self) -> None:
        """Avoid leaking a descriptor when a direct caller drops an unclosed plan."""
        self.close()
