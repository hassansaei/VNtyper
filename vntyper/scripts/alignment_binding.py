"""Lifetime-owned descriptor bindings for run-local alignment views."""

from __future__ import annotations

import logging
import os
import secrets
import stat
from contextlib import suppress
from pathlib import Path

logger = logging.getLogger(__name__)


class AlignmentBinding:
    """Hold one alignment inode and its exact run-local directory entry."""

    def __init__(self, input_path: str) -> None:
        """Open the supplied alignment without requiring a particular descriptor filesystem.

        Args:
            input_path: Alignment pathname to bind. A symlink is followed once here.

        Raises:
            RuntimeError: If the alignment cannot be opened as a regular file.
        """
        self._descriptor: int | None = None
        self._descriptor_identity: tuple[int, int] | None = None
        self._view_path: Path | None = None
        self._view_identity: tuple[int, int] | None = None
        self._view_kind: str | None = None
        self.input_path = input_path
        try:
            descriptor = os.open(input_path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
            self._descriptor = descriptor
            descriptor_stat = os.fstat(descriptor)
            if not stat.S_ISREG(descriptor_stat.st_mode):
                raise OSError("alignment is not a regular file")
            self._descriptor_identity = (descriptor_stat.st_dev, descriptor_stat.st_ino)
            self.view_target = f"/proc/{os.getpid()}/fd/{descriptor}"
        except OSError as error:
            self._close_descriptor()
            message = f"Unable to create a stable run-local alignment binding for {input_path}: {error}"
            logger.error(message)
            raise RuntimeError(message) from error

    @property
    def is_open(self) -> bool:
        """Return whether the bound descriptor is still available to consumers."""
        return self._descriptor is not None

    def _temporary_view_path(self, destination: Path) -> Path:
        for _ in range(100):
            temporary = destination.parent / f".{destination.name}.{secrets.token_hex(8)}.tmp"
            if not os.path.lexists(temporary):
                return temporary
        message = f"Unable to allocate a temporary alignment view beside {destination}"
        logger.error(message)
        raise RuntimeError(message)

    def _record_view(self, destination: Path, installed_stat: os.stat_result, kind: str) -> None:
        self._view_path = destination
        self._view_identity = (installed_stat.st_dev, installed_stat.st_ino)
        self._view_kind = kind

    def _install_proc_view(self, destination: Path) -> bool:
        try:
            target_stat = os.stat(self.view_target)
        except OSError:
            return False
        if (target_stat.st_dev, target_stat.st_ino) != self._descriptor_identity:
            return False
        temporary = self._temporary_view_path(destination)
        try:
            os.symlink(self.view_target, temporary)
            installed_stat = os.lstat(temporary)
            os.replace(temporary, destination)
        except OSError:
            with suppress(OSError):
                os.unlink(temporary)
            return False
        self._record_view(destination, installed_stat, "symlink")
        return True

    def _install_hardlink_view(self, destination: Path) -> None:
        temporary = self._temporary_view_path(destination)
        try:
            source_path = os.path.realpath(self.input_path)
            os.link(source_path, temporary, follow_symlinks=True)
            installed_stat = os.stat(temporary, follow_symlinks=False)
            if (
                not stat.S_ISREG(installed_stat.st_mode)
                or (
                    installed_stat.st_dev,
                    installed_stat.st_ino,
                )
                != self._descriptor_identity
            ):
                raise OSError("hardlink does not identify the already-open alignment")
            os.replace(temporary, destination)
        except OSError as error:
            with suppress(OSError):
                os.unlink(temporary)
            message = (
                f"Unable to create a stable run-local alignment binding for {self.input_path}: "
                "procfs descriptor binding is unavailable and a same-filesystem hardlink could not be created: "
                f"{error}"
            )
            logger.error(message)
            raise RuntimeError(message) from error
        self._record_view(destination, installed_stat, "hardlink")

    def install_view(self, view_path: str | Path) -> None:
        """Atomically install and own a stable run-local view of the opened inode.

        Args:
            view_path: Validated pipeline-owned pathname for the alignment view.

        Raises:
            RuntimeError: If neither a proc descriptor link nor a verified hardlink can be installed.
        """
        if self._view_path is not None:
            message = f"Alignment binding already owns an alignment view: {self._view_path}"
            logger.error(message)
            raise RuntimeError(message)
        if self._descriptor is None:
            raise RuntimeError("Cannot install an alignment view from a closed descriptor binding.")
        destination = Path(view_path)
        if self._install_proc_view(destination):
            return
        self._install_hardlink_view(destination)

    def _remove_owned_view(self) -> None:
        if self._view_path is None or self._view_identity is None or self._view_kind is None:
            return
        try:
            current_stat = os.lstat(self._view_path)
        except FileNotFoundError:
            self._view_path = None
            self._view_identity = None
            self._view_kind = None
            return
        except OSError as error:
            message = f"Unable to inspect owned alignment view {self._view_path} before descriptor release: {error}"
            logger.error(message)
            raise RuntimeError(message) from error
        current_identity = (current_stat.st_dev, current_stat.st_ino)
        expected_type = (
            stat.S_ISLNK(current_stat.st_mode) if self._view_kind == "symlink" else stat.S_ISREG(current_stat.st_mode)
        )
        target_matches = (
            self._view_kind != "symlink" or expected_type and os.readlink(self._view_path) == self.view_target
        )
        if current_identity != self._view_identity or not expected_type or not target_matches:
            message = (
                f"Refusing to remove alignment view because the owned alignment view was replaced: {self._view_path}"
            )
            logger.error(message)
            raise RuntimeError(message)
        try:
            os.unlink(self._view_path)
        except OSError as error:
            message = f"Unable to remove owned alignment view {self._view_path} before descriptor release: {error}"
            logger.error(message)
            raise RuntimeError(message) from error
        self._view_path = None
        self._view_identity = None
        self._view_kind = None

    def _close_descriptor(self) -> None:
        descriptor = self._descriptor
        if descriptor is None:
            return
        self._descriptor = None
        with suppress(OSError):
            os.close(descriptor)

    def close(self) -> None:
        """Remove the exact owned view, then release the descriptor exactly once."""
        self._remove_owned_view()
        self._close_descriptor()

    def __del__(self) -> None:
        """Finalize safe ownership, preserving the FD when cleanup refusal makes reuse unsafe."""
        try:
            self.close()
        except Exception as error:
            with suppress(Exception):
                logger.error(f"Preserving alignment descriptor because safe owned-view cleanup was refused: {error}")
