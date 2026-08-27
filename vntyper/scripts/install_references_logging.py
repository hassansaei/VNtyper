"""Additive logging lifecycle for reference installation.

The CLI owns the root logger's configuration.  This module only adds the installer's
transient file and, for direct module execution, the console handler that would
otherwise be absent.  Keeping this lifecycle outside ``install_references.py`` also
keeps inode-rotation decisions out of that already oversized I/O orchestrator.
"""

from __future__ import annotations

import logging
import shutil
import sys
from collections.abc import Callable
from pathlib import Path

logger = logging.getLogger(__name__)


class InstallLogHandler(logging.FileHandler):
    """The installer-owned file handler plus the state needed for safe finalisation."""

    destination_handlers: tuple[logging.FileHandler, ...]
    standalone_console_handler: logging.StreamHandler | None


def _normalise_log_path(path: str | Path) -> Path:
    """Return a stable absolute spelling for a possibly non-existent log path.

    Args:
        path: File path to normalise.

    Returns:
        Path: Absolute, symlink-resolved path without requiring it to exist.
    """
    return Path(path).expanduser().resolve(strict=False)


def file_handlers_for_path(handlers: list[logging.Handler], destination: Path) -> tuple[logging.FileHandler, ...]:
    """Select existing file handlers whose configured destination is ``destination``.

    Args:
        handlers: Root handlers present before the installer attaches its own.
        destination: Final installer-log path.

    Returns:
        tuple[logging.FileHandler, ...]: Matching handlers, preserving root order.
    """
    normalised_destination = _normalise_log_path(destination)
    return tuple(
        handler
        for handler in handlers
        if isinstance(handler, logging.FileHandler)
        and _normalise_log_path(handler.baseFilename) == normalised_destination
    )


def attach_install_log(output_dir: Path, log_file: Path) -> InstallLogHandler:
    """Attach one installer file handler without replacing CLI-owned handlers.

    When a CLI file handler already targets the final installer log, its pre-install
    records are copied into the transient file before the installer handler opens it.
    Later records are captured by the transient handler, so finalisation can replace
    the output directory without losing either half of the run.

    Args:
        output_dir: Final reference directory.
        log_file: Transient path that survives staged activation.

    Returns:
        InstallLogHandler: The attached installer-owned handler.
    """
    root = logging.getLogger()
    final_log = output_dir / "install_references.log"
    destination_handlers = file_handlers_for_path(root.handlers[:], final_log)
    for handler in destination_handlers:
        handler.flush()

    if destination_handlers and final_log.exists() and _normalise_log_path(log_file) != _normalise_log_path(final_log):
        # The final rename replaces the CLI-created destination with this transient
        # inode.  Preserve the operator's permission bits (and other copystat metadata)
        # rather than recreating them through the process umask.  `copy2` deliberately
        # does not copy owner/group: the transient remains owned by the process which
        # is performing this install, matching ordinary file-creation semantics.
        shutil.copy2(final_log, log_file)

    # Open the required file before mutating the standalone root configuration.  A
    # permissions/disk failure must leave a caller with exactly the handlers and level
    # it had on entry, not a half-configured console-only logger.
    file_handler = InstallLogHandler(log_file)
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))

    standalone_console_handler: logging.StreamHandler | None = None
    if not root.handlers:
        root.setLevel(logging.INFO)
        standalone_console_handler = logging.StreamHandler(sys.stdout)
        standalone_console_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        root.addHandler(standalone_console_handler)

    file_handler.destination_handlers = destination_handlers
    file_handler.standalone_console_handler = standalone_console_handler
    root.addHandler(file_handler)
    return file_handler


def _suspend_file_handler(handler: logging.FileHandler) -> bool:
    """Close a handler's current inode while retaining its configuration and identity.

    Args:
        handler: CLI-owned handler that targets the final installer log.

    Returns:
        bool: Whether an open stream must be restored after finalisation.
    """
    handler.acquire()
    try:
        if handler.stream is None:
            return False
        handler.flush()
        handler.stream.close()
        handler.stream = None
        return True
    finally:
        handler.release()


def _restore_file_handler(handler: logging.FileHandler, was_open: bool) -> None:
    """Reopen a suspended CLI handler on its configured, now-final path.

    Args:
        handler: CLI-owned handler retained on the root logger.
        was_open: Result returned by :func:`_suspend_file_handler`.
    """
    if not was_open:
        return
    handler.acquire()
    try:
        handler.stream = handler._open()  # noqa: SLF001 - stdlib FileHandler has no public reopen operation
    finally:
        handler.release()


def finish_install_log(file_handler: InstallLogHandler, finalize: Callable[[], None]) -> None:
    """Detach/close the installer handler, finalise, then restore CLI destinations.

    Nested ``finally`` blocks are intentional: close, rename, reopen and standalone
    cleanup must all be attempted even if an earlier cleanup operation fails.  The
    caller decides whether a resulting cleanup exception is authoritative or secondary
    to an install failure already unwinding.

    Args:
        file_handler: Handler returned by :func:`attach_install_log`.
        finalize: Move the completed transient log into the settled output directory.
    """
    root = logging.getLogger()
    suspended: list[tuple[logging.FileHandler, bool]] = []
    try:
        root.removeHandler(file_handler)
    finally:
        try:
            file_handler.close()
        finally:
            try:
                suspended.extend(
                    (handler, _suspend_file_handler(handler)) for handler in file_handler.destination_handlers
                )
            finally:
                try:
                    finalize()
                finally:
                    try:
                        for handler, was_open in suspended:
                            _restore_file_handler(handler, was_open)
                    finally:
                        console = file_handler.standalone_console_handler
                        if console is not None:
                            root.removeHandler(console)
                            console.close()
