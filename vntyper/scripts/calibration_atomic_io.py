"""No-clobber atomic directory installation for calibration artifacts."""

from __future__ import annotations

import ctypes
import errno
import logging
import os
import secrets
import shutil
import stat
import sys
from collections.abc import Callable
from pathlib import Path
from typing import NoReturn, Protocol, cast

logger = logging.getLogger(__name__)

_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_DIRECTORY = getattr(os, "O_DIRECTORY", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)
_RENAME_NOREPLACE = 1


class _RenameAt2(Protocol):
    argtypes: list[object]
    restype: object

    def __call__(self, old_dir: int, old_name: bytes, new_dir: int, new_name: bytes, flags: int) -> int: ...


def _load_renameat2() -> _RenameAt2 | None:
    if not sys.platform.startswith("linux"):
        return None
    try:
        function = ctypes.CDLL(None, use_errno=True).renameat2
    except (AttributeError, OSError):
        return None
    function.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    function.restype = ctypes.c_int
    return cast(_RenameAt2, function)


_renameat2 = _load_renameat2()


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _entry_exists(parent_descriptor: int, name: str) -> bool:
    try:
        os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return False
    return True


def _make_staging(output: Path, parent_descriptor: int) -> tuple[Path, tuple[int, int]]:
    for _ in range(100):
        name = f".{output.name}.{secrets.token_hex(8)}"
        try:
            os.mkdir(name, mode=0o700, dir_fd=parent_descriptor)
        except FileExistsError:
            continue
        staging = output.parent / name
        metadata = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
        return staging, (metadata.st_dev, metadata.st_ino)
    raise RuntimeError("calibration output could not allocate a unique staging directory")


def _is_owned_directory(path: Path, identity: tuple[int, int]) -> bool:
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        return False
    return stat.S_ISDIR(metadata.st_mode) and (metadata.st_dev, metadata.st_ino) == identity


def _remove_owned_staging(path: Path, identity: tuple[int, int]) -> None:
    if _is_owned_directory(path, identity):
        shutil.rmtree(path, ignore_errors=True)


def _activate(parent_descriptor: int, staging_name: str, output_name: str) -> None:
    if _renameat2 is None:
        raise RuntimeError("calibration atomic output requires Linux libc renameat2 with RENAME_NOREPLACE")
    ctypes.set_errno(0)
    result = _renameat2(
        parent_descriptor,
        os.fsencode(staging_name),
        parent_descriptor,
        os.fsencode(output_name),
        _RENAME_NOREPLACE,
    )
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
        _fail("calibration output already exists")
    if error_number in {errno.ENOSYS, errno.EINVAL, errno.ENOTSUP}:
        raise RuntimeError("calibration atomic output filesystem does not support renameat2 RENAME_NOREPLACE")
    raise OSError(error_number, os.strerror(error_number))


def atomic_output(output: Path, producer: Callable[[Path], bool]) -> bool:
    """Build and atomically install a new calibration artifact directory.

    Linux ``renameat2(RENAME_NOREPLACE)`` provides the no-clobber activation.
    The function fails closed on platforms or filesystems without that primitive.

    Args:
        output: New destination directory; any existing entry is refused.
        producer: Function that populates the private sibling staging directory
            and returns the completed scientific outcome status.

    Returns:
        The producer's strict boolean completed-outcome status.

    Raises:
        ValueError: If arguments, producer output, or destination are invalid.
        RuntimeError: If secure staging or no-clobber activation is unavailable.
        OSError: If staging or activation fails for another operating-system reason.
    """
    if not isinstance(output, Path) or output.name in {"", ".", ".."}:
        _fail("calibration output must be a non-root Path")
    if not callable(producer):
        _fail("calibration output producer must be callable")
    if not _NOFOLLOW:
        raise RuntimeError("calibration atomic output requires O_NOFOLLOW support")
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        parent_descriptor = os.open(output.parent, os.O_RDONLY | _DIRECTORY | _CLOEXEC | _NOFOLLOW)
    except OSError as error:
        raise ValueError("calibration output parent is unreadable or a symlink") from error
    staging: Path | None = None
    staging_identity: tuple[int, int] | None = None
    activated = False
    try:
        if not stat.S_ISDIR(os.fstat(parent_descriptor).st_mode):
            _fail("calibration output parent must be a directory")
        if _entry_exists(parent_descriptor, output.name):
            _fail("calibration output already exists")
        staging, staging_identity = _make_staging(output, parent_descriptor)
        successful = producer(staging)
        if not isinstance(successful, bool):
            _fail("calibration operation must return a completed-operation success value")
        if not _is_owned_directory(staging, staging_identity):
            raise RuntimeError("calibration staging directory changed during production")
        if not any(staging.iterdir()):
            _fail("calibration operation produced no artifacts")
        _activate(parent_descriptor, staging.name, output.name)
        activated = True
        return successful
    finally:
        os.close(parent_descriptor)
        if not activated and staging is not None and staging_identity is not None:
            _remove_owned_staging(staging, staging_identity)
