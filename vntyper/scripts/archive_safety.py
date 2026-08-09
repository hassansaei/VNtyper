"""Descriptor-anchored, fail-closed archives for VNtyper result trees."""

from __future__ import annotations

import logging
import os
import secrets
import shutil
import stat
import tarfile
import tempfile
import zipfile
from collections.abc import Callable, Iterable
from pathlib import Path

logger = logging.getLogger(__name__)

_ARCHIVE_SUFFIXES = {"zip": ".zip", "gztar": ".tar.gz"}
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
_FILE_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
_DirectoryWriter = Callable[[str, os.stat_result], None]
_FileWriter = Callable[[str, int, os.stat_result], None]


def _reject(message: str) -> None:
    logger.error(message)
    raise ValueError(message)


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _same_identity(left: os.stat_result, right: os.stat_result) -> bool:
    return (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)


def _archive_path(base_name: str | Path, archive_format: str) -> Path:
    try:
        suffix = _ARCHIVE_SUFFIXES[archive_format]
    except KeyError:
        _reject(f"Unsupported archive format: {archive_format}")
    return Path(f"{base_name}{suffix}")


def _validate_destination(
    destination: Path,
    root: Path | None,
    protected_paths: Iterable[str | Path],
    *,
    allow_existing_symlink: bool,
) -> None:
    destination_absolute = _absolute(destination)
    destination_resolved = destination_absolute.resolve(strict=False)
    if root is not None:
        root_absolute = _absolute(root)
        root_resolved = root_absolute.resolve(strict=False)
        if _is_within(destination_absolute, root_absolute) or _is_within(destination_resolved, root_resolved):
            _reject(f"Unsafe archive destination is inside its source tree: {destination}")

    for protected_path in protected_paths:
        protected = _absolute(protected_path)
        if (
            destination_absolute == protected
            or destination_resolved == protected.resolve(strict=False)
            or _same_file(destination, protected)
        ):
            _reject(f"Unsafe archive destination aliases protected input: {destination}")

    metadata = os.lstat(destination) if os.path.lexists(destination) else None
    if metadata is not None and stat.S_ISLNK(metadata.st_mode) and not allow_existing_symlink:
        _reject(f"Unsafe archive destination is a symbolic link: {destination}")
    if metadata is not None and not stat.S_ISREG(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode):
        _reject(f"Unsafe archive destination is not a regular file: {destination}")
    if metadata is not None and stat.S_ISREG(metadata.st_mode) and metadata.st_nlink > 1:
        _reject(f"Unsafe archive destination has multiple hard links: {destination}")


def clear_stale_archive(
    base_name: str | Path,
    archive_format: str,
    *,
    protected_paths: Iterable[str | Path] = (),
) -> None:
    """Remove one stale public archive entry without following a symlink.

    This is used by the web worker before launching a retry, so a failed
    subprocess cannot leave an older result downloadable as if it were new.
    Regular files with multiple links and lexical operator-input collisions are
    rejected; a symlink entry itself is safe to unlink because its target is
    never opened or removed.

    Args:
        base_name: Public archive path without its suffix.
        archive_format: ``zip`` or ``gztar``.
        protected_paths: Operator-owned paths which must not be unlinked.

    Raises:
        ValueError: If the destination is an unsafe regular-file alias.
        OSError: If the entry cannot be inspected or removed.
    """
    destination = _archive_path(base_name, archive_format)
    _validate_destination(destination, None, protected_paths, allow_existing_symlink=True)
    if not os.path.lexists(destination):
        return

    parent_descriptor = os.open(destination.parent, _DIRECTORY_FLAGS)
    quarantine_name = f".{destination.name}.stale-{secrets.token_hex(8)}"
    try:
        previous_metadata = os.stat(destination.name, dir_fd=parent_descriptor, follow_symlinks=False)
        os.rename(
            destination.name,
            quarantine_name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        quarantined_metadata = os.stat(quarantine_name, dir_fd=parent_descriptor, follow_symlinks=False)
        if not _same_identity(previous_metadata, quarantined_metadata):
            _reject(f"Archive destination changed during stale cleanup; retained quarantine '{quarantine_name}'.")
        os.unlink(quarantine_name, dir_fd=parent_descriptor)
    finally:
        os.close(parent_descriptor)


def _archive_directory(
    directory_descriptor: int,
    relative_directory: str,
    write_directory: _DirectoryWriter,
    write_file: _FileWriter,
) -> None:
    for entry_name in sorted(os.listdir(directory_descriptor)):
        relative_entry = f"{relative_directory}/{entry_name}" if relative_directory else entry_name
        metadata = os.stat(entry_name, dir_fd=directory_descriptor, follow_symlinks=False)
        if stat.S_ISLNK(metadata.st_mode):
            _reject(f"Refusing to archive unsafe symbolic link '{relative_entry}'.")
        if stat.S_ISDIR(metadata.st_mode):
            child_descriptor = os.open(entry_name, _DIRECTORY_FLAGS, dir_fd=directory_descriptor)
            try:
                opened_metadata = os.fstat(child_descriptor)
                if not stat.S_ISDIR(opened_metadata.st_mode) or not _same_identity(metadata, opened_metadata):
                    _reject(f"Refusing to archive unsupported filesystem entry '{relative_entry}'.")
                write_directory(relative_entry, opened_metadata)
                _archive_directory(child_descriptor, relative_entry, write_directory, write_file)
            finally:
                os.close(child_descriptor)
            continue
        if not stat.S_ISREG(metadata.st_mode):
            _reject(f"Refusing to archive unsupported filesystem entry '{relative_entry}'.")

        file_descriptor = os.open(entry_name, _FILE_FLAGS, dir_fd=directory_descriptor)
        try:
            opened_metadata = os.fstat(file_descriptor)
            if not stat.S_ISREG(opened_metadata.st_mode) or not _same_identity(metadata, opened_metadata):
                _reject(f"Refusing to archive unsupported filesystem entry '{relative_entry}'.")
            if opened_metadata.st_nlink != 1:
                _reject(f"Refusing to archive unsafe hard-linked file '{relative_entry}'.")
            write_file(relative_entry, file_descriptor, opened_metadata)
            if os.fstat(file_descriptor).st_nlink != 1:
                _reject(f"Refusing to archive unsafe hard-linked file '{relative_entry}'.")
        finally:
            os.close(file_descriptor)


def _write_zip(temporary_archive: Path, root_descriptor: int) -> None:
    with zipfile.ZipFile(temporary_archive, "x", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:

        def write_directory(relative: str, metadata: os.stat_result) -> None:
            info = zipfile.ZipInfo(f"{relative}/")
            info.external_attr = (metadata.st_mode & 0xFFFF) << 16
            archive.writestr(info, b"")

        def write_file(relative: str, descriptor: int, metadata: os.stat_result) -> None:
            info = zipfile.ZipInfo(relative)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (metadata.st_mode & 0xFFFF) << 16
            with os.fdopen(os.dup(descriptor), "rb") as source, archive.open(info, "w", force_zip64=True) as target:
                shutil.copyfileobj(source, target)

        _archive_directory(root_descriptor, "", write_directory, write_file)


def _write_tar(temporary_archive: Path, root_descriptor: int) -> None:
    with tarfile.open(temporary_archive, "x:gz", dereference=False) as archive:

        def populate_info(relative: str, metadata: os.stat_result) -> tarfile.TarInfo:
            info = tarfile.TarInfo(f"./{relative}")
            info.mode = stat.S_IMODE(metadata.st_mode)
            info.mtime = metadata.st_mtime
            info.uid = metadata.st_uid
            info.gid = metadata.st_gid
            return info

        def write_directory(relative: str, metadata: os.stat_result) -> None:
            info = populate_info(relative, metadata)
            info.type = tarfile.DIRTYPE
            archive.addfile(info)

        def write_file(relative: str, descriptor: int, metadata: os.stat_result) -> None:
            info = populate_info(relative, metadata)
            info.size = metadata.st_size
            with os.fdopen(os.dup(descriptor), "rb") as source:
                archive.addfile(info, source)

        _archive_directory(root_descriptor, "", write_directory, write_file)


def create_safe_archive(
    base_name: str | Path,
    archive_format: str,
    root_dir: str | Path,
    *,
    protected_paths: Iterable[str | Path] = (),
) -> str:
    """Atomically archive exclusively owned files through anchored descriptors.

    The root, every directory, and every file is opened with ``O_NOFOLLOW``.
    Source bytes are streamed from those already-open file descriptors, so a
    path replacement cannot redirect an archive read. Hard links, symlinks,
    and special files fail the whole attempt. The stale public destination is
    removed before work begins and only a complete temporary archive is
    installed.

    Args:
        base_name: Destination path without its format suffix.
        archive_format: ``zip`` or ``gztar``.
        root_dir: Result directory whose contents should be archived.
        protected_paths: Operator-owned inputs which the destination may not alias.

    Returns:
        The installed archive path.

    Raises:
        ValueError: If the format, destination, root, or a source entry is unsafe.
        OSError: If the filesystem changes concurrently or archive I/O fails.
    """
    archive_path = _archive_path(base_name, archive_format)
    root = Path(root_dir)
    root_metadata = os.lstat(root)
    if stat.S_ISLNK(root_metadata.st_mode):
        _reject("Refusing to archive unsafe symbolic link used as result root.")
    if not stat.S_ISDIR(root_metadata.st_mode):
        _reject(f"Archive root is not a directory: {root}")
    _validate_destination(archive_path, root, protected_paths, allow_existing_symlink=False)
    clear_stale_archive(base_name, archive_format, protected_paths=protected_paths)

    archive_base = Path(base_name)
    temporary_dir = Path(tempfile.mkdtemp(prefix=f".{archive_base.name}.archive-", dir=str(archive_base.parent)))
    temporary_archive = temporary_dir / f"payload{_ARCHIVE_SUFFIXES[archive_format]}"
    root_descriptor = -1
    primary_failure: BaseException | None = None
    try:
        root_descriptor = os.open(root, _DIRECTORY_FLAGS)
        opened_root = os.fstat(root_descriptor)
        if not stat.S_ISDIR(opened_root.st_mode) or not _same_identity(root_metadata, opened_root):
            _reject(f"Archive root is not a directory: {root}")
        if archive_format == "zip":
            _write_zip(temporary_archive, root_descriptor)
        else:
            _write_tar(temporary_archive, root_descriptor)
        os.replace(temporary_archive, archive_path)
    except BaseException as error:
        primary_failure = error
        raise
    finally:
        if root_descriptor >= 0:
            os.close(root_descriptor)
        try:
            shutil.rmtree(temporary_dir)
        except Exception as cleanup_error:
            if primary_failure is None:
                raise
            logger.error(
                f"Archive attempt failed and temporary cleanup also failed; retained partial at "
                f"{temporary_dir}: {cleanup_error}"
            )

    return str(archive_path)
