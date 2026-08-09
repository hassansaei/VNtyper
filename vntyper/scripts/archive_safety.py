"""Descriptor-anchored, fail-closed archives for VNtyper result trees."""

from __future__ import annotations

import logging
import os
import secrets
import shutil
import stat
import tarfile
import zipfile
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import BinaryIO

logger = logging.getLogger(__name__)

_ARCHIVE_SUFFIXES = {"zip": ".zip", "gztar": ".tar.gz"}
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
_FILE_FLAGS = os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
_DirectoryWriter = Callable[[str, os.stat_result], None]
_FileWriter = Callable[[str, int, os.stat_result], None]


def _reject(message: str) -> None:
    logger.error(message)
    raise ValueError(message)


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


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
    parent_descriptor: int,
) -> os.stat_result | None:
    destination_absolute = _absolute(destination)
    destination_resolved = destination_absolute.resolve(strict=False)
    if root is not None:
        root_absolute = _absolute(root)
        root_resolved = root_absolute.resolve(strict=False)
        if _is_within(destination_absolute, root_absolute) or _is_within(destination_resolved, root_resolved):
            _reject(f"Unsafe archive destination is inside its source tree: {destination}")

    try:
        metadata = os.stat(destination.name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        metadata = None

    for protected_path in protected_paths:
        protected = _absolute(protected_path)
        same_inode = False
        if metadata is not None:
            try:
                protected_metadata = os.stat(protected)
            except OSError:
                pass
            else:
                same_inode = _same_identity(metadata, protected_metadata)
        if (
            destination_absolute == protected
            or (
                destination_resolved == protected.resolve(strict=False)
                and not (allow_existing_symlink and metadata is not None and stat.S_ISLNK(metadata.st_mode))
            )
            or same_inode
        ):
            _reject(f"Unsafe archive destination aliases protected input: {destination}")

    if metadata is not None and stat.S_ISLNK(metadata.st_mode) and not allow_existing_symlink:
        _reject(f"Unsafe archive destination is a symbolic link: {destination}")
    if metadata is not None and not stat.S_ISREG(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode):
        _reject(f"Unsafe archive destination is not a regular file: {destination}")
    if metadata is not None and stat.S_ISREG(metadata.st_mode) and metadata.st_nlink > 1:
        _reject(f"Unsafe archive destination has multiple hard links: {destination}")

    return metadata


def _require_destination_identity(
    parent_descriptor: int,
    destination_name: str,
    expected: os.stat_result | None,
    phase: str,
) -> os.stat_result | None:
    """Require a public filename to retain its previously observed identity."""
    try:
        current = os.stat(destination_name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        current = None
    if expected is None and current is None:
        return None
    if expected is None or current is None or not _same_identity(expected, current):
        _reject(f"Archive destination changed {phase}: {destination_name}")
    return current


def _unlink_if_same(parent_descriptor: int, destination_name: str, expected: os.stat_result) -> bool:
    """Unlink a public name only while it still denotes the expected inode."""
    try:
        current = os.stat(destination_name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return False
    if not _same_identity(expected, current):
        return False
    os.unlink(destination_name, dir_fd=parent_descriptor)
    return True


def _clear_stale_at(
    parent_descriptor: int,
    destination_name: str,
    expected_metadata: os.stat_result | None,
) -> None:
    previous_metadata = _require_destination_identity(
        parent_descriptor,
        destination_name,
        expected_metadata,
        "after validation",
    )
    if previous_metadata is None:
        return
    quarantine_name = f".{destination_name}.stale-{secrets.token_hex(8)}"
    os.rename(
        destination_name,
        quarantine_name,
        src_dir_fd=parent_descriptor,
        dst_dir_fd=parent_descriptor,
    )
    quarantined_metadata = os.stat(quarantine_name, dir_fd=parent_descriptor, follow_symlinks=False)
    if not _same_identity(previous_metadata, quarantined_metadata):
        _reject(f"Archive destination changed during stale cleanup; retained quarantine '{quarantine_name}'.")
    os.unlink(quarantine_name, dir_fd=parent_descriptor)


def _quarantine_at(parent_descriptor: int, destination_name: str) -> str | None:
    """Move one regular public archive to an unserved name in its opened parent."""
    try:
        previous_metadata = os.stat(destination_name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(previous_metadata.st_mode) or previous_metadata.st_nlink != 1:
        _reject(f"Unsafe public archive cannot be quarantined: {destination_name}")
    quarantine_name = f".{destination_name}.failed-{secrets.token_hex(8)}"
    os.rename(
        destination_name,
        quarantine_name,
        src_dir_fd=parent_descriptor,
        dst_dir_fd=parent_descriptor,
    )
    quarantined_metadata = os.stat(quarantine_name, dir_fd=parent_descriptor, follow_symlinks=False)
    if not _same_identity(previous_metadata, quarantined_metadata):
        _reject(f"Public archive changed during quarantine; retained '{quarantine_name}'.")
    return quarantine_name


def _open_parent(destination: Path) -> tuple[int, os.stat_result]:
    expected_metadata = os.stat(destination.parent)
    descriptor = os.open(destination.parent, _DIRECTORY_FLAGS)
    opened_metadata = os.fstat(descriptor)
    if not _same_identity(expected_metadata, opened_metadata):
        os.close(descriptor)
        _reject(f"Unsafe archive parent changed while opening: {destination.parent}")
    return descriptor, opened_metadata


def _require_current_parent(destination: Path, expected_metadata: os.stat_result) -> None:
    try:
        current_metadata = os.stat(destination.parent)
    except OSError as error:
        raise OSError(f"Unsafe archive parent changed before install: {destination.parent}") from error
    if not _same_identity(expected_metadata, current_metadata):
        _reject(f"Unsafe archive parent changed before install: {destination.parent}")


def clear_stale_archive(
    base_name: str | Path,
    archive_format: str,
    *,
    protected_paths: Iterable[str | Path] = (),
) -> None:
    """Remove one stale public archive through a validated parent descriptor.

    Args:
        base_name: Public archive path without its suffix.
        archive_format: ``zip`` or ``gztar``.
        protected_paths: Operator-owned paths which must not be unlinked.

    Raises:
        ValueError: If the destination aliases protected state.
        OSError: If the parent or destination cannot be safely inspected.
    """
    destination = _archive_path(base_name, archive_format)
    parent_descriptor, parent_metadata = _open_parent(destination)
    primary_failure: BaseException | None = None
    try:
        destination_metadata = _validate_destination(
            destination,
            None,
            protected_paths,
            allow_existing_symlink=True,
            parent_descriptor=parent_descriptor,
        )
        _require_current_parent(destination, parent_metadata)
        _clear_stale_at(parent_descriptor, destination.name, destination_metadata)
        _require_current_parent(destination, parent_metadata)
    except BaseException as error:
        primary_failure = error
        raise
    finally:
        try:
            os.close(parent_descriptor)
        except Exception as cleanup_error:
            if primary_failure is None:
                raise
            logger.error(f"Stale archive cleanup failed and parent descriptor cleanup also failed: {cleanup_error}")


def quarantine_archive(
    base_name: str | Path,
    archive_format: str,
    *,
    protected_paths: Iterable[str | Path] = (),
) -> str | None:
    """Atomically hide a complete public archive while preserving its bytes.

    Args:
        base_name: Public archive path without its suffix.
        archive_format: ``zip`` or ``gztar``.
        protected_paths: Operator-owned paths which must not be renamed.

    Returns:
        The private quarantine path, or ``None`` when no public archive exists.

    Raises:
        ValueError: If the destination is not an exclusively owned regular file.
        OSError: If the parent changes or the rename cannot be completed safely.
    """
    destination = _archive_path(base_name, archive_format)
    parent_descriptor, parent_metadata = _open_parent(destination)
    primary_failure: BaseException | None = None
    try:
        _validate_destination(
            destination,
            None,
            protected_paths,
            allow_existing_symlink=False,
            parent_descriptor=parent_descriptor,
        )
        _require_current_parent(destination, parent_metadata)
        quarantine_name = _quarantine_at(parent_descriptor, destination.name)
        _require_current_parent(destination, parent_metadata)
        if quarantine_name is None:
            return None
        return str(destination.parent / quarantine_name)
    except BaseException as error:
        primary_failure = error
        raise
    finally:
        try:
            os.close(parent_descriptor)
        except Exception as cleanup_error:
            if primary_failure is None:
                raise
            logger.error(f"Archive quarantine failed and parent descriptor cleanup also failed: {cleanup_error}")


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


def _write_zip(temporary_archive: BinaryIO, root_descriptor: int) -> None:
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


def _write_tar(temporary_archive: BinaryIO, root_descriptor: int) -> None:
    with tarfile.open(fileobj=temporary_archive, mode="w:gz", dereference=False) as archive:

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
    parent_descriptor = -1
    temporary_name: str | None = None
    completed_metadata: os.stat_result | None = None
    try:
        parent_descriptor, parent_metadata = _open_parent(archive_path)
        destination_metadata = _validate_destination(
            archive_path,
            root,
            protected_paths,
            allow_existing_symlink=False,
            parent_descriptor=parent_descriptor,
        )
        _clear_stale_at(parent_descriptor, archive_path.name, destination_metadata)
        _require_current_parent(archive_path, parent_metadata)

        temporary_descriptor = -1
        for _ in range(100):
            temporary_name = f".{archive_path.name}.archive-{secrets.token_hex(8)}"
            try:
                temporary_descriptor = os.open(
                    temporary_name,
                    os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
                    0o600,
                    dir_fd=parent_descriptor,
                )
            except FileExistsError:
                continue
            break
        if temporary_descriptor < 0:
            raise OSError(f"Unable to allocate temporary archive beside {archive_path}")

        with os.fdopen(temporary_descriptor, "w+b") as temporary_archive:
            root_descriptor = os.open(root, _DIRECTORY_FLAGS)
            write_failure: BaseException | None = None
            try:
                opened_root = os.fstat(root_descriptor)
                if not stat.S_ISDIR(opened_root.st_mode) or not _same_identity(root_metadata, opened_root):
                    _reject(f"Archive root is not a directory: {root}")
                if archive_format == "zip":
                    _write_zip(temporary_archive, root_descriptor)
                else:
                    _write_tar(temporary_archive, root_descriptor)
                temporary_archive.flush()
                os.fsync(temporary_archive.fileno())
                completed_metadata = os.fstat(temporary_archive.fileno())
            except BaseException as error:
                write_failure = error
                raise
            finally:
                try:
                    os.close(root_descriptor)
                except Exception as cleanup_error:
                    if write_failure is None:
                        raise
                    logger.error(f"Archive write failed and source descriptor cleanup also failed: {cleanup_error}")

        _require_current_parent(archive_path, parent_metadata)
        if temporary_name is None:
            raise RuntimeError("Temporary archive name was lost before install.")
        os.link(
            temporary_name,
            archive_path.name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if completed_metadata is None:
            raise RuntimeError("Completed archive identity was lost before install.")
        try:
            _require_current_parent(archive_path, parent_metadata)
            _require_destination_identity(
                parent_descriptor,
                archive_path.name,
                completed_metadata,
                "during install",
            )
            os.unlink(temporary_name, dir_fd=parent_descriptor)
            temporary_name = None
            _require_current_parent(archive_path, parent_metadata)
            _require_destination_identity(
                parent_descriptor,
                archive_path.name,
                completed_metadata,
                "after install",
            )
        except Exception:
            try:
                _unlink_if_same(parent_descriptor, archive_path.name, completed_metadata)
            except Exception as rollback_error:
                logger.error(
                    f"Temporary archive cleanup failed and public archive rollback also failed: {rollback_error}"
                )
            raise
    finally:
        if temporary_name is not None and parent_descriptor >= 0:
            try:
                os.unlink(temporary_name, dir_fd=parent_descriptor)
            except FileNotFoundError:
                pass
            except Exception as cleanup_error:
                logger.error(
                    f"Archive attempt failed and temporary cleanup also failed; retained partial "
                    f"'{temporary_name}': {cleanup_error}"
                )
        if parent_descriptor >= 0:
            try:
                os.close(parent_descriptor)
            except Exception as cleanup_error:
                logger.error(f"Archive parent descriptor cleanup failed: {cleanup_error}")

    return str(archive_path)
