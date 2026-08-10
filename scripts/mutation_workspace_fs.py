from __future__ import annotations

import contextlib
import hashlib
import os
import stat
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

DIRECTORY_OPEN_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
PROC_SELF_FD = Path("/proc/self/fd")


@dataclass(frozen=True)
class RootCapability:
    """One open directory descriptor and the identity it pins."""

    path: Path
    descriptor: int
    device: int
    inode: int


def relative_parts(relative: str) -> tuple[str, ...]:
    """Return the components of one lexically safe workspace-relative path."""
    candidate = Path(relative)
    if (
        relative in {"", "."}
        or candidate.is_absolute()
        or candidate.as_posix() != relative
        or ".." in candidate.parts
        or ".git" in candidate.parts
    ):
        raise ValueError(f"unsafe workspace path: {relative}")
    return candidate.parts


def symlink_target_relative(relative: str, link_text: str) -> str:
    """Return a lexically confined target name for a copied relative symlink.

    Args:
        relative: Repository-relative path of the symlink.
        link_text: Exact relative link text to resolve lexically.

    Returns:
        The normalized repository-relative target path.

    Raises:
        RuntimeError: If the link is absolute or escapes the workspace.
    """
    if Path(link_text).is_absolute():
        raise RuntimeError(f"workspace symlink escapes workspace root: {relative}")
    target = Path(os.path.normpath(str(Path(relative).parent / link_text))).as_posix()
    try:
        relative_parts(target)
    except ValueError as exc:
        raise RuntimeError(f"workspace symlink escapes workspace root: {relative}") from exc
    return target


def open_root_capability(root: Path) -> RootCapability:
    """Open and pin one workspace root directory.

    Args:
        root: Exact root path to open without following its final component.

    Returns:
        An open capability recording the root's device and inode.

    Raises:
        RuntimeError: If the path is not the directory that was opened.
        OSError: If the directory cannot be opened or inspected.
    """
    descriptor = os.open(root, DIRECTORY_OPEN_FLAGS)
    try:
        root_stat = os.fstat(descriptor)
        capability = RootCapability(root, descriptor, root_stat.st_dev, root_stat.st_ino)
        assert_root_identity(capability)
    except BaseException:
        os.close(descriptor)
        raise
    return capability


def git_capability_path(capability: RootCapability) -> Path:
    """Return the inherited proc-fd path that names the pinned root inode.

    Args:
        capability: Open workspace root whose descriptor must be inherited.

    Returns:
        The proc-fd path naming the same directory device and inode.

    Raises:
        RuntimeError: If this environment cannot expose the opened directory to Git.
    """
    proc_path = PROC_SELF_FD / str(capability.descriptor)
    try:
        proc_stat = os.stat(proc_path)
        descriptor_stat = os.fstat(capability.descriptor)
    except OSError as exc:
        raise RuntimeError(f"proc/self/fd capability is unavailable: {proc_path}: {exc}") from exc
    if (
        not stat.S_ISDIR(proc_stat.st_mode)
        or proc_stat.st_dev != descriptor_stat.st_dev
        or proc_stat.st_ino != descriptor_stat.st_ino
    ):
        raise RuntimeError(f"proc/self/fd capability does not name the workspace root: {proc_path}")
    return proc_path


def remove_exact_root(root: Path, expected_identity: tuple[int, int]) -> None:
    """Remove one created root only while its parent entry has the expected identity.

    Args:
        root: Exact created directory entry to remove.
        expected_identity: Required device and inode for that entry.

    Raises:
        OSError: If anchored traversal or removal fails.
        RuntimeError: If the public entry no longer has the expected identity.
    """
    parent_fd = os.open(root.parent, DIRECTORY_OPEN_FLAGS)
    try:
        remove_entry_at(parent_fd, root.name, expected_identity=expected_identity)
    finally:
        os.close(parent_fd)


def assert_root_identity(capability: RootCapability) -> None:
    """Require the capability path to still name the pinned directory inode."""
    try:
        path_stat = os.stat(capability.path, follow_symlinks=False)
    except OSError as exc:
        raise RuntimeError(f"workspace root identity mismatch: {capability.path}") from exc
    if (
        not stat.S_ISDIR(path_stat.st_mode)
        or path_stat.st_dev != capability.device
        or path_stat.st_ino != capability.inode
    ):
        raise RuntimeError(f"workspace root identity mismatch: {capability.path}")


@contextlib.contextmanager
def using_root_capability(root: Path | RootCapability) -> Iterator[RootCapability]:
    """Borrow a pinned capability or temporarily open one for a direct helper call."""
    if isinstance(root, RootCapability):
        owned = False
        capability = root
    else:
        owned = True
        capability = open_root_capability(root)
    try:
        assert_root_identity(capability)
        yield capability
    finally:
        try:
            assert_root_identity(capability)
        finally:
            if owned:
                os.close(capability.descriptor)


@contextlib.contextmanager
def parent_directory_fd(root: Path | RootCapability, relative: str, *, create: bool) -> Iterator[tuple[int, str]]:
    """Open a relative path's parent without following directory symlinks."""
    parts = relative_parts(relative)
    with using_root_capability(root) as capability:
        descriptors = [os.dup(capability.descriptor)]
        try:
            for component in parts[:-1]:
                if create:
                    with contextlib.suppress(FileExistsError):
                        os.mkdir(component, dir_fd=descriptors[-1])
                descriptors.append(os.open(component, DIRECTORY_OPEN_FLAGS, dir_fd=descriptors[-1]))
            yield descriptors[-1], parts[-1]
        finally:
            for descriptor in reversed(descriptors):
                os.close(descriptor)


def remove_entry_at(parent_fd: int, name: str, *, expected_identity: tuple[int, int] | None = None) -> None:
    """Remove one directory entry without following its final component."""
    try:
        entry_stat = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    if expected_identity is not None and (entry_stat.st_dev, entry_stat.st_ino) != expected_identity:
        raise RuntimeError(f"workspace destination changed during removal: {name}")
    if not stat.S_ISDIR(entry_stat.st_mode):
        os.unlink(name, dir_fd=parent_fd)
        return

    child_fd = os.open(name, DIRECTORY_OPEN_FLAGS, dir_fd=parent_fd)
    try:
        opened_stat = os.fstat(child_fd)
        if (opened_stat.st_dev, opened_stat.st_ino) != (entry_stat.st_dev, entry_stat.st_ino):
            raise RuntimeError(f"workspace destination changed during removal: {name}")
        for child_name in os.listdir(child_fd):
            remove_entry_at(child_fd, child_name)
        current_stat = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (current_stat.st_dev, current_stat.st_ino) != (opened_stat.st_dev, opened_stat.st_ino):
            raise RuntimeError(f"workspace destination changed during removal: {name}")
        os.rmdir(name, dir_fd=parent_fd)
    finally:
        os.close(child_fd)


def remove_pinned_root_if_present(capability: RootCapability) -> None:
    """Remove the pinned root entry only when its exact inode remains at its path."""
    parent_fd = os.open(capability.path.parent, DIRECTORY_OPEN_FLAGS)
    try:
        try:
            path_stat = os.stat(capability.path.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            return
        if (
            not stat.S_ISDIR(path_stat.st_mode)
            or path_stat.st_dev != capability.device
            or path_stat.st_ino != capability.inode
        ):
            raise RuntimeError("workspace identity mismatch")
        remove_entry_at(
            parent_fd,
            capability.path.name,
            expected_identity=(capability.device, capability.inode),
        )
    finally:
        os.close(parent_fd)


def copy_regular_file_at(
    source_parent_fd: int,
    source_name: str,
    destination_parent_fd: int,
    destination_name: str,
) -> None:
    """Copy one regular file through anchored directory descriptors."""
    source_fd = os.open(source_name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=source_parent_fd)
    destination_fd: int | None = None
    try:
        source_stat = os.fstat(source_fd)
        if not stat.S_ISREG(source_stat.st_mode):
            raise RuntimeError(f"overlay source is not a regular file: {source_name}")
        remove_entry_at(destination_parent_fd, destination_name)
        destination_fd = os.open(
            destination_name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
            stat.S_IMODE(source_stat.st_mode),
            dir_fd=destination_parent_fd,
        )
        while chunk := os.read(source_fd, 1024 * 1024):
            offset = 0
            while offset < len(chunk):
                offset += os.write(destination_fd, chunk[offset:])
        os.fchmod(destination_fd, stat.S_IMODE(source_stat.st_mode))
        os.utime(destination_fd, ns=(source_stat.st_atime_ns, source_stat.st_mtime_ns))
        with contextlib.suppress(OSError):
            for attribute in os.listxattr(source_fd):
                with contextlib.suppress(OSError):
                    os.setxattr(destination_fd, attribute, os.getxattr(source_fd, attribute))
    except BaseException:
        if destination_fd is not None:
            os.close(destination_fd)
            destination_fd = None
            remove_entry_at(destination_parent_fd, destination_name)
        raise
    finally:
        if destination_fd is not None:
            os.close(destination_fd)
        os.close(source_fd)


def entry_stat(root: Path | RootCapability, relative: str) -> os.stat_result | None:
    """lstat one anchored relative entry, returning ``None`` when it is absent."""
    try:
        with parent_directory_fd(root, relative, create=False) as (parent_fd, name):
            try:
                return os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
            except FileNotFoundError:
                return None
    except FileNotFoundError:
        return None


def _sha256_fd(descriptor: int) -> str:
    """Return the SHA-256 digest read from a newly opened file descriptor."""
    digest = hashlib.sha256()
    while chunk := os.read(descriptor, 1024 * 1024):
        digest.update(chunk)
    return digest.hexdigest()


def captured_path_state(root: Path | RootCapability, relative: str) -> str | None:
    """Capture an anchored regular-file digest or exact symlink text."""
    try:
        with parent_directory_fd(root, relative, create=False) as (parent_fd, name):
            try:
                entry_status = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
            except FileNotFoundError:
                return None
            if stat.S_ISLNK(entry_status.st_mode):
                return f"symlink:{os.fsencode(os.readlink(name, dir_fd=parent_fd)).hex()}"
            if not stat.S_ISREG(entry_status.st_mode):
                raise RuntimeError(f"baseline path is not a regular file: {relative}")
            descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent_fd)
            try:
                if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                    raise RuntimeError(f"baseline path is not a regular file: {relative}")
                return f"sha256:{_sha256_fd(descriptor)}"
            finally:
                os.close(descriptor)
    except FileNotFoundError:
        return None
    except OSError as exc:
        raise RuntimeError(f"{relative}: baseline capture failed: {exc}") from exc
