from __future__ import annotations

import contextlib
import hashlib
import os
import shutil
import stat
import subprocess
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

_PORCELAIN_V1_STATUS_PAIRS = frozenset(
    {
        b" M",
        b" T",
        b" A",
        b" D",
        b" R",
        b" C",
        b"M ",
        b"MM",
        b"MT",
        b"MD",
        b"T ",
        b"TM",
        b"TT",
        b"TD",
        b"A ",
        b"AM",
        b"AT",
        b"AD",
        b"D ",
        b"R ",
        b"RM",
        b"RT",
        b"RD",
        b"C ",
        b"CM",
        b"CT",
        b"CD",
        b"DD",
        b"AU",
        b"UD",
        b"UA",
        b"DU",
        b"AA",
        b"UU",
        b"??",
        b"!!",
    }
)


@dataclass(frozen=True)
class OverlayChange:
    """One file operation needed to reproduce a working-tree change."""

    action: Literal["copy", "delete"]
    path: str


@dataclass(frozen=True)
class MutationWorkspace:
    """A detached mutation-sweep workspace and its restoration baseline."""

    real_root: Path
    sweep_root: Path
    head: str
    overlay_changes: tuple[OverlayChange, ...]
    baseline_manifest: tuple[OverlayChange, ...]
    baseline_status: bytes
    baseline_digests: dict[str, str]

    def target_path(self, relative: str) -> Path:
        """Return a validated target path in the disposable workspace.

        Args:
            relative: Repository-relative target name.

        Returns:
            The existing target path inside the disposable workspace.

        Raises:
            RuntimeError: If the target is missing or escapes the workspace.
        """
        try:
            return confined_path(self.sweep_root, relative, must_exist=True)
        except ValueError as exc:
            raise RuntimeError(str(exc)) from exc

    def verify_baseline(self) -> None:
        """Verify that the disposable workspace still matches its baseline.

        Raises:
            RuntimeError: If Git status or any captured path state has drifted.
        """
        for relative, expected in self.baseline_digests.items():
            path = self.sweep_root / relative
            if expected == "deleted":
                if os.path.lexists(path):
                    raise RuntimeError(f"{relative}: deletion mismatch")
                continue
            if expected.startswith("symlink:"):
                if not path.is_symlink() or os.fsencode(path.readlink()).hex() != expected.removeprefix("symlink:"):
                    raise RuntimeError(f"{relative}: symlink mismatch")
                continue
            if path.is_symlink() or not path.is_file() or _sha256(path) != expected.removeprefix("sha256:"):
                raise RuntimeError(f"{relative}: content mismatch")

        status_result = subprocess.run(
            ["git", "-C", str(self.sweep_root), "status", "--porcelain=v1", "-z", "--untracked-files=all"],
            cwd=self.real_root,
            capture_output=True,
            check=False,
        )
        if status_result.returncode != 0:
            raise RuntimeError(status_result.stderr.decode("utf-8", errors="replace").strip())
        try:
            current_manifest = parse_porcelain_z(status_result.stdout, self.sweep_root)
        except ValueError as exc:
            raise RuntimeError(f"malformed porcelain rename or status: {exc}") from exc
        if current_manifest != self.baseline_manifest:
            raise RuntimeError(
                f"workspace status mismatch: baseline={self.baseline_status!r}, current={status_result.stdout!r}"
            )


def _registered_worktree_roots(payload: bytes) -> tuple[Path, ...]:
    """Extract absolute worktree roots from Git's NUL-delimited porcelain output."""
    if not payload:
        return ()
    if not payload.endswith(b"\0\0"):
        raise RuntimeError("malformed worktree list: missing record terminator")
    roots: list[Path] = []
    for record in payload[:-2].split(b"\0\0"):
        fields = record.split(b"\0")
        if (
            not fields
            or not fields[0].startswith(b"worktree ")
            or any(field.startswith(b"worktree ") for field in fields[1:])
        ):
            raise RuntimeError("malformed worktree list: invalid record")
        root = Path(os.fsdecode(fields[0].removeprefix(b"worktree ")))
        if not root.is_absolute():
            raise RuntimeError(f"malformed worktree list: non-absolute root: {root}")
        roots.append(root.resolve(strict=False))
    return tuple(roots)


_DIRECTORY_OPEN_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)


def _relative_parts(relative: str) -> tuple[str, ...]:
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


@contextlib.contextmanager
def _parent_directory_fd(root: Path, relative: str, *, create: bool) -> Iterator[tuple[int, str]]:
    """Open a relative path's parent without following directory symlinks."""
    parts = _relative_parts(relative)
    descriptors = [os.open(root, _DIRECTORY_OPEN_FLAGS)]
    try:
        for component in parts[:-1]:
            if create:
                with contextlib.suppress(FileExistsError):
                    os.mkdir(component, dir_fd=descriptors[-1])
            descriptors.append(os.open(component, _DIRECTORY_OPEN_FLAGS, dir_fd=descriptors[-1]))
        yield descriptors[-1], parts[-1]
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _remove_entry_at(parent_fd: int, name: str) -> None:
    """Remove one directory entry without following its final component."""
    try:
        entry_stat = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    if not stat.S_ISDIR(entry_stat.st_mode):
        os.unlink(name, dir_fd=parent_fd)
        return

    child_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=parent_fd)
    try:
        opened_stat = os.fstat(child_fd)
        for child_name in os.listdir(child_fd):
            _remove_entry_at(child_fd, child_name)
        current_stat = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (current_stat.st_dev, current_stat.st_ino) != (opened_stat.st_dev, opened_stat.st_ino):
            raise RuntimeError(f"workspace destination changed during removal: {name}")
        os.rmdir(name, dir_fd=parent_fd)
    finally:
        os.close(child_fd)


def _copy_regular_file_at(
    source_parent_fd: int, source_name: str, destination_parent_fd: int, destination_name: str
) -> None:
    """Copy one regular file through anchored directory descriptors."""
    source_fd = os.open(source_name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=source_parent_fd)
    destination_fd: int | None = None
    try:
        source_stat = os.fstat(source_fd)
        if not stat.S_ISREG(source_stat.st_mode):
            raise RuntimeError(f"overlay source is not a regular file: {source_name}")
        _remove_entry_at(destination_parent_fd, destination_name)
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
            _remove_entry_at(destination_parent_fd, destination_name)
        raise
    finally:
        if destination_fd is not None:
            os.close(destination_fd)
        os.close(source_fd)


def _symlink_target_relative(relative: str, link_text: str) -> str:
    """Return a lexically confined target name for a copied relative symlink."""
    if Path(link_text).is_absolute():
        raise RuntimeError(f"workspace symlink escapes workspace root: {relative}")
    target = Path(os.path.normpath(str(Path(relative).parent / link_text))).as_posix()
    try:
        _relative_parts(target)
    except ValueError as exc:
        raise RuntimeError(f"workspace symlink escapes workspace root: {relative}") from exc
    return target


def _copy_overlay_entry(real_root: Path, sweep_root: Path, relative: str) -> None:
    """Copy a confined source to an anchored, no-follow destination."""
    try:
        with _parent_directory_fd(real_root, relative, create=False) as (source_parent_fd, source_name):
            source_stat = os.stat(source_name, dir_fd=source_parent_fd, follow_symlinks=False)
            with _parent_directory_fd(sweep_root, relative, create=True) as (destination_parent_fd, destination_name):
                if stat.S_ISLNK(source_stat.st_mode):
                    link_text = os.readlink(source_name, dir_fd=source_parent_fd)
                    target_relative = _symlink_target_relative(relative, link_text)
                    confined_path(real_root, target_relative, must_exist=True)
                    confined_path(sweep_root, target_relative, must_exist=True)
                    _remove_entry_at(destination_parent_fd, destination_name)
                    os.symlink(link_text, destination_name, dir_fd=destination_parent_fd)
                else:
                    _copy_regular_file_at(
                        source_parent_fd,
                        source_name,
                        destination_parent_fd,
                        destination_name,
                    )
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"unsafe workspace destination: {relative}: {exc}") from exc


def _delete_overlay_entry(sweep_root: Path, relative: str) -> None:
    """Delete a destination entry without following parent or final symlinks."""
    try:
        with _parent_directory_fd(sweep_root, relative, create=False) as (parent_fd, name):
            _remove_entry_at(parent_fd, name)
    except FileNotFoundError:
        return
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"unsafe workspace destination: {relative}: {exc}") from exc


def _sha256(path: Path) -> str:
    """Return the SHA-256 hex digest of one regular file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _capture_baseline_digests(
    sweep_root: Path,
    overlay_changes: tuple[OverlayChange, ...],
    targets: tuple[str, ...],
) -> dict[str, str]:
    """Capture file, symlink, and deletion state after applying the overlay."""
    captured: dict[str, str] = {}
    for change in overlay_changes:
        if change.action == "delete":
            captured[change.path] = "deleted"
            continue
        try:
            path = confined_path(sweep_root, change.path, must_exist=True)
        except ValueError as exc:
            raise RuntimeError(str(exc)) from exc
        if path.is_symlink():
            captured[change.path] = f"symlink:{os.fsencode(path.readlink()).hex()}"
        elif path.is_file():
            captured[change.path] = f"sha256:{_sha256(path)}"
        else:
            raise RuntimeError(f"baseline path is not a regular file: {change.path}")

    for relative in targets:
        try:
            target = confined_path(sweep_root, relative, must_exist=True)
        except ValueError as exc:
            raise RuntimeError(str(exc)) from exc
        if target.is_symlink() or not target.is_file():
            raise RuntimeError(f"mutation target is not a regular file: {relative}")
        captured[relative] = f"sha256:{_sha256(target)}"
    return captured


def _apply_overlay_changes(
    real_root: Path,
    sweep_root: Path,
    changes: tuple[OverlayChange, ...],
    excluded_outputs: tuple[Path, ...],
    registered_roots: tuple[Path, ...],
) -> None:
    """Apply validated working-tree operations to the disposable worktree."""
    resolved_outputs = tuple(
        (output if output.is_absolute() else real_root / output).resolve(strict=False) for output in excluded_outputs
    )
    nested_roots = tuple(
        root for root in registered_roots if root != real_root and root != sweep_root and root.is_relative_to(real_root)
    )
    for change in changes:
        try:
            source = confined_path(real_root, change.path, must_exist=change.action == "copy")
        except ValueError as exc:
            raise RuntimeError(str(exc)) from exc
        resolved_source = source.resolve(strict=False)
        if any(resolved_source == output or resolved_source.is_relative_to(output) for output in resolved_outputs):
            raise RuntimeError(f"excluded output cannot be overlaid: {change.path}")
        if any(resolved_source == root or resolved_source.is_relative_to(root) for root in nested_roots):
            raise RuntimeError(f"nested registered worktree cannot be overlaid: {change.path}")

        if change.action == "copy":
            _copy_overlay_entry(real_root, sweep_root, change.path)
        else:
            _delete_overlay_entry(sweep_root, change.path)


def detached_head_workspace(
    real_root: Path,
    targets: tuple[str, ...],
    excluded_outputs: tuple[Path, ...],
) -> contextlib.AbstractContextManager[MutationWorkspace]:
    """Create a detached mutation workspace and dispose it on context exit.

    Args:
        real_root: Source working tree whose current state is overlaid.
        targets: Repository-relative regular files selected for mutation.
        excluded_outputs: Source paths that must never be copied into the workspace.

    Returns:
        A context manager yielding the detached, overlaid mutation workspace.

    Raises:
        RuntimeError: If Git, overlay validation, baseline capture, or cleanup fails.
    """

    @contextlib.contextmanager
    def workspace_context() -> Iterator[MutationWorkspace]:
        resolved_real_root = real_root.resolve(strict=True)
        head_result = subprocess.run(
            ["git", "rev-parse", "--verify", "HEAD^{commit}"],
            cwd=resolved_real_root,
            capture_output=True,
            check=False,
        )
        if head_result.returncode != 0:
            raise RuntimeError(head_result.stderr.decode("utf-8", errors="replace").strip())
        head = head_result.stdout.decode("ascii").strip()

        sweep_root = Path(tempfile.mkdtemp(prefix="vntyper-mutation-")).resolve(strict=True)
        try:
            add_result = subprocess.run(
                ["git", "worktree", "add", "--detach", str(sweep_root), head],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
            if add_result.returncode != 0:
                raise RuntimeError(add_result.stderr.decode("utf-8", errors="replace").strip())
            disposable_head_result = subprocess.run(
                ["git", "-C", str(sweep_root), "rev-parse", "--verify", "HEAD^{commit}"],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
            if disposable_head_result.returncode != 0:
                raise RuntimeError(disposable_head_result.stderr.decode("utf-8", errors="replace").strip())
            disposable_head = disposable_head_result.stdout.decode("ascii").strip()
            if disposable_head != head:
                raise RuntimeError(f"disposable HEAD {disposable_head} differs from captured HEAD {head}")

            real_status_result = subprocess.run(
                ["git", "status", "--porcelain=v1", "-z", "--untracked-files=all"],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
            if real_status_result.returncode != 0:
                raise RuntimeError(real_status_result.stderr.decode("utf-8", errors="replace").strip())
            try:
                overlay_changes = parse_porcelain_z(real_status_result.stdout, resolved_real_root)
            except ValueError as exc:
                raise RuntimeError(f"malformed porcelain rename or status: {exc}") from exc

            worktree_list_result = subprocess.run(
                ["git", "worktree", "list", "--porcelain", "-z"],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
            if worktree_list_result.returncode != 0:
                raise RuntimeError(worktree_list_result.stderr.decode("utf-8", errors="replace").strip())
            registered_roots = _registered_worktree_roots(worktree_list_result.stdout)
            _apply_overlay_changes(
                resolved_real_root,
                sweep_root,
                overlay_changes,
                excluded_outputs,
                registered_roots,
            )

            sweep_status_result = subprocess.run(
                ["git", "-C", str(sweep_root), "status", "--porcelain=v1", "-z", "--untracked-files=all"],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
            if sweep_status_result.returncode != 0:
                raise RuntimeError(sweep_status_result.stderr.decode("utf-8", errors="replace").strip())
            try:
                baseline_manifest = parse_porcelain_z(sweep_status_result.stdout, sweep_root)
            except ValueError as exc:
                raise RuntimeError(f"malformed porcelain rename or status: {exc}") from exc
            baseline_digests = _capture_baseline_digests(sweep_root, overlay_changes, targets)

            yield MutationWorkspace(
                real_root=resolved_real_root,
                sweep_root=sweep_root,
                head=head,
                overlay_changes=overlay_changes,
                baseline_manifest=baseline_manifest,
                baseline_status=sweep_status_result.stdout,
                baseline_digests=baseline_digests,
            )
        finally:
            try:
                remove_result = subprocess.run(
                    ["git", "worktree", "remove", "--force", str(sweep_root)],
                    cwd=resolved_real_root,
                    capture_output=True,
                    check=False,
                )
            except (OSError, KeyboardInterrupt) as exc:
                raise RuntimeError(f"orphaned worktree: {sweep_root}: cleanup command failed: {exc}") from exc
            if remove_result.returncode != 0:
                stderr = remove_result.stderr.decode("utf-8", errors="replace").strip()
                raise RuntimeError(f"orphaned worktree: {sweep_root}: {stderr}")
            if sweep_root.exists():
                shutil.rmtree(sweep_root)

    return workspace_context()


def confined_path(root: Path, relative: str, *, must_exist: bool) -> Path:
    """Resolve a repository-relative path without allowing a worktree escape.

    Args:
        root: Root of the workspace containing the path.
        relative: Workspace-relative path to validate.
        must_exist: Whether the path must currently exist, including symlinks.

    Returns:
        The validated lexical path beneath the resolved workspace root.

    Raises:
        ValueError: If the name is unsafe, escapes the root, or must exist but does not.
    """
    candidate = Path(*_relative_parts(relative))
    resolved_root = root.resolve()
    lexical = resolved_root / candidate
    try:
        resolved_parent = lexical.parent.resolve(strict=must_exist)
    except FileNotFoundError:
        raise ValueError(f"workspace path does not exist: {relative}") from None
    if not resolved_parent.is_relative_to(resolved_root):
        raise ValueError(f"workspace path escapes workspace root: {relative}")
    if ".git" in resolved_parent.relative_to(resolved_root).parts:
        raise ValueError(f"unsafe workspace path: {relative}")
    if must_exist and not os.path.lexists(lexical):
        raise ValueError(f"workspace path does not exist: {relative}")
    if lexical.is_symlink():
        try:
            resolved = lexical.resolve(strict=True)
        except FileNotFoundError:
            raise ValueError(f"workspace path does not exist: {relative}") from None
        if not resolved.is_relative_to(resolved_root):
            raise ValueError(f"workspace path escapes workspace root: {relative}")
        if ".git" in resolved.relative_to(resolved_root).parts:
            raise ValueError(f"unsafe workspace path: {relative}")
    return lexical


def _affected_paths_from_porcelain_z(payload: bytes) -> tuple[str, ...]:
    """Decode paths affected by porcelain v1 ``-z`` records.

    Rename and copy records contribute both their current path and their original path.
    This function does not inspect the filesystem.

    Args:
        payload: Raw output from ``git status --porcelain=v1 -z``.

    Returns:
        Affected paths in record order, including rename and copy origins.

    Raises:
        ValueError: If the payload is malformed, ignored, or has ambiguous rename metadata.
    """
    if not payload:
        return ()
    if not payload.endswith(b"\0"):
        raise ValueError("malformed porcelain payload: missing NUL terminator")

    records = payload[:-1].split(b"\0")
    affected: list[str] = []
    index = 0
    while index < len(records):
        record = records[index]
        if len(record) < 4 or record[2:3] != b" " or not record[3:]:
            raise ValueError("malformed porcelain record")
        status = record[:2]
        if status == b"!!":
            raise ValueError("ignored porcelain record is not allowed")

        rename_fields = sum(field in b"RC" for field in status)
        if rename_fields > 1:
            raise ValueError("conflicting rename encoding")
        if status not in _PORCELAIN_V1_STATUS_PAIRS:
            raise ValueError(f"invalid porcelain status: {os.fsdecode(status)}")

        affected.append(os.fsdecode(record[3:]))
        if rename_fields == 1:
            index += 1
            if index >= len(records) or not records[index]:
                raise ValueError("rename record is missing original path")
            affected.append(os.fsdecode(records[index]))
        index += 1
    return tuple(affected)


def parse_porcelain_z(payload: bytes, real_root: Path) -> tuple[OverlayChange, ...]:
    """Resolve porcelain records into a safe snapshot of current filesystem state.

    Args:
        payload: Raw output from ``git status --porcelain=v1 -z``.
        real_root: Root of the working tree whose current files decide each action.

    Returns:
        One path-sorted copy or delete action per affected path.

    Raises:
        ValueError: If any porcelain record or affected path is unsafe.
    """
    affected_paths = _affected_paths_from_porcelain_z(payload)
    validated_paths = {path: confined_path(real_root, path, must_exist=False) for path in dict.fromkeys(affected_paths)}
    changes = tuple(
        OverlayChange("copy" if os.path.lexists(validated) else "delete", path)
        for path, validated in validated_paths.items()
    )
    return tuple(sorted(changes, key=lambda change: change.path))
