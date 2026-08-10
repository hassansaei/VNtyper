from __future__ import annotations

import contextlib
import os
import stat
import subprocess
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from mutation_workspace_fs import RootCapability as _RootCapability
from mutation_workspace_fs import assert_root_identity as _assert_root_identity
from mutation_workspace_fs import captured_path_state as _captured_path_state
from mutation_workspace_fs import copy_regular_file_at as _copy_regular_file_at
from mutation_workspace_fs import entry_stat as _entry_stat
from mutation_workspace_fs import open_root_capability as _open_root_capability
from mutation_workspace_fs import parent_directory_fd as _parent_directory_fd
from mutation_workspace_fs import relative_parts as _relative_parts
from mutation_workspace_fs import remove_entry_at as _remove_entry_at
from mutation_workspace_fs import remove_pinned_root_if_present as _remove_pinned_root_if_present
from mutation_workspace_fs import using_root_capability as _using_root_capability

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
    _root_capability: _RootCapability | None = field(default=None, repr=False, compare=False)

    def target_path(self, relative: str) -> Path:
        """Return a validated target path in the disposable workspace.

        Args:
            relative: Repository-relative target name.

        Returns:
            The existing target path inside the disposable workspace.

        Raises:
            RuntimeError: If the target is missing or escapes the workspace.
        """
        with _using_root_capability(self._root_capability or self.sweep_root):
            try:
                return confined_path(self.sweep_root, relative, must_exist=True)
            except ValueError as exc:
                raise RuntimeError(str(exc)) from exc

    def verify_baseline(self) -> None:
        """Verify that the disposable workspace still matches its baseline.

        Raises:
            RuntimeError: If Git status or any captured path state has drifted.
        """
        with _using_root_capability(self._root_capability or self.sweep_root) as capability:
            for relative, expected in self.baseline_digests.items():
                current = _captured_path_state(capability, relative)
                if expected == "deleted":
                    if current is not None:
                        raise RuntimeError(f"{relative}: deletion mismatch")
                    continue
                if expected.startswith("symlink:"):
                    if current != expected:
                        raise RuntimeError(f"{relative}: symlink mismatch")
                    continue
                if current != expected:
                    raise RuntimeError(f"{relative}: content mismatch")

            _assert_root_identity(capability)
            status_result = subprocess.run(
                ["git", "-C", str(self.sweep_root), "status", "--porcelain=v1", "-z", "--untracked-files=all"],
                cwd=self.real_root,
                capture_output=True,
                check=False,
            )
            _assert_root_identity(capability)
            if status_result.returncode != 0:
                raise RuntimeError(status_result.stderr.decode("utf-8", errors="replace").strip())
            try:
                current_manifest = _parse_porcelain_with_capability(status_result.stdout, capability)
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


@dataclass(frozen=True)
class _OverlayPlanEntry:
    """One fully validated overlay operation."""

    change: OverlayChange
    source_kind: Literal["delete", "regular", "symlink"]
    link_text: str | None = None
    link_target: str | None = None


def _inspect_overlay_source(real_root: Path, change: OverlayChange) -> _OverlayPlanEntry:
    """Inspect one already confined source without following its final component."""
    if change.action == "delete":
        return _OverlayPlanEntry(change, "delete")
    try:
        with _parent_directory_fd(real_root, change.path, create=False) as (parent_fd, name):
            source_stat = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
            if stat.S_ISREG(source_stat.st_mode):
                return _OverlayPlanEntry(change, "regular")
            if stat.S_ISLNK(source_stat.st_mode):
                link_text = os.readlink(name, dir_fd=parent_fd)
                link_target = _symlink_target_relative(change.path, link_text)
                confined_path(real_root, link_target, must_exist=True)
                return _OverlayPlanEntry(change, "symlink", link_text, link_target)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"unsafe overlay source: {change.path}: {exc}") from exc
    raise RuntimeError(f"overlay source is not a regular file or symlink: {change.path}")


def _validate_planned_symlink(
    entry: _OverlayPlanEntry,
    entries_by_path: dict[str, _OverlayPlanEntry],
    sweep_root: _RootCapability,
    visiting: frozenset[str] = frozenset(),
) -> None:
    """Require a planned symlink chain to end at a safe final file."""
    target = entry.link_target
    if target is None:
        raise RuntimeError(f"planned symlink target is missing: {entry.change.path}")
    if entry.change.path in visiting:
        raise RuntimeError(f"planned symlink target cycle: {entry.change.path}")
    planned_target = entries_by_path.get(target)
    if planned_target is None:
        if _entry_stat(sweep_root, target) is None:
            raise RuntimeError(f"planned symlink target is missing: {entry.change.path} -> {target}")
        try:
            confined_path(sweep_root.path, target, must_exist=True)
        except ValueError as exc:
            raise RuntimeError(f"planned symlink target is unsafe: {entry.change.path} -> {target}") from exc
        return
    if planned_target.source_kind == "delete":
        raise RuntimeError(f"planned symlink target is deleted: {entry.change.path} -> {target}")
    if planned_target.source_kind == "symlink":
        _validate_planned_symlink(planned_target, entries_by_path, sweep_root, visiting | {entry.change.path})


def _copy_overlay_entry(real_root: Path, sweep_root: _RootCapability, entry: _OverlayPlanEntry) -> None:
    """Copy a confined source to an anchored, no-follow destination."""
    relative = entry.change.path
    try:
        with _parent_directory_fd(real_root, relative, create=False) as (source_parent_fd, source_name):
            source_stat = os.stat(source_name, dir_fd=source_parent_fd, follow_symlinks=False)
            with _parent_directory_fd(sweep_root, relative, create=True) as (destination_parent_fd, destination_name):
                if entry.source_kind == "symlink":
                    if not stat.S_ISLNK(source_stat.st_mode):
                        raise RuntimeError(f"overlay source changed during copy: {relative}")
                    link_text = os.readlink(source_name, dir_fd=source_parent_fd)
                    if link_text != entry.link_text or entry.link_target is None:
                        raise RuntimeError(f"overlay source changed during copy: {relative}")
                    if _entry_stat(sweep_root, entry.link_target) is None:
                        raise RuntimeError(f"planned symlink target is missing: {relative} -> {entry.link_target}")
                    confined_path(sweep_root.path, entry.link_target, must_exist=True)
                    _remove_entry_at(destination_parent_fd, destination_name)
                    os.symlink(link_text, destination_name, dir_fd=destination_parent_fd)
                else:
                    if entry.source_kind != "regular":
                        raise RuntimeError(f"invalid overlay copy plan: {relative}")
                    _copy_regular_file_at(
                        source_parent_fd,
                        source_name,
                        destination_parent_fd,
                        destination_name,
                    )
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"unsafe workspace destination: {relative}: {exc}") from exc


def _delete_overlay_entry(sweep_root: Path | _RootCapability, relative: str) -> None:
    """Delete a destination entry without following parent or final symlinks."""
    try:
        with _parent_directory_fd(sweep_root, relative, create=False) as (parent_fd, name):
            _remove_entry_at(parent_fd, name)
    except FileNotFoundError:
        return
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"unsafe workspace destination: {relative}: {exc}") from exc


def _capture_baseline_digests(
    sweep_root: Path | _RootCapability,
    overlay_changes: tuple[OverlayChange, ...],
    targets: tuple[str, ...],
) -> dict[str, str]:
    """Capture file, symlink, and deletion state after applying the overlay."""
    captured: dict[str, str] = {}
    with _using_root_capability(sweep_root) as capability:
        for change in overlay_changes:
            if change.action == "delete":
                if _entry_stat(capability, change.path) is not None:
                    raise RuntimeError(f"{change.path}: deletion mismatch")
                captured[change.path] = "deleted"
                continue
            current = _captured_path_state(capability, change.path)
            if current is None:
                raise RuntimeError(f"workspace path does not exist: {change.path}")
            captured[change.path] = current

        for relative in targets:
            try:
                current = _captured_path_state(capability, relative)
            except RuntimeError as exc:
                if "baseline path is not a regular file" in str(exc):
                    raise RuntimeError(f"mutation target is not a regular file: {relative}") from exc
                raise
            if current is None:
                raise RuntimeError(f"workspace path does not exist: {relative}")
            if current.startswith("symlink:"):
                raise RuntimeError(f"mutation target is not a regular file: {relative}")
            captured[relative] = current
    return captured


def _parse_porcelain_with_capability(
    payload: bytes,
    capability: _RootCapability,
) -> tuple[OverlayChange, ...]:
    """Resolve porcelain actions from final state beneath a pinned root."""
    affected_paths = _affected_paths_from_porcelain_z(payload)
    validated = tuple(dict.fromkeys(affected_paths))
    for relative in validated:
        _relative_parts(relative)
    changes = tuple(
        OverlayChange("copy" if _entry_stat(capability, relative) is not None else "delete", relative)
        for relative in validated
    )
    return tuple(sorted(changes, key=lambda change: change.path))


def _apply_overlay_changes(
    real_root: Path,
    sweep_root: Path | _RootCapability,
    changes: tuple[OverlayChange, ...],
    excluded_outputs: tuple[Path, ...],
    registered_roots: tuple[Path, ...],
) -> None:
    """Apply validated working-tree operations to the disposable worktree."""
    resolved_outputs = tuple(
        (output if output.is_absolute() else real_root / output).resolve(strict=False) for output in excluded_outputs
    )
    nested_roots = tuple(
        root
        for root in registered_roots
        if root != real_root
        and root != (sweep_root.path if isinstance(sweep_root, _RootCapability) else sweep_root)
        and root.is_relative_to(real_root)
    )
    with _using_root_capability(sweep_root) as capability:
        planned: list[_OverlayPlanEntry] = []
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
            planned.append(_inspect_overlay_source(real_root, change))

        entries_by_path = {entry.change.path: entry for entry in planned}
        for entry in planned:
            if entry.source_kind == "symlink":
                _validate_planned_symlink(entry, entries_by_path, capability)

        for entry in planned:
            if entry.source_kind == "delete":
                _delete_overlay_entry(capability, entry.change.path)
            elif entry.source_kind == "regular":
                _copy_overlay_entry(real_root, capability, entry)

        applied_symlinks: set[str] = set()
        pending = [entry for entry in planned if entry.source_kind == "symlink"]
        while pending:
            ready = [
                entry
                for entry in pending
                if entry.link_target not in entries_by_path
                or entries_by_path[entry.link_target].source_kind == "regular"
                or entry.link_target in applied_symlinks
            ]
            if not ready:
                raise RuntimeError("planned symlink target cycle")
            for entry in ready:
                _copy_overlay_entry(real_root, capability, entry)
                applied_symlinks.add(entry.change.path)
                pending.remove(entry)


def _cleanup_detached_worktree(capability: _RootCapability, real_root: Path) -> None:
    """Deregister and remove only the exact pinned disposable worktree."""
    try:
        _assert_root_identity(capability)
    except RuntimeError as exc:
        raise RuntimeError(f"orphaned worktree: {capability.path}: workspace identity mismatch") from exc
    try:
        remove_result = subprocess.run(
            ["git", "worktree", "remove", "--force", str(capability.path)],
            cwd=real_root,
            capture_output=True,
            check=False,
        )
    except (OSError, KeyboardInterrupt) as exc:
        raise RuntimeError(f"orphaned worktree: {capability.path}: cleanup command failed: {exc}") from exc
    if remove_result.returncode != 0:
        stderr = remove_result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"orphaned worktree: {capability.path}: {stderr}")
    try:
        _remove_pinned_root_if_present(capability)
    except RuntimeError as exc:
        raise RuntimeError(f"orphaned worktree: {capability.path}: workspace identity mismatch") from exc
    except OSError as exc:
        raise RuntimeError(f"orphaned worktree: {capability.path}: direct cleanup failed: {exc}") from exc


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
        root_capability = _open_root_capability(sweep_root)
        try:
            try:
                add_result = subprocess.run(
                    ["git", "worktree", "add", "--detach", str(sweep_root), head],
                    cwd=resolved_real_root,
                    capture_output=True,
                    check=False,
                )
            except OSError as exc:
                raise RuntimeError(f"git worktree add failed: {exc}") from exc
            _assert_root_identity(root_capability)
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
                root_capability,
                overlay_changes,
                excluded_outputs,
                registered_roots,
            )

            _assert_root_identity(root_capability)
            sweep_status_result = subprocess.run(
                ["git", "-C", str(sweep_root), "status", "--porcelain=v1", "-z", "--untracked-files=all"],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
            _assert_root_identity(root_capability)
            if sweep_status_result.returncode != 0:
                raise RuntimeError(sweep_status_result.stderr.decode("utf-8", errors="replace").strip())
            try:
                baseline_manifest = _parse_porcelain_with_capability(sweep_status_result.stdout, root_capability)
            except ValueError as exc:
                raise RuntimeError(f"malformed porcelain rename or status: {exc}") from exc
            baseline_digests = _capture_baseline_digests(root_capability, overlay_changes, targets)

            yield MutationWorkspace(
                real_root=resolved_real_root,
                sweep_root=sweep_root,
                head=head,
                overlay_changes=overlay_changes,
                baseline_manifest=baseline_manifest,
                baseline_status=sweep_status_result.stdout,
                baseline_digests=baseline_digests,
                _root_capability=root_capability,
            )
        finally:
            try:
                _cleanup_detached_worktree(root_capability, resolved_real_root)
            finally:
                os.close(root_capability.descriptor)

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
