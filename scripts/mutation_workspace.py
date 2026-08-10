from __future__ import annotations

import contextlib
import hashlib
import os
import shutil
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
    if payload and not payload.endswith(b"\0"):
        raise RuntimeError("malformed worktree list: missing NUL terminator")
    roots: list[Path] = []
    for field in payload.split(b"\0"):
        if not field.startswith(b"worktree "):
            continue
        root = Path(os.fsdecode(field.removeprefix(b"worktree ")))
        if not root.is_absolute():
            raise RuntimeError(f"malformed worktree list: non-absolute root: {root}")
        roots.append(root.resolve(strict=False))
    return tuple(roots)


def _replace_with_overlay_copy(source: Path, destination: Path, sweep_root: Path) -> None:
    """Replace one disposable path with a confined regular file or symlink."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or destination.is_file():
        destination.unlink()
    elif destination.exists():
        shutil.rmtree(destination)

    if source.is_symlink():
        link_text = source.readlink()
        proposed_target = link_text if link_text.is_absolute() else destination.parent / link_text
        resolved_target = proposed_target.resolve(strict=False)
        resolved_sweep = sweep_root.resolve(strict=True)
        if not resolved_target.is_relative_to(resolved_sweep):
            raise RuntimeError(f"workspace symlink escapes workspace root: {destination.relative_to(resolved_sweep)}")
        if ".git" in resolved_target.relative_to(resolved_sweep).parts:
            raise RuntimeError(f"unsafe workspace path: {destination.relative_to(resolved_sweep)}")
        destination.symlink_to(link_text)
        return
    if not source.is_file():
        raise RuntimeError(f"overlay source is not a regular file: {source}")
    shutil.copy2(source, destination)


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
            destination = confined_path(sweep_root, change.path, must_exist=False)
        except ValueError as exc:
            raise RuntimeError(str(exc)) from exc
        resolved_source = source.resolve(strict=False)
        if any(resolved_source == output or resolved_source.is_relative_to(output) for output in resolved_outputs):
            raise RuntimeError(f"excluded output cannot be overlaid: {change.path}")
        if any(resolved_source == root or resolved_source.is_relative_to(root) for root in nested_roots):
            raise RuntimeError(f"nested registered worktree cannot be overlaid: {change.path}")

        if change.action == "copy":
            _replace_with_overlay_copy(source, destination, sweep_root)
        elif destination.is_symlink() or destination.is_file():
            destination.unlink()
        elif destination.exists():
            shutil.rmtree(destination)


def detached_head_workspace(
    real_root: Path,
    targets: tuple[str, ...],
    excluded_outputs: tuple[Path, ...],
) -> contextlib.AbstractContextManager[MutationWorkspace]:
    """Create a detached mutation workspace and dispose it on context exit."""

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
        add_result = subprocess.run(
            ["git", "worktree", "add", "--detach", str(sweep_root), head],
            cwd=resolved_real_root,
            capture_output=True,
            check=False,
        )
        try:
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
            remove_result = subprocess.run(
                ["git", "worktree", "remove", "--force", str(sweep_root)],
                cwd=resolved_real_root,
                capture_output=True,
                check=False,
            )
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
    candidate = Path(relative)
    if (
        relative in {"", "."}
        or candidate.is_absolute()
        or candidate.as_posix() != relative
        or ".." in candidate.parts
        or ".git" in candidate.parts
    ):
        raise ValueError(f"unsafe workspace path: {relative}")
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
