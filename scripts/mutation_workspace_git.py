from __future__ import annotations

import os
import subprocess
from collections.abc import Callable
from pathlib import Path

from mutation_workspace_fs import RootCapability, git_capability_path


def registered_worktree_roots(payload: bytes) -> tuple[Path, ...]:
    """Extract absolute worktree roots from Git's NUL-delimited porcelain output.

    Args:
        payload: Raw ``git worktree list --porcelain -z`` output.

    Returns:
        Absolute canonical worktree roots in Git's record order.

    Raises:
        RuntimeError: If record framing or a root field is malformed.
    """
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


GitRunner = Callable[..., subprocess.CompletedProcess[bytes]]


def run_git(
    command: list[str],
    cwd: Path,
    operation: str,
    capability: RootCapability | None,
    runner: GitRunner,
) -> subprocess.CompletedProcess[bytes]:
    """Run one Git boundary and normalize process-spawn failures.

    Args:
        command: Complete Git argument array.
        cwd: Repository root from which Git must run.
        operation: Stable diagnostic label for this boundary.
        capability: Optional disposable-root descriptor to inherit.
        runner: Subprocess-compatible runner used at the process boundary.

    Returns:
        The completed Git process without return-code enforcement.

    Raises:
        RuntimeError: If the capability is unavailable or Git cannot be spawned.
    """
    pass_fds = () if capability is None else (capability.descriptor,)
    if capability is not None:
        git_capability_path(capability)
    try:
        return runner(command, cwd=cwd, capture_output=True, check=False, pass_fds=pass_fds)
    except OSError as exc:
        raise RuntimeError(f"git {operation} failed: {exc}") from exc
