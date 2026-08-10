from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace
import mutation_workspace_fs

pytestmark = pytest.mark.unit


def _git(repo: Path, *arguments: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(["git", "-C", str(repo), *arguments], check=True, capture_output=True)


def _repository(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "capability-test@example.invalid")
    _git(repo, "config", "user.name", "Capability Test")
    (repo / "tracked.txt").write_text("tracked\n", encoding="utf-8")
    _git(repo, "add", "tracked.txt")
    _git(repo, "commit", "-qm", "base")
    return repo


def _replacement_directory(public_root: Path, displaced: Path) -> Path:
    public_root.rename(displaced)
    public_root.mkdir()
    sentinel = public_root / "sentinel.txt"
    sentinel.write_text("replacement\n", encoding="utf-8")
    return sentinel


def _replacement_symlink(public_root: Path, displaced: Path, external: Path) -> Path:
    public_root.rename(displaced)
    external.mkdir()
    sentinel = external / "sentinel.txt"
    sentinel.write_text("external replacement\n", encoding="utf-8")
    public_root.symlink_to(external, target_is_directory=True)
    return sentinel


def test_real_git_add_uses_the_opened_root_when_its_public_name_is_replaced(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = _repository(tmp_path)
    public_root = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    real_run = subprocess.run
    raced = False
    inherited_descriptor = False
    sentinel: Path | None = None

    def make_root(**_kwargs: object) -> str:
        public_root.mkdir()
        return str(public_root)

    def racing_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[bytes]:
        nonlocal inherited_descriptor, raced, sentinel
        if command[1:3] == ["worktree", "add"] and not raced:
            raced = True
            sentinel = _replacement_directory(public_root, displaced)
            pass_fds = kwargs.get("pass_fds", ())
            inherited_descriptor = bool(pass_fds) and command[-2] == f"/proc/self/fd/{pass_fds[0]}"
        return real_run(command, **kwargs)

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(mutation_workspace.subprocess, "run", racing_run)

    with (
        pytest.raises(RuntimeError, match="workspace identity mismatch"),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pytest.fail("a replaced public root must not be yielded")

    assert raced is True
    assert inherited_descriptor is True
    assert sentinel is not None
    assert sentinel.read_text(encoding="utf-8") == "replacement\n"
    assert not displaced.exists()
    assert str(displaced) not in _git(repo, "worktree", "list", "--porcelain").stdout.decode()


@pytest.mark.parametrize("replacement_kind", ["directory", "symlink"])
def test_real_git_cleanup_removes_the_opened_root_without_touching_a_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    replacement_kind: str,
) -> None:
    repo = _repository(tmp_path)
    public_root = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    external = tmp_path / "external"
    real_run = subprocess.run
    raced = False
    inherited_descriptor = False
    sentinel: Path | None = None

    def make_root(**_kwargs: object) -> str:
        public_root.mkdir()
        return str(public_root)

    def racing_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[bytes]:
        nonlocal inherited_descriptor, raced, sentinel
        if command[1:3] in (["worktree", "repair"], ["worktree", "remove"]) and not raced:
            raced = True
            if replacement_kind == "directory":
                sentinel = _replacement_directory(public_root, displaced)
            else:
                sentinel = _replacement_symlink(public_root, displaced, external)
            pass_fds = kwargs.get("pass_fds", ())
            inherited_descriptor = bool(pass_fds) and command[-1] == f"/proc/self/fd/{pass_fds[0]}"
        return real_run(command, **kwargs)

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(mutation_workspace.subprocess, "run", racing_run)

    with (
        pytest.raises(RuntimeError, match="workspace identity mismatch"),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pass

    assert raced is True
    assert inherited_descriptor is True
    assert sentinel is not None
    expected = "replacement\n" if replacement_kind == "directory" else "external replacement\n"
    assert sentinel.read_text(encoding="utf-8") == expected
    assert not displaced.exists()
    if replacement_kind == "symlink":
        assert public_root.is_symlink()
        assert os.readlink(public_root) == str(external)
    assert str(displaced) not in _git(repo, "worktree", "list", "--porcelain").stdout.decode()


def test_missing_proc_fd_capability_fails_closed_before_git_add(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = _repository(tmp_path)
    public_root = tmp_path / "sweep"
    calls: list[tuple[str, ...]] = []
    real_run = subprocess.run

    def make_root(**_kwargs: object) -> str:
        public_root.mkdir()
        return str(public_root)

    def recording_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[bytes]:
        calls.append(tuple(command))
        return real_run(command, **kwargs)

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(mutation_workspace.subprocess, "run", recording_run)
    monkeypatch.setattr(mutation_workspace_fs, "PROC_SELF_FD", tmp_path / "missing", raising=False)

    with (
        pytest.raises(RuntimeError, match="proc/self/fd capability is unavailable"),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pytest.fail("missing proc-fd support must not yield")

    assert not any(command[1:3] == ("worktree", "add") for command in calls)
    assert not public_root.exists()


@pytest.mark.parametrize(
    "failing_call",
    [
        ("git", "rev-parse", "--verify", "HEAD^{commit}"),
        ("git", "worktree", "add"),
        ("git", "-C", "CAPABILITY", "rev-parse"),
        ("git", "status"),
        ("git", "worktree", "list"),
        ("git", "-C", "CAPABILITY", "status"),
    ],
    ids=["head", "add", "disposable-head", "source-status", "worktree-list", "disposable-status"],
)
def test_every_setup_git_spawn_oserror_is_normalized(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failing_call: tuple[str, ...],
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    sweep = tmp_path / "sweep"
    head = "a" * 40

    def make_root(**_kwargs: object) -> str:
        sweep.mkdir()
        return str(sweep)

    def matches(command: list[str]) -> bool:
        normalized = tuple("CAPABILITY" if part.startswith("/proc/self/fd/") else part for part in command)
        return normalized[: len(failing_call)] == failing_call

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        if matches(command):
            raise OSError("spawn blocked")
        if command == ["git", "rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if command[1:3] == ["worktree", "add"]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command[-3:] == ["rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if command[1:3] == ["worktree", "repair"]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command[1:3] == ["worktree", "remove"]:
            if sweep.exists():
                shutil.rmtree(sweep)
            return subprocess.CompletedProcess(command, 0, b"", b"")
        return subprocess.CompletedProcess(command, 0, b"", b"")

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)

    with (
        pytest.raises(RuntimeError, match="git .* failed: spawn blocked"),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pytest.fail("a failed Git spawn must not yield")


def test_root_capability_open_failure_is_normalized_and_cleans_the_exact_empty_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    sweep = tmp_path / "sweep"
    head = "a" * 40

    def make_root(**_kwargs: object) -> str:
        sweep.mkdir()
        return str(sweep)

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        assert command == ["git", "rev-parse", "--verify", "HEAD^{commit}"]
        return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)
    monkeypatch.setattr(
        mutation_workspace, "_open_root_capability", lambda _root: (_ for _ in ()).throw(OSError("no fd"))
    )

    with (
        pytest.raises(RuntimeError, match="workspace capability open failed: no fd"),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pytest.fail("a workspace without a root capability must not yield")

    assert not sweep.exists()


def test_direct_cleanup_keyboard_interrupt_is_an_orphan_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    capability_root = tmp_path / "sweep"
    capability_root.mkdir()
    capability = mutation_workspace_fs.open_root_capability(capability_root)
    monkeypatch.setattr(
        mutation_workspace,
        "_remove_pinned_root_if_present",
        lambda _capability: (_ for _ in ()).throw(KeyboardInterrupt("interrupted")),
    )
    monkeypatch.setattr(
        mutation_workspace,
        "_run_git",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 0, b"", b""),
        raising=False,
    )

    try:
        with pytest.raises(RuntimeError, match="direct cleanup failed: interrupted"):
            mutation_workspace._cleanup_detached_worktree(capability, tmp_path)
        assert capability_root.exists()
    finally:
        os.close(capability.descriptor)


def test_cleanup_capability_lookup_keyboard_interrupt_is_an_orphan_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    sweep = tmp_path / "sweep"
    head = "a" * 40
    original_git_capability_path = mutation_workspace._git_capability_path
    capability_lookups = 0

    def make_root(**_kwargs: object) -> str:
        sweep.mkdir()
        return str(sweep)

    def interrupt_second_lookup(capability: mutation_workspace_fs.RootCapability) -> Path:
        nonlocal capability_lookups
        capability_lookups += 1
        if capability_lookups == 2:
            raise KeyboardInterrupt("interrupted")
        return original_git_capability_path(capability)

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        if command == ["git", "rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if command[-3:] == ["rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        return subprocess.CompletedProcess(command, 0, b"", b"")

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)
    monkeypatch.setattr(mutation_workspace, "_git_capability_path", interrupt_second_lookup)

    with (
        pytest.raises(
            RuntimeError,
            match=rf"orphaned worktree: {sweep}: cleanup command failed: interrupted",
        ),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pass

    assert capability_lookups == 2
    assert sweep.exists()
    shutil.rmtree(sweep)


@pytest.mark.parametrize(
    "cleanup_exception",
    [OSError("blocked"), RuntimeError("changed"), KeyboardInterrupt("interrupted")],
    ids=["oserror", "runtime-error", "keyboard-interrupt"],
)
def test_pre_add_cleanup_failure_is_an_orphan_runtime_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cleanup_exception: BaseException,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    sweep = tmp_path / "sweep"
    head = "a" * 40

    def make_root(**_kwargs: object) -> str:
        sweep.mkdir()
        return str(sweep)

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(
        mutation_workspace.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b""),
    )
    monkeypatch.setattr(
        mutation_workspace,
        "_git_capability_path",
        lambda _capability: (_ for _ in ()).throw(RuntimeError("proc blocked")),
    )
    monkeypatch.setattr(
        mutation_workspace,
        "_remove_pinned_root_if_present",
        lambda _capability: (_ for _ in ()).throw(cleanup_exception),
    )

    with (
        pytest.raises(
            RuntimeError,
            match=rf"orphaned workspace root: {sweep}: direct cleanup failed: {cleanup_exception}",
        ),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pytest.fail("a workspace without a Git capability must not yield")

    assert sweep.exists()
    shutil.rmtree(sweep)


def test_proc_fd_capability_rejects_a_path_to_a_different_inode(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "root"
    root.mkdir()
    other = tmp_path / "other"
    other.mkdir()
    capability = mutation_workspace_fs.open_root_capability(root)
    fake_proc = tmp_path / "proc-fd"
    fake_proc.mkdir()
    (fake_proc / str(capability.descriptor)).symlink_to(other, target_is_directory=True)
    monkeypatch.setattr(mutation_workspace_fs, "PROC_SELF_FD", fake_proc)

    try:
        with pytest.raises(RuntimeError, match="does not name the workspace root"):
            mutation_workspace_fs.git_capability_path(capability)
    finally:
        os.close(capability.descriptor)


def test_open_failure_cleanup_refuses_a_replacement_root_without_touching_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    sweep = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    head = "a" * 40
    sentinel: Path | None = None

    def make_root(**_kwargs: object) -> str:
        sweep.mkdir()
        return str(sweep)

    def replace_then_fail(_root: Path) -> mutation_workspace_fs.RootCapability:
        nonlocal sentinel
        sentinel = _replacement_directory(sweep, displaced)
        raise OSError("no fd")

    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", make_root)
    monkeypatch.setattr(
        mutation_workspace.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b""),
    )
    monkeypatch.setattr(mutation_workspace, "_open_root_capability", replace_then_fail)

    with (
        pytest.raises(RuntimeError, match="orphaned workspace root.*changed during removal"),
        mutation_workspace.detached_head_workspace(repo, (), ()),
    ):
        pytest.fail("a substituted root must not yield")

    assert sentinel is not None
    assert sentinel.read_text(encoding="utf-8") == "replacement\n"
    assert displaced.exists()


@pytest.mark.parametrize(
    ("repair_result", "expected"),
    [
        (OSError("repair spawn blocked"), "git worktree repair failed: repair spawn blocked"),
        (subprocess.CompletedProcess([], 1, b"", b"repair refused"), "repair refused"),
    ],
    ids=["spawn-oserror", "nonzero"],
)
def test_cleanup_repair_failure_retains_the_exact_opened_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    repair_result: OSError | subprocess.CompletedProcess[bytes],
    expected: str,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    root = tmp_path / "sweep"
    root.mkdir()
    capability = mutation_workspace_fs.open_root_capability(root)
    calls: list[tuple[str, ...]] = []

    def fake_run(command: list[str], **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        calls.append(tuple(command))
        if isinstance(repair_result, OSError):
            raise repair_result
        return repair_result

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)

    try:
        with pytest.raises(RuntimeError, match=expected):
            mutation_workspace._cleanup_detached_worktree(capability, repo)
        assert root.exists()
        assert len(calls) == 1
        assert calls[0][1:3] == ("worktree", "repair")
    finally:
        os.close(capability.descriptor)


def test_cleanup_fails_closed_if_proc_fd_disappears_after_workspace_creation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = _repository(tmp_path)
    missing_proc = tmp_path / "missing-proc"

    with (
        pytest.raises(RuntimeError, match="orphaned worktree.*proc/self/fd capability is unavailable"),
        mutation_workspace.detached_head_workspace(repo, (), ()) as workspace,
    ):
        monkeypatch.setattr(mutation_workspace_fs, "PROC_SELF_FD", missing_proc)

    assert workspace.sweep_root.exists()


def test_verify_baseline_git_spawn_oserror_is_normalized(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo = _repository(tmp_path)

    with (
        mutation_workspace.detached_head_workspace(repo, (), ()) as workspace,
        monkeypatch.context() as scoped_patch,
    ):
        scoped_patch.setattr(
            mutation_workspace.subprocess,
            "run",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("status spawn blocked")),
        )
        with pytest.raises(RuntimeError, match="git disposable status failed: status spawn blocked"):
            workspace.verify_baseline()
