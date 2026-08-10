from __future__ import annotations

import os
import shutil
import subprocess
import sys
from collections import Counter
from collections.abc import Callable
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace

pytestmark = pytest.mark.unit

GitCall = tuple[str, ...]


def _public_command(command: list[str]) -> tuple[str, ...]:
    return tuple(str(Path(part).resolve()) if part.startswith("/proc/self/fd/") else part for part in command)


def _install_identity_lifecycle(
    monkeypatch: pytest.MonkeyPatch,
    real: Path,
    sweep: Path,
    *,
    on_remove: Callable[[], None] | None = None,
) -> Counter[GitCall]:
    head = "a" * 40
    sweep.mkdir()
    seen: Counter[GitCall] = Counter()
    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", lambda **_kwargs: str(sweep))

    def fake_run(
        command: list[str],
        *,
        cwd: Path,
        **_kwargs: object,
    ) -> subprocess.CompletedProcess[bytes]:
        assert Path(cwd) == real
        call = _public_command(command)
        seen[call] += 1
        if command == ["git", "rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if call == ("git", "worktree", "add", "--detach", str(sweep), head):
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if call == ("git", "-C", str(sweep), "rev-parse", "--verify", "HEAD^{commit}"):
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if call in (
            ("git", "status", "--porcelain=v1", "-z", "--untracked-files=all"),
            ("git", "-C", str(sweep), "status", "--porcelain=v1", "-z", "--untracked-files=all"),
        ):
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command == ["git", "worktree", "list", "--porcelain", "-z"]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if call[1:3] == ("worktree", "repair"):
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if call[1:3] == ("worktree", "remove"):
            if on_remove is not None:
                on_remove()
            else:
                opened_root = Path(command[-1]).resolve()
                if opened_root.exists():
                    shutil.rmtree(opened_root)
            return subprocess.CompletedProcess(command, 0, b"", b"")
        raise AssertionError(f"unexpected subprocess call: {command}")

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)
    return seen


def _substitute_root(sweep: Path, displaced: Path, replacement_text: str = "replacement\n") -> Path:
    sweep.rename(displaced)
    sweep.mkdir()
    sentinel = sweep / "sentinel.txt"
    sentinel.write_text(replacement_text, encoding="utf-8")
    return sentinel


def test_cleanup_rejects_an_ordinary_directory_root_substitution_without_touching_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    seen = _install_identity_lifecycle(monkeypatch, real, sweep)

    with pytest.raises(RuntimeError) as exc_info, mutation_workspace.detached_head_workspace(real, (), ()):
        sentinel = _substitute_root(sweep, displaced)

    assert str(exc_info.value) == f"orphaned worktree: {sweep}: workspace identity mismatch"
    assert sentinel.read_text(encoding="utf-8") == "replacement\n"
    assert not displaced.exists()
    assert seen[("git", "worktree", "remove", "--force", str(displaced))] == 1


def test_verify_baseline_rejects_root_substitution_before_reading_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    replacement = tmp_path / "replacement"
    _install_identity_lifecycle(monkeypatch, real, sweep)

    with mutation_workspace.detached_head_workspace(real, (), ()) as workspace:
        sentinel = _substitute_root(sweep, displaced)
        with pytest.raises(RuntimeError, match="workspace root identity mismatch"):
            workspace.verify_baseline()
        sweep.rename(replacement)
        displaced.rename(sweep)

    assert (replacement / sentinel.name).read_text(encoding="utf-8") == "replacement\n"


def test_target_path_rejects_root_substitution_before_returning_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    replacement = tmp_path / "replacement"
    _install_identity_lifecycle(monkeypatch, real, sweep)

    with mutation_workspace.detached_head_workspace(real, (), ()) as workspace:
        original_target = sweep / "selected.py"
        original_target.write_text("original\n", encoding="utf-8")
        _substitute_root(sweep, displaced)
        (sweep / "selected.py").write_text("replacement\n", encoding="utf-8")
        with pytest.raises(RuntimeError, match="workspace root identity mismatch"):
            workspace.target_path("selected.py")
        sweep.rename(replacement)
        displaced.rename(sweep)


def test_post_git_cleanup_rejects_a_substituted_directory_without_deleting_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    sentinel: Path | None = None

    def substitute_on_remove() -> None:
        nonlocal sentinel
        sentinel = _substitute_root(sweep, displaced, "post-remove replacement\n")

    _install_identity_lifecycle(monkeypatch, real, sweep, on_remove=substitute_on_remove)

    with pytest.raises(RuntimeError) as exc_info, mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    assert str(exc_info.value) == f"orphaned worktree: {sweep}: workspace identity mismatch"
    assert sentinel is not None
    assert sentinel.read_text(encoding="utf-8") == "post-remove replacement\n"


def test_post_git_cleanup_rejects_a_dangling_symlink_without_deleting_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    displaced = tmp_path / "displaced"
    missing = tmp_path / "missing"

    def substitute_on_remove() -> None:
        sweep.rename(displaced)
        sweep.symlink_to(missing, target_is_directory=True)

    _install_identity_lifecycle(monkeypatch, real, sweep, on_remove=substitute_on_remove)

    with pytest.raises(RuntimeError) as exc_info, mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    assert str(exc_info.value) == f"orphaned worktree: {sweep}: workspace identity mismatch"
    assert sweep.is_symlink()
    assert os.readlink(sweep) == str(missing)


def test_post_git_direct_cleanup_error_is_an_exact_orphan_runtime_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    _install_identity_lifecycle(monkeypatch, real, sweep, on_remove=lambda: None)
    injected = False

    def fail_direct_cleanup(_capability: object) -> None:
        nonlocal injected
        injected = True
        raise OSError("blocked")

    monkeypatch.setattr(mutation_workspace, "_remove_pinned_root_if_present", fail_direct_cleanup)

    with pytest.raises(RuntimeError) as exc_info, mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    assert injected is True
    assert str(exc_info.value) == f"orphaned worktree: {sweep}: direct cleanup failed: blocked"
    assert sweep.exists()
