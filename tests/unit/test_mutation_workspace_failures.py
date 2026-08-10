from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace

pytestmark = pytest.mark.unit


def _public_command(command: list[str]) -> tuple[str, ...]:
    return tuple(str(Path(part).resolve()) if part.startswith("/proc/self/fd/") else part for part in command)


def _install_lifecycle_double(
    monkeypatch: pytest.MonkeyPatch,
    real: Path,
    sweep: Path,
    *,
    add_exception: OSError | KeyboardInterrupt | None = None,
    remove_exception: OSError | KeyboardInterrupt | None = None,
    remove_returncode: int = 0,
    remove_stderr: bytes = b"",
) -> list[tuple[str, ...]]:
    head = "a" * 40
    sweep.mkdir()
    calls: list[tuple[str, ...]] = []
    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", lambda **_kwargs: str(sweep))

    def fake_run(
        command: list[str],
        *,
        cwd: Path,
        **_kwargs: object,
    ) -> subprocess.CompletedProcess[bytes]:
        assert Path(cwd) == real
        call = _public_command(command)
        calls.append(call)
        if command[1:3] == ["worktree", "add"] and add_exception is not None:
            raise add_exception
        if call[1:3] == ("worktree", "remove"):
            if remove_exception is not None:
                raise remove_exception
            if remove_returncode == 0 and sweep.exists():
                shutil.rmtree(sweep)
            return subprocess.CompletedProcess(command, remove_returncode, b"", remove_stderr)
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
        if call == ("git", "worktree", "repair", str(sweep)):
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command == ["git", "worktree", "list", "--porcelain", "-z"]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        raise AssertionError(f"unexpected subprocess call: {command}")

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)
    return calls


@pytest.mark.parametrize(
    ("add_exception", "expected_type", "expected_message"),
    [
        (OSError("add exploded"), RuntimeError, "git worktree add failed: add exploded"),
        (KeyboardInterrupt("add interrupted"), KeyboardInterrupt, "add interrupted"),
    ],
    ids=["oserror", "keyboard-interrupt"],
)
def test_add_subprocess_exception_still_forces_exact_cleanup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    add_exception: OSError | KeyboardInterrupt,
    expected_type: type[BaseException],
    expected_message: str,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    calls = _install_lifecycle_double(monkeypatch, real, sweep, add_exception=add_exception)

    with (
        pytest.raises(expected_type, match=expected_message),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("a failed add must not yield")

    assert calls == [
        ("git", "rev-parse", "--verify", "HEAD^{commit}"),
        ("git", "worktree", "add", "--detach", str(sweep), "a" * 40),
        ("git", "worktree", "repair", str(sweep)),
        ("git", "worktree", "remove", "--force", str(sweep)),
    ]
    assert not sweep.exists()


def test_add_subprocess_exception_reports_orphan_when_cleanup_is_uncertain(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    calls = _install_lifecycle_double(
        monkeypatch,
        real,
        sweep,
        add_exception=OSError("add exploded"),
        remove_returncode=1,
        remove_stderr=b"not registered",
    )

    with (
        pytest.raises(RuntimeError) as exc_info,
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("an uncertain add must not yield")

    assert str(exc_info.value) == f"orphaned worktree: {sweep}: not registered"
    assert calls[-1] == ("git", "worktree", "remove", "--force", str(sweep))
    assert sweep.exists()


@pytest.mark.parametrize(
    "cleanup_exception",
    [OSError("remove exploded"), KeyboardInterrupt("remove interrupted")],
    ids=["oserror", "keyboard-interrupt"],
)
def test_cleanup_subprocess_exception_becomes_exact_orphan_runtime_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cleanup_exception: OSError | KeyboardInterrupt,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    calls = _install_lifecycle_double(monkeypatch, real, sweep, remove_exception=cleanup_exception)

    with pytest.raises(RuntimeError) as exc_info, mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    normalized = (
        f"git worktree remove failed: {cleanup_exception}"
        if isinstance(cleanup_exception, OSError)
        else str(cleanup_exception)
    )
    assert str(exc_info.value) == f"orphaned worktree: {sweep}: cleanup command failed: {normalized}"
    assert calls[-1] == ("git", "worktree", "remove", "--force", str(sweep))
    assert sweep.exists()


def test_parent_symlink_swap_during_overlay_never_touches_external_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    source = real / "parent/victim.txt"
    source.parent.mkdir(parents=True)
    source.write_text("overlay\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    (sweep / "parent").mkdir(parents=True)
    external = tmp_path / "external"
    external.mkdir()
    external_file = external / "victim.txt"
    external_file.write_text("sentinel\n", encoding="utf-8")
    original_open = mutation_workspace.os.open
    swapped = False

    def racing_open(
        path: str | bytes | os.PathLike[str] | os.PathLike[bytes],
        flags: int,
        mode: int = 0o777,
        *,
        dir_fd: int | None = None,
    ) -> int:
        nonlocal swapped
        if path == "parent" and dir_fd is not None and not swapped:
            opened_directory = Path(f"/proc/self/fd/{dir_fd}").resolve()
            if opened_directory == sweep:
                (sweep / "parent").rename(sweep / "original-parent")
                (sweep / "parent").symlink_to(external, target_is_directory=True)
                swapped = True
        return original_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(mutation_workspace.os, "open", racing_open)

    with pytest.raises(RuntimeError, match="unsafe workspace destination"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (mutation_workspace.OverlayChange("copy", "parent/victim.txt"),),
            (),
            (),
        )

    assert swapped is True
    assert external_file.read_text(encoding="utf-8") == "sentinel\n"


def test_overlay_delete_unlinks_external_final_symlink_without_following_it(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    external = tmp_path / "external.txt"
    external.write_text("sentinel\n", encoding="utf-8")
    destination = sweep / "deleted.txt"
    destination.symlink_to(external)

    mutation_workspace._apply_overlay_changes(
        real,
        sweep,
        (mutation_workspace.OverlayChange("delete", "deleted.txt"),),
        (),
        (),
    )

    assert not os.path.lexists(destination)
    assert external.read_text(encoding="utf-8") == "sentinel\n"


def test_overlay_copy_replaces_dangling_final_symlink_without_following_it(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    source = real / "changed.txt"
    source.write_text("overlay\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    missing_external = tmp_path / "missing-external.txt"
    destination = sweep / "changed.txt"
    destination.symlink_to(missing_external)

    mutation_workspace._apply_overlay_changes(
        real,
        sweep,
        (mutation_workspace.OverlayChange("copy", "changed.txt"),),
        (),
        (),
    )

    assert not destination.is_symlink()
    assert destination.read_text(encoding="utf-8") == "overlay\n"
    assert not missing_external.exists()


def _baseline_workspace(real: Path, sweep: Path) -> mutation_workspace.MutationWorkspace:
    return mutation_workspace.MutationWorkspace(
        real_root=real,
        sweep_root=sweep,
        head="a" * 40,
        overlay_changes=(),
        baseline_manifest=(),
        baseline_status=b"",
        baseline_digests={},
    )


@pytest.mark.parametrize(
    ("returncode", "stdout", "stderr", "expected"),
    [
        (1, b"", b"status failed", "status failed"),
        (0, b"R  renamed.py\0", b"", "malformed porcelain rename"),
        (0, b"?? changed.py\0", b"", "workspace status mismatch"),
    ],
    ids=["command-failure", "malformed-status", "semantic-mismatch"],
)
def test_verify_baseline_reaches_and_rejects_semantic_status_failures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    returncode: int,
    stdout: bytes,
    stderr: bytes,
    expected: str,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    (sweep / "changed.py").write_text("changed\n", encoding="utf-8")
    workspace = _baseline_workspace(real, sweep)
    monkeypatch.setattr(
        mutation_workspace.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], returncode, stdout, stderr),
    )

    with pytest.raises(RuntimeError, match=expected):
        workspace.verify_baseline()


@pytest.mark.parametrize("relative", ["../outside.py", "missing.py"])
def test_target_path_rejects_unsafe_or_missing_targets(tmp_path: Path, relative: str) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    workspace = _baseline_workspace(real, sweep)

    with pytest.raises(RuntimeError, match="unsafe workspace path|workspace path does not exist"):
        workspace.target_path(relative)


@pytest.mark.parametrize("kind", ["missing", "directory", "symlink"])
def test_baseline_capture_rejects_missing_or_nonregular_selected_targets(tmp_path: Path, kind: str) -> None:
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    relative = "selected.py"
    selected = sweep / relative
    if kind == "directory":
        selected.mkdir()
    elif kind == "symlink":
        target = sweep / "target.py"
        target.write_text("target\n", encoding="utf-8")
        selected.symlink_to(target)

    with pytest.raises(RuntimeError, match="workspace path does not exist|mutation target is not a regular file"):
        mutation_workspace._capture_baseline_digests(sweep, (), (relative,))


@pytest.mark.parametrize(
    "payload",
    [
        b"HEAD aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\0\0",
        b"worktree /tmp/first\0HEAD aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\0",
        (
            b"worktree /tmp/first\0HEAD aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\0"
            b"worktree /tmp/second\0HEAD bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\0\0"
        ),
    ],
    ids=["missing-worktree-header", "missing-record-separator", "embedded-worktree-header"],
)
def test_registered_worktree_roots_rejects_malformed_records(payload: bytes) -> None:
    with pytest.raises(RuntimeError, match="malformed worktree list"):
        mutation_workspace._registered_worktree_roots(payload)
