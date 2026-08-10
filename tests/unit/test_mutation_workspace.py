from __future__ import annotations

import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from unittest import mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace

pytestmark = pytest.mark.unit

GitCall = tuple[tuple[str, ...], Path]


def _record_probe(probes: list[object], path: object, result: bool) -> bool:
    probes.append(path)
    return result


def _install_workspace_git_double(
    monkeypatch: pytest.MonkeyPatch,
    real: Path,
    sweep: Path,
    *,
    disposable_head: str = "a" * 40,
    real_status: bytes = b"",
    sweep_status: bytes = b"",
    worktrees: bytes | None = None,
    remove_directory: bool = True,
) -> tuple[Counter[GitCall], dict[GitCall, tuple[int, bytes, bytes]]]:
    head = "a" * 40
    sweep.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", lambda **_kwargs: str(sweep))
    if worktrees is None:
        worktrees = (f"worktree {real}\0HEAD {head}\0\0worktree {sweep}\0HEAD {head}\0detached\0\0").encode()
    responses: dict[GitCall, tuple[int, bytes, bytes]] = {
        (("git", "rev-parse", "--verify", "HEAD^{commit}"), real): (0, f"{head}\n".encode(), b""),
        (("git", "worktree", "add", "--detach", str(sweep), head), real): (0, b"", b""),
        (("git", "-C", str(sweep), "rev-parse", "--verify", "HEAD^{commit}"), real): (
            0,
            f"{disposable_head}\n".encode(),
            b"",
        ),
        (("git", "status", "--porcelain=v1", "-z", "--untracked-files=all"), real): (
            0,
            real_status,
            b"",
        ),
        (("git", "worktree", "list", "--porcelain", "-z"), real): (0, worktrees, b""),
        (("git", "-C", str(sweep), "status", "--porcelain=v1", "-z", "--untracked-files=all"), real): (
            0,
            sweep_status,
            b"",
        ),
        (("git", "worktree", "remove", "--force", str(sweep)), real): (0, b"", b""),
    }
    seen: Counter[GitCall] = Counter()

    def fake_run(command: list[str], *, cwd: Path, **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        key = (tuple(command), Path(cwd))
        assert key in responses, f"unexpected subprocess call: {key}"
        seen[key] += 1
        returncode, stdout, stderr = responses[key]
        if remove_directory and command[1:3] == ["worktree", "remove"] and returncode == 0 and sweep.exists():
            shutil.rmtree(sweep)
        return subprocess.CompletedProcess(command, returncode, stdout, stderr)

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)
    return seen, responses


def test_workspace_adds_the_captured_head_detached_and_removes_exact_path(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    seen, responses = _install_workspace_git_double(monkeypatch, real, sweep)
    with mutation_workspace.detached_head_workspace(real, (), ()) as workspace:
        assert workspace.sweep_root == sweep
        assert workspace.head == "a" * 40
        assert workspace.overlay_changes == ()
        assert workspace.baseline_manifest == ()
    assert seen == Counter(dict.fromkeys(responses, 1))


def test_workspace_refuses_a_disposable_head_mismatch(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    seen, _responses = _install_workspace_git_double(monkeypatch, real, sweep, disposable_head="b" * 40)
    with (
        pytest.raises(RuntimeError, match=r"disposable HEAD.*b{40}.*captured HEAD.*a{40}"),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("wrong-HEAD workspace must not yield")
    assert sum(seen.values()) == 4


def test_workspace_refuses_malformed_porcelain_before_overlay(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    seen, _responses = _install_workspace_git_double(monkeypatch, real, sweep, real_status=b"R  new.py\0")
    with (
        pytest.raises(RuntimeError, match="malformed porcelain rename"),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("malformed status must not yield")
    assert sum(seen.values()) == 5


def test_overlay_operations_and_baseline_manifest_diverge_for_staged_revert(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    sweep = tmp_path / "sweep"
    for root in (real, sweep):
        (root / "tests/unit").mkdir(parents=True)
        (root / "docs").mkdir()
        (root / "tests/unit/test_changed.py").write_text("head\n", encoding="utf-8")
        (root / "tests/unit/test_reverted.py").write_text("head\n", encoding="utf-8")
        (root / "tests/unit/test_link.py").symlink_to("test_reverted.py")
    (real / "tests/unit/test_changed.py").write_text("changed\n", encoding="utf-8")
    (real / "tests/unit/test_new.py").write_text("new\n", encoding="utf-8")
    (sweep / "docs/removed.md").write_text("removed\n", encoding="utf-8")
    real_status = (
        b" M tests/unit/test_changed.py\0"
        b"?? tests/unit/test_new.py\0"
        b" D docs/removed.md\0"
        b"MM tests/unit/test_reverted.py\0"
        b"MM tests/unit/test_link.py\0"
    )
    sweep_status = b" M tests/unit/test_changed.py\0?? tests/unit/test_new.py\0 D docs/removed.md\0"
    _install_workspace_git_double(
        monkeypatch,
        real,
        sweep,
        real_status=real_status,
        sweep_status=sweep_status,
    )

    with mutation_workspace.detached_head_workspace(real, (), ()) as workspace:
        assert (workspace.sweep_root / "tests/unit/test_changed.py").read_text() == "changed\n"
        assert (workspace.sweep_root / "tests/unit/test_new.py").read_text() == "new\n"
        assert not (workspace.sweep_root / "docs/removed.md").exists()
        assert (workspace.sweep_root / "tests/unit/test_link.py").is_symlink()
        assert (workspace.sweep_root / "tests/unit/test_link.py").readlink() == Path("test_reverted.py")
        assert workspace.overlay_changes == (
            mutation_workspace.OverlayChange("delete", "docs/removed.md"),
            mutation_workspace.OverlayChange("copy", "tests/unit/test_changed.py"),
            mutation_workspace.OverlayChange("copy", "tests/unit/test_link.py"),
            mutation_workspace.OverlayChange("copy", "tests/unit/test_new.py"),
            mutation_workspace.OverlayChange("copy", "tests/unit/test_reverted.py"),
        )
        assert mutation_workspace.OverlayChange("copy", "tests/unit/test_reverted.py") in workspace.overlay_changes
        assert (
            mutation_workspace.OverlayChange("copy", "tests/unit/test_reverted.py") not in workspace.baseline_manifest
        )
        assert workspace.baseline_manifest == mutation_workspace.parse_porcelain_z(
            workspace.baseline_status,
            workspace.sweep_root,
        )
        workspace.verify_baseline()


@pytest.mark.parametrize("relative", [".git/config", "results/artifact.json", "summary.json"])
def test_excluded_overlay_paths_are_rejected(tmp_path: Path, monkeypatch, relative: str) -> None:
    real = tmp_path / "real"
    real.mkdir()
    source = real / relative
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text("excluded\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    _install_workspace_git_double(monkeypatch, real, sweep, real_status=f"?? {relative}\0".encode())
    excluded_outputs = (real / "results", real / "summary.json")

    with (
        pytest.raises(RuntimeError, match="unsafe workspace path|excluded output"),
        mutation_workspace.detached_head_workspace(real, (), excluded_outputs),
    ):
        pytest.fail("excluded overlay path must not yield")


def test_registered_nested_worktree_is_rejected(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    nested_file = real / "vendor/nested/changed.py"
    nested_file.parent.mkdir(parents=True)
    nested_file.write_text("nested\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    worktrees = (f"worktree {real}\0HEAD {'a' * 40}\0\0worktree {real / 'vendor/nested'}\0HEAD {'b' * 40}\0\0").encode()
    seen, _responses = _install_workspace_git_double(
        monkeypatch,
        real,
        sweep,
        real_status=b" M vendor/nested/changed.py\0",
        worktrees=worktrees,
    )

    with (
        pytest.raises(RuntimeError, match="nested registered worktree"),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("nested worktree overlay must not yield")

    list_call = (("git", "worktree", "list", "--porcelain", "-z"), real)
    assert seen[list_call] == 1


def test_registered_ancestor_worktree_does_not_exclude_real_root(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "primary/authorized"
    changed = real / "tests/unit/test_changed.py"
    changed.parent.mkdir(parents=True)
    changed.write_text("changed\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    worktrees = (f"worktree {real.parent}\0HEAD {'b' * 40}\0\0worktree {real}\0HEAD {'a' * 40}\0\0").encode()
    seen, _responses = _install_workspace_git_double(
        monkeypatch,
        real,
        sweep,
        real_status=b"?? tests/unit/test_changed.py\0",
        sweep_status=b"?? tests/unit/test_changed.py\0",
        worktrees=worktrees,
    )

    with mutation_workspace.detached_head_workspace(real, (), ()) as workspace:
        assert (workspace.sweep_root / "tests/unit/test_changed.py").read_text() == "changed\n"

    list_call = (("git", "worktree", "list", "--porcelain", "-z"), real)
    assert seen[list_call] == 1


def test_external_symlink_overlay_is_rejected(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    outside = tmp_path / "outside.py"
    outside.write_text("outside\n", encoding="utf-8")
    (real / "external.py").symlink_to(outside)
    sweep = tmp_path / "sweep"
    _install_workspace_git_double(monkeypatch, real, sweep, real_status=b"?? external.py\0")

    with (
        pytest.raises(RuntimeError, match="escapes workspace root"),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("external symlink overlay must not yield")


def _install_overlay_baseline_fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[Path, Path, tuple[str, ...]]:
    real = tmp_path / "real"
    sweep = tmp_path / "sweep"
    for root in (real, sweep):
        (root / "tests/unit").mkdir(parents=True)
        (root / "docs").mkdir()
        (root / "vntyper/scripts").mkdir(parents=True)
        (root / "tests/unit/test_changed.py").write_text("head\n", encoding="utf-8")
        (root / "tests/unit/target-a.py").write_text("a\n", encoding="utf-8")
        (root / "tests/unit/target-b.py").write_text("b\n", encoding="utf-8")
        (root / "tests/unit/test_link.py").symlink_to("target-a.py")
        (root / "vntyper/scripts/selected.py").write_text("selected\n", encoding="utf-8")
    (real / "tests/unit/test_changed.py").write_text("changed\n", encoding="utf-8")
    (sweep / "docs/removed.md").write_text("removed\n", encoding="utf-8")
    real_status = b" M tests/unit/test_changed.py\0MM tests/unit/test_link.py\0 D docs/removed.md\0"
    sweep_status = b" M tests/unit/test_changed.py\0 D docs/removed.md\0"
    _install_workspace_git_double(
        monkeypatch,
        real,
        sweep,
        real_status=real_status,
        sweep_status=sweep_status,
    )
    return real, sweep, ("vntyper/scripts/selected.py",)


@pytest.mark.parametrize(
    ("drift", "relative", "mismatch"),
    [
        ("content", "tests/unit/test_changed.py", "content mismatch"),
        ("symlink", "tests/unit/test_link.py", "symlink mismatch"),
        ("deletion", "docs/removed.md", "deletion mismatch"),
        ("target", "vntyper/scripts/selected.py", "content mismatch"),
    ],
)
def test_baseline_detects_content_symlink_and_deletion_drift(
    tmp_path: Path,
    monkeypatch,
    drift: str,
    relative: str,
    mismatch: str,
) -> None:
    real, _sweep, targets = _install_overlay_baseline_fixture(tmp_path, monkeypatch)

    with mutation_workspace.detached_head_workspace(real, targets, ()) as workspace:
        path = workspace.sweep_root / relative
        if drift in {"content", "target"}:
            path.write_text("drifted\n", encoding="utf-8")
        elif drift == "symlink":
            path.unlink()
            path.symlink_to("target-b.py")
        else:
            path.write_text("recreated\n", encoding="utf-8")

        with pytest.raises(RuntimeError, match=rf"{relative}.*{mismatch}"):
            workspace.verify_baseline()


def test_baseline_accepts_unchanged_post_overlay_state(tmp_path: Path, monkeypatch) -> None:
    real, _sweep, targets = _install_overlay_baseline_fixture(tmp_path, monkeypatch)

    with mutation_workspace.detached_head_workspace(real, targets, ()) as workspace:
        assert workspace.target_path(targets[0]) == workspace.sweep_root / targets[0]
        assert targets[0] in workspace.baseline_digests
        workspace.verify_baseline()


def test_cleanup_failure_retains_only_the_exact_orphan(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    seen, responses = _install_workspace_git_double(monkeypatch, real, sweep)
    remove_call = (("git", "worktree", "remove", "--force", str(sweep)), real)
    responses[remove_call] = (1, b"", b"busy")

    with pytest.raises(RuntimeError) as exc_info, mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    assert str(exc_info.value) == f"orphaned worktree: {sweep}: busy"
    assert sweep.exists()
    assert seen[remove_call] == 1
    forbidden = {"prune", "reset", "checkout", "clean"}
    assert not any(forbidden.intersection(command) for command, _cwd in seen)


def test_fresh_retry_never_reuses_an_orphan(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    first_sweep = tmp_path / "first-sweep"
    second_sweep = tmp_path / "second-sweep"
    first_sweep.mkdir()
    second_sweep.mkdir()
    head = "a" * 40
    monkeypatch.setattr(
        mutation_workspace.tempfile,
        "mkdtemp",
        mock.Mock(side_effect=[str(first_sweep), str(second_sweep)]),
    )
    responses: dict[GitCall, tuple[int, bytes, bytes]] = {
        (("git", "rev-parse", "--verify", "HEAD^{commit}"), real): (0, f"{head}\n".encode(), b""),
        (("git", "status", "--porcelain=v1", "-z", "--untracked-files=all"), real): (0, b"", b""),
        (("git", "worktree", "list", "--porcelain", "-z"), real): (0, b"", b""),
    }
    for sweep in (first_sweep, second_sweep):
        responses[(("git", "worktree", "add", "--detach", str(sweep), head), real)] = (0, b"", b"")
        responses[(("git", "-C", str(sweep), "rev-parse", "--verify", "HEAD^{commit}"), real)] = (
            0,
            f"{head}\n".encode(),
            b"",
        )
        responses[(("git", "-C", str(sweep), "status", "--porcelain=v1", "-z", "--untracked-files=all"), real)] = (
            0,
            b"",
            b"",
        )
        responses[(("git", "worktree", "remove", "--force", str(sweep)), real)] = (0, b"", b"")
    responses[(("git", "worktree", "remove", "--force", str(first_sweep)), real)] = (1, b"", b"busy")
    commands: list[tuple[str, ...]] = []

    def fake_run(command: list[str], *, cwd: Path, **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        key = (tuple(command), Path(cwd))
        assert key in responses, f"unexpected subprocess call: {key}"
        commands.append(tuple(command))
        returncode, stdout, stderr = responses[key]
        if command[1:3] == ["worktree", "remove"] and returncode == 0:
            Path(command[-1]).rmdir()
        return subprocess.CompletedProcess(command, returncode, stdout, stderr)

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)

    with (
        pytest.raises(RuntimeError, match=rf"orphaned worktree: {first_sweep}: busy"),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pass
    with mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    assert first_sweep.exists()
    add_commands = [command for command in commands if command[1:3] == ("worktree", "add")]
    assert add_commands == [
        ("git", "worktree", "add", "--detach", str(first_sweep), head),
        ("git", "worktree", "add", "--detach", str(second_sweep), head),
    ]


@pytest.mark.parametrize(
    ("stage", "expected_calls", "removes_registered_worktree"),
    [
        ("capture", 1, False),
        ("add", 3, True),
        ("disposable-head", 4, True),
        ("real-status", 5, True),
        ("worktree-list", 6, True),
        ("sweep-status", 7, True),
    ],
)
def test_lifecycle_git_failures_use_only_the_available_cleanup_path(
    tmp_path: Path,
    monkeypatch,
    stage: str,
    expected_calls: int,
    removes_registered_worktree: bool,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    seen, responses = _install_workspace_git_double(monkeypatch, real, sweep)
    commands = {
        "capture": ("git", "rev-parse", "--verify", "HEAD^{commit}"),
        "add": ("git", "worktree", "add", "--detach", str(sweep), "a" * 40),
        "disposable-head": ("git", "-C", str(sweep), "rev-parse", "--verify", "HEAD^{commit}"),
        "real-status": ("git", "status", "--porcelain=v1", "-z", "--untracked-files=all"),
        "worktree-list": ("git", "worktree", "list", "--porcelain", "-z"),
        "sweep-status": ("git", "-C", str(sweep), "status", "--porcelain=v1", "-z", "--untracked-files=all"),
    }
    responses[(commands[stage], real)] = (1, b"", b"boom")

    with pytest.raises(RuntimeError, match="boom"), mutation_workspace.detached_head_workspace(real, (), ()):
        pytest.fail("a failed Git lifecycle command must not yield")

    remove_call = (("git", "worktree", "remove", "--force", str(sweep)), real)
    assert sum(seen.values()) == expected_calls
    assert bool(seen[remove_call]) is removes_registered_worktree
    if stage == "add":
        assert not sweep.exists()


def test_successful_git_removal_deletes_only_a_remaining_exact_directory(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    _install_workspace_git_double(monkeypatch, real, sweep, remove_directory=False)

    with mutation_workspace.detached_head_workspace(real, (), ()):
        pass

    assert not sweep.exists()


def test_porcelain_z_maps_changes_renames_untracked_and_deletions(tmp_path: Path) -> None:
    payload = (
        b" M tests/unit/test_changed.py\0"
        b"?? tests/unit/test_new.py\0"
        b" D docs/removed.md\0"
        b"R  tests/unit/test_renamed.py\0tests/unit/test_old.py\0"
    )
    for relative in ("tests/unit/test_changed.py", "tests/unit/test_new.py", "tests/unit/test_renamed.py"):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("current\n", encoding="utf-8")
    assert mutation_workspace.parse_porcelain_z(payload, tmp_path) == (
        mutation_workspace.OverlayChange("delete", "docs/removed.md"),
        mutation_workspace.OverlayChange("copy", "tests/unit/test_changed.py"),
        mutation_workspace.OverlayChange("copy", "tests/unit/test_new.py"),
        mutation_workspace.OverlayChange("delete", "tests/unit/test_old.py"),
        mutation_workspace.OverlayChange("copy", "tests/unit/test_renamed.py"),
    )


@pytest.mark.parametrize("relative", ["", ".", "/absolute.py", "../escape.py", ".git/config", "a/.git/config"])
def test_confined_path_rejects_unsafe_names(tmp_path: Path, relative: str) -> None:
    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.confined_path(tmp_path, relative, must_exist=False)


def test_confined_path_rejects_a_symlink_escape(tmp_path: Path) -> None:
    outside = tmp_path.parent / "outside.py"
    outside.write_text("outside\n", encoding="utf-8")
    (tmp_path / "link.py").symlink_to(outside)
    with pytest.raises(ValueError, match="escapes workspace root"):
        mutation_workspace.confined_path(tmp_path, "link.py", must_exist=True)


def test_porcelain_rejects_unsafe_path_before_filesystem_probe(tmp_path: Path, monkeypatch) -> None:
    probes: list[object] = []
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: _record_probe(probes, path, False))
    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.parse_porcelain_z(b"?? ../outside.py\0", tmp_path)
    assert probes == []


def test_delete_plus_untracked_replacement_has_one_final_copy(tmp_path: Path) -> None:
    replacement = tmp_path / "tests/unit/replaced.py"
    replacement.parent.mkdir(parents=True)
    replacement.write_text("replacement\n", encoding="utf-8")
    payload = b"D  tests/unit/replaced.py\0?? tests/unit/replaced.py\0"
    assert mutation_workspace.parse_porcelain_z(payload, tmp_path) == (
        mutation_workspace.OverlayChange("copy", "tests/unit/replaced.py"),
    )


def test_deleted_path_has_one_final_delete(tmp_path: Path) -> None:
    payload = b" D tests/unit/deleted.py\0"
    assert mutation_workspace.parse_porcelain_z(payload, tmp_path) == (
        mutation_workspace.OverlayChange("delete", "tests/unit/deleted.py"),
    )


def test_duplicate_status_records_have_one_final_action(tmp_path: Path) -> None:
    changed = tmp_path / "tests/unit/changed.py"
    changed.parent.mkdir(parents=True)
    changed.write_text("changed\n", encoding="utf-8")
    payload = b" M tests/unit/changed.py\0M  tests/unit/changed.py\0"
    assert mutation_workspace.parse_porcelain_z(payload, tmp_path) == (
        mutation_workspace.OverlayChange("copy", "tests/unit/changed.py"),
    )


@pytest.mark.parametrize("rename_status", [b"R ", b" C"])
def test_affected_paths_expands_a_rename_without_filesystem_access(rename_status: bytes) -> None:
    payload = b" M current.py\0" + rename_status + b" renamed.py\0original.py\0"
    assert mutation_workspace._affected_paths_from_porcelain_z(payload) == (
        "current.py",
        "renamed.py",
        "original.py",
    )


def test_affected_paths_rejects_ignored_records() -> None:
    with pytest.raises(ValueError, match="ignored porcelain record"):
        mutation_workspace._affected_paths_from_porcelain_z(b"!! ignored.py\0")


def test_affected_paths_rejects_conflicting_rename_encoding() -> None:
    with pytest.raises(ValueError, match="conflicting rename encoding"):
        mutation_workspace._affected_paths_from_porcelain_z(b"RC renamed.py\0original.py\0")


@pytest.mark.parametrize("status", [b"XX", b"  ", b"?M", b"MR"])
def test_affected_paths_rejects_invalid_porcelain_status_pairs(status: bytes) -> None:
    payload = status + b" changed.py\0original.py\0"
    with pytest.raises(ValueError, match="invalid porcelain status"):
        mutation_workspace._affected_paths_from_porcelain_z(payload)


def test_affected_paths_rejects_a_rename_without_an_original_path() -> None:
    with pytest.raises(ValueError, match="missing original path"):
        mutation_workspace._affected_paths_from_porcelain_z(b"R  renamed.py\0")


def test_affected_paths_rejects_a_missing_nul_terminator() -> None:
    with pytest.raises(ValueError, match="missing NUL terminator"):
        mutation_workspace._affected_paths_from_porcelain_z(b" M changed.py")


@pytest.mark.parametrize("payload", [b"M changed.py\0", b" M \0"])
def test_affected_paths_rejects_malformed_records(payload: bytes) -> None:
    with pytest.raises(ValueError, match="malformed porcelain record"):
        mutation_workspace._affected_paths_from_porcelain_z(payload)


def test_empty_porcelain_payload_has_no_changes(tmp_path: Path) -> None:
    assert mutation_workspace.parse_porcelain_z(b"", tmp_path) == ()


def test_confined_path_rejects_a_missing_required_path(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="workspace path does not exist"):
        mutation_workspace.confined_path(tmp_path, "missing.py", must_exist=True)


def test_confined_path_rejects_a_parent_symlink_escape(tmp_path: Path) -> None:
    outside = tmp_path.parent / "outside"
    outside.mkdir()
    (tmp_path / "link").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="escapes workspace root"):
        mutation_workspace.confined_path(tmp_path, "link/child.py", must_exist=False)


@pytest.mark.parametrize(
    "payload",
    [
        b"?? safe.py\0?? ../outside.py\0",
        b"R  renamed.py\0../original.py\0",
    ],
)
def test_porcelain_validates_every_path_before_lexists(tmp_path: Path, monkeypatch, payload: bytes) -> None:
    (tmp_path / "safe.py").write_text("safe\n", encoding="utf-8")
    (tmp_path / "renamed.py").write_text("renamed\n", encoding="utf-8")
    probes: list[object] = []
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: _record_probe(probes, path, True))

    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.parse_porcelain_z(payload, tmp_path)

    assert probes == []


def test_confined_path_rejects_a_final_symlink_into_git_metadata(tmp_path: Path) -> None:
    metadata = tmp_path / ".git"
    metadata.mkdir()
    config = metadata / "config"
    config.write_text("metadata\n", encoding="utf-8")
    (tmp_path / "config-link").symlink_to(config)

    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.confined_path(tmp_path, "config-link", must_exist=True)


def test_confined_path_accepts_an_internal_symlink_outside_git_metadata(tmp_path: Path) -> None:
    target = tmp_path / "target.py"
    target.write_text("target\n", encoding="utf-8")
    link = tmp_path / "link.py"
    link.symlink_to(target)

    assert mutation_workspace.confined_path(tmp_path, "link.py", must_exist=True) == link


def test_confined_path_rejects_a_parent_symlink_into_git_metadata(tmp_path: Path) -> None:
    metadata = tmp_path / ".git"
    metadata.mkdir()
    (metadata / "config").write_text("metadata\n", encoding="utf-8")
    (tmp_path / "metadata-link").symlink_to(metadata, target_is_directory=True)

    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.confined_path(tmp_path, "metadata-link/config", must_exist=True)


def test_confined_path_reports_a_dangling_internal_symlink_as_a_value_error(tmp_path: Path) -> None:
    (tmp_path / "dangling.py").symlink_to(tmp_path / "missing.py")

    with pytest.raises(ValueError, match="workspace path does not exist"):
        mutation_workspace.confined_path(tmp_path, "dangling.py", must_exist=False)


def test_confined_path_reports_a_missing_required_parent_as_a_value_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="workspace path does not exist"):
        mutation_workspace.confined_path(tmp_path, "missing/child.py", must_exist=True)


def test_porcelain_rejects_duplicate_lexical_spellings_before_lexists(tmp_path: Path, monkeypatch) -> None:
    changed = tmp_path / "dir/changed.py"
    changed.parent.mkdir()
    changed.write_text("changed\n", encoding="utf-8")
    probes: list[object] = []
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: _record_probe(probes, path, True))

    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.parse_porcelain_z(b"?? dir/changed.py\0?? dir/./changed.py\0", tmp_path)

    assert probes == []
