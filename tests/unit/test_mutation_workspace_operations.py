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


def _git(repo: Path, *arguments: str) -> None:
    subprocess.run(["git", "-C", str(repo), *arguments], check=True, capture_output=True)


def test_real_git_overlay_accepts_a_link_before_its_later_sorted_untracked_target(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "probe@example.invalid")
    _git(repo, "config", "user.name", "Probe")
    (repo / "tracked.txt").write_text("tracked\n", encoding="utf-8")
    _git(repo, "add", "tracked.txt")
    _git(repo, "commit", "-qm", "baseline")
    (repo / "z_target").write_text("target\n", encoding="utf-8")
    (repo / "a_link").symlink_to("z_target")

    with mutation_workspace.detached_head_workspace(repo, ("z_target",), ()) as workspace:
        link = workspace.sweep_root / "a_link"
        target = workspace.sweep_root / "z_target"
        assert link.is_symlink()
        assert link.readlink() == Path("z_target")
        assert target.read_text(encoding="utf-8") == "target\n"
        workspace.verify_baseline()


def test_overlay_validates_every_source_before_writing_any_destination(tmp_path: Path) -> None:
    real = tmp_path / "real"
    sweep = tmp_path / "sweep"
    for root in (real, sweep):
        root.mkdir()
        (root / "a_safe.txt").write_text("head\n", encoding="utf-8")
    (real / "a_safe.txt").write_text("overlay\n", encoding="utf-8")
    excluded = real / "z_excluded.txt"
    excluded.write_text("excluded\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="excluded output"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (
                mutation_workspace.OverlayChange("copy", "a_safe.txt"),
                mutation_workspace.OverlayChange("copy", "z_excluded.txt"),
            ),
            (excluded,),
            (),
        )

    assert (sweep / "a_safe.txt").read_text(encoding="utf-8") == "head\n"


def test_recursive_destination_removal_does_not_follow_nested_external_symlink(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (real / "changed").write_text("overlay\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    nested = sweep / "changed/nested"
    nested.mkdir(parents=True)
    external = tmp_path / "external.txt"
    external.write_text("sentinel\n", encoding="utf-8")
    (nested / "external-link").symlink_to(external)

    mutation_workspace._apply_overlay_changes(
        real,
        sweep,
        (mutation_workspace.OverlayChange("copy", "changed"),),
        (),
        (),
    )

    assert (sweep / "changed").read_text(encoding="utf-8") == "overlay\n"
    assert external.read_text(encoding="utf-8") == "sentinel\n"


def test_recursive_destination_removal_detects_final_directory_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "root"
    entry = root / "entry"
    entry.mkdir(parents=True)
    (entry / "child.txt").write_text("child\n", encoding="utf-8")
    displaced = root / "displaced"
    original_stat = mutation_workspace.os.stat
    entry_stats = 0

    def racing_stat(
        path: str | bytes | os.PathLike[str] | os.PathLike[bytes],
        *,
        dir_fd: int | None = None,
        follow_symlinks: bool = True,
    ) -> os.stat_result:
        nonlocal entry_stats
        if path == "entry" and dir_fd is not None and not follow_symlinks:
            entry_stats += 1
            if entry_stats == 2:
                entry.rename(displaced)
                entry.mkdir()
                (entry / "replacement.txt").write_text("replacement\n", encoding="utf-8")
        return original_stat(path, dir_fd=dir_fd, follow_symlinks=follow_symlinks)

    monkeypatch.setattr(mutation_workspace.os, "stat", racing_stat)
    parent_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(RuntimeError, match="workspace destination changed during removal"):
            mutation_workspace._remove_entry_at(parent_fd, "entry")
    finally:
        os.close(parent_fd)

    assert (entry / "replacement.txt").read_text(encoding="utf-8") == "replacement\n"


def test_recursive_destination_removal_detects_replacement_before_opening_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "root"
    entry = root / "entry"
    entry.mkdir(parents=True)
    (entry / "original.txt").write_text("original\n", encoding="utf-8")
    displaced = root / "displaced"
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
        if path == "entry" and dir_fd is not None and not swapped:
            entry.rename(displaced)
            entry.mkdir()
            (entry / "replacement.txt").write_text("replacement\n", encoding="utf-8")
            swapped = True
        return original_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(mutation_workspace.os, "open", racing_open)
    parent_fd = original_open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(RuntimeError, match="workspace destination changed during removal"):
            mutation_workspace._remove_entry_at(parent_fd, "entry")
    finally:
        os.close(parent_fd)

    assert swapped is True
    assert (entry / "replacement.txt").read_text(encoding="utf-8") == "replacement\n"
    assert (displaced / "original.txt").read_text(encoding="utf-8") == "original\n"


def test_recursive_destination_removal_rejects_the_wrong_expected_inode(tmp_path: Path) -> None:
    root = tmp_path / "root"
    entry = root / "entry"
    entry.mkdir(parents=True)
    parent_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(RuntimeError, match="workspace destination changed during removal"):
            mutation_workspace._remove_entry_at(parent_fd, "entry", expected_identity=(0, 0))
    finally:
        os.close(parent_fd)

    assert entry.is_dir()


def test_regular_copy_rejects_a_nonregular_source(tmp_path: Path) -> None:
    source_root = tmp_path / "source"
    destination_root = tmp_path / "destination"
    source_root.mkdir()
    destination_root.mkdir()
    (source_root / "directory").mkdir()
    source_fd = os.open(source_root, os.O_RDONLY | os.O_DIRECTORY)
    destination_fd = os.open(destination_root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(RuntimeError, match="overlay source is not a regular file"):
            mutation_workspace._copy_regular_file_at(source_fd, "directory", destination_fd, "copied")
    finally:
        os.close(source_fd)
        os.close(destination_fd)

    assert not (destination_root / "copied").exists()


def test_partial_regular_copy_failure_removes_the_partial_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_root = tmp_path / "source"
    destination_root = tmp_path / "destination"
    source_root.mkdir()
    destination_root.mkdir()
    (source_root / "source.txt").write_text("source bytes\n", encoding="utf-8")
    source_fd = os.open(source_root, os.O_RDONLY | os.O_DIRECTORY)
    destination_fd = os.open(destination_root, os.O_RDONLY | os.O_DIRECTORY)
    monkeypatch.setattr(
        mutation_workspace.os, "write", lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("full"))
    )
    try:
        with pytest.raises(OSError, match="full"):
            mutation_workspace._copy_regular_file_at(source_fd, "source.txt", destination_fd, "copied.txt")
    finally:
        os.close(source_fd)
        os.close(destination_fd)

    assert not os.path.lexists(destination_root / "copied.txt")


@pytest.mark.parametrize("link_text", ["/external", "../../external"])
def test_symlink_plan_rejects_absolute_or_escaping_text(link_text: str) -> None:
    with pytest.raises(RuntimeError, match="workspace symlink escapes workspace root"):
        mutation_workspace._symlink_target_relative("nested/link", link_text)


def test_symlink_overlay_rejects_a_missing_unplanned_target(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    target = real / "target"
    target.write_text("target\n", encoding="utf-8")
    (real / "link").symlink_to("target")
    sweep = tmp_path / "sweep"
    sweep.mkdir()

    with pytest.raises(RuntimeError, match="workspace path does not exist|planned symlink target"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (mutation_workspace.OverlayChange("copy", "link"),),
            (),
            (),
        )


def test_symlink_overlay_rejects_an_external_committed_target(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (real / "target").write_text("target\n", encoding="utf-8")
    (real / "link").symlink_to("target")
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    external = tmp_path / "external"
    external.write_text("external\n", encoding="utf-8")
    (sweep / "target").symlink_to(external)

    with pytest.raises(RuntimeError, match="planned symlink target is unsafe"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (mutation_workspace.OverlayChange("copy", "link"),),
            (),
            (),
        )

    assert not os.path.lexists(sweep / "link")
    assert external.read_text(encoding="utf-8") == "external\n"


def test_overlay_delete_normalizes_an_os_error_without_losing_the_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    destination = sweep / "deleted.txt"
    destination.write_text("keep\n", encoding="utf-8")
    monkeypatch.setattr(
        mutation_workspace.os, "unlink", lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("denied"))
    )

    with pytest.raises(RuntimeError, match="unsafe workspace destination.*denied"):
        mutation_workspace._delete_overlay_entry(sweep, "deleted.txt")

    assert destination.read_text(encoding="utf-8") == "keep\n"


def test_baseline_capture_rejects_a_present_deletion_sentinel(tmp_path: Path) -> None:
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    (sweep / "deleted.txt").write_text("still present\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="deleted.txt.*deletion mismatch"):
        mutation_workspace._capture_baseline_digests(
            sweep,
            (mutation_workspace.OverlayChange("delete", "deleted.txt"),),
            (),
        )


def test_baseline_capture_rejects_a_missing_copied_path(tmp_path: Path) -> None:
    sweep = tmp_path / "sweep"
    sweep.mkdir()

    with pytest.raises(RuntimeError, match="workspace path does not exist: missing.txt"):
        mutation_workspace._capture_baseline_digests(
            sweep,
            (mutation_workspace.OverlayChange("copy", "missing.txt"),),
            (),
        )


def test_baseline_capture_normalizes_digest_read_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    (sweep / "changed.txt").write_text("changed\n", encoding="utf-8")
    monkeypatch.setattr(mutation_workspace.os, "read", lambda *_args: (_ for _ in ()).throw(OSError("unreadable")))

    with pytest.raises(RuntimeError, match="changed.txt.*baseline capture failed.*unreadable"):
        mutation_workspace._capture_baseline_digests(
            sweep,
            (mutation_workspace.OverlayChange("copy", "changed.txt"),),
            (),
        )


def test_overlay_source_validation_rejects_missing_copy_before_destination_creation(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    sweep.mkdir()

    with pytest.raises(RuntimeError, match="workspace path does not exist"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (mutation_workspace.OverlayChange("copy", "missing/changed.txt"),),
            (),
            (),
        )

    assert not (sweep / "missing").exists()


def test_overlay_source_validation_rejects_a_nonregular_copy(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (real / "directory").mkdir()
    sweep = tmp_path / "sweep"
    sweep.mkdir()

    with pytest.raises(RuntimeError, match="overlay source is not a regular file or symlink"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (mutation_workspace.OverlayChange("copy", "directory"),),
            (),
            (),
        )

    assert not (sweep / "directory").exists()


def test_overlay_source_validation_normalizes_an_os_error_before_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (real / "source.txt").write_text("source\n", encoding="utf-8")
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    original_stat = mutation_workspace.os.stat

    def failing_stat(
        path: str | bytes | os.PathLike[str] | os.PathLike[bytes],
        *,
        dir_fd: int | None = None,
        follow_symlinks: bool = True,
    ) -> os.stat_result:
        if path == "source.txt" and dir_fd is not None and not follow_symlinks:
            raise OSError("source raced")
        return original_stat(path, dir_fd=dir_fd, follow_symlinks=follow_symlinks)

    monkeypatch.setattr(mutation_workspace.os, "stat", failing_stat)

    with pytest.raises(RuntimeError, match="unsafe overlay source: source.txt: source raced"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (mutation_workspace.OverlayChange("copy", "source.txt"),),
            (),
            (),
        )

    assert not (sweep / "source.txt").exists()


def test_symlink_overlay_rejects_a_target_planned_for_deletion_before_writing(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (real / "target").write_text("target\n", encoding="utf-8")
    (real / "link").symlink_to("target")
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    (sweep / "target").write_text("committed\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="planned symlink target is deleted"):
        mutation_workspace._apply_overlay_changes(
            real,
            sweep,
            (
                mutation_workspace.OverlayChange("copy", "link"),
                mutation_workspace.OverlayChange("delete", "target"),
            ),
            (),
            (),
        )

    assert not os.path.lexists(sweep / "link")
    assert (sweep / "target").read_text(encoding="utf-8") == "committed\n"


def test_lifecycle_rejects_malformed_disposable_status(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    sweep.mkdir()
    head = "a" * 40
    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", lambda **_kwargs: str(sweep))

    def fake_run(
        command: list[str],
        *,
        cwd: Path,
        **_kwargs: object,
    ) -> subprocess.CompletedProcess[bytes]:
        command = [str(Path(part).resolve()) if part.startswith("/proc/self/fd/") else part for part in command]
        if command == ["git", "rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if command == ["git", "worktree", "add", "--detach", str(sweep), head]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command == ["git", "-C", str(sweep), "rev-parse", "--verify", "HEAD^{commit}"]:
            return subprocess.CompletedProcess(command, 0, f"{head}\n".encode(), b"")
        if command == ["git", "status", "--porcelain=v1", "-z", "--untracked-files=all"]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command == ["git", "worktree", "list", "--porcelain", "-z"]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command == ["git", "-C", str(sweep), "status", "--porcelain=v1", "-z", "--untracked-files=all"]:
            return subprocess.CompletedProcess(command, 0, b"R  renamed.py\0", b"")
        if command == ["git", "worktree", "repair", str(sweep)]:
            return subprocess.CompletedProcess(command, 0, b"", b"")
        if command == ["git", "worktree", "remove", "--force", str(sweep)]:
            shutil.rmtree(sweep)
            return subprocess.CompletedProcess(command, 0, b"", b"")
        raise AssertionError(f"unexpected subprocess call: {command} from {cwd}")

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)

    with (
        pytest.raises(RuntimeError, match="malformed porcelain rename"),
        mutation_workspace.detached_head_workspace(real, (), ()),
    ):
        pytest.fail("malformed disposable status must not yield")
