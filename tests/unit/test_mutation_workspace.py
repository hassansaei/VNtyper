from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace

pytestmark = pytest.mark.unit


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
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: probes.append(path) or False)
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
