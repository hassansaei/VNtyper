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
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: probes.append(path) or True)

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
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: probes.append(path) or True)

    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.parse_porcelain_z(b"?? dir/changed.py\0?? dir/./changed.py\0", tmp_path)

    assert probes == []
