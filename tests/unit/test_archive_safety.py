"""Fail-closed and atomic result archive creation."""

import logging
import os
import tarfile
import zipfile
from pathlib import Path

import pytest

from vntyper.scripts import archive_safety

pytestmark = pytest.mark.unit

PATIENT_BYTES = b"patient-alignment-bytes-that-must-never-enter-an-archive"


@pytest.mark.parametrize(
    ("archive_format", "suffix"),
    [pytest.param("zip", ".zip", id="zip"), pytest.param("gztar", ".tar.gz", id="gztar")],
)
def test_regular_result_trees_are_installed_atomically(tmp_path: Path, archive_format: str, suffix: str) -> None:
    """Both CLI-supported formats preserve ordinary files and nested paths.

    Args:
        tmp_path: Scratch directory standing in for a result tree.
        archive_format: The shutil-compatible format name.
        suffix: The archive suffix that format produces.
    """
    root = tmp_path / "results"
    nested = root / "nested"
    nested.mkdir(parents=True)
    (root / "result.txt").write_bytes(b"top-level")
    (nested / "details.tsv").write_bytes(b"nested-result")

    archive = archive_safety.create_safe_archive(tmp_path / "download", archive_format, root)

    assert archive == str(tmp_path / f"download{suffix}")
    if archive_format == "zip":
        with zipfile.ZipFile(archive) as package:
            assert package.read("result.txt") == b"top-level"
            assert package.read("nested/details.tsv") == b"nested-result"
    else:
        with tarfile.open(archive, "r:gz") as package:
            top_level = package.extractfile("./result.txt")
            nested_result = package.extractfile("./nested/details.tsv")
            assert top_level is not None and top_level.read() == b"top-level"
            assert nested_result is not None and nested_result.read() == b"nested-result"


@pytest.mark.parametrize(
    ("archive_format", "suffix"),
    [pytest.param("zip", ".zip", id="zip"), pytest.param("gztar", ".tar.gz", id="gztar")],
)
def test_archive_succeeds_after_owned_alignment_plan_close_removes_exact_view(
    tmp_path: Path, archive_format: str, suffix: str
) -> None:
    """H1's ``plan.close()`` boundary removes its view before either archive format."""
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    root = tmp_path / "results"
    root.mkdir()
    view = root / "alignment_view.bam"
    view.symlink_to(patient)
    (root / "result.txt").write_bytes(b"safe result")

    class OwnedPlan:
        def close(self) -> None:
            view.unlink()

    plan = OwnedPlan()
    plan.close()
    archive = Path(archive_safety.create_safe_archive(tmp_path / "download", archive_format, root))

    assert archive == tmp_path / f"download{suffix}"
    assert patient.read_bytes() == PATIENT_BYTES
    assert PATIENT_BYTES not in archive.read_bytes()


@pytest.mark.parametrize("link_target_exists", [True, False], ids=["live", "broken"])
def test_a_file_symlink_is_rejected_before_any_archive_is_installed(tmp_path: Path, link_target_exists: bool) -> None:
    """Live and broken aliases are unsafe regardless of target readability.

    Args:
        tmp_path: Scratch directory standing in for input and output trees.
        link_target_exists: Whether the symlink target exists.
    """
    patient = tmp_path / "patient.bam"
    if link_target_exists:
        patient.write_bytes(PATIENT_BYTES)
    root = tmp_path / "results"
    root.mkdir()
    (root / "alignment_view.bam").symlink_to(patient)

    with pytest.raises(ValueError, match="symbolic link.*alignment_view\\.bam"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()
    if link_target_exists:
        assert patient.read_bytes() == PATIENT_BYTES


def test_a_directory_symlink_is_rejected_without_reading_its_tree(tmp_path: Path) -> None:
    """A directory alias must not turn an external tree into archive members.

    Args:
        tmp_path: Scratch directory standing in for input and output trees.
    """
    patient_dir = tmp_path / "patient"
    patient_dir.mkdir()
    (patient_dir / "patient.cram").write_bytes(PATIENT_BYTES)
    root = tmp_path / "results"
    root.mkdir()
    (root / "external").symlink_to(patient_dir, target_is_directory=True)

    with pytest.raises(ValueError, match="symbolic link.*external"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()
    assert (patient_dir / "patient.cram").read_bytes() == PATIENT_BYTES


def test_a_hard_link_to_external_patient_bytes_is_rejected(tmp_path: Path) -> None:
    """A second directory entry is also an external alias, despite being regular.

    Args:
        tmp_path: Scratch directory standing in for input and output trees.
    """
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    root = tmp_path / "results"
    root.mkdir()
    os.link(patient, root / "alignment_view.bam")

    with pytest.raises(ValueError, match="hard-linked file.*alignment_view\\.bam"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()
    assert patient.read_bytes() == PATIENT_BYTES


def test_a_fifo_is_rejected_by_metadata_without_opening_it(tmp_path: Path) -> None:
    """Unsupported special entries must fail immediately rather than block on read.

    Args:
        tmp_path: Scratch directory standing in for a result tree.
    """
    root = tmp_path / "results"
    root.mkdir()
    os.mkfifo(root / "stream")

    with pytest.raises(ValueError, match="unsupported filesystem entry.*stream"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()


def test_a_symlinked_result_root_is_rejected(tmp_path: Path) -> None:
    """Validation begins at the root rather than only checking its descendants.

    Args:
        tmp_path: Scratch directory standing in for a result tree.
    """
    external = tmp_path / "external"
    external.mkdir()
    (external / "patient.bam").write_bytes(PATIENT_BYTES)
    root = tmp_path / "results"
    root.symlink_to(external, target_is_directory=True)

    with pytest.raises(ValueError, match="symbolic link.*result root"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()
    assert (external / "patient.bam").read_bytes() == PATIENT_BYTES


def test_an_archive_write_failure_removes_the_stale_destination_and_any_partial_temp(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failed retry must not leave an older archive available as current.

    Args:
        monkeypatch: Standard pytest fixture.
        tmp_path: Scratch directory standing in for a result tree.
    """
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"new result")
    destination = tmp_path / "download.zip"
    destination.write_bytes(b"previous complete archive")

    monkeypatch.setattr(archive_safety.os, "replace", lambda *args: (_ for _ in ()).throw(OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not destination.exists()
    assert not [path for path in tmp_path.iterdir() if path.name.startswith(".download.archive-")]


@pytest.mark.parametrize("collision", ["lexical", "resolved", "hardlink"], ids=str)
def test_archive_destination_cannot_alias_a_protected_operator_input(tmp_path: Path, collision: str) -> None:
    """Lexical, symlink-resolved, and inode aliases are rejected before replacement."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    destination = tmp_path / "download.zip"
    protected = patient
    if collision == "lexical":
        protected = destination
    elif collision == "resolved":
        destination.symlink_to(patient)
    else:
        os.link(patient, destination)

    with pytest.raises(ValueError, match="archive destination.*protected input"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root, protected_paths=(protected,))

    assert patient.read_bytes() == PATIENT_BYTES


def test_archive_destination_is_rejected_inside_its_source_tree(tmp_path: Path) -> None:
    """The temporary or public archive must never become a recursive source member."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")

    with pytest.raises(ValueError, match="archive destination.*source tree"):
        archive_safety.create_safe_archive(root / "download", "zip", root)

    assert not (root / "download.zip").exists()


def test_replacement_between_metadata_check_and_open_never_archives_external_bytes(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A concurrent path swap must hit O_NOFOLLOW instead of reopening the alias."""
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    root = tmp_path / "results"
    root.mkdir()
    result = root / "result.txt"
    result.write_bytes(b"safe result")
    original_stat = archive_safety.os.stat

    def replace_after_stat(path, *args, **kwargs):
        metadata = original_stat(path, *args, **kwargs)
        if path == result.name and kwargs.get("dir_fd") is not None:
            result.unlink()
            result.symlink_to(patient)
        return metadata

    monkeypatch.setattr(archive_safety.os, "stat", replace_after_stat)

    with pytest.raises((OSError, ValueError)):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    archive = tmp_path / "download.zip"
    assert not archive.exists() or PATIENT_BYTES not in archive.read_bytes()
    assert patient.read_bytes() == PATIENT_BYTES


def test_temporary_cleanup_failure_is_logged_without_masking_primary_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A retained partial is named while the archive failure remains the outcome."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    monkeypatch.setattr(archive_safety.os, "replace", lambda *args: (_ for _ in ()).throw(OSError("install failed")))
    monkeypatch.setattr(
        archive_safety.shutil, "rmtree", lambda *args, **kwargs: (_ for _ in ()).throw(OSError("cleanup failed"))
    )
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.archive_safety")

    with pytest.raises(OSError, match="install failed"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert "cleanup failed" in caplog.text
    assert ".download.archive-" in caplog.text


def test_unknown_formats_and_non_directory_roots_are_rejected(tmp_path: Path) -> None:
    """Caller mistakes fail before any archive or temporary directory is created.

    Args:
        tmp_path: Scratch directory standing in for a result tree.
    """
    result_file = tmp_path / "result.txt"
    result_file.write_bytes(b"result")

    with pytest.raises(ValueError, match="Unsupported archive format: rar"):
        archive_safety.create_safe_archive(tmp_path / "download", "rar", tmp_path)
    with pytest.raises(ValueError, match="not a directory"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", result_file)

    assert sorted(path.name for path in tmp_path.iterdir()) == ["result.txt"]
