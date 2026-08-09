"""Fail-closed and atomic result archive creation."""

import logging
import os
import tarfile
import tempfile
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

    monkeypatch.setattr(archive_safety.os, "link", lambda *args, **kwargs: (_ for _ in ()).throw(OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not destination.exists()
    assert not [path for path in tmp_path.iterdir() if path.name.startswith(".download")]


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


@pytest.mark.parametrize("unsafe_kind", ["symlink", "directory", "hardlink"], ids=str)
def test_unprotected_unsafe_archive_destinations_are_rejected(tmp_path: Path, unsafe_kind: str) -> None:
    """Destination entry safety does not depend on callers listing a protected path."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    destination = tmp_path / "download.zip"
    external = tmp_path / "external.bin"
    external.write_bytes(PATIENT_BYTES)
    if unsafe_kind == "symlink":
        destination.symlink_to(external)
    elif unsafe_kind == "directory":
        destination.mkdir()
    else:
        os.link(external, destination)

    with pytest.raises(ValueError, match=f"archive destination.*{'symbolic link|regular file|hard links'}"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert external.read_bytes() == PATIENT_BYTES


def test_stale_cleanup_preserves_a_protected_lexical_destination(tmp_path: Path) -> None:
    """The standalone retry cleanup applies the same protected-path decision."""
    destination = tmp_path / "download.zip"
    destination.write_bytes(PATIENT_BYTES)

    with pytest.raises(ValueError, match="archive destination.*protected input"):
        archive_safety.clear_stale_archive(tmp_path / "download", "zip", protected_paths=(destination,))

    assert destination.read_bytes() == PATIENT_BYTES


@pytest.mark.parametrize("initially_present", [False, True], ids=["appears", "is-replaced"])
def test_stale_cleanup_never_moves_a_destination_changed_after_validation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, initially_present: bool
) -> None:
    """A newly inserted filename belongs to the concurrent writer."""
    destination = tmp_path / "download.zip"
    if initially_present:
        destination.write_bytes(b"old stale result")
    original_validate = archive_safety._validate_destination

    def validate_then_insert(*args, **kwargs):
        result = original_validate(*args, **kwargs)
        if destination.exists():
            destination.unlink()
        destination.write_bytes(PATIENT_BYTES)
        return result

    monkeypatch.setattr(archive_safety, "_validate_destination", validate_then_insert)

    with pytest.raises(ValueError, match="changed after validation"):
        archive_safety.clear_stale_archive(tmp_path / "download", "zip")

    assert destination.read_bytes() == PATIENT_BYTES
    assert sorted(path.name for path in tmp_path.iterdir()) == ["download.zip"]


def test_archive_creation_never_cleans_a_destination_that_appears_after_validation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Initial cleanup refuses a new public entry instead of treating it as stale."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    destination = tmp_path / "download.zip"
    original_validate = archive_safety._validate_destination

    def validate_then_insert(*args, **kwargs):
        result = original_validate(*args, **kwargs)
        destination.write_bytes(PATIENT_BYTES)
        return result

    monkeypatch.setattr(archive_safety, "_validate_destination", validate_then_insert)

    with pytest.raises(ValueError, match="changed after validation"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert destination.read_bytes() == PATIENT_BYTES
    assert not [path for path in tmp_path.iterdir() if path.name.startswith(".download")]


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


def test_fifo_replacement_before_source_open_is_nonblocking_and_rejected(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A regular source replaced by a FIFO cannot hang before descriptor validation."""
    root = tmp_path / "results"
    root.mkdir()
    result = root / "result.txt"
    result.write_bytes(b"safe result")
    original_open = archive_safety.os.open
    replaced = False

    def replace_with_fifo_before_open(path, flags, *args, **kwargs):
        nonlocal replaced
        if str(path) == result.name and kwargs.get("dir_fd") is not None and not replaced:
            assert flags & os.O_NONBLOCK
            replaced = True
            result.unlink()
            os.mkfifo(result)
        return original_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(archive_safety.os, "open", replace_with_fifo_before_open)

    with pytest.raises(ValueError, match="unsupported filesystem entry 'result.txt'"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()


def test_temporary_cleanup_failure_is_logged_without_masking_primary_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A retained partial is named while the archive failure remains the outcome."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    monkeypatch.setattr(
        archive_safety.os, "link", lambda *args, **kwargs: (_ for _ in ()).throw(OSError("install failed"))
    )
    original_unlink = archive_safety.os.unlink

    def fail_temporary_cleanup(path, *args, **kwargs):
        if str(path).startswith(".download.zip.archive-"):
            raise OSError("cleanup failed")
        return original_unlink(path, *args, **kwargs)

    monkeypatch.setattr(archive_safety.os, "unlink", fail_temporary_cleanup)
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.archive_safety")

    with pytest.raises(OSError, match="install failed"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert "cleanup failed" in caplog.text
    assert ".download.zip.archive-" in caplog.text


def test_concurrent_archive_parent_replacement_never_touches_the_attacker_directory(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Validation and install must stay on one opened parent directory inode."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"safe result")
    public_parent = tmp_path / "public"
    public_parent.mkdir()
    original_parent = tmp_path / "original-public"
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    original_open = archive_safety.os.open
    original_mkdtemp = tempfile.mkdtemp
    replaced = False

    def replace_parent() -> None:
        nonlocal replaced
        if not replaced:
            replaced = True
            public_parent.rename(original_parent)
            public_parent.mkdir()
            os.link(patient, public_parent / "download.zip")

    def replace_after_parent_open(path, flags, *args, **kwargs):
        descriptor = original_open(path, flags, *args, **kwargs)
        if Path(path) == public_parent:
            replace_parent()
        return descriptor

    def replace_before_path_based_temp(*args, **kwargs):
        replace_parent()
        return original_mkdtemp(*args, **kwargs)

    monkeypatch.setattr(archive_safety.os, "open", replace_after_parent_open)
    monkeypatch.setattr(archive_safety, "tempfile", tempfile, raising=False)
    monkeypatch.setattr(tempfile, "mkdtemp", replace_before_path_based_temp)

    with pytest.raises((OSError, ValueError), match="archive parent.*changed"):
        archive_safety.create_safe_archive(public_parent / "download", "zip", root)

    attacker_destination = public_parent / "download.zip"
    assert os.path.samefile(attacker_destination, patient)
    assert attacker_destination.read_bytes() == PATIENT_BYTES
    assert patient.read_bytes() == PATIENT_BYTES
    assert not (original_parent / "download.zip").exists()


def test_stale_clear_rechecks_parent_after_the_anchored_operation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A parent swap during stale removal cannot be reported as a successful clear."""
    public_parent = tmp_path / "public"
    public_parent.mkdir()
    stale = public_parent / "download.zip"
    stale.write_bytes(b"stale archive")
    original_parent = tmp_path / "original-public"
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    original_clear = archive_safety._clear_stale_at

    def clear_then_swap(
        parent_descriptor: int, destination_name: str, expected_metadata: os.stat_result | None
    ) -> None:
        original_clear(parent_descriptor, destination_name, expected_metadata)
        public_parent.rename(original_parent)
        public_parent.mkdir()
        os.link(patient, public_parent / destination_name)

    monkeypatch.setattr(archive_safety, "_clear_stale_at", clear_then_swap)

    with pytest.raises((OSError, ValueError), match="archive parent.*changed"):
        archive_safety.clear_stale_archive(public_parent / "download", "zip")

    attacker_destination = public_parent / "download.zip"
    assert os.path.samefile(attacker_destination, patient)
    assert attacker_destination.read_bytes() == PATIENT_BYTES
    assert patient.read_bytes() == PATIENT_BYTES
    assert not (original_parent / "download.zip").exists()


def test_parent_swap_after_public_link_rolls_back_the_anchored_archive(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A post-link parent swap fails and removes the archive from the opened parent."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"safe result")
    public_parent = tmp_path / "public"
    public_parent.mkdir()
    original_parent = tmp_path / "original-public"
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    original_link = archive_safety.os.link

    def link_then_swap(source, target, *args, **kwargs):
        result = original_link(source, target, *args, **kwargs)
        if str(target) == "download.zip" and kwargs.get("dst_dir_fd") is not None:
            public_parent.rename(original_parent)
            public_parent.mkdir()
            original_link(patient, public_parent / "download.zip")
        return result

    monkeypatch.setattr(archive_safety.os, "link", link_then_swap)

    with pytest.raises((OSError, ValueError), match="archive parent.*changed"):
        archive_safety.create_safe_archive(public_parent / "download", "zip", root)

    attacker_destination = public_parent / "download.zip"
    assert os.path.samefile(attacker_destination, patient)
    assert attacker_destination.read_bytes() == PATIENT_BYTES
    assert patient.read_bytes() == PATIENT_BYTES
    assert not (original_parent / "download.zip").exists()
    assert not [path for path in original_parent.iterdir() if path.name.startswith(".download")]


def test_descriptor_cleanup_failure_occurs_before_public_archive_install(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """No cleanup that can fail remains after the archive becomes public."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    root_descriptor = -1
    original_open = archive_safety.os.open
    original_close = archive_safety.os.close

    def record_root_open(path, flags, *args, **kwargs):
        nonlocal root_descriptor
        descriptor = original_open(path, flags, *args, **kwargs)
        if Path(path) == root:
            root_descriptor = descriptor
        return descriptor

    def fail_after_root_close(descriptor: int) -> None:
        original_close(descriptor)
        if descriptor == root_descriptor:
            raise OSError("descriptor cleanup failed")

    monkeypatch.setattr(archive_safety.os, "open", record_root_open)
    monkeypatch.setattr(archive_safety.os, "close", fail_after_root_close)

    with pytest.raises(OSError, match="descriptor cleanup failed"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()


def test_post_install_temporary_unlink_failure_rolls_back_public_archive(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failure after the no-replace link removes the public archive again."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    original_unlink = archive_safety.os.unlink
    failed_once = False

    def fail_first_temporary_unlink(path, *args, **kwargs):
        nonlocal failed_once
        if str(path).startswith(".download.zip.archive-") and not failed_once:
            failed_once = True
            raise OSError("temporary unlink denied")
        return original_unlink(path, *args, **kwargs)

    monkeypatch.setattr(archive_safety.os, "unlink", fail_first_temporary_unlink)

    with pytest.raises(OSError, match="temporary unlink denied"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert not (tmp_path / "download.zip").exists()
    assert not [path for path in tmp_path.iterdir() if path.name.startswith(".download")]


def test_archive_return_refuses_a_public_name_replaced_after_install(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The returned filename must still name the exact completed archive inode."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    destination = tmp_path / "download.zip"
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    original_unlink = archive_safety.os.unlink
    replaced = False

    def replace_public_after_temporary_unlink(path, *args, **kwargs):
        nonlocal replaced
        result = original_unlink(path, *args, **kwargs)
        if str(path).startswith(".download.zip.archive-") and not replaced:
            replaced = True
            destination.unlink()
            os.link(patient, destination)
        return result

    monkeypatch.setattr(archive_safety.os, "unlink", replace_public_after_temporary_unlink)

    with pytest.raises(ValueError, match="changed after install"):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root)

    assert os.path.samefile(destination, patient)
    assert destination.read_bytes() == PATIENT_BYTES
    assert patient.read_bytes() == PATIENT_BYTES


def test_concurrent_protected_destination_insertion_is_never_replaced(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Atomic install must refuse a protected inode inserted after validation."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_bytes(b"result")
    patient = tmp_path / "patient.bam"
    patient.write_bytes(PATIENT_BYTES)
    destination = tmp_path / "download.zip"
    original_replace = archive_safety.os.replace
    original_link = archive_safety.os.link
    inserted = False

    def insert_before_replace(source, target, *args, **kwargs):
        nonlocal inserted
        if not inserted and str(target).endswith("download.zip"):
            inserted = True
            original_link(patient, destination)
        return original_replace(source, target, *args, **kwargs)

    def insert_before_noreplace_link(source, target, *args, **kwargs):
        nonlocal inserted
        if not inserted and str(target).endswith("download.zip"):
            inserted = True
            original_link(patient, destination)
        return original_link(source, target, *args, **kwargs)

    monkeypatch.setattr(archive_safety.os, "replace", insert_before_replace)
    monkeypatch.setattr(archive_safety.os, "link", insert_before_noreplace_link)

    with pytest.raises((FileExistsError, ValueError)):
        archive_safety.create_safe_archive(tmp_path / "download", "zip", root, protected_paths=(patient,))

    assert os.path.samefile(destination, patient)
    assert destination.read_bytes() == PATIENT_BYTES
    assert patient.read_bytes() == PATIENT_BYTES


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
