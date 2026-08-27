"""Installation must be all-or-nothing, and an archive must not write outside its root.

The installer this replaces skipped existing destinations without revalidating them
(`install_references.download_file`), unpacked archives with an unrestricted
`tar.extractall`, and overwrote config.json in place. A network blip therefore left a
half-populated reference tree that the next run treated as complete.
"""

import io
import logging
import stat
import tarfile
import zipfile

import pytest

from vntyper.scripts import reference_bundle
from vntyper.scripts.reference_bundle import safe_extract, staged_install, verify_sha256

pytestmark = pytest.mark.unit


def _tar_with(tmp_path, members):
    """Build a tar whose members are (arcname, payload) pairs.

    Uses addfile rather than add: `TarFile.gettarinfo(arcname="/etc/passwd")` silently
    strips the leading slash, so an `add`-based fixture cannot produce the absolute-path
    member the traversal test is supposed to reject.
    """
    archive = tmp_path / "bundle.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for arcname, payload in members:
            info = tarfile.TarInfo(arcname)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return archive


class TestVerifySha256:
    def test_a_matching_digest_passes_silently(self, tmp_path):
        import hashlib

        target = tmp_path / "f"
        target.write_bytes(b"reference")
        verify_sha256(target, hashlib.sha256(b"reference").hexdigest())

    def test_a_mismatch_names_the_file_and_both_digests(self, tmp_path):
        target = tmp_path / "chr1.hg38.fa"
        target.write_bytes(b"tampered")
        with pytest.raises(ValueError, match="chr1.hg38.fa"):
            verify_sha256(target, "0" * 64)


class TestSafeExtract:
    def test_an_ordinary_archive_extracts(self, tmp_path):
        archive = _tar_with(tmp_path, [("alignment/chr1.hg38.fa", b">chr1\nACGT\n")])
        destination = tmp_path / "out"
        members = safe_extract(archive, destination)
        assert members == ["alignment/chr1.hg38.fa"]
        assert (destination / "alignment/chr1.hg38.fa").read_bytes() == b">chr1\nACGT\n"

    @pytest.mark.parametrize("arcname", ["../escaped.fa", "/etc/passwd", "a/../../escaped.fa"])
    def test_a_member_escaping_the_root_is_rejected_by_name(self, tmp_path, arcname):
        archive = _tar_with(tmp_path, [(arcname, b"x")])
        with pytest.raises(ValueError, match="escap|outside|absolute"):
            safe_extract(archive, tmp_path / "out")

    def test_a_hard_link_escaping_the_root_is_rejected(self, tmp_path):
        """Hard-link targets are archive-root-relative, unlike symlink targets."""
        archive = tmp_path / "hard.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            info = tarfile.TarInfo("nested/evil")
            info.type = tarfile.LNKTYPE
            info.linkname = "../../outside.fa"
            tar.addfile(info)
        with pytest.raises(ValueError):
            safe_extract(archive, tmp_path / "out")

    def test_a_symlink_pointing_outside_the_root_is_rejected(self, tmp_path):
        archive = tmp_path / "link.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            info = tarfile.TarInfo("evil")
            info.type = tarfile.SYMTYPE
            info.linkname = "/etc/passwd"
            tar.addfile(info)
        with pytest.raises(ValueError):
            safe_extract(archive, tmp_path / "out")

    def test_a_link_resolving_inside_the_root_is_rejected_anyway(self, tmp_path):
        """Even a link whose target is unambiguously inside the root is refused.

        An earlier round asked for the opposite of this test, on the reasoning that a
        guard rejecting everything would pass every rejection test above. That reasoning
        is obsolete: resolving a link target before extraction is not safe at all (see
        `test_an_earlier_member_cannot_redirect_a_later_one`), so the guard now rejects
        links outright and this pins that. `test_an_ordinary_archive_extracts` is what
        proves the guard is not simply rejecting everything.
        """
        archive = tmp_path / "links.tar.gz"
        payload = b">chr1\nACGT\n"
        with tarfile.open(archive, "w:gz") as tar:
            regular = tarfile.TarInfo("alignment/chr1.hg38.fa")
            regular.size = len(payload)
            tar.addfile(regular, io.BytesIO(payload))

            symlink = tarfile.TarInfo("alignment/chr1.sym.fa")
            symlink.type = tarfile.SYMTYPE
            symlink.linkname = "chr1.hg38.fa"  # relative to the symlink's own directory
            tar.addfile(symlink)

        destination = tmp_path / "out"
        with pytest.raises(ValueError, match="is a link"):
            safe_extract(archive, destination)
        assert not (destination / "alignment/chr1.sym.fa").exists()

    def test_an_earlier_member_cannot_redirect_a_later_one(self, tmp_path, monkeypatch):
        """Regression: validating every member up front, then extracting, is not enough.

        Members in archive order:

        1. `pivot` - a symlink to `.`
        2. `escape` - a symlink to `pivot/..`
        3. `escape/owned` - a regular file

        A pre-extraction check on link targets clears all three: `destination/pivot/..`
        normalises lexically to `destination`, and the regular member has neither a `..`
        component nor an absolute path. Once extraction has created `pivot`, `escape`
        resolves to `destination.parent` and the regular member is written *through* it.

        `_EXTRACT_KWARGS` is emptied so this exercises the member loop rather than
        `tarfile`'s `filter="data"`. Without that, the test would pass on this
        interpreter for the wrong reason and would not have caught the bug - the real
        `_EXTRACT_KWARGS` is `{}` on a 3.10 older than 3.10.12, which `requires-python`
        allows.
        """
        import vntyper.scripts.reference_bundle as reference_bundle

        monkeypatch.setattr(reference_bundle, "_EXTRACT_KWARGS", {})

        archive = tmp_path / "ordering.tar.gz"
        payload = b"OWNED"
        with tarfile.open(archive, "w:gz") as tar:
            pivot = tarfile.TarInfo("pivot")
            pivot.type = tarfile.SYMTYPE
            pivot.linkname = "."
            tar.addfile(pivot)

            escape = tarfile.TarInfo("escape")
            escape.type = tarfile.SYMTYPE
            escape.linkname = "pivot/.."
            tar.addfile(escape)

            owned = tarfile.TarInfo("escape/owned")
            owned.size = len(payload)
            tar.addfile(owned, io.BytesIO(payload))

        destination = tmp_path / "root" / "out"
        with pytest.raises(ValueError, match="is a link"):
            safe_extract(archive, destination)

        assert not (destination.parent / "owned").exists(), "a member was written outside the destination"
        assert not (destination / "pivot").exists()

    def test_a_preexisting_destination_symlink_cannot_redirect_a_member(self, tmp_path, monkeypatch):
        """Regression (MAJOR 4, milestone-5 PR-2 review): the symlink does not have to be
        in the archive at all - one already sitting in the destination is enough.

        On Python 3.10.0-3.10.11 and 3.11.0-3.11.3, `tarfile`'s `filter="data"` is
        unavailable and this member loop is the *only* defence. It validates member
        names and types, but - before this fix - never resolved the destination path
        through a symlink that already existed there. `staged_install` seeds staging
        from the existing tree with `symlinks=True` (see its own docstring), so an
        `alignment` symlink pointing outside the reference root survives into a fresh
        install's staging directory, and an ordinary, fully-verified member such as
        `alignment/chr1.hg19.fa` would be written straight through it - landing outside
        `destination` even though every per-member check (no `..`, not absolute, not a
        link itself, a regular file) passes.

        `_EXTRACT_KWARGS` is emptied for the same reason as
        `test_an_earlier_member_cannot_redirect_a_later_one`: the real value is `{}` on
        an interpreter `requires-python` still allows, and the member loop has to hold
        the guarantee entirely on its own.
        """
        import vntyper.scripts.reference_bundle as reference_bundle

        monkeypatch.setattr(reference_bundle, "_EXTRACT_KWARGS", {})

        outside = tmp_path / "outside"
        outside.mkdir()
        destination = tmp_path / "out"
        destination.mkdir()
        (destination / "alignment").symlink_to(outside, target_is_directory=True)

        archive = tmp_path / "preexisting_symlink.tar.gz"
        payload = b">chr1\nACGT\n"
        with tarfile.open(archive, "w:gz") as tar:
            member = tarfile.TarInfo("alignment/chr1.hg19.fa")
            member.size = len(payload)
            tar.addfile(member, io.BytesIO(payload))

        with pytest.raises(ValueError, match="alignment/chr1.hg19.fa"):
            safe_extract(archive, destination)

        assert not (outside / "chr1.hg19.fa").exists(), "a member was written outside the destination"
        assert (destination / "alignment").is_symlink(), "the pre-existing symlink itself must be left alone"

    def test_the_destination_itself_cannot_be_a_symlink(self, tmp_path):
        outside = tmp_path / "outside"
        outside.mkdir()
        destination = tmp_path / "out"
        destination.symlink_to(outside, target_is_directory=True)
        archive = _tar_with(tmp_path, [("chr1.fa", b">chr1\nACGT\n")])

        with pytest.raises(ValueError, match="destination .* is a symlink"):
            safe_extract(archive, destination)

        assert not (outside / "chr1.fa").exists()

    def test_a_device_member_is_rejected(self, tmp_path):
        """Neither a regular file nor a directory, so it belongs in no reference bundle."""
        archive = tmp_path / "device.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            info = tarfile.TarInfo("dev/null")
            info.type = tarfile.CHRTYPE
            info.devmajor = 1
            info.devminor = 3
            tar.addfile(info)
        with pytest.raises(ValueError, match="not a regular file or directory"):
            safe_extract(archive, tmp_path / "out")

    def test_the_data_filter_refusing_the_archive_surfaces_as_a_value_error(self, tmp_path, monkeypatch):
        """`filter="data"` raises `tarfile.FilterError`, which is not a `ValueError`.

        Every caller is documented to see `ValueError`, and the member loop cannot be
        the only thing that honours that: on an interpreter carrying the filter, a
        refusal from `tarfile` itself must be translated or the contract holds on some
        Python versions and not others.
        """
        import vntyper.scripts.reference_bundle as reference_bundle

        if not reference_bundle._FILTER_ERRORS:
            pytest.skip("this interpreter predates tarfile's extraction filter")

        archive = _tar_with(tmp_path, [("alignment/chr1.hg38.fa", b">chr1\nACGT\n")])

        def refuse(self, *_args, **_kwargs):
            raise reference_bundle._FILTER_ERRORS[0](tarfile.TarInfo("alignment/chr1.hg38.fa"), "/elsewhere")

        monkeypatch.setattr(tarfile.TarFile, "extractall", refuse)
        with pytest.raises(ValueError, match="extraction filter refused"):
            safe_extract(archive, tmp_path / "out")


class TestSafeExtractZip:
    """ZIP extraction must enforce the same archive-root boundary as tar."""

    @staticmethod
    def _zip_with(tmp_path, name, payload=b"x"):
        archive = tmp_path / "seed.zip"
        with zipfile.ZipFile(archive, "w") as bundle:
            bundle.writestr(name, payload)
        return archive

    def test_a_member_with_a_parent_component_is_refused(self, tmp_path):
        archive = self._zip_with(tmp_path, "../escape.txt")

        with pytest.raises(ValueError, match="escapes the archive root"):
            reference_bundle.safe_extract_zip(archive, tmp_path / "dest")

        assert not (tmp_path / "escape.txt").exists()

    def test_an_absolute_member_is_refused(self, tmp_path):
        archive = self._zip_with(tmp_path, "/etc/owned")

        with pytest.raises(ValueError, match="absolute path"):
            reference_bundle.safe_extract_zip(archive, tmp_path / "dest")

    def test_a_link_member_is_refused(self, tmp_path):
        archive = tmp_path / "seed.zip"
        info = zipfile.ZipInfo("evil-link")
        info.external_attr = (stat.S_IFLNK | 0o777) << 16
        with zipfile.ZipFile(archive, "w") as bundle:
            bundle.writestr(info, "/etc")

        with pytest.raises(ValueError, match="is a link"):
            reference_bundle.safe_extract_zip(archive, tmp_path / "dest")

    def test_a_preexisting_destination_symlink_cannot_redirect_a_member(self, tmp_path):
        outside = tmp_path / "outside"
        outside.mkdir()
        destination = tmp_path / "dest"
        destination.mkdir()
        (destination / "data").symlink_to(outside, target_is_directory=True)
        archive = self._zip_with(tmp_path, "data/file.txt")

        with pytest.raises(ValueError, match="outside"):
            reference_bundle.safe_extract_zip(archive, destination)

        assert not (outside / "file.txt").exists()

    def test_the_destination_itself_cannot_be_a_symlink(self, tmp_path):
        outside = tmp_path / "outside"
        outside.mkdir()
        destination = tmp_path / "dest"
        destination.symlink_to(outside, target_is_directory=True)
        archive = self._zip_with(tmp_path, "file.txt")

        with pytest.raises(ValueError, match="destination .* is a symlink"):
            reference_bundle.safe_extract_zip(archive, destination)

        assert not (outside / "file.txt").exists()

    def test_a_clean_archive_extracts_and_returns_its_file_members(self, tmp_path):
        archive = tmp_path / "seed.zip"
        with zipfile.ZipFile(archive, "w") as bundle:
            bundle.writestr("a.txt", "a")
            bundle.writestr("sub/", b"")
            bundle.writestr("sub/b.txt", "b")

        members = reference_bundle.safe_extract_zip(archive, tmp_path / "dest")

        assert members == ["a.txt", "sub/b.txt"]
        assert (tmp_path / "dest" / "sub" / "b.txt").read_text(encoding="utf-8") == "b"


class TestStagedInstall:
    def test_an_existing_installation_is_carried_into_staging(self, tmp_path):
        """Installing hg38 after hg19 must not delete hg19, or the retained README."""
        target = tmp_path / "refs"
        (target / "alignment").mkdir(parents=True)
        (target / "alignment" / "chr1.hg19.fa").write_text("hg19")
        (target / "README.md").write_text("kept")
        with staged_install(target) as staging:
            (staging / "alignment" / "chr1.hg38.fa").write_text("hg38")
        assert (target / "alignment" / "chr1.hg19.fa").read_text() == "hg19"
        assert (target / "alignment" / "chr1.hg38.fa").read_text() == "hg38"
        assert (target / "README.md").read_text() == "kept"

    def test_a_seeding_failure_leaves_no_staging_directory_behind(self, tmp_path, monkeypatch):
        """Seeding only reads `target`, but a mid-copy failure - disk full, permission
        error, an interrupted copy of a large reference tree - must still hit the
        cleanup path, or every failed retry leaks another staging directory."""
        import vntyper.scripts.reference_bundle as reference_bundle

        target = tmp_path / "refs"
        target.mkdir()
        (target / "chr1.fa").write_text("original")

        def failing_copytree(*_args, **_kwargs):
            raise OSError("no space left on device")

        monkeypatch.setattr(reference_bundle.shutil, "copytree", failing_copytree)
        with pytest.raises(OSError, match="no space left on device"), staged_install(target):
            pass

        assert list(tmp_path.glob(".refs.staging.*")) == []
        assert (target / "chr1.fa").read_text() == "original"

    def test_a_successful_stage_is_activated_atomically(self, tmp_path):
        target = tmp_path / "refs"
        with staged_install(target) as staging:
            (staging / "chr1.fa").write_text(">chr1")
            assert not target.exists(), "target must not appear until activation"
        assert (target / "chr1.fa").read_text() == ">chr1"

    def test_a_failure_leaves_no_partial_tree_behind(self, tmp_path):
        target = tmp_path / "refs"
        with pytest.raises(RuntimeError), staged_install(target) as staging:
            (staging / "half.fa").write_text("x")
            raise RuntimeError("network died mid-extract")
        assert not target.exists()

    def test_a_previous_tree_is_preserved_when_it_cannot_be_restored(self, tmp_path, monkeypatch):
        """Losing the only surviving copy is worse than leaving a stray directory.

        Activation is the one window in which the previous tree is not at `target`.
        When activation fails there and the rename back fails too, the moved-aside
        tree is the only copy that exists, so it stays on disk under a named path and
        is reported - never deleted.
        """
        import pathlib

        target = tmp_path / "refs"
        target.mkdir()
        (target / "chr1.fa").write_text("original")
        real_rename = pathlib.Path.rename

        def block_renames_into_target(self, other):
            if pathlib.Path(other) == target:
                raise OSError("rename into target blocked")
            return real_rename(self, other)

        monkeypatch.setattr(pathlib.Path, "rename", block_renames_into_target)
        with pytest.raises(OSError, match="blocked"), staged_install(target) as staging:
            (staging / "chr1.fa").write_text("replacement")

        preserved = list(tmp_path.glob(".refs.previous.*"))
        assert preserved, "the previous tree must survive under a named path"
        assert (preserved[0] / "chr1.fa").read_text() == "original"
        assert not list(tmp_path.glob(".refs.staging.*")), "staging must still be cleaned up"

    def test_a_failed_quarantine_does_not_report_a_preservation_that_never_happened(
        self, tmp_path, monkeypatch, caplog
    ):
        """If `target.rename(previous)` itself is the rename that fails, the reserved
        `previous` name was already vacated by `rmdir()` and nothing was ever moved
        there - `target` still holds the tree safely. The operator report must not
        send them hunting for data that was never displaced.
        """
        import pathlib

        target = tmp_path / "refs"
        target.mkdir()
        (target / "chr1.fa").write_text("original")
        real_rename = pathlib.Path.rename

        def fail_quarantine(self, other):
            if pathlib.Path(self) == target:
                raise OSError("quarantine blocked")
            return real_rename(self, other)

        monkeypatch.setattr(pathlib.Path, "rename", fail_quarantine)
        with (
            caplog.at_level(logging.ERROR, logger="vntyper.scripts.reference_bundle"),
            pytest.raises(OSError, match="quarantine blocked"),
            staged_install(target) as staging,
        ):
            (staging / "chr1.fa").write_text("replacement")

        assert (target / "chr1.fa").read_text() == "original"
        assert list(tmp_path.glob(".refs.previous.*")) == []
        assert "preserved at" not in caplog.text, (
            "nothing was actually moved aside; the operator must not be told to go find it"
        )

    def test_a_failure_leaves_a_previous_installation_untouched(self, tmp_path):
        target = tmp_path / "refs"
        target.mkdir()
        (target / "chr1.fa").write_text("original")
        with pytest.raises(RuntimeError), staged_install(target) as staging:
            (staging / "chr1.fa").write_text("replacement")
            raise RuntimeError("boom")
        assert (target / "chr1.fa").read_text() == "original"
        assert not list(tmp_path.glob(".refs.previous.*")), (
            "a failure before activation must not move the live tree aside at all"
        )
