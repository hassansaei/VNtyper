"""Installation must be all-or-nothing, and an archive must not write outside its root.

The installer this replaces skipped existing destinations without revalidating them
(`install_references.download_file`), unpacked archives with an unrestricted
`tar.extractall`, and overwrote config.json in place. A network blip therefore left a
half-populated reference tree that the next run treated as complete.
"""

import io
import logging
import tarfile

import pytest

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
        safe_extract(archive, destination)
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

    def test_a_link_resolving_inside_the_root_is_extracted(self, tmp_path):
        """A legitimate in-root link must survive extraction, not just be rejected.

        The four tests above only prove the guard rejects bad members; none of them
        prove it lets a good one through, and a guard that rejects everything passes
        every one of them.
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

            hardlink = tarfile.TarInfo("alignment/chr1.hard.fa")
            hardlink.type = tarfile.LNKTYPE
            hardlink.linkname = "alignment/chr1.hg38.fa"  # relative to the archive root
            tar.addfile(hardlink)

        destination = tmp_path / "out"
        safe_extract(archive, destination)

        assert (destination / "alignment/chr1.hg38.fa").read_bytes() == payload
        assert (destination / "alignment/chr1.sym.fa").read_bytes() == payload
        assert (destination / "alignment/chr1.hard.fa").read_bytes() == payload


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
