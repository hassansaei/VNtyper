"""Resolving an alignment's existing BAI, under either name it can carry (#210).

`vntyper/scripts/alignment_index.py` answers one question - given an alignment, is there
already an index this pipeline can read, and where - because the alternative is to build a
second index that nothing else knows about. A BAI has two spellings, `<file>.bai` and
`<stem>.bai`; the pipeline used to reconstruct only the first, so the `sample.bai` that the
upload endpoint and the worker both deliberately accept was invisible to it.

**It is not htslib's resolution order and does not claim to be.** htslib tries CSI, in both
spellings, before BAI. The returned path goes straight into
`extract_unmapped_from_offset.get_last_chunk_end`, which parses the BAI container itself
and rejects anything else, so ignoring a CSI is the correct behaviour rather than a missing
feature - and is pinned below in both directions: the resolution ignores it, and the reader
downstream would reject it.

These tests moved here with `resolve_bam_index` when it came out of
`fastq_bam_processing.py`. What that resolution is *for* - no index is ever written into
the read-only input directory - stays in `test_input_tree_is_never_written.py`, which
drives the whole stage; these pin the rule itself.
"""

from pathlib import Path

import pytest

from vntyper.scripts.alignment_index import resolve_any_index, resolve_bam_index

pytestmark = pytest.mark.unit


def test_cram_index_uses_appended_spelling_before_stem_spelling(tmp_path: Path) -> None:
    """CRAM index resolution returns the first existing CRAI candidate."""
    cram = tmp_path / "sample.cram"
    cram.write_bytes(b"CRAM\x02")
    appended = tmp_path / "sample.cram.crai"
    appended.write_bytes(b"CRAI")
    stem = tmp_path / "sample.crai"
    stem.write_bytes(b"CRAI")

    assert resolve_any_index(cram, "cram") == str(appended)


def test_cram_index_uses_stem_spelling_when_appended_is_missing(tmp_path: Path) -> None:
    """A CRAM's stem-based CRAI spelling is found when it is the only candidate."""
    cram = tmp_path / "sample.cram"
    cram.write_bytes(b"CRAM\x02")
    stem = tmp_path / "sample.crai"
    stem.write_bytes(b"CRAI")

    assert resolve_any_index(cram, "cram") == str(stem)


def test_bai_beside_cram_is_not_a_cram_index(tmp_path: Path) -> None:
    """A BAI next to a CRAM does not satisfy CRAI resolution."""
    cram = tmp_path / "sample.cram"
    cram.write_bytes(b"CRAM\x02")
    (tmp_path / "sample.bai").write_bytes(b"BAI\x01")

    assert resolve_any_index(cram, "cram") is None


def test_bam_csi_is_found_by_any_index_but_not_bai_index(tmp_path: Path) -> None:
    """General resolution accepts CSI while the BAI-only reader does not."""
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    csi = tmp_path / "sample.bam.csi"
    csi.write_bytes(b"CSI\x01")

    assert resolve_any_index(bam, "bam") == str(csi)
    assert resolve_bam_index(bam) is None


def test_no_supported_index_resolves_to_none(tmp_path: Path) -> None:
    """No existing candidate is reported as unavailable."""
    cram = tmp_path / "sample.cram"
    cram.write_bytes(b"CRAM\x02")

    assert resolve_any_index(cram, "cram") is None


def test_an_alternate_index_name_is_found_and_no_second_index_is_built(tmp_path: Path) -> None:
    """#210: the upload path accepts ``sample.bai``; htslib resolves it, so must we.

    Args:
        tmp_path: Pytest temporary directory.
    """
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    alternate = tmp_path / "sample.bai"
    alternate.write_bytes(b"BAI\x01")

    assert resolve_bam_index(bam) == str(alternate)


def test_the_conventional_index_name_wins_over_the_alternate(tmp_path: Path) -> None:
    """htslib tries ``<file>.bai`` first.

    Args:
        tmp_path: Pytest temporary directory.
    """
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    (tmp_path / "sample.bam.bai").write_bytes(b"BAI\x01")
    (tmp_path / "sample.bai").write_bytes(b"BAI\x01")

    assert resolve_bam_index(bam) == str(tmp_path / "sample.bam.bai")


def test_no_index_at_all_resolves_to_none(tmp_path: Path) -> None:
    """Nothing to reuse is reported as such, so the caller builds one.

    Args:
        tmp_path: Pytest temporary directory.
    """
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")

    assert resolve_bam_index(bam) is None


def test_a_csi_index_is_ignored_and_a_bai_is_built_instead(tmp_path: Path) -> None:
    """The BAI-only limitation, pinned so it is deliberate rather than accidental.

    htslib resolves CSI before BAI, in both the appended (``sample.bam.csi``) and the
    substituted (``sample.csi``) spelling, so this function is **not** htslib-equivalent
    and must not claim to be. Reusing a CSI here would be a defect, not a feature: the
    only consumer of the returned path is
    ``extract_unmapped_from_offset.get_last_chunk_end``, which parses the BAI container
    directly and rejects anything whose first four bytes are not ``BAI\\x01``. Returning
    a CSI would turn a working run into a ``ValueError`` mid-stage; returning None makes
    the caller build the BAI it needs, into the run's output directory.

    This test **passes before and after** the docstring correction that accompanies it -
    the behaviour was already right and is what is being pinned. What was wrong, and what
    ``test_the_docstring_does_not_claim_htslib_parity`` covers, is the claim above it.

    Args:
        tmp_path: Pytest temporary directory.
    """
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    (tmp_path / "sample.bam.csi").write_bytes(b"CSI\x01")
    (tmp_path / "sample.csi").write_bytes(b"CSI\x01")

    assert resolve_bam_index(bam) is None


def test_the_bai_reader_downstream_really_does_reject_a_csi(tmp_path: Path) -> None:
    """Why BAI-only is the correct behaviour rather than a gap, stated as a test.

    Args:
        tmp_path: Pytest temporary directory.
    """
    from vntyper.scripts.extract_unmapped_from_offset import get_last_chunk_end

    csi = tmp_path / "sample.bam.csi"
    csi.write_bytes(b"CSI\x01" + b"\x00" * 64)

    with pytest.raises(ValueError):
        get_last_chunk_end(str(csi))


def test_the_docstring_does_not_claim_htslib_parity() -> None:
    """The limitation belongs where the function is read, not only in its tests.

    The docstring opened "Find an existing BAM index the way htslib itself does", and a
    reader taking that at face value would expect a CSI to be honoured. It is not, on
    purpose. A claim of parity is the kind of thing that gets acted on - by widening the
    candidate list to match it - so the text is pinned here rather than left to drift.
    """
    doc = resolve_bam_index.__doc__ or ""

    assert "the way htslib itself does" not in doc, "the parity claim is false and must not come back"
    assert "CSI" in doc, "the docstring must say which index formats are deliberately not resolved"
