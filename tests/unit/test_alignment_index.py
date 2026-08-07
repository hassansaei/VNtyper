"""Resolving an alignment's existing index the way htslib does (#210).

`vntyper/scripts/alignment_index.py` answers one question - given an alignment, is there
already an index for it, and where - and it has to answer it the way htslib does, because
the pipeline's alternative is to build a second index that nothing else knows about.
htslib tries `<file>.bai` and then `<stem>.bai`; the pipeline used to reconstruct only the
first, so the `sample.bai` that the upload endpoint and the worker both deliberately
accept was invisible to it.

These tests moved here with `resolve_bam_index` when it came out of
`fastq_bam_processing.py`. What that resolution is *for* - no index is ever written into
the read-only input directory - stays in `test_input_tree_is_never_written.py`, which
drives the whole stage; these three pin the rule itself.
"""

from pathlib import Path

import pytest

from vntyper.scripts.alignment_index import resolve_bam_index

pytestmark = pytest.mark.unit


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
