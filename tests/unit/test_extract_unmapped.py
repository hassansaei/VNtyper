"""Unit tests for extract_unmapped_from_offset.py.

get_last_chunk_end walks a BAI index by hand to find where the mapped chunks
end, and the unmapped scan starts there. If it returns too high an offset,
reads are silently skipped and the genotype is computed on a fraction of the
input -- with no error anywhere.
"""

import io
import struct
from unittest.mock import MagicMock

import pysam
import pytest

from vntyper.scripts.extract_unmapped_from_offset import (
    extract_unmapped_reads_from_offset,
    get_last_chunk_end,
    read_uint32,
    read_uint64,
)

pytestmark = pytest.mark.unit


def _bai(references):
    """Build a minimal BAI. `references` is a list of lists of (beg, end) chunks.

    Matches, byte for byte, the order `get_last_chunk_end` walks: a 4-byte
    magic (unparsed, just skipped), then n_ref (uint32), then per reference:
    n_bins (uint32); for each bin, a bin number (uint32, read but unused), a
    chunk count (uint32) and that many (chunk_beg, chunk_end) uint64 pairs;
    then n_intv (uint32) and n_intv * 8 bytes of linear index that the real
    parser skips with `bai.seek` rather than reading.

    All chunks for one reference are written into a single synthetic bin
    here (n_bins is 0 or 1). That is a simplification versus a real BAI,
    which spreads chunks across many bins, but it is invisible to
    `get_last_chunk_end`: it maxes `chunk_end` over every chunk in every
    bin, so how those chunks are grouped into bins does not affect the
    result.
    """
    out = [b"BAI\x01", struct.pack("<I", len(references))]
    for chunks in references:
        out.append(struct.pack("<I", 1 if chunks else 0))  # n_bins
        if chunks:
            out.append(struct.pack("<I", 4681))  # bin number, unused by the parser
            out.append(struct.pack("<I", len(chunks)))  # n_chunks
            for beg, end in chunks:
                out.append(struct.pack("<QQ", beg, end))
        out.append(struct.pack("<I", 2))  # n_intv: a 2-entry linear index
        out.append(struct.pack("<QQ", 0, 0))  # the linear index itself, skipped over
    return b"".join(out)


# ---------------------------------------------------------------------------
# read_uint32 / read_uint64
# ---------------------------------------------------------------------------


def test_read_uint32_is_little_endian_and_unsigned():
    assert read_uint32(io.BytesIO(b"\x01\x00\x00\x00")) == 1
    assert read_uint32(io.BytesIO(b"\xff\xff\xff\xff")) == 4294967295


def test_read_uint64_is_little_endian_and_unsigned():
    assert read_uint64(io.BytesIO(b"\x01" + b"\x00" * 7)) == 1
    assert read_uint64(io.BytesIO(b"\xff" * 8)) == 18446744073709551615


# ---------------------------------------------------------------------------
# get_last_chunk_end
# ---------------------------------------------------------------------------


def test_the_maximum_chunk_end_is_returned_not_the_last_one(tmp_path):
    """Chunks are not sorted by end offset; taking the last would truncate."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900), (100, 500)]]))
    assert get_last_chunk_end(str(path)) == 900


def test_the_maximum_is_taken_across_every_reference(tmp_path):
    """The linear index of reference 1 must be skipped correctly to reach 2."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 100)], [(0, 7000)], [(0, 300)]]))
    assert get_last_chunk_end(str(path)) == 7000


def test_a_reference_with_no_bins_contributes_nothing(tmp_path):
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[], [(0, 42)]]))
    assert get_last_chunk_end(str(path)) == 42


def test_an_index_with_no_chunks_at_all_returns_zero(tmp_path):
    """max_vo starts at 0; an empty index must not raise."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[]]))
    assert get_last_chunk_end(str(path)) == 0


def test_a_truncated_index_does_not_silently_return_a_low_offset_today(tmp_path):
    """CHARACTERISATION of a live fail-open. Do not "fix" this here.

    `read_uint32`/`read_uint64` do `int.from_bytes(f.read(n), ...)`. At EOF,
    `f.read(n)` returns `b""` and `int.from_bytes(b"", ...)` is 0 rather than
    raising. Truncate a BAI to just its 4-byte magic and `n_ref =
    read_uint32(bai)` silently reads 0 instead of the real reference count,
    so `get_last_chunk_end` believes the index describes zero references and
    returns 0 -- even though the untruncated index (built from the same
    `_bai([[(0, 900)]])` call) records a chunk ending at offset 900.

    The caller (`extract_unmapped_reads_from_offset`) then seeks the input
    BAM to offset 0 and scans the entire file rather than raising for a
    corrupt/truncated index -- reads are not silently dropped in *this*
    direction, but a real-world truncation partway through a large index
    (not just after the magic) would return some offset lower than the true
    one, silently dropping whatever mapped-chunk data existed past the
    truncation point, with no error anywhere. This test pins today's
    behaviour so it cannot drift further by accident; whether a truncated
    BAI *should* raise is a separate product decision, not made here (see
    issue-truncated-bai-fail-open.md).
    """
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:4])
    assert get_last_chunk_end(str(path)) == 0


# ---------------------------------------------------------------------------
# extract_unmapped_reads_from_offset
# ---------------------------------------------------------------------------


def _mock_alignment_files(monkeypatch, reads):
    """Wire `pysam.AlignmentFile` mocks matching the module's real call shape.

    `extract_unmapped_reads_from_offset` calls `pysam.AlignmentFile` twice in
    sequence: once for the input BAM ("rb"), once for the output BAM ("wb",
    opened with `header=`, not `template=`). Both are used as context
    managers. The input is iterated directly with a plain `for read in
    inbam:` -- it does not call `.fetch()`.
    """
    inbam = MagicMock(name="inbam")
    outbam = MagicMock(name="outbam")
    inbam.__enter__.return_value = inbam
    inbam.__exit__.return_value = False
    outbam.__enter__.return_value = outbam
    outbam.__exit__.return_value = False
    inbam.__iter__.return_value = iter(reads)
    inbam.header = MagicMock(name="header")
    monkeypatch.setattr(pysam, "AlignmentFile", MagicMock(side_effect=[inbam, outbam]))
    return inbam, outbam


def test_the_scan_seeks_to_the_offset_the_index_reported(tmp_path, monkeypatch):
    """The whole point of the module: start where the mapped chunks ended."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 7000)]]))
    inbam, _outbam = _mock_alignment_files(monkeypatch, reads=[])

    extract_unmapped_reads_from_offset(str(tmp_path / "in.bam"), str(path), str(tmp_path / "out.bam"))

    inbam.seek.assert_called_once_with(7000)


def test_only_unmapped_reads_are_written(tmp_path, monkeypatch):
    """A mapped read reaching the output would double-count it downstream."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[]]))
    mapped = MagicMock(is_unmapped=False)
    unmapped = MagicMock(is_unmapped=True)
    _inbam, outbam = _mock_alignment_files(monkeypatch, reads=[mapped, unmapped])

    extract_unmapped_reads_from_offset(str(tmp_path / "in.bam"), str(path), str(tmp_path / "out.bam"))

    outbam.write.assert_called_once_with(unmapped)


def test_the_output_file_is_opened_with_the_input_files_header(tmp_path, monkeypatch):
    """`header=inbam.header`, not `template=` -- get the kwarg name right."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[]]))
    inbam, _outbam = _mock_alignment_files(monkeypatch, reads=[])

    extract_unmapped_reads_from_offset(str(tmp_path / "in.bam"), str(path), str(tmp_path / "out.bam"))

    _args, kwargs = pysam.AlignmentFile.call_args_list[1]
    assert kwargs["header"] is inbam.header
    assert "template" not in kwargs


def test_both_files_are_closed_when_the_scan_raises(tmp_path, monkeypatch):
    """`with` must cover the write side too, or a failed run leaks a handle
    and leaves a truncated BAM that looks complete."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[]]))

    def _raising_reads():
        yield MagicMock(is_unmapped=False)
        raise OSError("disk gone")

    inbam, outbam = _mock_alignment_files(monkeypatch, reads=[])
    inbam.__iter__.return_value = _raising_reads()

    with pytest.raises(OSError):
        extract_unmapped_reads_from_offset(str(tmp_path / "in.bam"), str(path), str(tmp_path / "out.bam"))

    inbam.__exit__.assert_called_once()
    outbam.__exit__.assert_called_once()
