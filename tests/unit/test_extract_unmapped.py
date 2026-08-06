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


def test_read_uint32_raises_on_a_short_read():
    """SPECIFICATION (issue-truncated-bai-fail-open.md, decision 1): `f.read(4)`
    returning fewer than 4 bytes at EOF must raise, not silently zero-fill via
    `int.from_bytes(b"", ...) == 0`. The message names both the expected and the
    actually-found byte count so a truncated index is diagnosable from the log."""
    with pytest.raises(ValueError, match=r"expected 4 bytes.*found 2"):
        read_uint32(io.BytesIO(b"\x01\x00"))


def test_read_uint32_raises_on_a_completely_empty_read():
    with pytest.raises(ValueError, match=r"expected 4 bytes.*found 0"):
        read_uint32(io.BytesIO(b""))


def test_read_uint64_raises_on_a_short_read():
    with pytest.raises(ValueError, match=r"expected 8 bytes.*found 3"):
        read_uint64(io.BytesIO(b"\x01\x00\x00"))


def test_read_uint64_raises_on_a_completely_empty_read():
    with pytest.raises(ValueError, match=r"expected 8 bytes.*found 0"):
        read_uint64(io.BytesIO(b""))


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


def test_a_truncated_index_raises_instead_of_silently_returning_a_low_offset(tmp_path):
    """SPECIFICATION (issue-truncated-bai-fail-open.md, decision 1). Upgraded
    from the characterisation test `..._today` that pinned the fail-open.

    `read_uint32`/`read_uint64` used to do `int.from_bytes(f.read(n), ...)`
    with no length check. At EOF, `f.read(n)` returns `b""` and
    `int.from_bytes(b"", ...)` is 0 rather than raising. Truncating a BAI to
    just its 4-byte magic made `n_ref = read_uint32(bai)` silently read 0
    instead of the real reference count, so `get_last_chunk_end` believed the
    index described zero references and returned 0 -- even though the
    untruncated index (built from the same `_bai([[(0, 900)]])` call) records
    a chunk ending at offset 900.

    The module's own docstring for `extract_unmapped_reads_from_offset`
    already promised `IOError`/`ValueError` for "if the BAI file is
    invalid" -- a truncated index is invalid, so it must raise, matching that
    contract instead of contradicting it.
    """
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:4])
    with pytest.raises(ValueError, match=r"[Tt]runcated BAI.*n_ref"):
        get_last_chunk_end(str(path))


# Byte offset at which each field starts in `_bai([[(0, 900)]])` (a single
# reference, single bin, single chunk, 2-entry linear index): magic(4) +
# n_ref(4) + n_bins(4) + bin_number(4) + n_chunks(4) + chunk_beg(8) +
# chunk_end(8) + n_intv(4) + linear_index(16) = 56 bytes total. Truncating to
# exactly the start offset of a field means that field's read finds 0 bytes;
# truncating one byte further in means it finds 1.
_FIELD_START_OFFSETS = [
    (4, "n_ref"),
    (8, "n_bins"),
    (12, "bin number"),
    (16, "n_chunks"),
    (20, "chunk_beg"),
    (28, "chunk_end"),
    (36, "n_intv"),
]


@pytest.mark.parametrize("truncate_at,field_name", _FIELD_START_OFFSETS)
def test_a_truncation_at_any_fixed_width_field_raises_and_names_that_field(tmp_path, truncate_at, field_name):
    """Every fixed-width read `get_last_chunk_end` performs -- not only the
    `n_ref` one the first truncation test happens to exercise -- must raise
    on a short read, and the message must name the field that came up short
    so a truncated BAI is diagnosable from the log alone."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:truncate_at])
    with pytest.raises(ValueError, match=rf"[Tt]runcated BAI.*{field_name}.*found 0"):
        get_last_chunk_end(str(path))


def test_a_truncation_mid_field_reports_the_partial_byte_count(tmp_path):
    """Not every short read lands on a field boundary: a truncation one byte
    into `n_ref` must still raise, and report 1 byte found, not 0."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:5])
    with pytest.raises(ValueError, match=r"[Tt]runcated BAI.*n_ref.*found 1"):
        get_last_chunk_end(str(path))


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
