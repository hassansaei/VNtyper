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
    BAI_MAGIC,
    extract_unmapped_reads_from_offset,
    get_last_chunk_end,
    read_exact,
    read_uint32,
    read_uint64,
)

pytestmark = pytest.mark.unit


def _bai(references, magic=b"BAI\x01"):
    """Build a minimal BAI. `references` is a list of lists of (beg, end) chunks.

    Matches, byte for byte, the order `get_last_chunk_end` walks: a 4-byte
    magic, which the parser reads and compares against `BAI\\x01`, then n_ref
    (uint32), then per reference: n_bins (uint32); for each bin, a bin number
    (uint32, read but unused), a chunk count (uint32) and that many
    (chunk_beg, chunk_end) uint64 pairs; then n_intv (uint32) and n_intv * 8
    bytes of linear index, which the parser steps over with `bai.seek` after
    checking those bytes are actually present in the file.

    All chunks for one reference are written into a single synthetic bin
    here (n_bins is 0 or 1). That is a simplification versus a real BAI,
    which spreads chunks across many bins, but it is invisible to
    `get_last_chunk_end`: it maxes `chunk_end` over every chunk in every
    bin, so how those chunks are grouped into bins does not affect the
    result.

    Args:
        references: One list of `(chunk_beg, chunk_end)` pairs per reference.
        magic: The four signature bytes to write. Overridden only by the
            tests that check the signature is validated rather than skipped.

    Returns:
        bytes: The encoded index.
    """
    out = [magic, struct.pack("<I", len(references))]
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
# read_exact
# ---------------------------------------------------------------------------


def test_the_magic_constant_is_the_bai_signature():
    """Pinned as a literal so a typo in the constant cannot pass unnoticed.

    Confirmed against a real `samtools index` output, whose first four bytes
    are these.
    """
    assert BAI_MAGIC == b"BAI\x01"


@pytest.mark.parametrize("size", [1, 4, 8, 16])
def test_read_exact_returns_exactly_the_bytes_asked_for(size):
    """The bytes come back unchanged, and only `size` of them are consumed."""
    stream = io.BytesIO(bytes(range(32)))
    assert read_exact(stream, size, "field") == bytes(range(size))
    assert stream.tell() == size


def test_read_exact_raises_naming_the_field_and_both_byte_counts():
    """The one behaviour every caller depends on: a short read is not silently short.

    `f.read(n)` returns what it has rather than raising, so without this check
    a truncated field reaches the caller as a shorter `bytes` -- which
    `int.from_bytes` then zero-fills into a plausible number.
    """
    with pytest.raises(ValueError, match=r"expected 8 bytes for magic, found 3"):
        read_exact(io.BytesIO(b"BAI"), 8, "magic")


def test_read_exact_raises_on_an_empty_stream():
    with pytest.raises(ValueError, match=r"expected 4 bytes for magic, found 0"):
        read_exact(io.BytesIO(b""), 4, "magic")


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
# get_last_chunk_end: the two structural checks the docstring promised but the
# parser did not perform -- the signature, and the length of the linear index.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "magic", [b"CSI\x01", b"BAI\x02", b"\x00\x00\x00\x00", b"CRAI"], ids=["csi", "v2", "zeros", "crai"]
)
def test_an_index_whose_signature_is_not_the_bai_magic_is_rejected(tmp_path, magic):
    """The four magic bytes must be read and compared, not skipped.

    `get_last_chunk_end` consumed the signature with a bare `bai.read(4)` whose
    result was thrown away, so the only thing that decided whether a file was
    accepted was whether the bytes *after* the signature happened to parse. The
    concrete case this makes safe: a file with a wrong signature and a
    parseable `n_ref = 0` returned a plausible offset of `0` -- indistinguishable
    from the value a legitimately empty index returns -- instead of raising.

    The `.csi` case is the realistic one. A CSI index is a different format
    carrying the same information, it lives beside a BAM under a different
    extension, and its first four bytes are `CSI\\x01`; handed one, the parser
    would previously have walked its bgzf-compressed body as if it were a BAI.

    Args:
        tmp_path: Pytest temporary directory.
        magic: Four signature bytes that are not `BAI\\x01`.
    """
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([], magic=magic))
    with pytest.raises(ValueError, match=r"[Ii]nvalid BAI index.*magic"):
        get_last_chunk_end(str(path))


def test_the_valid_signature_is_still_accepted(tmp_path):
    """The guard above must reject the wrong magic without rejecting the right one."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]]))
    assert get_last_chunk_end(str(path)) == 900


def test_a_file_too_short_to_hold_a_signature_names_the_magic_not_n_ref(tmp_path):
    """A file shorter than four bytes runs out inside the signature itself.

    Reported against the field that actually came up short. Skipping the magic
    with an unchecked `read(4)` pushed the first detected shortfall onto
    `n_ref`, which sends whoever reads the log looking at the wrong offset.
    """
    path = tmp_path / "x.bai"
    path.write_bytes(b"BA")
    with pytest.raises(ValueError, match=r"[Tt]runcated BAI index.*magic.*found 2"):
        get_last_chunk_end(str(path))


def test_a_missing_linear_index_tail_raises_instead_of_seeking_past_the_end(tmp_path):
    """The declared linear-index length must be checked against the real file length.

    The `n_intv * 8` linear-index body is stepped over with `bai.seek(...,
    SEEK_CUR)`, and seeking past the end of a file succeeds on every platform --
    so the seek returning normally proves nothing about whether those bytes
    exist. For every reference but the last, the next `read_uint32` lands at EOF
    and the truncation is caught anyway; for the *final* reference nothing
    follows, so the walk simply ended and `max_vo` was returned from an index
    that stops mid-structure.

    Truncating `_bai([[(0, 900)]])` to 40 bytes leaves every fixed-width field
    intact -- including `n_intv = 2` -- and removes all 16 bytes of the linear
    index it declares.
    """
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:40])
    with pytest.raises(ValueError, match=r"[Tt]runcated BAI index.*linear index.*found 0"):
        get_last_chunk_end(str(path))


def test_a_partially_present_linear_index_reports_how_many_bytes_were_there(tmp_path):
    """Half a linear index is as invalid as none, and the message says how much was found."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:48])
    with pytest.raises(ValueError, match=r"[Tt]runcated BAI index.*linear index.*found 8"):
        get_last_chunk_end(str(path))


def test_the_optional_trailing_unplaced_read_count_is_tolerated(tmp_path):
    """The length check is a lower bound, not an equality: real BAIs carry a tail.

    `samtools index` writes an optional trailing `n_no_coor` (uint64, the count
    of reads with no coordinate) after the last reference. `get_last_chunk_end`
    never reads it, and the linear-index length check must not turn its presence
    into a validation failure -- an exact-length check would reject every index
    samtools produces.
    """
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]]) + struct.pack("<Q", 17))
    assert get_last_chunk_end(str(path)) == 900


def test_a_large_declared_linear_index_is_rejected_without_allocating_it(tmp_path):
    """A corrupt `n_intv` must fail fast rather than trying to buy the memory.

    Consuming the linear index with `f.read(n_intv * 8)` would honour a
    corrupt count literally: `n_intv = 2**32 - 1` asks for 34 GB. The check
    compares offsets against the file's real length instead, so it costs
    nothing and still raises.
    """
    path = tmp_path / "x.bai"
    body = _bai([[(0, 900)]])
    path.write_bytes(body[:36] + struct.pack("<I", 0xFFFFFFFF))
    with pytest.raises(ValueError, match=r"[Tt]runcated BAI index.*linear index"):
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
