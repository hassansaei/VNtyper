"""Unit tests for the VCF usability contract on Kestrel's output (#223).

A Kestrel run that exits 0 but writes a VCF nothing can parse used to become a confident
``Negative``. The severe case is not the zero-byte file named in the issue title -- that at
least corresponds to a genuinely empty result. It is a VCF that **carries real indel records
and has lost its ``#CHROM`` header**: ``filter_vcf`` copies data lines through whether or not
a header is present, so the derived files inherit the missing header,
``read_vcf_without_comments`` returns two empty frames, and ``process_kestrel_output`` renders
two empty frames as the ``Negative`` placeholder. Exit 0, a normal-looking report, and **zero
ERROR records**.

The contract is **ordered**, not merely present. ``read_vcf_without_comments`` sets its
columns from the ``#CHROM`` line and collects a data line only once that line has been seen,
so records arriving before the header are discarded in silence -- a file whose header is late
passes a presence check and still parses to nothing.
``test_variant_parsing.py::test_data_line_before_the_header_is_dropped_not_leaked_in`` pins
that dropping as correct behaviour of the parser, which is exactly why the *file* has to be
rejected instead.

A header with no records is **usable**, and that distinction is the point of the whole change:
"Kestrel ran and found nothing" is a real negative and must keep working, while "Kestrel
produced nothing readable" must not render identically to it.
"""

import pytest

from vntyper.scripts.file_processing import describe_unusable_vcf

pytestmark = pytest.mark.unit

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
RECORD = "chr1\t155160000\t.\tC\tCG\t.\t.\tDP=100\n"


def test_a_well_formed_vcf_is_usable(tmp_path):
    """Meta lines, a header, then records -- the ordinary case."""
    path = tmp_path / "output.vcf"
    path.write_text(META + HEADER + RECORD, encoding="utf-8")

    assert describe_unusable_vcf(path) is None


def test_a_header_with_no_records_is_usable(tmp_path):
    """A genuinely empty Kestrel result.

    This is the legitimate ``Negative`` and it must keep working. A check that rejected an
    empty result would convert every true negative into an aborted run, which is a strictly
    worse failure than the one being fixed.
    """
    path = tmp_path / "output.vcf"
    path.write_text(META + HEADER, encoding="utf-8")

    assert describe_unusable_vcf(path) is None


def test_a_header_with_no_meta_lines_is_usable(tmp_path):
    """The contract is about the header, not about what precedes it."""
    path = tmp_path / "output.vcf"
    path.write_text(HEADER + RECORD, encoding="utf-8")

    assert describe_unusable_vcf(path) is None


def test_a_zero_byte_vcf_names_the_missing_header(tmp_path):
    """The mild case from the issue title: nothing was written at all."""
    path = tmp_path / "output.vcf"
    path.write_text("", encoding="utf-8")

    assert describe_unusable_vcf(path) == "it has no #CHROM header line"


def test_a_meta_only_vcf_names_the_missing_header(tmp_path):
    """Kestrel wrote its preamble and then stopped."""
    path = tmp_path / "output.vcf"
    path.write_text(META + "##source=Kestrel\n", encoding="utf-8")

    assert describe_unusable_vcf(path) == "it has no #CHROM header line"


def test_a_headerless_vcf_carrying_records_says_how_many_are_discarded(tmp_path):
    """The severe case: real indel records, silently dropped, reported as a negative.

    The count is in the message because it is the measure of what was lost. "No header" reads
    like an empty file; "no header, so all 2 data lines are discarded" reads like the data
    loss it is.
    """
    path = tmp_path / "output.vcf"
    path.write_text(META + RECORD + RECORD, encoding="utf-8")

    reason = describe_unusable_vcf(path)

    assert reason is not None
    assert "no #CHROM header line" in reason
    assert "2 data line(s)" in reason


def test_records_before_a_late_header_are_reported(tmp_path):
    """An ordered contract, not a presence one.

    This file *has* a valid ``#CHROM`` line, so a presence check accepts it -- while the
    record before that line is still dropped by the parser and can leave both derived frames
    empty. If this test fails, the check has been weakened to "a #CHROM line exists somewhere".
    """
    path = tmp_path / "output.vcf"
    path.write_text(META + RECORD + HEADER + RECORD, encoding="utf-8")

    reason = describe_unusable_vcf(path)

    assert reason is not None
    assert "1 data line(s) appear before the #CHROM header" in reason


def test_blank_lines_before_the_header_are_not_records(tmp_path):
    """Whitespace is not data. Counting it would reject well-formed files."""
    path = tmp_path / "output.vcf"
    path.write_text(META + "\n\n" + HEADER + RECORD, encoding="utf-8")

    assert describe_unusable_vcf(path) is None


def test_an_absent_path_is_reported_rather_than_raising(tmp_path):
    """The caller decides what to do about it, so this reports rather than raises.

    ``read_vcf_without_comments`` turns a missing file into an empty frame, which two steps
    later is the ``Negative`` placeholder -- so "cannot be read" has to be distinguishable
    from "read fine, found nothing".
    """
    reason = describe_unusable_vcf(tmp_path / "absent.vcf")

    assert reason is not None
    assert "could not be read" in reason


def test_a_directory_is_reported_rather_than_raising(tmp_path):
    """``open()`` on a directory raises ``IsADirectoryError``, an ``OSError`` subclass."""
    reason = describe_unusable_vcf(tmp_path)

    assert reason is not None
    assert "could not be read" in reason


def test_a_truncated_chrom_line_is_not_a_header(tmp_path):
    """``read_vcf_without_comments`` accepts *any* line beginning ``#CHROM`` as its column
    definition, without checking the columns. A bare ``#CHROM`` therefore sets a one-element
    column list, collects no records, and produces an empty frame -- reaching the ``Negative``
    placeholder through a file that a presence-and-order check would have called usable.
    """
    path = tmp_path / "output.vcf"
    path.write_text(META + "#CHROM\n" + RECORD, encoding="utf-8")

    reason = describe_unusable_vcf(path)

    assert reason is not None
    assert "truncated" in reason
    assert "1 column(s)" in reason


def test_a_header_missing_its_last_mandatory_column_is_refused(tmp_path):
    """Seven of the eight mandatory columns is still truncated."""
    path = tmp_path / "output.vcf"
    path.write_text(META + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n", encoding="utf-8")

    reason = describe_unusable_vcf(path)

    assert reason is not None
    assert "truncated" in reason


def test_a_header_with_format_and_sample_columns_is_usable(tmp_path):
    """The eight mandatory columns are a prefix, not the whole line. Kestrel writes a sample
    column, and rejecting it would refuse every real VCF the pipeline produces."""
    path = tmp_path / "output.vcf"
    path.write_text(META + HEADER.rstrip("\n") + "\tFORMAT\tsample\n" + RECORD, encoding="utf-8")

    assert describe_unusable_vcf(path) is None


def test_a_header_whose_columns_are_misspelled_is_refused(tmp_path):
    """Order and spelling both matter: the parser names DataFrame columns from this line, and
    every downstream consumer looks them up by exact name."""
    path = tmp_path / "output.vcf"
    path.write_text(META + "#CHROM\tPOSITION\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")

    assert describe_unusable_vcf(path) is not None
