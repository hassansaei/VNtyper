"""Unit tests for file_processing.py -- the indel split that feeds Kestrel.

`filter_vcf` decides what counts as an indel and `filter_indel_vcf` decides
which side of the insertion/deletion split it lands on. Both were untested.
A row dropped here never reaches motif processing, scoring or the report.
"""

import pytest

from vntyper.scripts.file_processing import filter_indel_vcf, filter_vcf

pytestmark = pytest.mark.unit

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(tmp_path, name, *rows):
    path = tmp_path / name
    path.write_text(HEADER + "".join(f"X-5\t{i}\t.\t{ref}\t{alt}\t.\tPASS\t.\n" for i, (ref, alt) in enumerate(rows)))
    return path


def _rows(path):
    return [line for line in path.read_text().splitlines() if not line.startswith(("##", "#CHROM"))]


def test_the_header_is_carried_through_verbatim(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CC"))), str(out))
    assert out.read_text().startswith("##fileformat=VCFv4.2\n#CHROM")


def test_a_one_base_to_many_insertion_is_kept(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGGCA"))), str(out))
    assert len(_rows(out)) == 1


def test_a_many_to_one_base_deletion_is_kept(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("CGGCA", "C"))), str(out))
    assert len(_rows(out)) == 1


def test_a_single_base_substitution_is_dropped(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "G"))), str(out))
    assert _rows(out) == []


def test_an_indel_where_neither_allele_is_one_base_is_dropped_today(tmp_path):
    """CHARACTERISATION of a real gap, not a specification.

    file_processing.py:28-30 keeps a row only when EXACTLY ONE of REF/ALT has
    length snv_length (1). REF="AC", ALT="ACGGG" is a genuine 3-base insertion
    and is silently discarded before Kestrel post-processing ever sees it.
    Whether Kestrel can emit such a record is a domain question; pinned here so
    the behaviour cannot drift, and filed rather than fixed.
    """
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("AC", "ACGGG"))), str(out))
    assert _rows(out) == []


def test_a_line_with_too_few_columns_raises(tmp_path):
    """file_processing.py:27 unpacks six fields; a short line is not tolerated."""
    path = tmp_path / "short.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tC\n")
    with pytest.raises(ValueError):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_insertions_and_deletions_are_written_to_separate_files(tmp_path):
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    filter_indel_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGG"), ("CGG", "C"))), str(ins), str(dele))
    assert len(_rows(ins)) == 1
    assert len(_rows(dele)) == 1


def test_both_output_files_get_the_header(tmp_path):
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    filter_indel_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGG"))), str(ins), str(dele))
    assert dele.read_text().startswith("##fileformat")
    assert _rows(dele) == []


def test_equal_length_alleles_are_classified_as_a_deletion(tmp_path):
    """file_processing.py:60-63: the insertion test is strictly `alt > ref`.

    The docstring is explicit that a deletion is REF longer than *or equal
    to* ALT, so a same-length pair (which filter_vcf would normally have
    dropped as a substitution before this function ever sees it) falls to
    the `else` branch. This pins the documented `>` boundary directly, since
    none of the other rows in this file exercise ref_len == alt_len here.
    """
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    path = tmp_path / "indel.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tC\tG\t.\tPASS\t.\n")
    filter_indel_vcf(str(path), str(ins), str(dele))
    assert _rows(ins) == []
    assert len(_rows(dele)) == 1


def test_a_multi_base_ref_insertion_is_classified_as_a_deletion_today(tmp_path):
    """CHARACTERISATION of a real mis-classification, not a specification.

    file_processing.py:60-63 tests only `len(ref) == 1 and len(alt) > 1` for
    the insertion side; the `else` is an unconditional catch-all. REF="AC",
    ALT="ACGGG" is an insertion and lands in the DELETION file, where the
    deletion frame rule (3n+2) is applied to it instead of the insertion rule
    (3n+1) -- the two are not interchangeable (#181). It cannot arrive today
    because `filter_vcf` drops it first (see above), so the two gaps mask each
    other. Pinned and filed; not fixed.
    """
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    path = tmp_path / "indel.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tAC\tACGGG\t.\tPASS\t.\n")
    filter_indel_vcf(str(path), str(ins), str(dele))
    assert _rows(ins) == []
    assert len(_rows(dele)) == 1
