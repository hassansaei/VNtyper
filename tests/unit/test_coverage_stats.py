"""
Unit tests for :mod:`vntyper.scripts.coverage_stats`.

This module holds the pure part of VNTR coverage calculation, extracted from
``fastq_bam_processing.calculate_vntr_coverage`` so it can be tested without
samtools, a BAM file or a subprocess.

**Contract C1 - the TSV schema is frozen.** ``generate_report.py`` reads
``region_length``, ``uncovered_bases`` and ``percent_uncovered`` out of this file
by exact key, with ``0`` defaults, and ``tests/helpers.validate_coverage_output``
reads all eight. A renamed column does not fail - it silently reports zero
coverage in the HTML report. The header row asserted here is byte-identical to
the one the pre-extraction code wrote.

The interesting arithmetic is ``percent_uncovered``. ``samtools depth`` emits a
row **only for positions with at least one read**, so the number of rows is the
covered-base count and everything else in the region is uncovered. The region
length comes from parsing the region string, not from the depth file - which is
why a malformed region silently produces ``0`` and is pinned below.
"""

import logging

import pytest

from vntyper.scripts.coverage_stats import (
    COVERAGE_COLUMNS,
    format_coverage_summary,
    parse_region_length,
    read_depth_values,
    summarise_coverage,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit

COVERAGE_STATS_LOGGER = "vntyper.scripts.coverage_stats"


def _write_depth(tmp_path, rows, name="depth.txt"):
    """Write a synthetic ``samtools depth`` file: contig, position, depth per line."""
    path = tmp_path / name
    path.write_text("".join(f"{contig}\t{pos}\t{depth}\n" for contig, pos, depth in rows))
    return path


# ---------------------------------------------------------------------------
# C1 - the frozen TSV schema
# ---------------------------------------------------------------------------


def test_the_coverage_columns_are_the_frozen_lowercase_schema():
    """
    Contract C1, stated as a literal.

    These names are lowercase and `generate_report.py` looks them up with
    ``.get(key, 0)``. Renaming one does not raise anywhere - it makes the report
    display 0 bp uncovered on a sample with no coverage at all.
    """
    assert COVERAGE_COLUMNS == (
        "mean",
        "median",
        "stdev",
        "min",
        "max",
        "region_length",
        "uncovered_bases",
        "percent_uncovered",
    )


def test_the_tsv_header_and_row_are_byte_identical_to_the_pre_extraction_output():
    """
    The exact bytes ``calculate_vntr_coverage`` wrote before the extraction.

    Captured by running the pre-refactor function over a three-position depth
    file spanning ``chr1:155160500-155162000``.
    """
    stats = summarise_coverage([10, 20, 30], total_region_length=1501)

    assert format_coverage_summary(stats) == (
        "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\n"
        "20.00\t20.00\t10.00\t10\t30\t1501\t1498\t99.80\n"
    )


def test_the_returned_dictionary_keeps_the_pre_extraction_types():
    """
    The dict feeds ``summary.py`` as well as the TSV, so the value types matter.

    ``statistics.median`` returns an ``int`` for an odd-length integer sample and
    a ``float`` for an even-length one; ``stdev`` is the integer ``0`` for a
    single observation. Coercing any of these to float would be a schema change
    in everything but name.
    """
    stats = summarise_coverage([10, 20, 30], total_region_length=1501)

    assert stats == {
        "mean": 20.0,
        "median": 20,
        "stdev": 10.0,
        "min": 10,
        "max": 30,
        "region_length": 1501,
        "uncovered_bases": 1498,
        "percent_uncovered": pytest.approx(99.80013324450367),
    }
    assert isinstance(stats["median"], int)
    assert isinstance(stats["min"], int)
    assert isinstance(stats["max"], int)
    assert isinstance(stats["region_length"], int)
    assert isinstance(stats["uncovered_bases"], int)


def test_every_frozen_column_is_present_in_the_returned_dictionary():
    """The dict keys and the TSV columns are the same set, so they cannot drift."""
    stats = summarise_coverage([5], total_region_length=10)

    assert set(stats) == set(COVERAGE_COLUMNS)


# ---------------------------------------------------------------------------
# The arithmetic
# ---------------------------------------------------------------------------


def test_full_coverage_reports_no_uncovered_bases():
    """Every position in the region has a depth row, so nothing is uncovered."""
    stats = summarise_coverage([10, 10, 10, 10], total_region_length=4)

    assert stats["uncovered_bases"] == 0
    assert stats["percent_uncovered"] == 0.0
    assert stats["mean"] == 10.0


def test_partial_coverage_counts_the_missing_positions_as_uncovered():
    """
    The core `samtools depth` semantics: absent rows are zero-coverage positions.

    Ten covered positions inside a 100 bp region is 90% uncovered, and the mean is
    the mean of the **covered** positions only - not of the region. That is the
    pre-existing definition and it is preserved deliberately: changing it would
    move every coverage number in every historical report.
    """
    stats = summarise_coverage([20] * 10, total_region_length=100)

    assert stats["uncovered_bases"] == 90
    assert stats["percent_uncovered"] == 90.0
    assert stats["mean"] == 20.0, "the mean is over covered positions, not over the region"


def test_a_single_covered_position_reports_a_zero_standard_deviation():
    """``statistics.stdev`` needs two observations; one must not raise."""
    stats = summarise_coverage([42], total_region_length=1)

    assert stats["stdev"] == 0
    assert stats["mean"] == 42.0
    assert stats["median"] == 42


def test_zero_coverage_everywhere_is_a_hard_error_not_a_zero_mean():
    """
    An empty depth file means the region matched nothing - the silent-false-negative case.

    A wrong contig name, a region off the end of the chromosome or a BAM with no
    reads all land here. Returning ``mean=0`` would let the pipeline carry on and
    report a confident negative genotype, so this raises instead.
    """
    with pytest.raises(RuntimeError, match="No coverage data found"):
        summarise_coverage([], total_region_length=1500)


def test_a_zero_length_region_reports_zero_percent_rather_than_dividing_by_zero():
    """A region string that could not be parsed yields length 0; that must not crash."""
    stats = summarise_coverage([10, 20], total_region_length=0)

    assert stats["percent_uncovered"] == 0
    assert stats["region_length"] == 0
    assert stats["uncovered_bases"] == -2, "preserved: covered bases still subtract from an unknown length"


def test_the_median_of_an_even_sample_is_the_midpoint():
    stats = summarise_coverage([10, 20, 30, 40], total_region_length=4)

    assert stats["median"] == 25.0
    assert stats["mean"] == 25.0


# ---------------------------------------------------------------------------
# Region parsing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("region", "expected"),
    [
        ("chr1:155160500-155162000", 1501),
        ("1:155160500-155162000", 1501),
        ("NC_000001.10:155160500-155162000", 1501),
        ("chr1:100-100", 1),
    ],
)
def test_the_region_length_is_inclusive_of_both_endpoints(region, expected):
    """``end - start + 1`` - the same arithmetic as before the extraction."""
    assert parse_region_length(region) == expected


@pytest.mark.parametrize(
    "region",
    ["chr1", "chr1:155160500", "chr1:a-b", "", "chr1:1-2-3", "chr1:155160500:155162000"],
)
def test_an_unparseable_region_yields_zero_and_warns(region, caplog):
    """
    A malformed region degrades to length 0 rather than raising.

    That is the pre-existing behaviour and it is preserved, but it is worth
    knowing what it costs: ``percent_uncovered`` then reports ``0`` for a sample
    whose region could not be read at all. The warning is the only signal, so it
    must be at WARNING - ``pytest.ini`` sets ``log_level = DEBUG``, which is why
    this asserts on ``levelno`` rather than on ``caplog.text``.
    """
    with caplog.at_level(logging.INFO, logger=COVERAGE_STATS_LOGGER):
        assert parse_region_length(region) == 0

    warnings = [record for record in caplog.records if record.levelno == logging.WARNING]
    assert warnings, f"a malformed region {region!r} must warn; silence here hides a zero denominator"


# ---------------------------------------------------------------------------
# Depth file reading
# ---------------------------------------------------------------------------


def test_a_synthetic_depth_file_is_read_as_its_third_column(tmp_path):
    """``samtools depth`` writes contig, 1-based position, depth."""
    path = _write_depth(tmp_path, [("chr1", 155160500, 12), ("chr1", 155160501, 0), ("chr1", 155160502, 7)])

    assert read_depth_values(path) == [12, 0, 7]


def test_blank_lines_in_a_depth_file_are_skipped(tmp_path):
    """A trailing newline must not become a parse error."""
    path = tmp_path / "depth.txt"
    path.write_text("chr1\t1\t5\n\nchr1\t2\t6\n\n")

    assert read_depth_values(path) == [5, 6]


def test_an_empty_depth_file_reads_as_no_values(tmp_path):
    """Empty is the region-matched-nothing case; ``summarise_coverage`` turns it into an error."""
    path = tmp_path / "depth.txt"
    path.write_text("")

    assert read_depth_values(path) == []


def test_end_to_end_over_a_synthetic_depth_file(tmp_path):
    """
    The whole pure path: depth file plus region string to the exact TSV bytes.

    Twelve covered positions out of a 20 bp region, depths 1..12.
    """
    path = _write_depth(tmp_path, [("chr1", 100 + i, i + 1) for i in range(12)])

    values = read_depth_values(path)
    length = parse_region_length("chr1:100-119")
    stats = summarise_coverage(values, total_region_length=length)

    assert length == 20
    assert stats["mean"] == pytest.approx(6.5)
    assert stats["median"] == 6.5
    assert stats["min"] == 1
    assert stats["max"] == 12
    assert stats["uncovered_bases"] == 8
    assert stats["percent_uncovered"] == pytest.approx(40.0)
    assert format_coverage_summary(stats) == (
        "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\n"
        "6.50\t6.50\t3.61\t1\t12\t20\t8\t40.00\n"
    )
