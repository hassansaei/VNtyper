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

**Every statistic is over the region, not over the covered positions (#171).**
``samtools depth`` is now called with ``-a``, so it emits a row for every position
in the region including the zero-depth ones, and uncovered positions that a legacy
file omits are restored as zeros before anything is computed. Ten positions covered
at 20x inside a 100 bp region is therefore ``mean = 2.0`` and ``min = 0``, not
``mean = 20.0``. ``uncovered_bases`` counts zero-depth positions rather than
deriving them by subtraction, which is what keeps it correct under ``-a``.

The region length comes from parsing the region string, not from the depth file -
which is why a malformed region silently produces ``0`` and is pinned below.
"""

import logging

import pytest

from vntyper.scripts.coverage_stats import (
    COVERAGE_COLUMNS,
    COVERAGE_METRIC_COLUMNS,
    MAX_REGION_SPAN_BASES,
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
        "coverage_qc",
    )


def test_the_tsv_schema_is_the_metrics_plus_the_verdict():
    """The split is the point: `summarise_coverage` measures, the caller judges."""
    assert COVERAGE_COLUMNS == COVERAGE_METRIC_COLUMNS + ("coverage_qc",)


def test_the_tsv_header_and_row_are_the_exact_region_wide_bytes():
    """
    The exact bytes the TSV carries for three covered positions in a 1501 bp region.

    The header is unchanged from the pre-extraction output. The values are not:
    they were ``20.00 20.00 10.00 10 30 1501 1498 99.80`` until #171, which is the
    depth over the three positions that had reads. Recomputed by hand from the
    region-wide definition, over ``[10, 20, 30] + [0] * 1498``:

    * ``mean``   = 60 / 1501            = 0.0399733... -> ``0.04``
    * ``median`` = element 750 of 1501 sorted values, which is one of the 1498
      zeros                             = 0            -> ``0.00``
    * ``stdev``  = sqrt((1400 - 1501 * mean**2) / 1500) = 0.965263... -> ``0.97``
    * ``min``    = 0, because 1498 positions of the region carry no reads
    * ``percent_uncovered`` = 1498 / 1501 * 100 = 99.80013... -> ``99.80``,
      numerically identical to before: the metric that was already right stays right.

    ``coverage_qc`` is the caller's verdict, written verbatim - it is not in
    ``_TWO_DECIMAL_COLUMNS``. This sample is below the shipped mean threshold of 100
    *and* above the uncovered threshold of 50, so ``FAIL`` is what
    ``calculate_vntr_coverage`` would merge in (#172).
    """
    stats = summarise_coverage([10, 20, 30], total_region_length=1501)
    stats["coverage_qc"] = "FAIL"

    assert format_coverage_summary(stats) == (
        "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n"
        "0.04\t0.00\t0.97\t0\t30\t1501\t1498\t99.80\tFAIL\n"
    )


def test_formatting_a_summary_with_no_verdict_raises_rather_than_writing_eight_columns():
    """The contract that keeps a caller from writing a summary with no verdict in it.

    ``summarise_coverage`` deliberately does not produce ``coverage_qc``; if a caller
    forgets to merge it, the TSV must not silently lose its ninth column (#172).
    """
    stats = summarise_coverage([10, 20, 30], total_region_length=1501)

    with pytest.raises(KeyError, match="coverage_qc"):
        format_coverage_summary(stats)


def test_the_returned_dictionary_keeps_the_pre_extraction_types():
    """
    The dict feeds ``summary.py`` as well as the TSV, so the value types matter.

    ``statistics.median`` returns an ``int`` for an odd-length integer sample and
    a ``float`` for an even-length one; ``stdev`` is the integer ``0`` for a
    single observation. Coercing any of these to float would be a schema change
    in everything but name.

    The values moved with #171 - the base set is the 1501 bp region, not the three
    covered positions - but the types did not. ``median`` is still an ``int``
    because 1501 is odd; it is now the int ``0`` rather than the int ``20``.
    """
    stats = summarise_coverage([10, 20, 30], total_region_length=1501)

    assert stats == {
        "mean": pytest.approx(60 / 1501),
        "median": 0,
        "stdev": pytest.approx(0.9652639013706887),
        "min": 0,
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


def test_every_measured_column_is_present_in_the_returned_dictionary():
    """The dict keys and the measured columns are the same set, so they cannot drift.

    Compared against :data:`COVERAGE_METRIC_COLUMNS` rather than the full TSV schema
    since #172: ``coverage_qc`` is a verdict, and ``summarise_coverage`` holds no
    thresholds, so it cannot produce one. ``calculate_vntr_coverage`` fills it in.
    """
    stats = summarise_coverage([5], total_region_length=10)

    assert set(stats) == set(COVERAGE_METRIC_COLUMNS)


# ---------------------------------------------------------------------------
# The arithmetic
# ---------------------------------------------------------------------------


def test_full_coverage_reports_no_uncovered_bases():
    """
    Every position in the region carries reads, so nothing is uncovered.

    The one shape where the region-wide definition and the covered-position one
    agree: no zeros are restored, so #171 left these three values untouched.
    """
    stats = summarise_coverage([10, 10, 10, 10], total_region_length=4)

    assert stats["uncovered_bases"] == 0
    assert stats["percent_uncovered"] == 0.0
    assert stats["mean"] == 10.0


def test_the_mean_is_over_the_region_not_over_covered_positions():
    """#171: ten positions covered at 20x inside a 100 bp region is a mean of 2.0.

    The old definition divided by the covered-base count and reported 20.0, which is
    the depth *where there were reads*, not the depth of the region.
    """
    stats = summarise_coverage([20] * 10, total_region_length=100)

    assert stats["mean"] == 2.0
    assert stats["uncovered_bases"] == 90
    assert stats["percent_uncovered"] == 90.0


def test_uncovered_bases_counts_zero_rows_when_samtools_emitted_every_position():
    """#171: with `-a` the file already carries the zeros, so subtraction reports 0.

    This is the regression that adding `-a` alone would have introduced: under the old
    formula `region_length - len(values)` is `4 - 4 = 0` here, and the one metric that
    was correct becomes the one that is always wrong.
    """
    stats = summarise_coverage([10, 0, 0, 10], total_region_length=4)

    assert stats["uncovered_bases"] == 2
    assert stats["percent_uncovered"] == 50.0
    assert stats["mean"] == 5.0
    assert stats["min"] == 0, "the true minimum of a partly uncovered region is zero"


def test_a_legacy_depth_file_without_a_is_padded_to_the_region():
    """A file with no zero rows still yields the region-wide statistics."""
    stats = summarise_coverage([10, 20, 30], total_region_length=6)

    assert stats["uncovered_bases"] == 3
    assert stats["mean"] == 10.0
    assert stats["min"] == 0
    assert stats["max"] == 30


@pytest.mark.parametrize(
    ("values", "length"),
    [([20] * 10, 100), ([10, 20, 30], 1501), ([5] * 999, 1000), ([7], 2)],
)
def test_the_closed_form_identity_reconciles_old_and_new_means(values, length):
    """#171's free regression test: `mean_new == mean_old * (1 - pct_old / 100)`.

    With `S = sum(v)`, `c = len(v)`, `T = length`: `(S/c) * (c/T) == S/T`. It lets a
    historical report be reconciled with a re-run without new fixtures. Guarded on
    `T > 0`, because an unparseable region forces `percent_uncovered = 0`.
    """
    assert length > 0
    old_mean = sum(values) / len(values)
    old_percent = (length - len(values)) / length * 100

    stats = summarise_coverage(values, total_region_length=length)

    assert stats["mean"] == pytest.approx(old_mean * (1 - old_percent / 100))


def test_a_region_shorter_than_the_depth_file_is_refused():
    """#224. Was ``..._does_not_invent_negative_padding``, which asserted the four observed
    values simply stand and used only non-zero depths -- so it never exercised
    ``percent_uncovered`` and the impossible value it produces was masked entirely.

    Under ``samtools depth -a`` the row count is the region span, so a longer file means the
    depth table and the region string genuinely disagree and there is no correct statistic to
    return. Padding is for a SHORT file; a long one is refused.
    """
    with pytest.raises(ValueError, match="4 rows for a region of 2 bases"):
        summarise_coverage([10, 20, 30, 40], total_region_length=2)


def test_an_overlong_depth_list_is_refused_rather_than_reporting_over_100_percent():
    """#224: the value that made this worth fixing. Four uncovered positions in a two-base
    region reported ``percent_uncovered: 200.0`` -- a number no consumer should ever see."""
    with pytest.raises(ValueError, match="4 rows for a region of 2 bases"):
        summarise_coverage([0, 0, 0, 0], total_region_length=2)


def test_a_reversed_region_is_refused_rather_than_reporting_a_negative_length():
    """#224: ``chr1:2000-1000`` returned -999, and -999 reached the report under a PASS
    verdict, because ``percent_uncovered`` short-circuits to 0 whenever the length is <= 0
    and 0 is below every configured threshold."""
    with pytest.raises(ValueError, match="end 1000 is before start 2000"):
        parse_region_length("chr1:2000-1000")


def test_a_region_starting_below_one_is_refused():
    """Coordinates are 1-based; a 0 or negative start is not a region."""
    with pytest.raises(ValueError, match="start 0 is below 1"):
        parse_region_length("chr1:0-500")


def test_a_region_wider_than_the_bound_is_refused_naming_the_configured_value():
    """#224: the bound exists because disk, not memory, is the first failure mode.

    ``samtools depth -a -r chr1:1-250000000`` emitted roughly 48 million rows and exhausted
    30 GB before Python allocated anything. The message must name the configured region so an
    operator can find it in their config file.
    """
    with pytest.raises(ValueError, match=r"chr1:1-250000000"):
        parse_region_length("chr1:1-250000000")


def test_the_span_bound_is_inclusive():
    """Exactly at the bound is accepted; one base wider is not."""
    assert parse_region_length(f"chr1:1-{MAX_REGION_SPAN_BASES}") == MAX_REGION_SPAN_BASES

    with pytest.raises(ValueError, match="above the"):
        parse_region_length(f"chr1:1-{MAX_REGION_SPAN_BASES + 1}")


def test_an_invalid_region_raises_rather_than_being_swallowed_into_zero(caplog):
    """The validation must sit OUTSIDE ``parse_region_length``'s parsing ``try``.

    That block catches ``(ValueError, IndexError)`` and returns 0, so a validation raise
    placed inside it would be degraded to a silent zero-length region -- which is the exact
    failure mode being fixed rather than a fallback. If this test ever reports 0 instead of
    raising, the guard has been moved back inside the ``try``.
    """
    with caplog.at_level(logging.ERROR, logger=COVERAGE_STATS_LOGGER), pytest.raises(ValueError):
        parse_region_length("chr1:5000-1000")

    errors = [record for record in caplog.records if record.levelno == logging.ERROR]
    assert errors, "an invalid region must be logged at ERROR, not silently degraded"


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
    assert stats["uncovered_bases"] == 0, (
        "#171: a negative uncovered count was nonsense. `absent` is floored at 0, so an "
        "unparseable region reports 0 uncovered rather than -2. The WARNING from "
        "parse_region_length is still the only signal that the region was unreadable."
    )


def test_the_median_of_an_even_sample_is_the_midpoint():
    """A fully covered four-base region: no zeros are restored, so #171 moved nothing here."""
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

    Twelve covered positions out of a 20 bp region, depths 1..12. Recomputed by hand
    from #171's region-wide definition, over ``[1..12] + [0] * 8``:

    * ``mean``   = 78 / 20                                = 3.9      -> ``3.90``
    * ``median`` = mean of sorted elements 9 and 10, which are 2 and 3
      (the first eight are the restored zeros)             = 2.5      -> ``2.50``
    * ``stdev``  = sqrt((650 - 20 * 3.9**2) / 19) = sqrt(18.2) = 4.266145... -> ``4.27``
    * ``min``    = 0; ``max`` = 12
    * ``percent_uncovered`` = 8 / 20 * 100 = 40.0, unchanged.

    Before #171 this read ``6.50 6.50 3.61 1 12 20 8 40.00`` - the mean, median and
    stdev of the twelve positions that had reads. The ninth field is #172's verdict:
    ``FAIL`` on the mean alone here, since 40% uncovered is inside the shipped 50%.
    """
    path = _write_depth(tmp_path, [("chr1", 100 + i, i + 1) for i in range(12)])

    values = read_depth_values(path)
    length = parse_region_length("chr1:100-119")
    stats = summarise_coverage(values, total_region_length=length)
    stats["coverage_qc"] = "FAIL"

    assert length == 20
    assert stats["mean"] == pytest.approx(3.9)
    assert stats["median"] == 2.5
    assert stats["min"] == 0
    assert stats["max"] == 12
    assert stats["uncovered_bases"] == 8
    assert stats["percent_uncovered"] == pytest.approx(40.0)
    assert format_coverage_summary(stats) == (
        "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n"
        "3.90\t2.50\t4.27\t0\t12\t20\t8\t40.00\tFAIL\n"
    )
