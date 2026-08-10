# tests/unit/test_variant_parsing.py

"""
Unit tests for variant_parsing.py.

Covers both public functions:
  - filter_by_alt_values_and_finalize(): 'GG' alt requires a minimum Depth_Score;
    exclude_alts removed from final DataFrame; left/right columns retained.
  - read_vcf_without_comments(): parses a (possibly gzipped) VCF, skipping '##'
    meta lines, taking '#CHROM' as the header, and collecting only the lines that
    come after a header has been seen.

Mutation-killing coverage for read_vcf_without_comments()
-----------------------------------------------------------
Three surviving mutants (docs/development/mutation-testing.md) live in this
function's two compound boolean conditions:

    line 59: elif not line.startswith("##") and header:
    line 68: if data and header:

For line 59, let P = "line.startswith('##')" (line looks like a meta/comment
line) and H = "header is truthy" (the '#CHROM' line has already been seen).
The correct condition is `(not P) and H`. Truth table against the two mutants
(M1 = `not` deleted -> `P and H`; M2 = `and` -> `or` -> `(not P) or H`):

    P     H     correct=(not P) and H   M1=(P and H)   M2=((not P) or H)
    False False        False                False            True   <- M2 differs (test: pre-header data line)
    False True         True                 False            True   <- M1 differs
    True  False        False                False            False
    True  True         False                True             True   <- M1 and M2 both differ (test: post-header comment)

Row (False, False) is the case that isolates M2 alone: a data-shaped line
appearing *before* any '#CHROM' header must be dropped under `and`, but would
leak through under `or`. Row (True, True) isolates M1: a '##' comment line
appearing *after* the header is seen must still be skipped under the real
`not`, but would be captured as data if `not` is deleted.

For line 68, let D = "data is non-empty" and H = "header is truthy". The
correct condition is `D and H`; the mutant is `D or H`. Because line 59's real
(unmutated) `and` guarantees data can only be populated once header is
already truthy, the only reachable case that distinguishes `D and H` from
`D or H` is D=False, H=True: a file containing only the '#CHROM' header line
and no data rows. Real code returns a columnless `pd.DataFrame()`; the `or`
mutant returns `pd.DataFrame([], columns=header)`, which has the header's
columns even though it is still row-empty.
"""

import gzip
import logging

import pandas as pd
import pytest

from vntyper.scripts.variant_parsing import (
    filter_by_alt_values_and_finalize,
    read_vcf_without_comments,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


@pytest.fixture
def kestrel_config_mock():
    """
    Provide a mock kestrel_config dict specifically for ALT filtering tests.

    Real config may contain more fields, but we only need 'alt_filtering' here.
    """
    return {
        "alt_filtering": {
            "gg_alt_value": "GG",
            "gg_depth_score_threshold": 0.02,
            "exclude_alts": ["BAD_ALT", "ZZZ"],
        }
    }


def test_filter_by_alt_values_empty_df(kestrel_config_mock):
    """Test that an empty DataFrame simply returns empty without error."""
    df = pd.DataFrame()
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)
    assert out.empty, "Empty input should yield empty output."


def test_filter_by_alt_values_missing_columns(kestrel_config_mock):
    """Test that missing 'ALT' or 'Depth_Score' columns raises KeyError."""
    df = pd.DataFrame(
        {
            "ALT": ["GG", "ABC"],
            # 'Depth_Score' is missing here
        }
    )

    with pytest.raises(KeyError) as exc_info:
        filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    assert "Missing required columns" in str(exc_info.value), (
        "Expected KeyError due to missing required 'Depth_Score' column."
    )


def test_filter_by_alt_values_gg_filter_below_threshold(kestrel_config_mock):
    """
    Test that GG variants below the depth score threshold are correctly flagged.

    When ALT='GG' and Depth_Score < threshold (0.02), the row is marked
    with alt_filter_pass=False but retained in the DataFrame.
    This follows the flag-and-defer pattern where filtering is done downstream.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "ALT": ["GG", "GG", "XYZ"],
            "Depth_Score": [
                0.019,  # Below threshold (0.02) for GG => should fail
                0.02,  # At threshold for GG => should pass
                0.5,  # Non-GG ALT => unaffected by GG filter
            ],
        }
    )

    # Act
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # Assert - all rows are retained
    assert len(out) == 3, "All 3 rows should be retained in the DataFrame"

    # Assert - alt_filter_pass column exists
    assert "alt_filter_pass" in out.columns, "Should have 'alt_filter_pass' column"

    # Assert - row 0: GG with Depth_Score=0.019 should fail (below threshold)
    assert not out.loc[0, "alt_filter_pass"], (
        "Row 0 (GG with Depth_Score=0.019) should have alt_filter_pass=False (below threshold)"
    )

    # Assert - row 1: GG with Depth_Score=0.02 should pass (at threshold)
    assert out.loc[1, "alt_filter_pass"], (
        "Row 1 (GG with Depth_Score=0.02) should have alt_filter_pass=True (meets threshold)"
    )

    # Assert - row 2: XYZ should pass (not affected by GG-specific filter)
    assert out.loc[2, "alt_filter_pass"], "Row 2 (XYZ) should have alt_filter_pass=True (non-GG ALT not affected)"

    # Assert - verify the count of passing rows
    pass_count = out["alt_filter_pass"].sum()
    assert pass_count == 2, "Should have exactly 2 rows with alt_filter_pass=True"


def test_filter_by_alt_values_exclude_alts(kestrel_config_mock):
    """
    Test that ALTs in 'exclude_alts' list are correctly flagged.

    ALTs in the exclude_alts list (e.g., 'BAD_ALT', 'ZZZ') are marked
    with alt_filter_pass=False but retained in the DataFrame.
    This follows the flag-and-defer pattern where filtering is done downstream.
    """
    # Arrange
    # Config has exclude_alts = ["BAD_ALT", "ZZZ"]
    df = pd.DataFrame(
        {
            "ALT": ["GG", "BAD_ALT", "OK_ALT", "ZZZ", "ANOTHER"],
            "Depth_Score": [0.5, 0.1, 0.3, 0.2, 0.6],  # just some placeholder scores
        }
    )

    # Act
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # Assert - all rows are retained
    assert len(out) == 5, "All 5 rows should be retained in the DataFrame"

    # Assert - alt_filter_pass column exists
    assert "alt_filter_pass" in out.columns, "Should have 'alt_filter_pass' column"

    # Assert - check each row's filter status
    # Row 0: GG with good Depth_Score => should pass
    assert out.loc[0, "alt_filter_pass"], "Row 0 (GG with Depth_Score=0.5) should have alt_filter_pass=True"

    # Row 1: BAD_ALT is in exclude_alts => should fail
    assert not out.loc[1, "alt_filter_pass"], "Row 1 (BAD_ALT) should have alt_filter_pass=False (in exclude_alts)"

    # Row 2: OK_ALT is not in exclude_alts => should pass
    assert out.loc[2, "alt_filter_pass"], "Row 2 (OK_ALT) should have alt_filter_pass=True (not in exclude_alts)"

    # Row 3: ZZZ is in exclude_alts => should fail
    assert not out.loc[3, "alt_filter_pass"], "Row 3 (ZZZ) should have alt_filter_pass=False (in exclude_alts)"

    # Row 4: ANOTHER is not in exclude_alts => should pass
    assert out.loc[4, "alt_filter_pass"], "Row 4 (ANOTHER) should have alt_filter_pass=True (not in exclude_alts)"

    # Assert - verify the count: 3 should pass (GG, OK_ALT, ANOTHER)
    pass_count = out["alt_filter_pass"].sum()
    assert pass_count == 3, "Should have exactly 3 rows with alt_filter_pass=True (GG, OK_ALT, ANOTHER)"

    # Assert - verify that excluded ALTs are still in the DataFrame
    all_alts = out["ALT"].tolist()
    assert "BAD_ALT" in all_alts, "'BAD_ALT' should still be present in DataFrame"
    assert "ZZZ" in all_alts, "'ZZZ' should still be present in DataFrame"


def test_filter_by_alt_values_drop_left_right(kestrel_config_mock):
    """
    Test that intermediate columns like 'left' and 'right' are retained for debugging.

    The refactored implementation keeps all columns including intermediate ones
    like 'left' and 'right' for debugging purposes. These columns may be dropped
    in a later finalization step downstream.
    """
    # Arrange
    df = pd.DataFrame(
        {
            "ALT": ["GG", "ABC"],
            "Depth_Score": [0.05, 0.02],
            "left": ["some_left_data", "some_left_data"],
            "right": ["some_right_data", "some_right_data"],
        }
    )

    # Act
    out = filter_by_alt_values_and_finalize(df, kestrel_config_mock)

    # Assert - both rows are retained
    assert len(out) == 2, "Expected 2 rows total."

    # Assert - alt_filter_pass column exists
    assert "alt_filter_pass" in out.columns, "Should have 'alt_filter_pass' column"

    # Assert - both rows pass the filter
    # Row 0: GG with Depth_Score=0.05 (above threshold) => should pass
    # Row 1: ABC (not in exclude_alts) => should pass
    assert out.loc[0, "alt_filter_pass"], "Row 0 (GG) should pass the filter"
    assert out.loc[1, "alt_filter_pass"], "Row 1 (ABC) should pass the filter"

    # Assert - 'left' and 'right' columns are RETAINED for debugging
    # This is a change from the old behavior where they were dropped
    assert "left" in out.columns and "right" in out.columns, (
        "Intermediate columns 'left' and 'right' should be retained for debugging"
    )

    # Assert - verify the data in left/right columns is unchanged
    assert all(out["left"] == "some_left_data"), "'left' column data should be preserved"
    assert all(out["right"] == "some_right_data"), "'right' column data should be preserved"


# ---------------------------------------------------------------------------
# read_vcf_without_comments()
# ---------------------------------------------------------------------------
#
# HEADER_LINE below has 5 columns; every data-shaped line used in these tests
# also has exactly 5 tab-separated fields, so a line that leaks through a
# mutated filter still produces a row pandas can construct (no accidental
# shape-mismatch crash masking the real assertion).

HEADER_LINE = "#CHROM\tPOS\tID\tREF\tALT"


def _write_vcf(path, lines, gzipped):
    """Write `lines` (no trailing newlines) to `path`, one per line, plain or gzip."""
    content = "\n".join(lines) + "\n"
    if gzipped:
        with gzip.open(path, "wt") as f:
            f.write(content)
    else:
        with open(path, "w") as f:
            f.write(content)


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "gzip"])
def test_read_vcf_without_comments_parses_meta_header_and_data_rows(tmp_path, gzipped):
    """Baseline: '##' meta lines are skipped, '#CHROM' becomes the header, data rows are kept.

    Covers the open_func selection branch (plain vs. gzip) and the normal,
    fully-populated path through both compound conditions. This also happens
    to kill the line 59 `not`-deleted mutant on its own: under that mutant no
    ordinary data row (which never starts with '##') could ever be appended,
    so the two data rows below would vanish entirely.
    """
    vcf_path = tmp_path / ("sample.vcf.gz" if gzipped else "sample.vcf")
    _write_vcf(
        vcf_path,
        [
            "##fileformat=VCFv4.2",
            "##source=test",
            HEADER_LINE,
            "chr1\t100\trs1\tA\tT",
            "chr1\t200\trs2\tG\tC",
        ],
        gzipped,
    )

    out = read_vcf_without_comments(str(vcf_path))

    assert list(out.columns) == ["#CHROM", "POS", "ID", "REF", "ALT"]
    assert out.shape == (2, 5)
    assert out.loc[0, "ID"] == "rs1"
    assert out.loc[0, "POS"] == "100"
    assert out.loc[1, "ID"] == "rs2"
    assert out.loc[1, "ALT"] == "C"


def test_read_vcf_without_comments_missing_file_returns_empty_dataframe(tmp_path, caplog):
    """A nonexistent path hits the FileNotFoundError branch and returns an empty DataFrame."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.variant_parsing")
    missing_path = tmp_path / "does_not_exist.vcf"

    out = read_vcf_without_comments(str(missing_path))

    assert out.empty
    assert list(out.columns) == []
    assert any("VCF file not found" in r.message and str(missing_path) in r.message for r in caplog.records), (
        "Expected a 'VCF file not found' error log naming the missing path."
    )


def test_read_vcf_without_comments_unreadable_path_returns_empty_dataframe(tmp_path, caplog):
    """A path that raises something other than FileNotFoundError (e.g. IsADirectoryError)
    hits the generic `except Exception` branch and still returns an empty DataFrame
    rather than propagating.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.variant_parsing")
    # Opening a directory with open(..., "rt") raises IsADirectoryError, a subclass of
    # OSError but NOT of FileNotFoundError, so this exercises the second except clause.
    directory_path = tmp_path

    out = read_vcf_without_comments(str(directory_path))

    assert out.empty
    assert list(out.columns) == []
    records = [record for record in caplog.records if "Error reading VCF file" in record.getMessage()]
    assert len(records) == 1
    assert records[0].levelno == logging.ERROR


def test_read_vcf_without_comments_no_header_line_returns_empty_dataframe(tmp_path):
    """A file with only '##' meta lines and no '#CHROM' line never sets `header`,
    so no line is ever collected as data and the result is empty.
    """
    vcf_path = tmp_path / "no_header.vcf"
    _write_vcf(vcf_path, ["##fileformat=VCFv4.2", "##source=test"], gzipped=False)

    out = read_vcf_without_comments(str(vcf_path))

    assert out.empty
    assert list(out.columns) == []


def test_data_line_before_the_header_is_dropped_not_leaked_in(tmp_path):
    """Kills the line 59 `elif not line.startswith("##") and header:` `and` -> `or` mutant.

    Truth-table row P=False, H=False (see module docstring): a line that does not
    start with '##' but arrives *before* the '#CHROM' header is seen must be
    dropped under the real `and` (header is falsy, so the whole condition is
    falsy). Under the `or` mutant, `not line.startswith("##")` alone is True,
    so the pre-header line would be appended to `data` regardless of `header`,
    producing an extra row. No existing test put a data-shaped line before the
    header, so nothing previously distinguished `and` from `or` here.
    """
    vcf_path = tmp_path / "leading_data_line.vcf"
    _write_vcf(
        vcf_path,
        [
            "chr0\t1\tpre_header\tA\tT",  # not '##', but header is still None here
            HEADER_LINE,
            "chr1\t100\trs1\tA\tT",
        ],
        gzipped=False,
    )

    out = read_vcf_without_comments(str(vcf_path))

    # Under the real `and`, only the post-header row survives.
    assert len(out) == 1, "The pre-header line must not be captured as a data row."
    assert out.iloc[0]["ID"] == "rs1"
    assert "pre_header" not in out["ID"].tolist()


def test_meta_comment_after_the_header_is_still_skipped_not_captured(tmp_path):
    """Kills the line 59 `elif not line.startswith("##") and header:` `not`-deleted mutant.

    Truth-table row P=True, H=True (see module docstring): a '##' comment line
    appearing *after* the header is seen must still be skipped, because the
    real condition requires `not line.startswith("##")`. If `not` is deleted,
    the condition becomes `line.startswith("##") and header`, which is True
    for this line (and False for every subsequent ordinary data row, since
    those never start with '##') -- so the mutant would capture the comment
    line as the sole data row and drop the real one.
    """
    vcf_path = tmp_path / "late_comment.vcf"
    _write_vcf(
        vcf_path,
        [
            HEADER_LINE,
            "##late_comment\tX\tY\tZ\tW",  # 5 fields, matches the header's column count
            "chr1\t100\trs1\tA\tT",
        ],
        gzipped=False,
    )

    out = read_vcf_without_comments(str(vcf_path))

    assert len(out) == 1, "Exactly the one real data row should survive."
    assert out.iloc[0]["ID"] == "rs1", "The '##' comment row must not be captured as data."
    assert out.iloc[0]["#CHROM"] == "chr1", "The '##' comment row must not be captured as data."


def test_header_only_no_data_rows_returns_a_columnless_empty_dataframe(tmp_path, caplog):
    """Kills the line 68 `if data and header:` `and` -> `or` mutant.

    Truth-table case D=False, H=True (see module docstring): a file with a
    '#CHROM' header line and no data rows leaves `data == []` (falsy) and
    `header` truthy. The real `and` short-circuits to False, so the function
    returns a bare `pd.DataFrame()` with no columns. Under the `or` mutant,
    `header` alone makes the condition True, so it would instead return
    `pd.DataFrame([], columns=header)` -- still row-empty, but carrying the
    header's 5 columns. `out.empty` is True either way, so the assertion has
    to inspect the columns, not just emptiness.
    """
    vcf_path = tmp_path / "header_only.vcf"
    _write_vcf(vcf_path, [HEADER_LINE], gzipped=False)

    with caplog.at_level(logging.ERROR, logger="vntyper.scripts.variant_parsing"):
        out = read_vcf_without_comments(str(vcf_path))

    assert out.empty
    assert list(out.columns) == [], (
        "A header with no data rows must produce a columnless empty DataFrame, not one carrying the header's columns."
    )
    assert not [record for record in caplog.records if record.levelno >= logging.ERROR]
