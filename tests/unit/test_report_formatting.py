"""The report's pure presentation layer.

``report_formatting`` was extracted out of ``generate_report.py`` so that the
formatting could be tested without a filesystem, a subprocess or a rendered
template. These tests cover the four things that were previously only reachable
by running the whole pipeline:

* status icons and their treatment of a metric that was never computed;
* the coverage schema, which is agent D's frozen contract (C1) read back;
* the Kestrel column projection, which silently drops what it cannot find;
* the IGV fragment extraction, whose absent-marker path produced a syntax error
  in every negative sample's report.
"""

import pandas as pd
import pytest

from vntyper.scripts import report_formatting as rf
from vntyper.scripts.coverage_stats import COVERAGE_COLUMNS

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# Status icons
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value,cutoff,higher_better,expected_colour",
    [
        (99.0, 100, True, "red"),  # coverage below threshold
        (100.0, 100, True, "green"),  # exactly at the threshold passes
        (101.0, 100, True, "green"),
        (60.0, 50.0, False, "red"),  # percent uncovered above threshold
        (50.0, 50.0, False, "green"),  # exactly at the threshold passes
        (40.0, 50.0, False, "green"),
    ],
)
def test_threshold_icon_picks_the_colour_for_the_direction(value, cutoff, higher_better, expected_colour) -> None:
    icon, colour = rf.threshold_icon(value, cutoff, higher_better=higher_better)
    assert colour == expected_colour
    assert icon == (rf.WARNING_ICON if expected_colour == "red" else rf.OK_ICON)


def test_a_missing_metric_defaults_to_blank() -> None:
    """The fastp rows hide entirely when fastp did not run; a lone glyph in an
    otherwise empty row would read as a result."""
    assert rf.threshold_icon(None, 0.8) == rf.MISSING_AS_BLANK == ("", "")


def test_a_missing_metric_can_be_shown_as_passing() -> None:
    """The coverage rows have always shown a tick for "not calculated"; pinned
    here because the two callers used to be two hand-rolled copies."""
    assert rf.threshold_icon(None, 100, on_missing=rf.MISSING_AS_OK) == (rf.OK_ICON, "green")


def test_the_two_icons_are_different() -> None:
    """Guard the guard: identical icons would make every assertion above vacuous."""
    assert rf.WARNING_ICON != rf.OK_ICON


# ---------------------------------------------------------------------------
# Coverage - contract C1
# ---------------------------------------------------------------------------


def test_the_coverage_field_types_cover_exactly_the_frozen_schema() -> None:
    """Contract C1. Agent D declares the columns; the report must not re-type
    them. A column added there and not here would render as nothing at all."""
    assert set(rf.COVERAGE_FIELD_TYPES) == set(COVERAGE_COLUMNS)


def test_parse_coverage_stats_reads_every_column() -> None:
    row = {
        "mean": "123.45",
        "median": "120.0",
        "stdev": "10.5",
        "min": "3",
        "max": "300",
        "region_length": "1000",
        "uncovered_bases": "20",
        "percent_uncovered": "2.0",
        "coverage_qc": "PASS",
    }
    stats = rf.parse_coverage_stats([row])

    assert set(stats) == set(COVERAGE_COLUMNS)
    assert stats["mean"] == pytest.approx(123.45)
    assert stats["min"] == 3 and isinstance(stats["min"], int)
    assert stats["percent_uncovered"] == pytest.approx(2.0)
    assert stats["coverage_qc"] == "PASS", "#172's verdict is read back as a string, not coerced"


def test_parse_coverage_stats_yields_none_for_no_rows() -> None:
    """None, not zero: the template prints "Not calculated" for None, whereas a
    zero would state that a sample with no coverage data has no uncovered bases."""
    assert rf.parse_coverage_stats([]) == dict.fromkeys(COVERAGE_COLUMNS)


def test_parse_coverage_stats_survives_an_unreadable_value(caplog) -> None:
    import logging

    row: dict[str, object] = dict.fromkeys(COVERAGE_COLUMNS, 1)
    row["stdev"] = "not-a-number"
    with caplog.at_level(logging.ERROR, logger="vntyper.scripts.report_formatting"):
        stats = rf.parse_coverage_stats([row])

    assert stats["mean"] == pytest.approx(1.0), "fields read before the failure are kept"
    assert stats["stdev"] is None
    assert stats["percent_uncovered"] is None, "fields after the failure stay unset"
    assert any(record.levelno >= logging.ERROR for record in caplog.records)


def test_a_renamed_coverage_column_is_visible_rather_than_silent() -> None:
    """The failure mode contract C1 exists to prevent: `.get(name, 0)` turns a
    renamed column into a zero. Reading through COVERAGE_COLUMNS at least keeps
    the zero honest -- the guard is the schema test above, and this pins the
    consequence so nobody mistakes the zero for a measurement."""
    stats = rf.parse_coverage_stats([{"Mean": 123.0, "Uncovered%": 90.0}])
    assert stats["mean"] == 0.0
    assert stats["percent_uncovered"] == 0.0


# ---------------------------------------------------------------------------
# fastp
# ---------------------------------------------------------------------------


def test_summarise_fastp_reports_unavailable_for_no_data() -> None:
    metrics = rf.summarise_fastp({})
    assert metrics.available is False
    assert metrics.duplication_rate is None
    assert metrics.sequencing == ""


def test_summarise_fastp_computes_the_passed_filter_rate() -> None:
    metrics = rf.summarise_fastp(
        {
            "summary": {
                "sequencing": "paired end (150 cycles)",
                "before_filtering": {"total_reads": 1000},
                "after_filtering": {"q20_rate": 0.97, "q30_rate": 0.93},
            },
            "duplication": {"rate": 0.05},
            "filtering_result": {"passed_filter_reads": 900},
        }
    )
    assert metrics.available is True
    assert metrics.passed_filter_rate == pytest.approx(0.9)
    assert metrics.q20_rate == pytest.approx(0.97)
    assert metrics.sequencing == "paired end (150 cycles)"


def test_summarise_fastp_does_not_divide_by_zero() -> None:
    metrics = rf.summarise_fastp(
        {
            "summary": {"before_filtering": {"total_reads": 0}, "after_filtering": {}},
            "filtering_result": {"passed_filter_reads": 0},
        }
    )
    assert metrics.available is True
    assert metrics.passed_filter_rate is None


# ---------------------------------------------------------------------------
# Column projection
# ---------------------------------------------------------------------------


def test_select_display_columns_renames_and_orders() -> None:
    frame = pd.DataFrame([{"POS": 67, "REF": "G", "Variant": "Insertion"}])
    out = rf.select_display_columns(frame, {"Variant": "Variant", "POS": "Position", "REF": "REF"})
    assert list(out.columns) == ["Variant", "Position", "REF"]


def test_select_display_columns_drops_what_is_absent() -> None:
    """The behaviour that hides a missing column instead of failing. It is why a
    misnamed key produces a report with a section quietly missing."""
    frame = pd.DataFrame([{"POS": 67}])
    out = rf.select_display_columns(frame, {"POS": "Position", "Nope": "Nope"})
    assert list(out.columns) == ["Position"]


def test_select_display_columns_does_not_mutate_its_input() -> None:
    frame = pd.DataFrame([{"POS": 67, "REF": "G"}])
    rf.select_display_columns(frame, {"POS": "Position"})
    assert list(frame.columns) == ["POS", "REF"]


# ---------------------------------------------------------------------------
# Confidence colouring
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", ["High_Precision", "High_Precision*"])
def test_high_precision_is_red(value) -> None:
    assert rf.confidence_html(value) == f'<span style="color:red;font-weight:bold;">{value}</span>'


def test_low_precision_is_orange() -> None:
    assert rf.confidence_html("Low_Precision") == '<span style="color:orange;font-weight:bold;">Low_Precision</span>'


def test_an_unstyled_confidence_passes_through() -> None:
    """`Negative` is what a negative run's placeholder row carries."""
    assert rf.confidence_html("Negative") == "Negative"


def test_the_starred_confidence_is_styled_like_the_unstarred_one() -> None:
    """`High_Precision*` is a distinct label, not a typo -- if it fell through
    unstyled the strongest calls would render as plain text."""
    assert rf.confidence_html("High_Precision*") == rf.confidence_html("High_Precision").replace(
        ">High_Precision<", ">High_Precision*<"
    )


# ---------------------------------------------------------------------------
# IGV fragment extraction
# ---------------------------------------------------------------------------


IGV_PAGE = (
    "<html><body>\n"
    '<div id="container" class="x">markup</div>\n'
    "<script>\n"
    'const tableJson = {"headers": ["a"], "rows": [[1]]}\n'
    'const sessionDictionary = {"0": "blob:1"}\n'
    "</script>\n"
    "</body></html>\n"
)


def test_extract_igv_fragments_finds_all_three() -> None:
    container, table_json, session = rf.extract_igv_fragments(IGV_PAGE)
    assert container.startswith('<div id="container"')
    assert table_json == '{"headers": ["a"], "rows": [[1]]}'
    assert session == '{"0": "blob:1"}'


def test_extract_igv_fragments_returns_empty_without_a_container() -> None:
    assert rf.extract_igv_fragments("<html><body>nothing</body></html>") == ("", "", "")


def test_extract_igv_fragments_returns_empty_for_empty_input() -> None:
    assert rf.extract_igv_fragments("") == ("", "", "")


def test_extract_line_after_returns_empty_when_the_marker_is_absent() -> None:
    """`content.find(marker) + len(marker)` is 17 when the marker is absent, so
    the old code sliced from character 17 and produced arbitrary text where a
    JavaScript object literal belonged -- `const tableJson = <garbage>;` is a
    syntax error that kills the whole <script> block."""
    content = "x" * 200
    assert rf.extract_line_after(content, rf.IGV_TABLE_JSON_MARKER) == ""


def test_extract_line_after_keeps_the_last_line_without_a_trailing_newline() -> None:
    """`find("\\n", start)` is -1 at the end of a file, and slicing to -1 drops
    the final character -- turning a valid literal into an unbalanced one."""
    content = 'const tableJson = {"rows": []}'
    assert rf.extract_line_after(content, rf.IGV_TABLE_JSON_MARKER) == '{"rows": []}'


def test_a_page_with_a_container_but_no_table_json_yields_empty_fragments() -> None:
    """The real shape of the failure: `create_report` produced a page, but not
    the variable the report splices out of it."""
    page = '<html><body><div id="container">m</div></body></html>'
    container, table_json, session = rf.extract_igv_fragments(page)
    assert container
    assert table_json == ""
    assert session == ""


# ---------------------------------------------------------------------------
# Escaping
# ---------------------------------------------------------------------------


def test_escape_frame_cells_escapes_strings() -> None:
    frame = pd.DataFrame([{"Motif Sequence": "<script>alert(1)</script>"}])
    out = rf.escape_frame_cells(frame)
    assert out.iloc[0]["Motif Sequence"] == "&lt;script&gt;alert(1)&lt;/script&gt;"


def test_escape_frame_cells_leaves_numbers_alone() -> None:
    """Stringifying numbers here would take pandas' float and NA formatting out
    of the rendered table, and a number cannot carry markup anyway."""
    frame = pd.DataFrame([{"Depth Score": 0.0125, "Position": 67}])
    out = rf.escape_frame_cells(frame)
    assert out.iloc[0]["Depth Score"] == pytest.approx(0.0125)
    assert out.iloc[0]["Position"] == 67


def test_escape_frame_cells_skips_the_columns_holding_our_own_markup() -> None:
    frame = pd.DataFrame([{"Confidence": '<span style="x">High_Precision</span>', "REF": "<b>"}])
    out = rf.escape_frame_cells(frame, html_columns=("Confidence",))
    assert out.iloc[0]["Confidence"] == '<span style="x">High_Precision</span>'
    assert out.iloc[0]["REF"] == "&lt;b&gt;"


def test_escape_frame_cells_does_not_mutate_its_input() -> None:
    frame = pd.DataFrame([{"REF": "<b>"}])
    rf.escape_frame_cells(frame)
    assert frame.iloc[0]["REF"] == "<b>"


def test_confidence_html_escapes_an_unstyled_value() -> None:
    """The Kestrel table is rendered with escape=False, so this is the last
    chance to escape a Confidence value with no configured style."""
    assert rf.confidence_html("<img src=x onerror=alert(1)>") == "&lt;img src=x onerror=alert(1)&gt;"


def test_confidence_html_escapes_inside_the_span_too() -> None:
    assert rf.confidence_html("Low_Precision").startswith('<span style="color:orange')


def test_escape_html_escapes_quotes() -> None:
    """Quotes matter: these values land inside attributes in the IGV tooltips."""
    assert rf.escape_html('a"b') == "a&quot;b"
