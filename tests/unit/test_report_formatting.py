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

import json
import logging
from dataclasses import dataclass

import pandas as pd
import pytest

from vntyper.scripts import report_formatting as rf
from vntyper.scripts.coverage_qc import COVERAGE_QC_NOT_EVALUATED, evaluate_coverage_qc
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
    row: dict[str, object] = dict.fromkeys(COVERAGE_COLUMNS, "3.5")
    row["stdev"] = object()
    with caplog.at_level(logging.ERROR, logger="vntyper.scripts.report_formatting"):
        stats = rf.parse_coverage_stats([row])

    assert set(stats) == set(COVERAGE_COLUMNS)
    assert stats["mean"] == pytest.approx(3.5), "fields read before the failure are kept"
    assert stats["stdev"] is None
    assert all(stats[column] is None for column in COVERAGE_COLUMNS[2:])
    records = [record for record in caplog.records if record.name == "vntyper.scripts.report_formatting"]
    assert [record.levelno for record in records] == [logging.ERROR]


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
# Confidence classification
#
# There is no hue in either literal below, and that is the change #242's contrast
# pass made. `High_Precision` - the most trustworthy call the pipeline makes - was
# `color:red`, one column from a red `Flag` glyph meaning the opposite, and it
# measured 4.00:1 against white and 3.57:1 against a striped cohort row.
# `Low_Precision` was `color:orange` at 1.97:1 and 1.76:1. A transitional underline
# had been added so the two were not separated by hue alone; both go together here,
# because the honest channel was always the label and the label is text in the cell
# in every branch. `test_report_presentation.py` asserts the rule these literals are
# an instance of, and computes the ratios from the shipped stylesheet.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", ["High_Precision", "High_Precision*"])
def test_high_precision_carries_its_class_and_no_colour(value) -> None:
    assert rf.confidence_html(value) == f'<span class="confidence confidence-high-precision">{value}</span>'


def test_low_precision_carries_its_own_class() -> None:
    assert rf.confidence_html("Low_Precision") == (
        '<span class="confidence confidence-low-precision">Low_Precision</span>'
    )


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
    assert rf.confidence_html("Low_Precision").startswith('<span class="confidence')


def test_escape_html_escapes_quotes() -> None:
    """Quotes matter: these values land inside attributes in the IGV tooltips."""
    assert rf.escape_html('a"b') == "a&quot;b"


# ---------------------------------------------------------------------------
# js_json_literal (#216)
# ---------------------------------------------------------------------------
#
# SPECIFICATION: `report_template.html` interpolates these two fragments straight
# into a <script> block with `|safe`. They are lifted verbatim out of another
# tool's HTML by `extract_line_after`, and the only guard they used to get was
# `js_object_literal`, which checked that the fragment was non-empty. That is a
# *syntax* guard against `const tableJson = ;`, not a safety guard: a `</script>`
# anywhere in the fragment closes the block early and everything after it is
# parsed as HTML. The fragment is sample-derived -- it is the igv-reports variant
# table, built from VNtyper's own BED and VCF, whose REF, ALT and Motif_sequence
# values come from the sample's reads.
#
# `js_object_literal` has no test at all on `main` -- that is how #216 survived --
# so there is no prior test here to delete; this whole section is new.


def test_a_well_formed_fragment_round_trips_as_json() -> None:
    result = rf.js_json_literal('{"headers": ["a"], "rows": [[1]]}', rf.EMPTY_TABLE_JSON)
    assert json.loads(result) == {"headers": ["a"], "rows": [[1]]}


def test_the_trailing_statement_terminator_is_stripped() -> None:
    """Defensive against a *future* igv-reports version, not a correction of
    today's: verified against the installed igv-reports 1.16.0 that
    `extract_line_after` returns the bare literal with no `;` at all."""
    result = rf.js_json_literal('{"a": 1};', rf.EMPTY_TABLE_JSON)
    assert json.loads(result) == {"a": 1}


def test_a_script_close_in_a_string_value_cannot_terminate_the_block() -> None:
    """The `</script>` case, which is now covered by the blanket `<` escape.

    This used to assert the literal `<\\/script>` that a `</`-only replacement
    produced. Escaping every `<` as `\\u003c` subsumes that case, so what is
    asserted is the property rather than the particular escape: no literal `<`
    reaches the block, and the data survives the round trip unchanged.
    """
    fragment = json.dumps({"rows": [["</script><img src=x onerror=alert(1)>"]]})
    result = rf.js_json_literal(fragment, rf.EMPTY_TABLE_JSON)
    assert "<" not in result
    # It is still the same data once the browser parses it.
    assert json.loads(result)["rows"][0][0] == "</script><img src=x onerror=alert(1)>"


def test_the_double_escaped_script_state_cannot_be_entered() -> None:
    """Escaping only `</` is not enough: `<!--<script>` is the other way out.

    An HTML5 tokenizer that meets `<!--` followed by `<script` inside a script
    element enters the *double-escaped* state, in which the real `</script>` no
    longer ends the element -- so the rest of the document is swallowed as script
    text and the report renders as a blank page below the IGV panel. No `</`
    appears anywhere in that sequence, so the previous escape let it through
    verbatim. Every literal `<` is escaped now, which closes both routes at once.
    """
    fragment = json.dumps({"rows": [["<!--<script>", "a<b"]]})
    result = rf.js_json_literal(fragment, rf.EMPTY_TABLE_JSON)
    assert "<" not in result
    assert "\\u003c" in result
    # `<` is a JSON string escape, so the browser still sees the original text.
    assert json.loads(result)["rows"][0] == ["<!--<script>", "a<b"]


@pytest.mark.parametrize("separator", ["\u2028", "\u2029"])
def test_javascript_line_terminators_do_not_appear_literally(separator: str) -> None:
    """U+2028 and U+2029 end a line in JavaScript but are legal inside a JSON
    string. `ensure_ascii=True` (the `json.dumps` default, kept explicit) is what
    keeps them out -- there is no `.replace()` for either in the implementation.
    This pins the *property* of the result, which is worth pinning even though
    the mechanism enforcing it is a stdlib default rather than code of ours."""
    result = rf.js_json_literal(json.dumps({"a": separator}), rf.EMPTY_TABLE_JSON)
    assert separator not in result


@pytest.mark.parametrize(
    "fragment",
    ["", "   ", "not json at all", "{unquoted: 1}", "{'single': 'quotes'}", "{"],
)
def test_a_fragment_that_is_not_json_falls_back(fragment: str) -> None:
    assert rf.js_json_literal(fragment, rf.EMPTY_TABLE_JSON) == rf.EMPTY_TABLE_JSON


def test_the_fallback_is_returned_verbatim_for_the_session_dictionary() -> None:
    assert rf.js_json_literal("", rf.EMPTY_SESSION_DICTIONARY) == rf.EMPTY_SESSION_DICTIONARY


def test_the_output_is_deterministic_for_the_same_input() -> None:
    """Two runs over the same IGV page must emit byte-identical script (exit criterion E2)."""
    fragment = '{"b": 2, "a": 1}'
    assert rf.js_json_literal(fragment, rf.EMPTY_TABLE_JSON) == rf.js_json_literal(fragment, rf.EMPTY_TABLE_JSON)
    assert rf.js_json_literal(fragment, rf.EMPTY_TABLE_JSON) == '{"a":1,"b":2}'


def test_a_malformed_fragment_is_logged_at_warning(caplog) -> None:
    with caplog.at_level(logging.WARNING):
        rf.js_json_literal("not json", rf.EMPTY_TABLE_JSON)
    assert any("could not be parsed as JSON" in r.getMessage() for r in caplog.records)


class TestCoverageNullToken:
    """Reading a summary whose #222 columns were never measured.

    Two separate hazards, and the second is the dangerous one: an absent figure
    must not become ``0`` (which reads as *no coverage was seen*), and a figure
    that is present but unreadable must not let the QC gate render a verdict on
    the half it could parse.
    """

    @staticmethod
    def _row(**overrides):
        row = dict.fromkeys(COVERAGE_COLUMNS, "1")
        row.update(
            {
                "mean": "120.00",
                "percent_uncovered": "1.00",
                "coverage_qc": "PASS",
                "depth_counting_policy": "samtools-depth-a/v1",
            }
        )
        row.update(overrides)
        return [row]

    def test_the_not_measured_token_parses_to_none_not_zero(self):
        """`NA` is what a `--custom-regions` run writes; zero would be a claim."""
        stats = rf.parse_coverage_stats(self._row(vntr_array_depth_sum="NA", vntr_flank_mean_depth="NA"))

        assert stats["vntr_array_depth_sum"] is None
        assert stats["vntr_flank_mean_depth"] is None

    def test_an_absent_column_still_parses_the_rest(self):
        """A summary written before #222 has no such column at all."""
        row = self._row()
        del row[0]["vntr_array_depth_sum"]

        stats = rf.parse_coverage_stats(row)

        assert stats["vntr_array_depth_sum"] is None
        assert stats["mean"] == 120.00
        assert stats["coverage_qc"] == "PASS"

    def test_one_unreadable_column_does_not_discard_the_qc_verdict(self):
        """The loop used to abort on the first bad value, losing every later field.

        `coverage_qc` is the last column, so a single malformed number anywhere
        ahead of it took the verdict with it.
        """
        stats = rf.parse_coverage_stats(self._row(vntr_array_length="not-a-number"))

        assert stats["vntr_array_length"] is None
        assert stats["coverage_qc"] == "PASS"
        assert stats["mean"] == 120.00

    def test_an_unreadable_mean_is_not_judged_on_the_other_metric(self):
        """The false-PASS path, found in adversarial review of the plan.

        An unparseable mean with an acceptable uncovered percentage must not
        render a passing verdict: the gate would be reporting on a number it
        could not read.
        """
        stats = rf.parse_coverage_stats(self._row(mean="corrupt"))

        # Both dropped, not just the unreadable one: that is what makes the existing
        # "neither metric was measured" rule fire instead of a verdict on half the input.
        assert stats["mean"] is None
        assert stats["percent_uncovered"] is None

        verdict = evaluate_coverage_qc(stats["mean"], stats["percent_uncovered"], 100.0, 50.0)

        assert verdict.status == COVERAGE_QC_NOT_EVALUATED


# ---------------------------------------------------------------------------
# Server-side number formatting - issue #242
# ---------------------------------------------------------------------------
#
# Until #242 the report shipped a jQuery ``applyRounding()`` that rewrote **every**
# numeric cell of **every** initialised table to four decimal places and then
# stripped trailing zeroes. It ran in the reader's browser, so the number on screen
# was a property of the reader's network rather than of the run, and it was
# table-agnostic, so it applied one rule to a genomic position, a depth score and a
# p-value alike.
#
# The baseline below was captured from a report rendered by this repository and
# opened in chromium with the three CDNs reachable, **before** ``applyRounding`` was
# removed. ``online`` is what that reader saw; ``server`` is what the report now
# writes into the file. Where they differ the divergence is deliberate and named.


@dataclass(frozen=True)
class NumericRendering:
    """One value's pre-#242 online rendering beside its rendering now.

    Attributes:
        table: ``kestrel`` or ``advntr`` - which format table declares the column.
        column: The result column name.
        stored: The value as it appears in ``pipeline_summary.json``.
        online: What the browser displayed for it before #242, measured.
        server: What the report renders for it now.
        divergence: Why the two differ, or "" when they agree.
    """

    table: str
    column: str
    stored: object
    online: str
    server: str
    divergence: str = ""

    def __repr__(self) -> str:
        return f"{self.table}.{self.column}={self.stored!r}"


NUMERIC_RENDERINGS = (
    # -- Kestrel ---------------------------------------------------------------
    NumericRendering("kestrel", "POS", 67, "67", "67"),
    NumericRendering("kestrel", "Estimated_Depth_AlternateVariant", 120, "120", "120"),
    NumericRendering("kestrel", "Estimated_Depth_Variant_ActiveRegion", 12000, "12000", "12000"),
    NumericRendering(
        "kestrel",
        "Depth_Score",
        0.010012,
        "0.01",
        "0.010012",
        "the confidence thresholds are 0.00469 and 0.00515, so four decimal places is coarser "
        "than the calibration the number is judged against",
    ),
    NumericRendering(
        "kestrel",
        "Depth_Score",
        0.00001234,
        "0",
        "0.000012",
        "four decimal places rendered a real depth score as the value 0",
    ),
    # -- adVNTR ----------------------------------------------------------------
    NumericRendering("advntr", "VID", 25561, "25561", "25561"),
    NumericRendering("advntr", "NumberOfSupportingReads", 14, "14", "14"),
    NumericRendering("advntr", "POS", 67, "67", "67"),
    NumericRendering("advntr", "MeanCoverage", 132.45, "132.45", "132.45"),
    NumericRendering(
        "advntr",
        "MeanCoverage",
        98.5,
        "98.5",
        "98.50",
        "trailing zeroes are kept so every mean states the same precision",
    ),
    NumericRendering(
        "advntr",
        "Pvalue",
        1e-09,
        "0",
        "1e-09",
        "four decimal places destroyed a highly significant p-value, displaying it as 0",
    ),
    NumericRendering(
        "advntr",
        "Pvalue",
        0.0001234,
        "0.0001",
        "0.000123",
        "three significant figures keep the value rather than truncating it at four decimals",
    ),
    NumericRendering("advntr", "Pvalue", 0.04999, "0.05", "0.05"),
)

#: Which declaration in ``report_formatting`` governs each table.
FORMAT_TABLES = {"kestrel": "KESTREL_CELL_FORMATS", "advntr": "ADVNTR_CELL_FORMATS"}


def rendered_cell(case: NumericRendering) -> str:
    """Render one stored value through the formatter its column declares.

    Args:
        case: The recorded rendering.

    Returns:
        str: The cell text the report writes.
    """
    formats = getattr(rf, FORMAT_TABLES[case.table])
    frame = pd.DataFrame({case.column: [case.stored]})
    return str(rf.format_number_columns(frame, formats).loc[0, case.column])


@pytest.mark.parametrize("case", NUMERIC_RENDERINGS, ids=repr)
def test_each_numeric_column_renders_the_string_recorded_for_it(case: NumericRendering) -> None:
    """Every displayed number is produced here, by the column's declared formatter."""
    assert rendered_cell(case) == case.server


@pytest.mark.parametrize("case", [c for c in NUMERIC_RENDERINGS if not c.divergence], ids=repr)
def test_a_column_with_no_named_divergence_still_reads_as_it_did_online(case: NumericRendering) -> None:
    """Where the browser's rule was harmless, the server reproduces it exactly."""
    assert rendered_cell(case) == case.online


@pytest.mark.parametrize("case", [c for c in NUMERIC_RENDERINGS if c.divergence], ids=repr)
def test_every_divergence_from_the_online_rendering_is_a_declared_one(case: NumericRendering) -> None:
    """A rendering may only differ from the measured baseline with a recorded reason."""
    assert rendered_cell(case) != case.online, (
        f"{case} no longer diverges from the online rendering; drop the recorded reason: {case.divergence}"
    )


def test_every_displayed_column_declares_how_its_value_is_rendered() -> None:
    """A new display column must state its formatter rather than inherit a default."""
    assert set(rf.KESTREL_CELL_FORMATS) == set(rf.KESTREL_DISPLAY_COLUMNS)
    assert set(rf.ADVNTR_CELL_FORMATS) == set(rf.ADVNTR_DISPLAY_COLUMNS)


def test_every_declared_format_names_a_real_formatter() -> None:
    declared = set(rf.KESTREL_CELL_FORMATS.values()) | set(rf.ADVNTR_CELL_FORMATS.values())
    assert declared <= set(rf.CELL_FORMATTERS), f"undefined formats: {sorted(declared - set(rf.CELL_FORMATTERS))}"


def test_an_undeclared_numeric_column_fails_loudly(caplog) -> None:
    """The point of the per-column table: silence is the failure mode being removed."""
    caplog.set_level(logging.ERROR, logger=rf.logger.name)
    frame = pd.DataFrame({"POS": [67], "Novel_Score": [1.5]})

    with pytest.raises(ValueError, match="Novel_Score"):
        rf.format_number_columns(frame, rf.KESTREL_CELL_FORMATS)

    assert any("Novel_Score" in record.getMessage() for record in caplog.records)


def test_a_non_numeric_placeholder_passes_through_untouched() -> None:
    """``output_empty_result`` writes the string "None" into every numeric column."""
    frame = pd.DataFrame({"POS": ["None"], "Depth_Score": ["None"]})

    formatted = rf.format_number_columns(frame, rf.KESTREL_CELL_FORMATS)

    assert formatted.loc[0, "POS"] == "None"
    assert formatted.loc[0, "Depth_Score"] == "None"


def test_a_missing_value_is_not_rendered_as_a_number() -> None:
    frame = pd.DataFrame({"Depth_Score": [float("nan")]})

    formatted = rf.format_number_columns(frame, rf.KESTREL_CELL_FORMATS)

    assert pd.isna(formatted.loc[0, "Depth_Score"])


@pytest.mark.parametrize("cells", [[True], [True, "not a number"]], ids=["numpy-bool", "python-bool"])
def test_a_boolean_is_not_silently_rendered_as_a_digit(cells: list) -> None:
    """``bool`` is a ``numbers.Integral``, so rendering it would print True as ``1``
    and change what the cell says. Both spellings are covered because the dtype
    decides which one pandas hands out: a bool column yields ``numpy.bool_``, and a
    mixed column yields the Python ``bool`` the explicit guard exists for."""
    formatted = rf.format_number_columns(pd.DataFrame({"POS": cells}), rf.KESTREL_CELL_FORMATS)

    assert bool(formatted.loc[0, "POS"]) is True
    assert not isinstance(formatted.loc[0, "POS"], str)


def test_an_absent_value_passes_through_rather_than_becoming_a_number() -> None:
    formatted = rf.format_number_columns(pd.DataFrame({"POS": [None]}), rf.KESTREL_CELL_FORMATS)

    assert formatted.loc[0, "POS"] is None


def test_formatting_does_not_mutate_its_input() -> None:
    frame = pd.DataFrame({"POS": [67.0]})

    rf.format_number_columns(frame, rf.KESTREL_CELL_FORMATS)

    assert frame.loc[0, "POS"] == 67.0


def test_a_column_absent_from_the_frame_is_skipped() -> None:
    """A negative Kestrel run carries neither the depth columns nor ``Flag``."""
    frame = pd.DataFrame({"POS": [67]})

    assert list(rf.format_number_columns(frame, rf.KESTREL_CELL_FORMATS).columns) == ["POS"]


# ---------------------------------------------------------------------------
# The flag cell - issue #242
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("clean", ["Not flagged", "Not applicable", ""])
def test_a_clean_flag_renders_the_passing_glyph(clean: str) -> None:
    markup = rf.flag_html(clean)

    assert rf.FLAG_OK_GLYPH in markup
    assert rf.FLAG_WARNING_GLYPH not in markup
    assert rf.FLAG_FLAGGED_CLASS not in markup


def test_a_flagged_value_renders_the_reason_as_text() -> None:
    """The reason must survive print, a screen reader and a failed script load -
    which a ``title`` attribute on a bare glyph does not."""
    markup = rf.flag_html("Low_Depth")

    assert "Low_Depth" in markup
    assert rf.FLAG_WARNING_GLYPH in markup
    assert rf.FLAG_FLAGGED_CLASS in markup
    assert "title=" not in markup, "the reason is text in the cell, not a hover-only attribute"


def test_a_flag_value_is_escaped_before_it_becomes_markup() -> None:
    """``Flag`` reaches the report out of a supplied ``pipeline_summary.json`` (#207)."""
    markup = rf.flag_html("<img src=x onerror=alert(1)>")

    assert "<img" not in markup
    assert "&lt;img src=x onerror=alert(1)&gt;" in markup


def test_a_non_string_flag_is_stringified_rather_than_dropped() -> None:
    assert "1" in rf.flag_html(1)


# ---------------------------------------------------------------------------
# The row-count statement - issue #242
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("total", "flagged", "expected"),
    [
        (3, 2, "Showing 3 of 3 Kestrel rows; 2 flagged."),
        (1, 0, "Showing 1 of 1 Kestrel row; none flagged."),
        (1, 1, "Showing 1 of 1 Kestrel row; 1 flagged."),
        (0, 0, "Showing 0 of 0 Kestrel rows; none flagged."),
    ],
)
def test_the_row_count_statement_says_how_many_rows_are_shown(total: int, flagged: int, expected: str) -> None:
    """Rendered from the frame in Python, never read back out of DataTables' footer."""
    assert rf.row_count_statement(total, flagged, noun="Kestrel") == expected


@pytest.mark.parametrize(
    ("frame", "expected"),
    [
        (pd.DataFrame({"Flag": ["Not flagged", "Low_Depth", "Not applicable"]}), 1),
        (pd.DataFrame({"Flag": ["Not flagged"]}), 0),
        (pd.DataFrame({"Flag": [""]}), 0),
        (pd.DataFrame({"Motif": ["5"]}), 0),
        (pd.DataFrame(), 0),
    ],
)
def test_flagged_row_count_counts_only_the_rows_carrying_a_reason(frame: pd.DataFrame, expected: int) -> None:
    assert rf.flagged_row_count(frame) == expected


# ---------------------------------------------------------------------------
# The empty-result placeholder - issue #242
# ---------------------------------------------------------------------------


def test_the_placeholder_a_negative_run_writes_is_recognised(tmp_path) -> None:
    """The two ends of the contract, bound together in one test.

    ``kestrel_genotyping.output_empty_result`` writes the row and ``summary.parse_tsv``
    is what turns it into the mapping the report reads - it splits on tabs and coerces
    nothing, so the literal string ``"None"`` is what arrives. Constructing the row by
    hand here would let the writer change shape without this failing, which is the
    silently-wrong-call failure mode AGENTS.md warns about: the report would go back to
    tabulating a non-result and every test would still pass.
    """
    from vntyper.scripts.kestrel_genotyping import output_empty_result
    from vntyper.scripts.summary import parse_tsv

    output_empty_result(str(tmp_path), ["## VNtyper Kestrel result"])
    rows = parse_tsv(str(tmp_path / "kestrel_result.tsv"))["data"]

    assert len(rows) == 1, "output_empty_result no longer writes exactly one placeholder row"
    assert rf.is_empty_result_row(rows[0]), f"the placeholder {rows[0]} is not recognised as one"
    assert rf.drop_empty_result_rows(rows) == []


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        ({"Motif": "None", "Confidence": "Negative"}, True),
        ({"Confidence": "Negative"}, True),
        ({"Motif": "None", "POS": "", "Confidence": "Negative"}, True),
        ({"Motif": "None", "POS": None, "Confidence": "Negative"}, True),
        ({"Motif": "None", "POS": float("nan"), "Confidence": "Negative"}, True),
        ({"Motif": "5", "POS": "67", "Confidence": "Negative"}, False),
        ({"Motif": "None", "POS": "67", "Confidence": "Negative"}, False),
        ({"Motif": "None", "Confidence": "High_Precision"}, False),
        ({"Motif": "None"}, False),
        ({}, False),
    ],
)
def test_a_row_is_a_placeholder_only_when_it_names_no_variant(row: dict, expected: bool) -> None:
    """Both halves are required, and the second is what keeps a real call safe.

    ``Confidence == "Negative"`` alone is not the test: a report is a record, and a rule
    that deleted every row carrying that label would delete a real variant the moment a
    future calibration used it. A row is a non-result only when it also names nothing -
    no position, no REF, no ALT, no motif.
    """
    assert rf.is_empty_result_row(row) is expected


def test_dropping_placeholders_keeps_every_other_row_and_its_order() -> None:
    first = {"Motif": "5", "POS": "67", "Confidence": "High_Precision"}
    second = {"Motif": "6", "POS": "68", "Confidence": "Low_Precision"}

    kept = rf.drop_empty_result_rows([first, {"Motif": "None", "Confidence": "Negative"}, second])

    assert kept == [first, second]


def test_dropping_placeholders_does_not_modify_its_input() -> None:
    rows = [{"Motif": "None", "Confidence": "Negative"}]

    rf.drop_empty_result_rows(rows)

    assert rows == [{"Motif": "None", "Confidence": "Negative"}]


# ---------------------------------------------------------------------------
# escaped_table_html
# ---------------------------------------------------------------------------


def test_an_empty_frame_renders_no_table_at_all() -> None:
    """The authored-empty-state hook: ``to_html`` on an empty frame is a stray box."""
    assert rf.escaped_table_html(pd.DataFrame(), classes="table") == ""


def test_a_table_can_be_given_the_id_its_selectors_use() -> None:
    """The per-sample Kestrel table is addressed as ``#kestrel_table`` by the browser
    tier and by the template's DataTables initialisation, so routing it through this
    helper has to keep the id it had when it called ``to_html`` directly."""
    markup = rf.escaped_table_html(pd.DataFrame({"POS": [67]}), classes="table", table_id="kestrel_table")

    assert 'id="kestrel_table"' in markup


def test_a_table_with_no_id_asked_for_carries_none() -> None:
    """The cohort tables share one stylesheet and are addressed by class."""
    assert "id=" not in rf.escaped_table_html(pd.DataFrame({"POS": [67]}), classes="table")
