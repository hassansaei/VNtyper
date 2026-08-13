"""
report_formatting.py

Module Purpose:
---------------
The pure presentation half of the HTML report: turning values into the strings
the Jinja2 template interpolates. No filesystem, no subprocess, no pandas I/O --
every function here takes values and returns values, which is what makes the
report's formatting testable at all.

It was extracted from ``generate_report.py`` (861 LOC, 4% covered) under
AGENTS.md rule 3. What stayed behind is the part that cannot be called without a
filesystem: loading ``pipeline_summary.json``, shelling out to ``create_report``
and writing the rendered HTML.

Four groups live here:

* **Status icons.** A value, a cutoff and a direction produce a coloured glyph
  and a colour name. The report shows eight of them and they all shared one
  ``higher_better`` flag plus two hand-rolled copies with a different treatment
  of a missing value; :func:`threshold_icon` takes that treatment as an argument
  instead.
* **Table formatting.** Column selection and renaming for the Kestrel table, the
  colour-coding of the ``Confidence`` column, the ``Flag`` cell, and the
  per-column number formatting the browser used to do (#242).
* **Row counts.** The visible/total statement printed beside each results table,
  computed from the frame rather than read back out of DataTables' footer.
* **IGV fragment extraction.** ``create_report`` writes a standalone HTML page;
  the report splices three pieces out of it. The splicing is string work and is
  pure, so it lives here; opening the file does not.

Functions:
    threshold_icon: Value + cutoff to (icon, colour)
    select_display_columns: Project and rename a results frame for display
    confidence_html: Colour-code one ``Confidence`` value
    flag_html: One ``Flag`` value to a glyph plus its reason in words
    format_number_columns: Render a frame's numbers with the formatter each column declares
    flagged_row_count: How many rows of a results frame carry a flag
    row_count_statement: The visible/total sentence shown beside a results table
    extract_line_after: Text following a marker, up to the next newline
    extract_igv_fragments: IGV page text to (container, tableJson, sessionDictionary)
"""

from __future__ import annotations

import html
import json
import logging
import numbers
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any

import pandas as pd

# Contract C1: the coverage TSV field names are declared once, by the module that
# produces them. Re-typing the strings here is how the report silently loses
# coverage when a column is renamed - `.get(name, 0)` raises nothing.
from vntyper.scripts.coverage_stats import _BUILD_COMPARABLE_COLUMNS, COVERAGE_COLUMNS, COVERAGE_NULL_TOKEN

logger = logging.getLogger(__name__)

#: Red warning triangle, shown when a metric fails its threshold.
WARNING_ICON = '<span style="color:red;font-weight:bold;">&#9888;</span>'

#: Green tick, shown when a metric passes its threshold.
OK_ICON = '<span style="color:green;font-weight:bold;">&#10004;</span>'

#: What a metric with no value at all renders as. The coverage rows show a tick
#: (the report has always treated "not calculated" as "not a problem"); the fastp
#: rows show nothing, because the whole fastp section is hidden when fastp did not
#: run and a lone glyph in an otherwise empty row reads as a result.
MISSING_AS_OK: tuple[str, str] = (OK_ICON, "green")
MISSING_AS_BLANK: tuple[str, str] = ("", "")

#: Kestrel result column -> report heading, in display order. Columns absent from
#: the results frame are dropped, so this is a superset of any one run's output:
#: a negative run's ``kestrel_result.tsv`` carries neither ``Flag`` nor the depth
#: columns.
#:
#: The motif key is ``Motif``, not ``Motifs`` (contract C3). Both columns exist and
#: both are load-bearing upstream -- the shipped Kestrel flagging expressions
#: ``eval`` against singular ``Motif`` while duplicate ordering uses plural
#: ``Motifs``, and a missing name there evaluates to False rather than raising
#: (AGENTS.md trap 3), so neither may be renamed. ``Motifs`` is the raw
#: ``left-right`` pair Kestrel emits; ``Motif`` is the motif the variant was
#: annotated onto, which is what ``cohort_summary.py`` shows and what the heading
#: has always claimed to be.
KESTREL_DISPLAY_COLUMNS: dict[str, str] = {
    "Motif": "Motif",
    "Variant": "Variant",
    "POS": "Position",
    "REF": "REF",
    "ALT": "ALT",
    "Motif_sequence": "Motif Sequence",
    "Estimated_Depth_AlternateVariant": "Depth (Variant)",
    "Estimated_Depth_Variant_ActiveRegion": "Depth (Region)",
    "Depth_Score": "Depth Score",
    "Confidence": "Confidence",
    "Flag": "Flag",
}

#: adVNTR result columns, in display order. Unlike Kestrel these are not renamed.
ADVNTR_DISPLAY_COLUMNS: tuple[str, ...] = (
    "VID",
    "Variant",
    "NumberOfSupportingReads",
    "MeanCoverage",
    "Pvalue",
    "RU",
    "POS",
    "REF",
    "ALT",
    "Flag",
)

#: How a displayed value is turned into text. Every column of both results tables
#: names one of these, so adding a display column without deciding how its value is
#: rendered fails in :func:`format_number_columns` rather than silently inheriting
#: whatever pandas would have printed.
#:
#: Until #242 there was no such decision to make: the template shipped a jQuery
#: ``applyRounding()`` that rewrote every numeric cell of every initialised table to
#: four decimal places and stripped trailing zeroes. It ran in the reader's browser,
#: so the figure on screen depended on whether three CDNs resolved, and it was
#: table-agnostic, so one rule governed a genomic position, a depth score and a
#: p-value alike.
FORMAT_TEXT = "text"
FORMAT_INTEGER = "integer"
FORMAT_TWO_DECIMALS = "two-decimals"
FORMAT_SIX_DECIMALS = "six-decimals"
FORMAT_THREE_SIGNIFICANT = "three-significant"

#: Kestrel result column -> how its value is rendered. Keys are exactly
#: :data:`KESTREL_DISPLAY_COLUMNS`'s.
#:
#: ``Depth_Score`` is six decimal places rather than the browser's four: the
#: confidence calibration it is judged against is stated to five
#: (``kestrel_config.json``: 0.00469 and 0.00515), and four decimal places printed a
#: score of 1.234e-05 as the value ``0``.
KESTREL_CELL_FORMATS: dict[str, str] = {
    "Motif": FORMAT_TEXT,
    "Variant": FORMAT_TEXT,
    "POS": FORMAT_INTEGER,
    "REF": FORMAT_TEXT,
    "ALT": FORMAT_TEXT,
    "Motif_sequence": FORMAT_TEXT,
    "Estimated_Depth_AlternateVariant": FORMAT_INTEGER,
    "Estimated_Depth_Variant_ActiveRegion": FORMAT_INTEGER,
    "Depth_Score": FORMAT_SIX_DECIMALS,
    "Confidence": FORMAT_TEXT,
    "Flag": FORMAT_TEXT,
}

#: The same decisions keyed by the heading each Kestrel column is renamed to, which
#: is what the display frame carries by the time it is formatted.
KESTREL_DISPLAY_CELL_FORMATS: dict[str, str] = {
    KESTREL_DISPLAY_COLUMNS[column]: cell_format for column, cell_format in KESTREL_CELL_FORMATS.items()
}

#: adVNTR result column -> how its value is rendered. Keys are exactly
#: :data:`ADVNTR_DISPLAY_COLUMNS`'s.
#:
#: ``Pvalue`` is three significant figures rather than four decimal places, because
#: ``parseFloat((1e-9).toFixed(4)).toString()`` is ``"0"`` -- the online report
#: displayed a highly significant p-value as zero. Significant figures keep the
#: magnitude of a value whose whole meaning is its magnitude.
ADVNTR_CELL_FORMATS: dict[str, str] = {
    "VID": FORMAT_INTEGER,
    "Variant": FORMAT_TEXT,
    "NumberOfSupportingReads": FORMAT_INTEGER,
    "MeanCoverage": FORMAT_TWO_DECIMALS,
    "Pvalue": FORMAT_THREE_SIGNIFICANT,
    "RU": FORMAT_TEXT,
    "POS": FORMAT_INTEGER,
    "REF": FORMAT_TEXT,
    "ALT": FORMAT_TEXT,
    "Flag": FORMAT_TEXT,
}

#: ``Flag`` values that mean "nothing to report". These are the three the shipped
#: client-side code treated as clean, kept verbatim so the glyph a reader sees does
#: not change meaning with this release.
FLAG_CLEAN_VALUES: tuple[str, ...] = ("Not flagged", "Not applicable", "")

#: Heavy check mark, shown beside a clean ``Flag``.
FLAG_OK_GLYPH = "&#10004;"

#: Heavy multiplication X, shown beside a flagged one.
FLAG_WARNING_GLYPH = "&#10006;"

#: Classes on the ``Flag`` cell's wrapper. ``FLAG_FLAGGED_CLASS`` is what the
#: "Highlight flagged values" switch styles; it is emphasis only, and no rule
#: anywhere may use it to remove a row (#242).
FLAG_FLAGGED_CLASS = "flag-flagged"
FLAG_CLEAN_CLASS = "flag-clean"

#: Inline style per ``Confidence`` value. A value not listed renders unstyled.
CONFIDENCE_STYLES: dict[str, str] = {
    "Low_Precision": "color:orange;font-weight:bold;",
    "High_Precision": "color:red;font-weight:bold;",
    "High_Precision*": "color:red;font-weight:bold;",
}

#: How each coverage field is coerced for display. The keys must be exactly
#: :data:`~vntyper.scripts.coverage_stats.COVERAGE_COLUMNS`; a contract test
#: enforces that, so adding a coverage column fails here rather than in the HTML.
COVERAGE_FIELD_TYPES: dict[str, type] = {
    "mean": float,
    "median": float,
    "stdev": float,
    "min": int,
    "max": int,
    "region_length": int,
    "uncovered_bases": int,
    "percent_uncovered": float,
    # #222's build-comparable columns. `depth_counting_policy` is a token, not a number:
    # the array sum is only comparable across builds under the policy it names.
    "vntr_array_length": int,
    "vntr_array_depth_sum": int,
    "vntr_array_depth_sum_per_unit_length": float,
    "depth_sum_reference_length": int,
    "vntr_flank_bases": int,
    "vntr_flank_mean_depth": float,
    "depth_counting_policy": str,
    # The QC verdict is a label, not a measurement: `str` keeps "PASS"/"FAIL" intact
    # and keeps it out of the two-decimal rendering the numeric fields get (#172).
    "coverage_qc": str,
}

#: Markers delimiting the three fragments spliced out of the IGV report page.
IGV_CONTAINER_MARKER = '<div id="container"'
IGV_BODY_END_MARKER = "</body>"
IGV_TABLE_JSON_MARKER = "const tableJson = "
IGV_SESSION_DICTIONARY_MARKER = "const sessionDictionary = "

#: Valid JavaScript literals for a report with no IGV panel. The template
#: interpolates the extracted fragments directly into a ``<script>`` block, so an
#: empty one is a syntax error rather than an empty table.
EMPTY_TABLE_JSON = '{"headers": [], "rows": []}'
EMPTY_SESSION_DICTIONARY = "{}"


def threshold_icon(
    value: float | None,
    cutoff: float,
    *,
    higher_better: bool = True,
    on_missing: tuple[str, str] = MISSING_AS_BLANK,
) -> tuple[str, str]:
    """Return the status glyph and colour for one metric.

    Args:
        value: The measured value, or None when the metric was not computed.
        cutoff: The threshold to compare against.
        higher_better: True when a value below ``cutoff`` is the problem; False
            when a value above it is.
        on_missing: What to return when ``value`` is None.

    Returns:
        tuple[str, str]: An HTML glyph and a CSS colour name.
    """
    if value is None:
        logger.debug("threshold_icon called with value=None; returning %s.", on_missing)
        return on_missing
    failed = value < cutoff if higher_better else value > cutoff
    logger.debug(
        "Value %s against cutoff %s (higher_better=%s): %s.",
        value,
        cutoff,
        higher_better,
        "failed" if failed else "passed",
    )
    return (WARNING_ICON, "red") if failed else (OK_ICON, "green")


def select_display_columns(df: pd.DataFrame, columns: dict[str, str]) -> pd.DataFrame:
    """Project a results frame onto its display columns and rename them.

    Columns absent from ``df`` are skipped rather than raising, because the
    Kestrel result schema differs between a positive and a negative run.

    Args:
        df: The results frame.
        columns: Source column name -> display heading, in display order.

    Returns:
        pd.DataFrame: A new frame carrying only the columns that were present,
        renamed and in the order ``columns`` declares.
    """
    present = [column for column in columns if column in df.columns]
    missing = [column for column in columns if column not in df.columns]
    if missing:
        logger.debug("Display columns absent from the results frame: %s", missing)
    return df[present].rename(columns={column: columns[column] for column in present})


def confidence_html(value: Any) -> str:
    """Colour-code one ``Confidence`` value for display.

    The label is escaped either way. This is the one cell of the Kestrel table
    that legitimately carries markup, which is why the table is rendered with
    ``escape=False``; a value with no configured style therefore used to reach
    the HTML untouched.

    Args:
        value: The confidence label.

    Returns:
        str: The escaped label, wrapped in a coloured span when it has a
        configured style.
    """
    key = value if isinstance(value, str) else str(value)
    text = escape_html(key)
    style = CONFIDENCE_STYLES.get(key)
    if style is None:
        return text
    return f'<span style="{style}">{text}</span>'


def flag_html(value: Any) -> str:
    """Render one ``Flag`` value as a glyph followed by its reason in words.

    The reason used to be visible only as a Bootstrap tooltip built in the browser:
    ``updateFlagColumn`` replaced the cell with a bare tick or cross and put the
    reason in a ``title`` attribute. That made it invisible in print, invisible to a
    screen reader once Bootstrap moves ``title`` into ``data-original-title``, and
    absent entirely when the script did not run (#242). Rendering it here means the
    cell says why with no script at all.

    The value is sample-derived -- ``vntyper report`` and ``vntyper cohort`` both
    consume a supplied ``pipeline_summary.json`` (#207) -- so it is escaped before it
    becomes markup, and the column it goes into must be named in ``html_columns``
    when the table is rendered so that :func:`escape_frame_cells` leaves it alone
    and escapes everything else.

    Args:
        value: The ``Flag`` cell value.

    Returns:
        str: A span carrying the glyph and the escaped reason.
    """
    key = value if isinstance(value, str) else str(value)
    clean = key in FLAG_CLEAN_VALUES
    glyph = FLAG_OK_GLYPH if clean else FLAG_WARNING_GLYPH
    css_class = FLAG_CLEAN_CLASS if clean else FLAG_FLAGGED_CLASS
    colour = "green" if clean else "red"
    return (
        f'<span class="{css_class}">'
        f'<span class="flag-glyph" style="color:{colour};font-weight:bold;" aria-hidden="true">{glyph}</span> '
        f'<span class="flag-reason">{escape_html(key)}</span>'
        f"</span>"
    )


def _numeric(value: Any) -> float | None:
    """Coerce one cell to a float, or None when it is not a number.

    ``numbers.Real`` rather than ``(int, float)``: pandas hands out numpy scalars,
    and ``numpy.int64`` is not a Python ``int`` on this platform. ``bool`` is
    excluded explicitly because it *is* a ``numbers.Integral`` and rendering True as
    ``1`` would be a silent type change.

    Args:
        value: The cell value.

    Returns:
        float | None: The number, or None for a missing or non-numeric value.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    if isinstance(value, numbers.Real):
        return None if pd.isna(value) else float(value)
    return None


def _as_text(value: Any) -> Any:
    """Leave a value exactly as it is."""
    return value


def _as_integer(value: Any) -> Any:
    """Render a number with no decimal part; pass anything else through."""
    number = _numeric(value)
    return value if number is None else str(int(round(number)))


def _as_two_decimals(value: Any) -> Any:
    """Render a number to two decimal places, trailing zeroes kept."""
    number = _numeric(value)
    return value if number is None else f"{number:.2f}"


def _as_six_decimals(value: Any) -> Any:
    """Render a number to six decimal places, trailing zeroes kept."""
    number = _numeric(value)
    return value if number is None else f"{number:.6f}"


def _as_three_significant(value: Any) -> Any:
    """Render a number to three significant figures, in scientific form when small."""
    number = _numeric(value)
    return value if number is None else f"{number:.3g}"


#: Format name -> the function that applies it.
CELL_FORMATTERS: dict[str, Callable[[Any], Any]] = {
    FORMAT_TEXT: _as_text,
    FORMAT_INTEGER: _as_integer,
    FORMAT_TWO_DECIMALS: _as_two_decimals,
    FORMAT_SIX_DECIMALS: _as_six_decimals,
    FORMAT_THREE_SIGNIFICANT: _as_three_significant,
}


def format_number_columns(df: pd.DataFrame, formats: Mapping[str, str]) -> pd.DataFrame:
    """Render every cell of a results frame with the formatter its column declares.

    Args:
        df: The projected display frame. Not modified.
        formats: Column -> format name. Must cover every column of ``df``.

    Returns:
        pd.DataFrame: A new frame whose formatted columns hold display strings. A
        value that is not a number - the string ``"None"`` a negative Kestrel run
        writes into every numeric column, or a missing value - passes through
        untouched.

    Raises:
        ValueError: If ``df`` carries a column ``formats`` says nothing about. This
            is the point of declaring the table per column: a new column that
            inherited a default would be rendered at whatever precision pandas chose
            for the rest of its column, which is how ``40.0`` came out as ``40.00``
            and ``1e-9`` as ``1.000000e-09`` in the same table.
    """
    undeclared = [str(column) for column in df.columns if column not in formats]
    if undeclared:
        msg = f"No display format declared for report column(s): {', '.join(sorted(undeclared))}"
        logger.error(msg)
        raise ValueError(msg)

    formatted = df.copy()
    for column, name in formats.items():
        if column in formatted.columns:
            formatted[column] = formatted[column].map(CELL_FORMATTERS[name])
    return formatted


def flagged_row_count(df: pd.DataFrame) -> int:
    """Count the rows of a results frame whose ``Flag`` states a reason.

    Args:
        df: A results frame, formatted or not, possibly without a ``Flag`` column -
            a Kestrel run with no flagging rules configured produces none
            (AGENTS.md trap 4).

    Returns:
        int: The number of flagged rows; 0 when the column is absent.
    """
    if df.empty or "Flag" not in df.columns:
        return 0
    return sum(
        1
        for value in df["Flag"]
        if not (isinstance(value, float) and pd.isna(value))
        and (value if isinstance(value, str) else str(value)) not in FLAG_CLEAN_VALUES
    )


def row_count_statement(total: int, flagged: int, *, noun: str) -> str:
    """State how many rows of a results table the reader is being shown.

    Every row the pipeline wrote is rendered, so the visible count and the total are
    the same number - and printing both is the point: a reader who has been handed a
    report with rows silently removed has no way to know, and this sentence is the
    one that says nothing was withheld. It is computed here rather than read out of
    DataTables' "Showing 1 to 3 of 3 entries" footer, which only exists when the
    reader's browser could reach a CDN (#242).

    Args:
        total: Rows in the frame.
        flagged: How many of them carry a flag.
        noun: The algorithm's name, e.g. ``Kestrel``.

    Returns:
        str: The sentence to print beside the table.
    """
    rows = "row" if total == 1 else "rows"
    count = str(flagged) if flagged else "none"
    return f"Showing {total} of {total} {noun} {rows}; {count} flagged."


def escape_html(value: Any) -> str:
    """HTML-escape one value, including quotes.

    Args:
        value: Any value; non-strings are stringified first.

    Returns:
        str: The escaped text.
    """
    return html.escape(value if isinstance(value, str) else str(value), quote=True)


def escape_frame_cells(df: pd.DataFrame, html_columns: tuple[str, ...] = ()) -> pd.DataFrame:
    """Escape every string cell of a frame that is about to be rendered raw.

    Only ``str`` cells are touched. Numbers and NaN cannot carry markup, and
    stringifying them here would take pandas' own float and NA formatting out of
    the rendered table.

    Args:
        df: The frame to escape.
        html_columns: Columns already holding markup we constructed ourselves,
            left alone.

    Returns:
        pd.DataFrame: A new frame; the input is not modified.
    """
    escaped = df.copy()
    for column in escaped.columns:
        if column in html_columns:
            continue
        escaped[column] = escaped[column].map(
            lambda cell: html.escape(cell, quote=True) if isinstance(cell, str) else cell
        )
    return escaped


def escaped_table_html(df: pd.DataFrame, classes: str, html_columns: tuple[str, ...] = ()) -> str:
    """Render a frame as an HTML table with every sample-derived cell escaped.

    ``DataFrame.to_html(escape=False)`` is needed whenever *any* column holds markup
    VNtyper built, and it disables escaping for the whole table - so the columns holding
    a sample's own strings go out as HTML too. This pairs it with
    :func:`escape_frame_cells`, which escapes everything except the columns named as
    ours, so the exemption is stated per column instead of applying to all of them.

    Args:
        df: The frame to render. Not modified.
        classes: The CSS classes for the ``<table>`` element.
        html_columns: Columns already holding markup this codebase constructed.
            Anything not named here is escaped.

    Returns:
        str: The table markup, or "" for an empty frame - ``to_html`` on one produces a
            headerless table that renders as a stray empty box.
    """
    if df.empty:
        return ""
    return escape_frame_cells(df, html_columns=html_columns).to_html(
        classes=classes,
        index=False,
        escape=False,
    )


def parse_coverage_stats(data: list[dict[str, Any]]) -> dict[str, Any]:
    """Coerce the first coverage row into the values the report renders.

    Args:
        data: The ``Coverage Calculation`` step's rows. Only the first is used.

    Returns:
        dict[str, Any]: One entry per
        :data:`~vntyper.scripts.coverage_stats.COVERAGE_COLUMNS`. A field is None when
        there is no coverage row, when the column is absent, when it carries
        :data:`~vntyper.scripts.coverage_stats.COVERAGE_NULL_TOKEN`, or when its value
        cannot be coerced.

    Note:
        Coercion is per field. It used to abort the whole row on the first failure, which
        took every later column with it - and ``coverage_qc`` is the last one, so a single
        malformed number discarded the QC verdict (#222).

        A missing column reads as ``None``, never ``0``. Zero is a measurement: it says the
        region was looked at and no reads were found. A summary written before a column
        existed, or by a run that could not compute it, has said nothing at all.
    """
    stats: dict[str, Any] = dict.fromkeys(COVERAGE_COLUMNS)
    if not data or not isinstance(data, list):
        return stats

    row = data[0]
    for name in COVERAGE_COLUMNS:
        raw = row.get(name)
        if name in _BUILD_COMPARABLE_COLUMNS:
            # #222's columns are isolated from the originals in both directions. Absent or
            # not-measured reads as None, never 0 - zero would say the region was looked at
            # and found empty - and an unreadable one is logged and skipped rather than
            # aborting the row, so a column added in 2026 cannot cost a reader the eight
            # statistics and the QC verdict that summaries have always carried.
            if raw is None or raw in (COVERAGE_NULL_TOKEN, ""):
                continue
            try:
                stats[name] = COVERAGE_FIELD_TYPES[name](raw)
            except (TypeError, ValueError) as e:
                logger.error("Error parsing VNTR coverage statistic %s=%r: %s", name, raw, e)
            continue

        # The original columns keep their pre-#222 behaviour exactly: absent coerces from 0,
        # and the first unreadable value stops the row so the fields after it stay None
        # rather than fabricating zeroes. Changing that would change what every existing
        # consumer reads out of a summary it already has.
        try:
            stats[name] = COVERAGE_FIELD_TYPES[name](raw if raw is not None else 0)
        except Exception as e:
            logger.error("Error parsing VNTR coverage statistics: %s", e)
            return stats

    logger.debug("All VNTR coverage statistics extracted successfully: %s", stats)
    return stats


@dataclass(frozen=True)
class FastpMetrics:
    """The four fastp rates the report shows, plus the sequencing setup line.

    Attributes:
        available: Whether fastp output was found at all. The report hides the
            whole section when it was not.
        duplication_rate: Read duplication rate, or None.
        q20_rate: Post-filtering Q20 rate, or None.
        q30_rate: Post-filtering Q30 rate, or None.
        passed_filter_rate: Reads passing the filter over reads before it, or
            None when no reads were seen.
        sequencing: The free-text sequencing setup, e.g. "paired end (150 cycles)".
    """

    available: bool
    duplication_rate: float | None
    q20_rate: float | None
    q30_rate: float | None
    passed_filter_rate: float | None
    sequencing: str


def summarise_fastp(fastp_data: dict[str, Any]) -> FastpMetrics:
    """Reduce parsed fastp JSON to the metrics the report shows.

    Args:
        fastp_data: The parsed ``output.json`` fastp wrote, or ``{}``.

    Returns:
        FastpMetrics: The metrics, with ``available`` False for empty input.
    """
    if not fastp_data:
        return FastpMetrics(False, None, None, None, None, "")

    summary = fastp_data.get("summary", {})
    after_filtering = summary.get("after_filtering", {})
    before_filtering = summary.get("before_filtering", {})
    filtering_result = fastp_data.get("filtering_result", {})

    total_reads_before = before_filtering.get("total_reads", 1)
    passed_filter_reads = filtering_result.get("passed_filter_reads", 0)
    if total_reads_before > 0:
        passed_filter_rate = passed_filter_reads / total_reads_before
        logger.debug("Passed filter rate calculated: %.2f", passed_filter_rate)
    else:
        passed_filter_rate = None
        logger.debug("Total reads before filtering is zero; passed filter rate set to None.")

    return FastpMetrics(
        available=True,
        duplication_rate=fastp_data.get("duplication", {}).get("rate", None),
        q20_rate=after_filtering.get("q20_rate", None),
        q30_rate=after_filtering.get("q30_rate", None),
        passed_filter_rate=passed_filter_rate,
        sequencing=summary.get("sequencing", ""),
    )


def extract_line_after(content: str, marker: str) -> str:
    """Return the text between ``marker`` and the end of its line.

    Both edge cases are load-bearing and both used to be wrong:

    * an absent marker made ``find`` return -1, and ``-1 + len(marker)`` is a
      valid index, so the old code sliced from character 17 and returned
      arbitrary page text where a JavaScript object literal belonged;
    * a marker on the last line with no trailing newline made ``find("\\n", ...)``
      return -1, and slicing to -1 dropped the final character.

    Args:
        content: The document to search.
        marker: The literal that precedes the wanted text.

    Returns:
        str: The stripped remainder of the marker's line, or "" when the marker
        is absent.
    """
    start = content.find(marker)
    if start == -1:
        logger.debug("Marker %r not found.", marker)
        return ""
    start += len(marker)
    end = content.find("\n", start)
    if end == -1:
        end = len(content)
    return content[start:end].strip()


def js_json_literal(fragment: str, fallback: str) -> str:
    """Re-serialise an extracted fragment as a literal that is safe in a ``<script>``.

    ``report_template.html`` interpolates the return value directly into a script
    block as ``const tableJson = {{ table_json|safe }};``. The fragment reaching
    here was lifted verbatim out of the igv-reports page by
    :func:`extract_line_after` and is sample-derived, so it is re-parsed and
    re-emitted rather than trusted: what the template receives is always the
    output of :func:`json.dumps`, never the extracted text.

    A single trailing ``;`` is stripped before parsing. This is defensive against
    a *future* igv-reports version, not a correction of today's: verified against
    the installed igv-reports 1.16.0, ``templates/variant_template.html:155-156``
    is ``const tableJson = "@TABLE_JSON@"`` with no terminator, and
    ``report.py:178-183`` substitutes the placeholder including its quotes, so the
    fragment reaching here never carries a trailing ``;`` today.

    ``json.dumps`` is called with ``ensure_ascii=True`` (the stdlib default, kept
    explicit here because it is load-bearing): it escapes every non-ASCII
    codepoint, which includes U+2028 and U+2029 -- line terminators to a
    JavaScript parser that are legal inside a JSON string. That leaves ``<`` as
    the one remaining script-context hazard, and **every** literal ``<`` is
    escaped rather than only the ``</`` of a closing tag. Escaping ``</`` alone
    stops a direct ``</script>`` but not the HTML5 tokenizer's *double-escaped*
    script state: ``<!--`` followed by ``<script`` inside a script element -- a
    sequence containing no ``</`` at all -- puts the parser into a state where the
    real ``</script>`` no longer terminates the element, and the remainder of the
    document is consumed as script text. Escaping ``<`` subsumes both routes.

    ``\\u003c`` is a JSON string escape, so what the browser parses back is the
    original ``<`` and the data reaching the page is unchanged. Only string values
    can contain a ``<`` -- JSON's structural characters do not -- so the
    replacement cannot touch the literal's syntax.

    Keys are sorted and separators are minimised so that two runs over the same
    IGV page emit byte-identical script.

    Args:
        fragment: The extracted literal, possibly empty, possibly not JSON.
        fallback: A valid literal to use when the fragment cannot be parsed. An
            empty fragment would otherwise produce ``const tableJson = ;`` -- a
            syntax error that takes the whole script block down, and with it the
            variant table, the flag toggles and the coverage switch.

    Returns:
        str: A JSON literal that parses as JavaScript and cannot escape the
        script block.
    """
    candidate = fragment.strip().removesuffix(";").strip()
    if not candidate:
        return fallback
    try:
        value = json.loads(candidate)
    except ValueError as e:
        logger.warning(f"IGV fragment could not be parsed as JSON and was discarded: {e}")
        return fallback

    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return encoded.replace("<", "\\u003c")


def extract_igv_fragments(content: str) -> tuple[str, str, str]:
    """Split the three fragments the report needs out of an IGV report page.

    Args:
        content: The full text of the page ``create_report`` wrote.

    Returns:
        tuple[str, str, str]: The ``#container`` markup, the ``tableJson``
        literal and the ``sessionDictionary`` literal. All three are empty when
        the container cannot be located.
    """
    igv_start = content.find(IGV_CONTAINER_MARKER)
    igv_end = content.find(IGV_BODY_END_MARKER)

    if igv_start == -1 or igv_end == -1:
        logger.error("Failed to extract IGV content from report.")
        return "", "", ""

    igv_content = content[igv_start:igv_end].strip()
    table_json = extract_line_after(content, IGV_TABLE_JSON_MARKER)
    session_dictionary = extract_line_after(content, IGV_SESSION_DICTIONARY_MARKER)

    logger.info("Successfully extracted IGV content, tableJson, and sessionDictionary.")
    return igv_content, table_json, session_dictionary
