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

Three groups live here:

* **Status icons.** A value, a cutoff and a direction produce a coloured glyph
  and a colour name. The report shows eight of them and they all shared one
  ``higher_better`` flag plus two hand-rolled copies with a different treatment
  of a missing value; :func:`threshold_icon` takes that treatment as an argument
  instead.
* **Table formatting.** Column selection and renaming for the Kestrel table, and
  the colour-coding of the ``Confidence`` column.
* **IGV fragment extraction.** ``create_report`` writes a standalone HTML page;
  the report splices three pieces out of it. The splicing is string work and is
  pure, so it lives here; opening the file does not.

Functions:
    threshold_icon: Value + cutoff to (icon, colour)
    select_display_columns: Project and rename a results frame for display
    confidence_html: Colour-code one ``Confidence`` value
    extract_line_after: Text following a marker, up to the next newline
    extract_igv_fragments: IGV page text to (container, tableJson, sessionDictionary)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import pandas as pd

# Contract C1: the coverage TSV field names are declared once, by the module that
# produces them. Re-typing the strings here is how the report silently loses
# coverage when a column is renamed - `.get(name, 0)` raises nothing.
from vntyper.scripts.coverage_stats import COVERAGE_COLUMNS

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

    Args:
        value: The confidence label.

    Returns:
        str: The label wrapped in a coloured span, or the label unchanged when
        it has no configured style.
    """
    style = CONFIDENCE_STYLES.get(value if isinstance(value, str) else str(value))
    if style is None:
        return value
    return f'<span style="{style}">{value}</span>'


def parse_coverage_stats(data: list[dict[str, Any]]) -> dict[str, Any]:
    """Coerce the first coverage row into the values the report renders.

    Args:
        data: The ``Coverage Calculation`` step's rows. Only the first is used.

    Returns:
        dict[str, Any]: One entry per
        :data:`~vntyper.scripts.coverage_stats.COVERAGE_COLUMNS`. A field is None
        when there is no coverage row, and fields after a coercion failure stay
        None so a partially unreadable row does not fabricate zeroes.
    """
    stats: dict[str, Any] = dict.fromkeys(COVERAGE_COLUMNS)
    if not data or not isinstance(data, list):
        return stats

    row = data[0]
    try:
        for name in COVERAGE_COLUMNS:
            stats[name] = COVERAGE_FIELD_TYPES[name](row.get(name, 0))
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


def js_object_literal(fragment: str, fallback: str) -> str:
    """Return ``fragment`` if it can stand as a JavaScript literal, else ``fallback``.

    The template interpolates these straight into a ``<script>`` block as
    ``const tableJson = {{ table_json|safe }};``. An empty fragment produces
    ``const tableJson = ;`` -- a syntax error that takes the whole script block
    down, and with it the variant table, the flag toggles and the coverage
    switch, on every sample with no IGV report.

    Args:
        fragment: The extracted literal, possibly empty.
        fallback: A valid literal to use instead.

    Returns:
        str: Something that parses.
    """
    return fragment if fragment.strip() else fallback


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
