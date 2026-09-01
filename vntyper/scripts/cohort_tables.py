"""
vntyper/scripts/cohort_tables.py

Module Purpose:
---------------
Build the three HTML tables of the cohort report: the Kestrel results, the adVNTR
results, and the per-sample statistics.

Each one is a column selection, an optional piece of markup this codebase constructs,
and an escaping decision. The escaping decision is the reason this module exists as a
module: the cohort report is what the web service hands back for a whole cohort to open
in a browser, and every cell of every table is a value read out of a sample's own
``pipeline_summary.json`` - the sample name (or its pseudonym), the motif, the flag, the
assembly, the version. None of that is markup VNtyper built, so none of it may reach the
page as markup (#190).

Exactly one cell in the three tables is markup VNtyper built: the ``Confidence`` span
:func:`confidence_span` produces. It is named as an ``html_columns`` exemption at the one
call site that needs it rather than by turning escaping off for a whole table, which is
what ``to_html(escape=False)`` used to do here.

Extracted from ``cohort_summary.py`` in Task 22 of the #181-#197 follow-ups.
``tests/unit/test_cohort_tables.py`` covers the seam and
``tests/unit/test_cohort_summary_escaping.py`` covers the rendered page.
"""

import logging
from typing import Any

import pandas as pd

from vntyper.scripts.cohort_profiles import PROFILE_EXPORT_COLUMNS
from vntyper.scripts.molecular_identity_presentation import IDENTITY_COLUMNS
from vntyper.scripts.report_formatting import confidence_html, escaped_table_html, is_empty_result_row

logger = logging.getLogger(__name__)

#: CSS classes every table in the cohort report carries. Named once so the three tables
#: cannot drift apart, and so the renderer call sites read as what they are.
#:
#: All three are the cohort report's own now (#242). ``table-bordered``, ``hover``,
#: ``compact``, ``order-column`` and ``table-sm`` came from Bootstrap and DataTables and
#: went with them; ``table`` is the hook the shared token layer's cell rules select on,
#: ``table-striped`` is drawn by a rule in the cohort template itself, and ``sortable``
#: is what the cohort's script looks for when it turns column headings into buttons.
#:
#: The stripe stays here and stays absent from the per-sample report, deliberately: a
#: cohort table is a triage scan across many samples where reading along a wide row is
#: the task, and its flagged rows are filtered rather than tinted, so a stripe competes
#: with nothing. In a per-sample report it competes with the flagged-value highlight,
#: which is the one row treatment there that carries meaning.
TABLE_CLASSES = "table table-striped sortable"

#: Kestrel result columns, in display order. ``Sample`` first; columns absent from the
#: frame are skipped, because a negative run's results carry neither ``Flag`` nor the
#: depth columns. Columns present in the frame but absent from this list - the working
#: columns the sample-level reduction writes, most of all - are dropped.
KESTREL_DISPLAY_COLUMNS: tuple[str, ...] = (
    "Sample",
    *PROFILE_EXPORT_COLUMNS,
    "Motif",
    "Variant",
    "POS",
    "REF",
    "ALT",
    "Motif_sequence",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "Depth_Score",
    "Confidence",
    "Flag",
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    *IDENTITY_COLUMNS,
)

#: adVNTR result columns, in display order.
ADVNTR_DISPLAY_COLUMNS: tuple[str, ...] = (
    "Sample",
    *PROFILE_EXPORT_COLUMNS,
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
    "Evidence_Disposition",
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    *IDENTITY_COLUMNS,
)

#: The one escaping exemption in the whole cohort report: the column
#: :func:`confidence_span` rewrites into a span. It is stated as a column name rather
#: than as a table-wide `escape=False`, so widening it is a visible edit.
KESTREL_HTML_COLUMNS: tuple[str, ...] = ("Confidence",)


def confidence_span(value: Any) -> Any:
    """
    Wrap a confidence label in its class, escaping the label itself.

    The class and the escaping are :func:`report_formatting.confidence_html`'s, so the
    two reports cannot label the same value differently. The colour that used to be
    inlined here is gone: ``orange`` measured 1.76:1 against this table's own striped
    rows and 1.97:1 against a plain one, and ``red`` - on ``High_Precision``, the most
    trustworthy call there is - measured 3.57:1 and 4.00:1. See
    :data:`report_formatting.CONFIDENCE_CLASSES` for why neither is replaced by another
    hue.

    Args:
        value: The ``Confidence`` cell. Any type; only ``str`` is styled.

    Returns:
        The cell as markup - a classed ``<span>`` for a recognised label, the escaped
        text for anything else, and non-strings unchanged so pandas keeps formatting
        numbers and NA itself. That last case is the one behavioural difference from
        :func:`report_formatting.confidence_html`, which stringifies.
    """
    if not isinstance(value, str):
        return value
    return confidence_html(value)


def _display_cell(value: Any) -> Any:
    """Replace a missing scalar with the blank used by cohort HTML tables.

    Args:
        value: One cell from a cohort presentation frame.

    Returns:
        The empty string for pandas missing scalars; every real value unchanged.
    """
    if pd.api.types.is_scalar(value) and bool(pd.isna(value)):
        return ""
    return value


def _normalize_display_cells(df: pd.DataFrame) -> pd.DataFrame:
    """Build a presentation copy in which missing cells are blank.

    Args:
        df: A cohort table frame. Not modified.

    Returns:
        A copy suitable for HTML rendering. Real strings, numbers and booleans retain
        their values and types.
    """
    normalized = df.copy()
    for column in normalized.columns:
        normalized[column] = normalized[column].map(_display_cell)
    return normalized


def _is_placeholder_row(row: pd.Series) -> bool:
    """Recognise Kestrel's empty-result row after the cohort adds ``Sample``.

    Args:
        row: One row from the aggregated Kestrel frame.

    Returns:
        True only for the negative empty-result placeholder.
    """
    result = {key: _display_cell(value) for key, value in row.items() if key != "Sample"}
    return is_empty_result_row(result)


def kestrel_table_html(kestrel_df: pd.DataFrame) -> str:
    """Render the cohort's Kestrel results table.

    The ``Confidence`` column is rewritten into a classed span on a **copy**, so the
    frame the caller goes on to export as CSV/TSV/JSON stays plain text. That column is
    then the table's single escaping exemption: the exemption is for the span, not for
    whatever lands in the column - an unrecognised value falls through unstyled and is
    still a sample's own string - so :func:`confidence_span` escapes the text it renders
    in every branch.

    Every other column is a sample's own string - the sample name most obviously, but
    the motif, the flag and the alleles too - so it is escaped.

    Args:
        kestrel_df (pandas.DataFrame): The aggregated Kestrel rows. Not modified.

    Returns:
        str: The table markup, or "" when there are no results.
    """
    # Create a separate copy for HTML formatting so that machine-readable outputs remain plain.
    kestrel_df_html = kestrel_df.copy()
    if not kestrel_df_html.empty and "Confidence" in kestrel_df_html.columns:
        placeholder = kestrel_df_html.apply(_is_placeholder_row, axis=1)
        if bool(placeholder.any()):
            logger.info(
                "Suppressed %d Kestrel empty-result placeholder row(s) from the cohort table.",
                int(placeholder.sum()),
            )
            kestrel_df_html = kestrel_df_html.loc[~placeholder]
    if "Confidence" in kestrel_df_html.columns:
        kestrel_df_html["Confidence"] = kestrel_df_html["Confidence"].apply(confidence_span)

    # Reorder Kestrel DataFrame columns: place Sample first then the remaining columns.
    kestrel_columns = [col for col in KESTREL_DISPLAY_COLUMNS if col in kestrel_df_html.columns]
    display_df = _normalize_display_cells(kestrel_df_html[kestrel_columns])
    return escaped_table_html(display_df, TABLE_CLASSES, html_columns=KESTREL_HTML_COLUMNS)


def advntr_table_html(advntr_df: pd.DataFrame) -> str:
    """Render the cohort's adVNTR results table.

    Nothing constructs markup for this table, so it has no exemption at all.

    Args:
        advntr_df (pandas.DataFrame): The aggregated adVNTR rows. Not modified.

    Returns:
        str: The table markup, or "" when there are no results.
    """
    # Reorder advntr DataFrame columns: ensure Sample is first.
    advntr_columns = [col for col in ADVNTR_DISPLAY_COLUMNS if col in advntr_df.columns]
    display_df = _normalize_display_cells(advntr_df[advntr_columns])
    return escaped_table_html(display_df, TABLE_CLASSES)


def stats_table_html(additional_stats_df: pd.DataFrame) -> str:
    """
    Render the per-sample statistics table for the cohort report.

    Every value in it is read out of a sample's ``pipeline_summary.json`` - the sample
    name (or its pseudonym, through the web service), the assembly, the VNtyper version,
    the pipeline description - so none of it is markup this codebase built and none of it
    is exempt from escaping.

    Args:
        additional_stats_df (pandas.DataFrame): One row per sample, ``Sample`` first.

    Returns:
        str: The table markup, or "" when there are no statistics to show.
    """
    return escaped_table_html(_normalize_display_cells(additional_stats_df), TABLE_CLASSES)


def additional_stats_frame(additional_stats_list: list[dict[str, Any]]) -> pd.DataFrame:
    """Assemble the per-sample statistics rows into the frame the table renders.

    The ``coverage`` value each sample contributes is a nested mapping, which is spread
    into one ``cov_``-prefixed column per metric. A sample whose pipeline ran no
    ``Coverage Calculation`` step contributes an empty mapping and gets NA in those
    columns rather than removing them for the whole cohort.

    Args:
        additional_stats_list: One mapping per sample, as
            :func:`~vntyper.scripts.cohort_inputs.load_pipeline_summary_for_sample`
            returns it, each already carrying its ``Sample`` key.

    Returns:
        pd.DataFrame: One row per sample, ``Sample`` first.
    """
    additional_stats_df = pd.DataFrame(additional_stats_list)
    # For coverage, flatten the dict (if available)
    if "coverage" in additional_stats_df.columns:
        coverage_df = additional_stats_df["coverage"].apply(pd.Series)
        coverage_df = coverage_df.add_prefix("cov_")
        additional_stats_df = pd.concat([additional_stats_df.drop(columns=["coverage"]), coverage_df], axis=1)
    # Reorder columns to place "Sample" first if it exists
    if "Sample" in additional_stats_df.columns:
        cols = ["Sample"] + [col for col in additional_stats_df.columns if col != "Sample"]
        additional_stats_df = additional_stats_df[cols]
    return additional_stats_df
