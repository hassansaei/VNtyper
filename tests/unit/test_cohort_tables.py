"""The three tables of the cohort report: which columns, in which order, escaped how.

This is where #190 was fixed and where it can be un-fixed. The cohort report is the
artefact the web service hands back for a whole cohort to open in a browser, and every
cell of every one of these tables is a value read out of a sample's own
`pipeline_summary.json` - the sample name, the motif, the assembly, the version. Exactly
one cell in the three tables is markup VNtyper built: the `Confidence` colour span. That
one column is named in `html_columns`; everything else is escaped.

`tests/unit/test_cohort_summary_escaping.py` proves the property end to end, through a
rendered page. This file proves it at the seam, where the exemption list is written down,
and adds the column contract the escaping suite does not test: which columns appear, in
what order, and what happens to the ones that do not.

Characterisation throughout, with one exception noted in its own docstring.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from vntyper.scripts.cohort_inputs import parse_pipeline_summary
from vntyper.scripts.cohort_tables import (
    ADVNTR_DISPLAY_COLUMNS,
    KESTREL_DISPLAY_COLUMNS,
    KESTREL_HTML_COLUMNS,
    TABLE_CLASSES,
    additional_stats_frame,
    advntr_table_html,
    confidence_span,
    kestrel_table_html,
    stats_table_html,
)

pytestmark = pytest.mark.unit

INJECTION = "<script>alert(1)</script>"
ESCAPED = "&lt;script&gt;alert(1)&lt;/script&gt;"


def _probe(column: str) -> str:
    """Build the injection probe written into one column.

    One probe per column rather than one shared string, so a failure names the column
    that reached the page as markup instead of only saying that something did.

    Args:
        column: The column the probe is written into. Column names here are word
            characters only, so nothing in the name needs escaping itself.

    Returns:
        str: Markup that executes if it reaches the page unescaped.
    """
    return f"<script>{column}</script>"


def _escaped_probe(column: str) -> str:
    """The form :func:`_probe` must appear in once the table has escaped it.

    Written out rather than computed with ``html.escape``, so the assertion does not
    reduce to running the escaper under test against itself.

    Args:
        column: The column the probe was written into.

    Returns:
        str: The escaped probe.
    """
    return f"&lt;script&gt;{column}&lt;/script&gt;"


# ---------------------------------------------------------------------------
# The one constructed fragment
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label,colour",
    [
        ("Low_Precision", "orange"),
        ("High_Precision", "red"),
        ("High_Precision*", "red"),
    ],
)
def test_each_recognised_confidence_label_gets_its_colour(label: str, colour: str) -> None:
    assert confidence_span(label) == f'<span style="color:{colour};font-weight:bold;">{label}</span>'


def test_an_unrecognised_label_is_returned_as_escaped_text_with_no_span() -> None:
    """The exemption is for the span this module builds, not for the column.

    A `Confidence` value matching no colour falls through unstyled, and it is still a
    sample's own string, so it has to leave here already escaped - the renderer will not
    escape it a second time.
    """
    assert confidence_span(INJECTION) == ESCAPED


def test_quotes_are_escaped_too() -> None:
    """The value lands inside a `<td>`, but escaping quotes costs nothing and removes
    the question of whether a future caller puts it in an attribute."""
    assert confidence_span('a"b') == "a&quot;b"


@pytest.mark.parametrize("value", [5, 0.5, None, float("nan")])
def test_a_non_string_is_returned_unchanged(value: object) -> None:
    """Deliberately not stringified: pandas formats numbers and NA itself when it
    renders the table, and stringifying here would take that formatting away.

    This is the one behavioural difference from `report_formatting.confidence_html`,
    which does stringify. The two are not interchangeable.

    Identity, not equality. The assertion used to be
    `result is value or (isinstance(result, float) and isinstance(value, float))`, whose
    second disjunct is satisfied by *any* float for the two float parameters - so
    returning `0.0`, or `float(value) + 1`, passed it. `is` is also what makes the NaN
    parameter mean anything: `float("nan") == float("nan")` is False, so an equality
    assertion could not have been written for it at all.
    """
    assert confidence_span(value) is value


# ---------------------------------------------------------------------------
# The Kestrel table
# ---------------------------------------------------------------------------


def test_the_kestrel_display_columns_are_the_recorded_seventeen() -> None:
    assert KESTREL_DISPLAY_COLUMNS == (
        "Sample",
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
        "Nomenclature_Kestrel",
        "Nomenclature_adVNTR",
        "Ambiguity_Interval",
        "Repeat_Form",
        "Nomenclature_Note",
    )


def test_the_kestrel_table_renders_its_columns_in_the_declared_order() -> None:
    frame = pd.DataFrame([{"Flag": "Not flagged", "Sample": "s1", "Motif": "5"}])

    headings = _headings(kestrel_table_html(frame))

    assert headings == ["Sample", "Motif", "Flag"]


def test_a_column_the_frame_does_not_have_is_skipped_rather_than_raising() -> None:
    """A negative run's `kestrel_result.tsv` carries neither `Flag` nor the depth
    columns, so the display list is a superset of any one run's output."""
    frame = pd.DataFrame([{"Sample": "s1"}])

    assert _headings(kestrel_table_html(frame)) == ["Sample"]


def test_a_column_that_is_not_a_display_column_is_dropped() -> None:
    frame = pd.DataFrame([{"Sample": "s1", "Del": "x", "__row_result": "High_Precision"}])

    assert _headings(kestrel_table_html(frame)) == ["Sample"]


def test_the_confidence_column_becomes_a_colour_span() -> None:
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": "High_Precision"}])

    assert '<span style="color:red;font-weight:bold;">High_Precision</span>' in kestrel_table_html(frame)


def test_confidence_is_the_only_column_exempt_from_escaping() -> None:
    """The exemption list itself, asserted whole rather than by membership.

    `kestrel_table_html` ends in `to_html(escape=False)`, so any column named here
    reaches the page as markup. Widening the list is therefore a security change and has
    to fail a test.

    This test previously built a frame of seven columns and counted seven escapes, and
    that could not fail: `KESTREL_DISPLAY_COLUMNS` has twelve members, so
    `KESTREL_HTML_COLUMNS = ("Confidence", "POS")` passed it - the frame had no `POS` -
    and passed the oracle fingerprint and the whole unit tier with it. Asserting the
    tuple is what closes that; the test below is what makes the *consequence* of
    widening it visible as well.
    """
    assert KESTREL_HTML_COLUMNS == ("Confidence",)


def test_every_kestrel_display_column_but_confidence_is_escaped() -> None:
    """Every column of the table except the one constructed fragment, probed.

    The probes are derived from `KESTREL_DISPLAY_COLUMNS`, so a column added to the
    table is covered the day it is added rather than the day someone remembers to extend
    a hand-written list. The one exclusion is spelled `"Confidence"` as a literal and
    **not** taken from `KESTREL_HTML_COLUMNS`: reading the exemption list here would make
    the test follow whatever that list says, which is exactly how the previous version
    could not see `POS` being exempted. Written this way, exempting a second column fails
    this test as well as the one above - including if it is done by passing
    `html_columns` at the call site without touching the constant.
    """
    assert "Confidence" in KESTREL_DISPLAY_COLUMNS  # else the exclusion below excludes nothing
    escapable = [column for column in KESTREL_DISPLAY_COLUMNS if column != "Confidence"]
    frame = pd.DataFrame([{**{column: _probe(column) for column in escapable}, "Confidence": "High_Precision"}])

    html = kestrel_table_html(frame)

    for column in escapable:
        assert _probe(column) not in html, f"{column} reached the page as markup"
        assert _escaped_probe(column) in html, f"{column} did not reach the page at all"


def test_the_frame_handed_to_the_kestrel_table_is_not_modified() -> None:
    """The colour span is written onto a copy; the caller still exports plain text.

    Unlike the sample-level reduction, which does mutate its input, this half of the
    original code took a copy and the copy is preserved.
    """
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": "High_Precision"}])

    kestrel_table_html(frame)

    assert frame["Confidence"].tolist() == ["High_Precision"]


def test_an_empty_kestrel_frame_renders_as_nothing_at_all() -> None:
    """`to_html` on an empty frame produces a headerless table that shows as a stray
    empty box."""
    assert kestrel_table_html(pd.DataFrame()) == ""


def test_the_kestrel_table_carries_the_shared_css_classes() -> None:
    frame = pd.DataFrame([{"Sample": "s1"}])

    assert f'class="dataframe {TABLE_CLASSES}"' in kestrel_table_html(frame)


# ---------------------------------------------------------------------------
# The adVNTR table
# ---------------------------------------------------------------------------


def test_the_advntr_display_columns_are_the_recorded_sixteen() -> None:
    assert ADVNTR_DISPLAY_COLUMNS == (
        "Sample",
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
        "Nomenclature",
        "Nomenclature_Tier",
        "Nomenclature_Kestrel",
        "Nomenclature_adVNTR",
        "Ambiguity_Interval",
        "Repeat_Form",
        "Nomenclature_Note",
    )


def test_the_advntr_table_renders_its_columns_in_the_declared_order() -> None:
    frame = pd.DataFrame([{"Flag": "Not flagged", "Sample": "s1", "VID": "25561"}])

    assert _headings(advntr_table_html(frame)) == ["Sample", "VID", "Flag"]


def test_the_advntr_table_has_no_escaping_exemption_at_all() -> None:
    """Nothing constructs markup for this table, so no column is exempt - including
    `Confidence`, were one ever to appear in an adVNTR frame.

    `advntr_table_html` takes `escaped_table_html`'s default `html_columns=()`, so there
    is no exemption list to pin the way there is for Kestrel; what pins it is that
    *every* display column is probed. This used to cover three of the eleven columns and
    count three escapes, which is the same shape of hole the Kestrel test had: adding an
    exemption for `POS` or `Pvalue` passed it.
    """
    frame = pd.DataFrame([{column: _probe(column) for column in ADVNTR_DISPLAY_COLUMNS}])

    html = advntr_table_html(frame)

    for column in ADVNTR_DISPLAY_COLUMNS:
        assert _probe(column) not in html, f"{column} reached the page as markup"
        assert _escaped_probe(column) in html, f"{column} did not reach the page at all"


def test_an_empty_advntr_frame_renders_as_nothing_at_all() -> None:
    assert advntr_table_html(pd.DataFrame()) == ""


# ---------------------------------------------------------------------------
# The per-sample statistics table
# ---------------------------------------------------------------------------


def _stats_row_of_probes() -> dict[str, Any]:
    """Build one statistics row carrying a probe in every column the cohort produces.

    The keys are taken from `parse_pipeline_summary({})`, which is the function that
    decides what a statistics row contains, so a statistic added there is probed without
    this file being touched. `coverage` is the one key that is not a column: it is a
    nested mapping `additional_stats_frame` spreads into `cov_`-prefixed columns.

    Returns:
        dict[str, Any]: One row, as `load_pipeline_summary_for_sample` would return it
        with its `Sample` key added.
    """
    _, _, defaults = parse_pipeline_summary({})
    row: dict[str, Any] = {"Sample": _probe("Sample")}
    row.update({key: _probe(key) for key in defaults if key != "coverage"})
    row["coverage"] = {"mean": _probe("cov_mean")}
    return row


def test_the_statistics_table_escapes_every_column() -> None:
    """Assembly, version and pipeline are read out of each sample's summary.

    This table has no display list and no exemption list - `stats_table_html` renders
    whatever frame it is handed with `escaped_table_html`'s default `html_columns=()` -
    so the probe set is derived from the frame's own columns, and the frame is built the
    way a cohort run builds it. That is what makes the assertion total rather than a
    sample of three columns, which is what it used to be.
    """
    frame = additional_stats_frame([_stats_row_of_probes()])
    assert "cov_mean" in frame.columns  # the coverage flattening path was exercised

    html = stats_table_html(frame)

    for column in frame.columns:
        assert _probe(column) not in html, f"{column} reached the page as markup"
        assert _escaped_probe(column) in html, f"{column} did not reach the page at all"


def test_an_empty_statistics_table_renders_as_nothing_at_all() -> None:
    assert stats_table_html(pd.DataFrame()) == ""


def test_the_statistics_frame_puts_sample_first() -> None:
    stats = [{"version": "2.0.6", "Sample": "s1", "assembly": "hg38"}]

    frame = additional_stats_frame(stats)

    assert list(frame.columns) == ["Sample", "version", "assembly"]


def test_the_statistics_frame_flattens_the_coverage_dictionary() -> None:
    """Coverage arrives as a nested mapping and is spread into `cov_`-prefixed columns."""
    stats = [{"Sample": "s1", "coverage": {"mean": "31.2", "median": "30"}}]

    frame = additional_stats_frame(stats)

    assert list(frame.columns) == ["Sample", "cov_mean", "cov_median"]
    assert frame.loc[0, "cov_mean"] == "31.2"


def test_a_sample_with_no_coverage_still_gets_the_coverage_columns() -> None:
    """One sample without a `Coverage Calculation` step must not drop the columns for
    the samples that have one; the missing cells come back as NA."""
    stats = [{"Sample": "s1", "coverage": {"mean": "31.2"}}, {"Sample": "s2", "coverage": {}}]

    frame = additional_stats_frame(stats)

    assert list(frame.columns) == ["Sample", "cov_mean"]
    assert pd.isna(frame.loc[1, "cov_mean"])


def test_a_statistics_list_with_no_coverage_key_is_left_alone() -> None:
    stats = [{"Sample": "s1", "version": "2.0.6"}]

    assert list(additional_stats_frame(stats).columns) == ["Sample", "version"]


def test_an_empty_statistics_list_produces_an_empty_frame() -> None:
    assert additional_stats_frame([]).empty


def _headings(html: str) -> list[str]:
    """Pull the table headings out of rendered markup.

    Args:
        html: A rendered table.

    Returns:
        list[str]: The heading text, in document order.
    """
    import re

    return re.findall(r"<th>(.*?)</th>", html)
