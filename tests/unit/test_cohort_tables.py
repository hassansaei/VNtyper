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

import pandas as pd
import pytest

from vntyper.scripts.cohort_tables import (
    ADVNTR_DISPLAY_COLUMNS,
    KESTREL_DISPLAY_COLUMNS,
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
    """
    result = confidence_span(value)

    assert result is value or (isinstance(result, float) and isinstance(value, float))


# ---------------------------------------------------------------------------
# The Kestrel table
# ---------------------------------------------------------------------------


def test_the_kestrel_display_columns_are_the_recorded_twelve() -> None:
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
    frame = pd.DataFrame(
        [
            {
                "Sample": INJECTION,
                "Motif": INJECTION,
                "Variant": INJECTION,
                "REF": INJECTION,
                "ALT": INJECTION,
                "Motif_sequence": INJECTION,
                "Flag": INJECTION,
                "Confidence": "High_Precision",
            }
        ]
    )

    html = kestrel_table_html(frame)

    assert INJECTION not in html
    assert html.count(ESCAPED) == 7


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


def test_the_advntr_display_columns_are_the_recorded_eleven() -> None:
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
    )


def test_the_advntr_table_renders_its_columns_in_the_declared_order() -> None:
    frame = pd.DataFrame([{"Flag": "Not flagged", "Sample": "s1", "VID": "25561"}])

    assert _headings(advntr_table_html(frame)) == ["Sample", "VID", "Flag"]


def test_the_advntr_table_has_no_escaping_exemption_at_all() -> None:
    """Nothing constructs markup for this table, so no column is exempt - including
    `Confidence`, were one ever to appear in an adVNTR frame."""
    frame = pd.DataFrame([{"Sample": INJECTION, "VID": INJECTION, "Flag": INJECTION}])

    html = advntr_table_html(frame)

    assert INJECTION not in html
    assert html.count(ESCAPED) == 3


def test_an_empty_advntr_frame_renders_as_nothing_at_all() -> None:
    assert advntr_table_html(pd.DataFrame()) == ""


# ---------------------------------------------------------------------------
# The per-sample statistics table
# ---------------------------------------------------------------------------


def test_the_statistics_table_escapes_every_column() -> None:
    """Assembly, version and pipeline are read out of each sample's summary."""
    frame = pd.DataFrame([{"Sample": INJECTION, "assembly": INJECTION, "version": "2.0.6"}])

    html = stats_table_html(frame)

    assert INJECTION not in html
    assert html.count(ESCAPED) == 2


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
