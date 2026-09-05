"""
The cohort report must not render sample-derived values as markup.

The per-sample report was brought to autoescaping with per-cell escaping; the cohort
report was not, and it is the one the web service hands back for a whole cohort to open
in a browser. Three ``to_html(escape=False)`` calls plus ``| safe`` in the template meant
every cell of every table - sample name, motif, flag, assembly, version - reached the
page as HTML.

A sample name is caller-chosen. In the CLI it is the directory name; through the web
service it is the pseudonym recorded against an uploaded job. Neither is markup VNtyper
built, so neither may be rendered as markup.

What *is* markup VNtyper built, and stays that way, is listed in one place here so the
list cannot quietly grow: the `Confidence` column's colour span, the two Plotly figures,
and the pandas table scaffolding itself.
"""

import html as html_module
import re
from pathlib import Path

import pandas as pd
import pytest

from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.cohort_tables import ADVNTR_DISPLAY_COLUMNS, KESTREL_DISPLAY_COLUMNS
from vntyper.scripts.nomenclature_presentation import KESTREL_BAM_SEMANTICS

pytestmark = pytest.mark.unit

#: The canonical probe. If it appears unescaped anywhere in the page, it executes.
INJECTION = "<script>alert(1)</script>"

#: The escaped form the page must carry instead.
ESCAPED = "&lt;script&gt;alert(1)&lt;/script&gt;"


def _normalized_visible_text(document: str) -> str:
    """Return whitespace-normalized text from rendered cohort markup."""
    return " ".join(html_module.unescape(re.sub(r"<[^>]+>", " ", document)).split())


def _render(
    tmp_path: Path,
    kestrel_df=None,
    advntr_df=None,
    additional_stats_html="",
    advntr_evidence_provenance=None,
) -> str:
    """Render the cohort report into ``tmp_path`` and return the HTML."""
    cohort_summary.generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=pd.DataFrame() if kestrel_df is None else kestrel_df,
        advntr_df=pd.DataFrame() if advntr_df is None else advntr_df,
        summary_file="cohort_summary.html",
        config=load_config(None),
        additional_stats_html=additional_stats_html,
        advntr_evidence_provenance=advntr_evidence_provenance,
    )
    return (tmp_path / "cohort_summary.html").read_text(encoding="utf-8")


def test_recorded_advntr_evidence_provenance_is_visible_and_escaped(tmp_path: Path) -> None:
    from vntyper.modules.advntr.artifact_evidence import ASSERTION

    digest = "8bb68bd5fba539feee6feb240f113aaa24fc65b5a1e55776c58cea83db5654b0"

    document = _render(
        tmp_path,
        advntr_evidence_provenance=[
            {"sample": INJECTION, "revision": digest, "assertion": ASSERTION},
            {"sample": "hostile", "revision": digest, "assertion": INJECTION},
            {"sample": "legacy", "revision": "artifact-evidence revision not recorded", "assertion": ""},
        ],
    )

    assert digest in document
    assert ASSERTION in document
    assert "artifact-evidence revision not recorded" in document
    assert INJECTION not in document
    assert document.count(ESCAPED) == 2


def test_a_sample_name_is_not_rendered_as_markup_in_the_kestrel_table(tmp_path) -> None:
    """The Kestrel table's `Sample` column comes from the sample directory or a pseudonym."""
    frame = pd.DataFrame([{"Sample": INJECTION, "Motif": "5", "Confidence": "High_Precision"}])

    html = _render(tmp_path, kestrel_df=frame)

    assert INJECTION not in html
    assert ESCAPED in html


def test_a_sample_name_is_not_rendered_as_markup_in_the_advntr_table(tmp_path) -> None:
    """The adVNTR table never had any markup of its own, so nothing in it is exempt."""
    frame = pd.DataFrame([{"Sample": INJECTION, "VID": "25561", "Flag": "Not flagged"}])

    html = _render(tmp_path, advntr_df=frame)

    assert INJECTION not in html
    assert ESCAPED in html


@pytest.mark.parametrize("column", [c for c in KESTREL_DISPLAY_COLUMNS if c != "Confidence"])
def test_no_kestrel_column_is_rendered_as_markup(tmp_path, column) -> None:
    """Every column of the table is sample-derived, not only the sample name.

    The parametrization is derived from `KESTREL_DISPLAY_COLUMNS` rather than listed by
    hand. The hand-written list it replaces named six of the twelve columns and omitted
    `POS` and the three depth columns, so exempting any of those from escaping did not
    fail here either - the same hole `test_cohort_tables.py` had. `Confidence` is
    excluded as a literal, not by reading `KESTREL_HTML_COLUMNS`, so widening that
    constant makes this parametrization grow rather than shrink.
    """
    frame = pd.DataFrame([{"Sample": "s1", column: INJECTION}])

    html = _render(tmp_path, kestrel_df=frame)

    assert INJECTION not in html
    assert ESCAPED in html


@pytest.mark.parametrize("column", ADVNTR_DISPLAY_COLUMNS)
def test_no_advntr_column_is_rendered_as_markup(tmp_path, column) -> None:
    """The adVNTR table has no exemption at all, so every one of its columns is probed."""
    frame = pd.DataFrame([{"Sample": "s1", "VID": "25561", column: INJECTION}])

    html = _render(tmp_path, advntr_df=frame)

    assert INJECTION not in html
    assert ESCAPED in html


def test_the_additional_statistics_table_is_escaped(tmp_path) -> None:
    """Assembly, version and pipeline strings are read out of each sample's summary."""
    stats = pd.DataFrame([{"Sample": INJECTION, "assembly": INJECTION, "version": "2.0.5"}])

    html = _render(tmp_path, additional_stats_html=cohort_summary.stats_table_html(stats))

    assert INJECTION not in html
    assert ESCAPED in html


@pytest.mark.parametrize("token", ["thin-haplotype-record-support", "low-haplotype-record-support"])
def test_each_kestrel_bam_token_shows_the_cohort_semantics_note(tmp_path, token) -> None:
    """Either BAM-specific support flag is sufficient to explain the evidence unit."""
    frame = pd.DataFrame([{"Sample": "s1", "Nomenclature_Tier": "B", "Nomenclature_Flags": token}])

    html = _render(tmp_path, kestrel_df=frame)
    visible_text = _normalized_visible_text(html)

    assert visible_text.count(KESTREL_BAM_SEMANTICS) == 1
    assert f"<dt>{token}" in html


def test_cohort_semantics_note_is_absent_without_a_kestrel_bam_token(tmp_path) -> None:
    """A general nomenclature key must not imply that BAM rescue supplied evidence."""
    frame = pd.DataFrame([{"Sample": "s1", "Nomenclature_Tier": "B", "Nomenclature_Flags": "known-variant"}])

    html = _render(tmp_path, kestrel_df=frame)

    assert "<h2>Reading key</h2>" in html
    assert KESTREL_BAM_SEMANTICS not in _normalized_visible_text(html)


def test_cohort_reading_key_uses_the_required_report_anchor(tmp_path) -> None:
    """The key follows the complete missing-adVNTR block and precedes statistics."""
    frame = pd.DataFrame(
        [
            {
                "Sample": "s1",
                "Nomenclature_Tier": "B",
                "Nomenclature_Flags": "low-haplotype-record-support",
            }
        ]
    )
    stats = cohort_summary.stats_table_html(pd.DataFrame([{"Sample": "s1", "version": "2.0.5"}]))

    html = _render(tmp_path, kestrel_df=frame, additional_stats_html=stats)

    assert html.index("contributed no adVNTR result") < html.index("<h2>Reading key</h2>")
    assert html.index("<h2>Reading key</h2>") < html.index("<h2>Additional Statistics</h2>")


def test_hostile_nomenclature_legend_values_are_autoescaped(tmp_path, monkeypatch) -> None:
    """The context seam passes legend strings to Jinja autoescaping, never through ``safe``."""
    legend = [
        {"term": INJECTION, "label": INJECTION, "meaning": INJECTION},
        {
            "term": "low-haplotype-record-support",
            "label": "",
            "meaning": "Resolved haplotype-record support is low.",
        },
    ]
    monkeypatch.setattr(cohort_summary, "nomenclature_legend", lambda *_frames: legend)

    html = _render(tmp_path)

    assert INJECTION not in html
    assert html.count(ESCAPED) == 3
    assert _normalized_visible_text(html).count(KESTREL_BAM_SEMANTICS) == 1


# ---------------------------------------------------------------------------
# What must keep working
# ---------------------------------------------------------------------------


def test_the_confidence_span_is_still_markup(tmp_path) -> None:
    """
    The one genuinely constructed fragment in these tables, preserved deliberately.

    `Confidence` is rewritten to a classed `<span>` before rendering. It is built here
    from a value the config supplies, so it is the single exemption, and it is passed as
    an explicit `html_columns` entry rather than by turning escaping off wholesale.

    The span carries a class and no colour since #242's contrast pass: `color:red` on
    `High_Precision` measured 3.57:1 against this table's striped rows. See
    `report_formatting.CONFIDENCE_CLASSES`.
    """
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": "High_Precision"}])

    html = _render(tmp_path, kestrel_df=frame)

    assert '<span class="confidence confidence-high-precision">High_Precision</span>' in html


def test_a_low_precision_call_keeps_its_own_class(tmp_path) -> None:
    """Both branches of the rule render, not just the one the first test hits."""
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": "Low_Precision"}])

    html = _render(tmp_path, kestrel_df=frame)

    assert '<span class="confidence confidence-low-precision">Low_Precision</span>' in html


def test_a_confidence_value_that_is_not_ours_is_still_escaped(tmp_path) -> None:
    """
    The exemption is for the span this module builds, not for whatever lands in the column.

    An unrecognised `Confidence` passes through the colour rule unchanged, so it reaches
    the page as the cell's own text and has to be escaped like any other value.
    """
    frame = pd.DataFrame([{"Sample": "s1", "Confidence": INJECTION}])

    html = _render(tmp_path, kestrel_df=frame)

    assert INJECTION not in html


def test_the_tables_still_render_as_tables(tmp_path) -> None:
    """
    Autoescaping applies to the *values*, never to the scaffolding.

    pandas builds the table markup and the plots are Plotly's own HTML; escaping those
    would replace the report with a page of visible angle brackets.
    """
    frame = pd.DataFrame([{"Sample": "s1", "Motif": "5", "Confidence": "High_Precision"}])

    html = _render(tmp_path, kestrel_df=frame)

    assert "<table" in html
    assert "<th>Sample</th>" in html
    assert ">s1</td>" in html
    assert html.lstrip().startswith("<!DOCTYPE html>")


def test_an_empty_cohort_still_renders(tmp_path) -> None:
    """The no-results page is a usable answer; a traceback is not."""
    html = _render(tmp_path)

    assert html.lstrip().startswith("<!DOCTYPE html>")


def test_the_call_frequency_table_is_escaped(tmp_path) -> None:
    """Call frequency table in the rendered cohort report escapes malicious input."""
    frame = pd.DataFrame(
        [
            {
                "Sample": INJECTION,
                "Molecular_Identity": INJECTION,
                "Molecular_Identity_Status": "unique",
                "Motifs": "5-5",
                "POS": "100",
                "REF": "A",
                "ALT": INJECTION,
                "Variant": INJECTION,
                "Depth_Score": "15.0",
            }
        ]
    )
    html = _render(tmp_path, kestrel_df=frame)

    assert INJECTION not in html
    assert ESCAPED in html

