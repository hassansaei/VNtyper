"""The per-sample HTML report, rendered end to end from a pipeline summary.

``generate_summary_report`` was 4% covered and only reachable by running the
whole pipeline, which is why four separate defects lived in it undetected. These
tests build a ``pipeline_summary.json`` in ``tmp_path``, render the **real**
shipped template, and assert on the HTML that comes out.

Nothing here mocks the template or the report configuration: a test that asserts
on a context dict cannot see a template that ignores the value, and the template
is where two of the four defects were.

IGV generation is never triggered -- ``bed_file`` is left unset, so
``create_report`` is never invoked and the tier stays pure.
"""

import json
import logging
import re
from pathlib import Path
from unittest.mock import Mock

import pytest

import vntyper
from vntyper.cli import load_config
from vntyper.scripts import generate_report, summary_steps
from vntyper.scripts.generate_report import generate_summary_report

pytestmark = pytest.mark.unit

TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"

#: A Kestrel row as it reaches the report: the pipeline emits both the raw motif
#: pair (`Motifs`) and the annotated motif (`Motif`), and both are load-bearing
#: upstream (contract C3 / AGENTS.md trap 3).
KESTREL_ROW = {
    "Motifs": "X-5",
    "Motif": "5",
    "Variant": "Insertion",
    "POS": 67,
    "REF": "G",
    "ALT": "GG",
    "Motif_sequence": "GGCCACCACCCTG",
    "Estimated_Depth_AlternateVariant": 120,
    "Estimated_Depth_Variant_ActiveRegion": 12000,
    "Depth_Score": 0.01,
    "Confidence": "High_Precision",
    "Flag": "Not flagged",
}

COVERAGE_ROW = {
    "mean": 250.0,
    "median": 248.0,
    "stdev": 12.5,
    "min": 100,
    "max": 400,
    "region_length": 1000,
    "uncovered_bases": 5,
    "percent_uncovered": 0.5,
}


def write_summary(output_dir: Path, *steps: dict, **top_level) -> Path:
    """Write a ``pipeline_summary.json`` with the given steps.

    Args:
        output_dir: Directory to write into.
        *steps: Step mappings, already shaped.
        **top_level: Extra top-level keys, e.g. ``input_files``.

    Returns:
        Path: The file written.
    """
    payload = {"version": "9.9.9", "input_files": {}, "steps": list(steps), **top_level}
    path = output_dir / "pipeline_summary.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def tabular_step(name: str, rows: list[dict]) -> dict:
    """Shape a tsv/csv-derived step the way ``summary.py`` records it."""
    return {"step": name, "parsed_result": {"comments": [], "data": rows}}


def render(output_dir: Path, config=None, **kwargs) -> str:
    """Render the report into ``output_dir`` and return the HTML."""
    generate_summary_report(
        output_dir=str(output_dir),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=kwargs.pop("log_file", None),
        config=config if config is not None else load_config(None),
        **kwargs,
    )
    return (output_dir / "summary_report.html").read_text(encoding="utf-8")


@pytest.fixture
def positive_summary(tmp_path):
    """A Kestrel-positive, adVNTR-absent run with good coverage."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        input_files={"bam": "sample.bam"},
    )
    return tmp_path


# ---------------------------------------------------------------------------
# The report renders at all
# ---------------------------------------------------------------------------


def test_a_report_is_written(positive_summary) -> None:
    html = render(positive_summary)
    assert html.lstrip().startswith("<!DOCTYPE html>")
    assert "Summary Report" in html


def test_the_config_is_required(tmp_path) -> None:
    with pytest.raises(ValueError, match="Config dictionary must be provided"):
        generate_summary_report(
            output_dir=str(tmp_path),
            template_dir=str(TEMPLATE_DIR),
            report_file="r.html",
            log_file=None,
            config=None,
        )


def test_a_missing_pipeline_summary_still_renders(tmp_path) -> None:
    """A report with nothing in it is a usable diagnostic; a traceback is not."""
    html = render(tmp_path)
    assert "Summary Report" in html
    assert "Not calculated" in html


# ---------------------------------------------------------------------------
# Coverage - contract C1 read through to the HTML
# ---------------------------------------------------------------------------


def test_every_coverage_statistic_reaches_the_html(positive_summary) -> None:
    html = render(positive_summary)
    for value in ("250.0", "248.0", "12.5", "100", "400", "1000", "0.5"):
        assert value in html, f"coverage value {value} is missing from the report"


def test_low_coverage_is_flagged_red(tmp_path) -> None:
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [{**COVERAGE_ROW, "mean": 3.0}]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
    )
    html = render(tmp_path)
    assert "color:red" in html


def test_coverage_that_was_never_calculated_says_so(tmp_path) -> None:
    write_summary(tmp_path, tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]))
    assert "Not calculated" in render(tmp_path)


# ---------------------------------------------------------------------------
# The coverage QC verdict - #172
# ---------------------------------------------------------------------------


def _coverage_qc_cell(html: str) -> str:
    """The value rendered in the report's Coverage QC row."""
    match = re.search(r"<td>Coverage QC</td>\s*<td>(.*?)</td>", html, re.DOTALL)
    assert match, "the report has no Coverage QC row"
    return match.group(1).strip()


def test_a_well_covered_sample_reports_a_passing_coverage_qc(positive_summary) -> None:
    """#172: 250x mean and 0.5% uncovered clears both shipped thresholds.

    ``COVERAGE_ROW`` deliberately carries no ``coverage_qc`` key - it is a summary
    written before this change - so this also pins that ``vntyper report`` recomputes
    the verdict rather than reading a stored one.
    """
    assert _coverage_qc_cell(render(positive_summary)) == "PASS"


def test_a_patchy_vntr_reports_a_failing_coverage_qc(tmp_path) -> None:
    """The headline case: an acceptable mean with most of the VNTR uncovered.

    Before #172 the report showed a red icon beside the uncovered fraction and no
    verdict, and the screening sentence still read as quality-passing.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [{**COVERAGE_ROW, "percent_uncovered": 80.0}]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
    )

    assert _coverage_qc_cell(render(tmp_path)) == "FAIL"


def test_the_reported_verdict_agrees_with_the_figure_it_prints_at_the_boundary(tmp_path) -> None:
    """Adversarial review A1, end to end.

    ``pipeline_summary.json`` carries what the TSV carries: a raw mean of 99.999 was
    serialised as ``100.00``. A report evaluating that displayed figure must print
    ``PASS`` beside it - printing ``FAIL`` next to a mean of 100.00 and a threshold of
    100 is the contradiction the rounding contract exists to prevent.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [{**COVERAGE_ROW, "mean": "100.00"}]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
    )
    html = render(tmp_path)

    assert "100.0" in html, "the figure the verdict was computed from must be the one displayed"
    assert _coverage_qc_cell(html) == "PASS"


# ---------------------------------------------------------------------------
# The Kestrel table
# ---------------------------------------------------------------------------


def test_the_kestrel_table_carries_the_display_headings(positive_summary) -> None:
    html = render(positive_summary)
    for heading in ("Position", "Depth (Variant)", "Depth (Region)", "Depth Score", "Confidence", "Flag"):
        assert f"<th>{heading}</th>" in html


def test_the_motif_column_reaches_the_report(positive_summary) -> None:
    """Defect W2. The pipeline emits an annotated `Motif` **and** a raw `Motifs`
    pair; the report keyed on `Motifs` and renamed it to the heading "Motif", so
    a positive run showed the raw pair under a heading claiming to be the
    annotated motif, and a negative run -- whose placeholder row has no `Motifs`
    column at all -- lost the column entirely. `cohort_summary.py` keys on
    `Motif`, so the two reports disagreed about the same sample."""
    html = render(positive_summary)
    assert "<th>Motif</th>" in html, "the Motif column is missing from the Kestrel table"


def test_the_motif_column_shows_the_annotated_motif_not_the_raw_pair(positive_summary) -> None:
    """`Motifs` is the raw left-right pair Kestrel emits ("X-5"); `Motif` is the
    motif the variant was assigned to ("5"). Both exist because both are
    load-bearing upstream (AGENTS.md trap 3) -- the report shows the annotated
    one, which is what `cohort_summary.py` shows."""
    html = render(positive_summary)
    assert ">5</td>" in html
    assert ">X-5</td>" not in html


def test_a_negative_run_still_has_a_motif_column(tmp_path) -> None:
    """The negative placeholder row carries `Motif` and no `Motifs` at all, so
    before the fix the column was absent from every negative report too."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [{"Motif": "None", "Confidence": "Negative"}]),
    )
    assert "<th>Motif</th>" in render(tmp_path)


def test_the_kestrel_display_columns_key_on_the_annotated_motif() -> None:
    """Contract C3, pinned at the declaration: the fix is on the report side.
    `motif_processing.py` keeps both names and neither may be renamed -- the
    Kestrel flagging rules `eval` against `Motif` and the duplicate ordering
    against `Motifs`, and a missing name evaluates to False rather than raising."""
    from vntyper.scripts.report_formatting import KESTREL_DISPLAY_COLUMNS

    assert KESTREL_DISPLAY_COLUMNS["Motif"] == "Motif"
    assert "Motifs" not in KESTREL_DISPLAY_COLUMNS


def test_the_confidence_column_is_colour_coded(positive_summary) -> None:
    html = render(positive_summary)
    assert '<span style="color:red;font-weight:bold;">High_Precision</span>' in html


def test_a_negative_run_renders_its_placeholder_row(tmp_path) -> None:
    """`output_empty_result` writes a `Motif` column and no `Motifs` at all --
    the shape that made the missing-column defect invisible to a positive run."""
    negative_row = {
        "Motif": "None",
        "Variant": "None",
        "POS": "None",
        "REF": "None",
        "ALT": "None",
        "Motif_sequence": "None",
        "Estimated_Depth_AlternateVariant": "None",
        "Estimated_Depth_Variant_ActiveRegion": "None",
        "Depth_Score": "None",
        "Confidence": "Negative",
    }
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [negative_row]),
    )
    html = render(tmp_path)
    assert "Negative" in html


def test_kestrel_conversion_failure_preserves_both_frames(monkeypatch, caplog) -> None:
    """A formatting conversion failure keeps matching evidence and escaped display data."""
    monkeypatch.setattr(generate_report.pd, "to_numeric", Mock(side_effect=ValueError("bad depth")))
    caplog.set_level(logging.WARNING, logger=generate_report.logger.name)
    caplog.clear()

    display_frame, matching_frame = generate_report.build_kestrel_frames(
        [{**KESTREL_ROW, "Motif_sequence": "<untrusted>"}]
    )

    assert len(display_frame) == len(matching_frame) == 1
    assert matching_frame.loc[0, "Confidence"] == "High_Precision"
    assert matching_frame.loc[0, "Motif Sequence"] == "<untrusted>"
    assert display_frame.loc[0, "Confidence"] == '<span style="color:red;font-weight:bold;">High_Precision</span>'
    assert display_frame.loc[0, "Motif Sequence"] == "&lt;untrusted&gt;"
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [
        (logging.WARNING, "Could not convert 'Depth Score' to numeric: bad depth")
    ]


# ---------------------------------------------------------------------------
# The screening summary box - defect W3
# ---------------------------------------------------------------------------


def summary_box_classes(html: str) -> str:
    """Return the class attribute of the screening summary box."""
    import re

    match = re.search(r'<p class="(summary-box[^"]*)"', html)
    assert match, "the screening summary box is missing from the report"
    return match.group(1)


def test_a_negative_screening_is_not_styled_as_a_finding(tmp_path) -> None:
    """Defect W3. The template decided emphasis with
    `{% if 'negative' not in summary_text %}summary-positive{% endif %}`, and
    not one of the configured messages contains that word -- only the fallback
    default does. So "No variant detected by either genotyping method" was
    rendered in the same bold style as a confirmed pathogenic frameshift."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, []),
        tabular_step(summary_steps.STEP_ADVNTR, []),
    )
    html = render(tmp_path)
    assert "No variant detected by either genotyping method" in html
    assert "summary-positive" not in summary_box_classes(html)


def test_a_positive_screening_is_styled_as_a_finding(positive_summary) -> None:
    html = render(positive_summary)
    assert "summary-positive" in summary_box_classes(html)


def test_an_advntr_only_finding_is_styled_as_a_finding(tmp_path) -> None:
    """Kestrel negative, adVNTR positive: still a finding, and the configured
    message for it contains no giveaway word either way."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, []),
        tabular_step(summary_steps.STEP_ADVNTR, [{"VID": "25561", "Flag": "Not flagged"}]),
    )
    assert "summary-positive" in summary_box_classes(render(tmp_path))


def test_a_run_with_no_results_at_all_is_not_styled_as_a_finding(tmp_path) -> None:
    """No pipeline summary: Kestrel negative, adVNTR never run. Its configured
    message is "No variant detected." -- which also lacks the giveaway word, so
    an empty run rendered as a finding too."""
    assert "summary-positive" not in summary_box_classes(render(tmp_path))


def test_the_template_no_longer_decides_emphasis_from_the_message_text() -> None:
    """Pinned at the template, because the substring test is the kind of thing
    that gets reintroduced by someone reading the rendered output."""
    template = (TEMPLATE_DIR / "report_template.html").read_text(encoding="utf-8")
    assert "summary-positive" in template, "the class vanished; this assertion would be vacuous"
    assert "'negative' not in summary_text" not in template
    assert "summary_is_positive" in template


# ---------------------------------------------------------------------------
# adVNTR
# ---------------------------------------------------------------------------


def test_an_absent_advntr_step_says_it_was_not_performed(positive_summary) -> None:
    assert "adVNTR genotyping was not performed." in render(positive_summary)


def test_an_advntr_step_with_no_rows_says_nothing_was_found(tmp_path) -> None:
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        tabular_step(summary_steps.STEP_ADVNTR, []),
    )
    assert "No pathogenic variants identified by adVNTR." in render(tmp_path)


def test_advntr_rows_are_tabulated(tmp_path) -> None:
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        tabular_step(
            summary_steps.STEP_ADVNTR,
            [{"VID": "25561", "Variant": "I22_2_G_LEN1", "NumberOfSupportingReads": 9, "Flag": "Not flagged"}],
        ),
    )
    html = render(tmp_path)
    assert "25561" in html
    assert "<th>NumberOfSupportingReads</th>" in html


# ---------------------------------------------------------------------------
# Cross-match
# ---------------------------------------------------------------------------


def _cross_match_paragraph(html: str) -> str:
    """Return the `<p>` that carries the cross-match sentence, class attribute and all."""
    match = re.search(r'<p class="summary-box[^"]*">\s*(?:At least one match|No matches)[^<]*</p>', html)
    assert match, "the cross-match paragraph should be in the report"
    return match.group(0)


def test_a_cross_match_hit_is_reported(tmp_path) -> None:
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        tabular_step(summary_steps.STEP_ADVNTR, [{"VID": "25561", "Flag": "Not flagged"}]),
        tabular_step(summary_steps.STEP_CROSS_MATCH, [{"Match": "Yes"}]),
    )
    html = render(tmp_path)

    assert "At least one match was found" in html
    assert "summary-positive" in _cross_match_paragraph(html)


def test_a_cross_match_miss_is_not_styled_as_a_hit(tmp_path) -> None:
    """
    Emphasis must come from the computed state, never from searching the sentence.

    The template asked whether the message contained "match" - and *both* fixed sentences
    do, so "No matches were found between Kestrel and adVNTR results." rendered in the
    positive style. This is the identical defect already fixed for the screening summary,
    in a second place, and the previous test missed it by asserting only on the text.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        tabular_step(summary_steps.STEP_ADVNTR, [{"VID": "25561", "Flag": "Not flagged"}]),
        tabular_step(summary_steps.STEP_CROSS_MATCH, [{"Match": "No"}]),
    )
    html = render(tmp_path)

    assert "No matches were found" in html
    assert "summary-positive" not in _cross_match_paragraph(html)


@pytest.mark.parametrize(
    ("rows", "expected"),
    [
        ([{"Match": "Yes"}], True),
        ([{"Match": "No"}], False),
        ([{"Match": "No"}, {"Match": "Yes"}], True),
        ([{"Match": "no"}], False),
        ([{}], False),
        ([], False),
    ],
)
def test_the_cross_match_state_is_computed_from_the_rows(rows, expected) -> None:
    """One row matching is a hit; the sentence is a consequence of the state, not its source."""
    summary = {"steps": [tabular_step(summary_steps.STEP_CROSS_MATCH, rows)]}

    _message, is_positive = generate_report.build_cross_match_summary(summary)

    assert is_positive is expected


def test_a_missing_cross_match_step_is_neither_positive_nor_worded() -> None:
    """No step, no section - and the flag must not default to the emphasised state."""
    assert generate_report.build_cross_match_summary({"steps": []}) == ("", False)


def test_no_cross_match_step_means_no_cross_match_section(positive_summary) -> None:
    assert "Cross-Match Summary" not in render(positive_summary)


# ---------------------------------------------------------------------------
# The IGV script block - defect C6
# ---------------------------------------------------------------------------


def test_a_report_without_igv_still_declares_valid_javascript(positive_summary) -> None:
    """Defect C6. `content.find(marker) + len(marker)` is 17 for an absent
    marker, so the extraction returned arbitrary page text and the report
    contained `const tableJson = <garbage>;` -- a syntax error that takes the
    whole <script> block down, and with it the variant table, the flag toggle
    and the coverage switch. That happened on every sample with no IGV report."""
    html = render(positive_summary)
    assert 'const tableJson = {"headers": [], "rows": []};' in html
    assert "const sessionDictionary = {};" in html
    assert "const tableJson = ;" not in html


def test_the_igv_fragments_are_used_when_a_report_exists(positive_summary, monkeypatch) -> None:
    """The other side of the same fix: a real IGV page still reaches the HTML.

    #216: the fragments are re-serialised through ``js_json_literal`` rather than
    passed through verbatim, so what lands in the page is ``json.dumps``' compact,
    key-sorted output -- not the extracted text's own spacing.
    """
    from vntyper.scripts import generate_report

    monkeypatch.setattr(
        generate_report,
        "extract_igv_content",
        lambda path: ('<div id="container">panel</div>', '{"headers": ["a"], "rows": []}', '{"0": "blob:x"}'),
    )
    monkeypatch.setattr(generate_report, "run_igv_report", lambda *a, **k: None)

    bed = positive_summary / "output.bed"
    bed.write_text("chr1\t1\t2\n", encoding="utf-8")
    (positive_summary / "igv_report.html").write_text("ignored", encoding="utf-8")

    html = render(positive_summary, bed_file=str(bed))
    assert 'const tableJson = {"headers":["a"],"rows":[]};' in html
    assert 'const sessionDictionary = {"0":"blob:x"};' in html


def test_igv_extraction_failure_returns_empty_fragment(monkeypatch, caplog) -> None:
    """An unreadable optional IGV page preserves its exact three-part fallback."""
    monkeypatch.setattr("builtins.open", Mock(side_effect=OSError("unreadable IGV")))
    caplog.set_level(logging.ERROR, logger=generate_report.logger.name)
    caplog.clear()

    assert generate_report.extract_igv_content("igv_report.html") == ("", "", "")
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [
        (logging.ERROR, "Unexpected error reading IGV report: unreadable IGV")
    ]


# ---------------------------------------------------------------------------
# The BAM header block
# ---------------------------------------------------------------------------


def test_the_bam_header_warning_is_shown(tmp_path) -> None:
    """`BAM Header Parsing` records a flat object, not {"data": [...]}."""
    write_summary(
        tmp_path,
        {
            "step": summary_steps.STEP_BAM_HEADER,
            "parsed_result": {
                "warning": "Declared assembly disagrees with chr1 length",
                "alignment_pipeline": "bwa mem",
                "assembly_text": "GRCh38",
                "assembly_contig": "chr1",
            },
        },
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
    )
    html = render(tmp_path)
    assert "Declared assembly disagrees with chr1 length" in html
    assert "bwa mem" in html


# ---------------------------------------------------------------------------
# The pipeline log
# ---------------------------------------------------------------------------


def test_the_pipeline_log_is_embedded(positive_summary) -> None:
    log_file = positive_summary / "pipeline.log"
    log_file.write_text("a distinctive log line\n", encoding="utf-8")
    assert "a distinctive log line" in render(positive_summary, log_file=str(log_file))


def test_no_log_file_says_so_rather_than_failing(positive_summary) -> None:
    assert "No pipeline log file was provided." in render(positive_summary, log_file=None)


def test_fastp_failure_returns_empty_mapping(monkeypatch, caplog) -> None:
    """Unreadable optional fastp metrics are reported as an empty mapping."""
    monkeypatch.setattr(generate_report.os.path, "exists", lambda _path: True)
    monkeypatch.setattr("builtins.open", Mock(side_effect=OSError("unreadable fastp")))
    caplog.set_level(logging.ERROR, logger=generate_report.logger.name)
    caplog.clear()

    assert generate_report.load_fastp_output("output.json") == {}
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [
        (logging.ERROR, "Failed to load or parse fastp output: unreadable fastp")
    ]


def test_pipeline_log_failure_returns_failure_message(monkeypatch, caplog) -> None:
    """A log read failure differs from an absent log and remains visible to the user."""
    monkeypatch.setattr(generate_report.os.path, "exists", lambda _path: True)
    monkeypatch.setattr("builtins.open", Mock(side_effect=OSError("unreadable log")))
    caplog.set_level(logging.ERROR, logger=generate_report.logger.name)
    caplog.clear()

    assert generate_report.load_pipeline_log("pipeline.log") == "Failed to load pipeline log."
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [
        (logging.ERROR, "Failed to read pipeline log file: unreadable log")
    ]


def test_pipeline_summary_failure_returns_empty_mapping(monkeypatch, caplog) -> None:
    """An unreadable pipeline summary preserves report rendering's empty state."""
    monkeypatch.setattr(generate_report.os.path, "exists", lambda _path: True)
    monkeypatch.setattr("builtins.open", Mock(side_effect=OSError("unreadable summary")))
    caplog.set_level(logging.ERROR, logger=generate_report.logger.name)
    caplog.clear()

    assert generate_report.load_pipeline_summary("pipeline_summary.json") == {}
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [
        (logging.ERROR, "Failed to load pipeline summary: unreadable summary")
    ]


# ---------------------------------------------------------------------------
# Escaping - defect C5
# ---------------------------------------------------------------------------

PAYLOAD = "<script>alert(1)</script>"
ESCAPED = "&lt;script&gt;alert(1)&lt;/script&gt;"


def test_an_input_file_name_is_escaped(tmp_path) -> None:
    """The report is a file people forward. Nothing that reaches it from a
    sample -- a file name, a BAM header, a motif sequence, a log line -- was
    escaped, and the Kestrel table was rendered with `escape=False`."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        input_files={"bam": f"{PAYLOAD}.bam"},
    )
    html = render(tmp_path)
    assert PAYLOAD not in html
    assert ESCAPED in html


def test_a_bam_header_field_is_escaped(tmp_path) -> None:
    """The header block is parsed straight out of the BAM's @PG and @SQ lines."""
    write_summary(
        tmp_path,
        {"step": summary_steps.STEP_BAM_HEADER, "parsed_result": {"alignment_pipeline": PAYLOAD, "warning": PAYLOAD}},
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
    )
    html = render(tmp_path)
    assert PAYLOAD not in html


def test_a_kestrel_cell_is_escaped(tmp_path) -> None:
    """The Kestrel table is rendered with `escape=False` so the colour-coded
    Confidence span survives; every other cell must be escaped by hand."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, "Motif_sequence": PAYLOAD}]),
    )
    html = render(tmp_path)
    assert PAYLOAD not in html
    assert ESCAPED in html


def test_a_malicious_flag_in_a_stored_summary_is_server_escaped(tmp_path) -> None:
    """#207's real trust boundary: the stored artefact, not the rule engine.

    On a normal pipeline run `Flag` is safe by construction -- `flagging.py` sets
    it from the *keys* of `flagging_rules`, config-declared identifiers, never
    from a DataFrame value. The untrusted path is the one #207 names: `Flag`
    values come from a `pipeline_summary.json` that ``vntyper report`` and
    ``vntyper cohort`` both consume as a supplied artefact -- a client-uploaded
    ZIP, or an output directory received from elsewhere -- so a hand-edited or
    otherwise adversarial summary is the thing this test plants.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, "Flag": PAYLOAD}]),
    )
    html = render(tmp_path)
    assert PAYLOAD not in html
    assert ESCAPED in html


def test_an_unstyled_confidence_value_is_escaped(tmp_path) -> None:
    """The Confidence column is the one cell that legitimately carries markup.
    A value with no configured style used to pass through untouched."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, "Confidence": PAYLOAD}]),
    )
    assert PAYLOAD not in render(tmp_path)


def test_the_pipeline_log_is_escaped(positive_summary) -> None:
    """`<pre>` does not stop a `</pre><script>` in a log line."""
    log_file = positive_summary / "pipeline.log"
    log_file.write_text(f"</pre>{PAYLOAD}\n", encoding="utf-8")
    html = render(positive_summary, log_file=str(log_file))
    assert PAYLOAD not in html
    assert "&lt;/pre&gt;" in html


def test_escaping_does_not_neuter_the_status_icons(positive_summary) -> None:
    """The icons are pre-built HTML fragments we construct ourselves. Turning
    autoescaping on without marking them would print the span markup as text."""
    html = render(positive_summary)
    assert '<span style="color:green;font-weight:bold;">&#10004;</span>' in html
    assert "&lt;span style=" not in html


def test_escaping_does_not_neuter_the_screening_message(positive_summary) -> None:
    """Configured messages carry `<br>` line breaks."""
    html = render(positive_summary)
    assert "<br>" in html


def test_escaping_does_not_neuter_the_results_tables(positive_summary) -> None:
    html = render(positive_summary)
    assert 'id="kestrel_table"' in html
    assert '<span style="color:red;font-weight:bold;">High_Precision</span>' in html


def test_escaping_does_not_neuter_the_igv_script_block(positive_summary) -> None:
    html = render(positive_summary)
    assert 'const tableJson = {"headers": [], "rows": []};' in html


#: Every context value the template is allowed to interpolate unescaped, and why.
SAFE_BY_DESIGN = {
    "kestrel_highlight": "pandas table; its cells are escaped by escape_frame_cells",
    "advntr_highlight": "pandas table rendered with escape=True, or a fixed <p>",
    "summary_text": "a configured clinical message carrying <br> line breaks",
    "cross_match_message": "one of two fixed sentences built in generate_report",
    "table_json": "a JavaScript literal spliced out of the IGV report",
    "session_dictionary": "a JavaScript literal spliced out of the IGV report",
    "mean_vntr_coverage_icon": "an HTML fragment built by report_formatting",
    "percent_vntr_uncovered_icon": "an HTML fragment built by report_formatting",
    "duplication_rate_icon": "an HTML fragment built by report_formatting",
    "q20_icon": "an HTML fragment built by report_formatting",
    "q30_icon": "an HTML fragment built by report_formatting",
    "passed_filter_icon": "an HTML fragment built by report_formatting",
}


def test_nothing_new_is_marked_safe_in_the_template() -> None:
    """Autoescaping only protects what is not marked `|safe`, and marking a value
    safe is a one-word edit. This is the list, and adding to it is deliberate."""
    import re

    template = (TEMPLATE_DIR / "report_template.html").read_text(encoding="utf-8")
    marked = set(re.findall(r"\{\{\s*([A-Za-z_][\w.]*)\s*\|\s*safe\s*\}\}", template))
    assert marked, "found no |safe expressions; this assertion would be vacuous"
    assert marked <= set(SAFE_BY_DESIGN), f"newly unescaped values: {sorted(marked - set(SAFE_BY_DESIGN))}"


def test_the_template_environment_autoescapes(positive_summary, monkeypatch) -> None:
    """The environment production *builds*, not the source that builds it.

    This used to assert that ``generate_summary_report``'s source contained the substring
    ``autoescape=`` -- which ``autoescape=False`` satisfies, so the assertion could not
    fail for the reason it existed. Capture the real ``Environment`` instead and ask it
    what it does with an HTML template name.
    """
    from jinja2 import Environment

    from vntyper.scripts import generate_report

    captured: list[Environment] = []
    real_environment = generate_report.Environment

    def spy(*args, **kwargs):
        env = real_environment(*args, **kwargs)
        captured.append(env)
        return env

    monkeypatch.setattr(generate_report, "Environment", spy)
    render(positive_summary)

    assert captured, "generate_summary_report built no Jinja2 environment; this test would be vacuous"
    env = captured[0]
    # `select_autoescape` returns a per-template-name policy; jinja2 types the attribute
    # as `bool`, so widen it before asking whether it is the callable form.
    policy: object = env.autoescape
    assert policy, "autoescaping is off; every |safe fragment is now the only thing escaped"
    assert callable(policy), "autoescape is not a select_autoescape policy"
    assert policy("report_template.html") is True, "an .html template must be autoescaped"
    assert env.from_string("{{ v }}").render(v=PAYLOAD) == ESCAPED


# ---------------------------------------------------------------------------
# The version and input files come from the summary, not from the caller
# ---------------------------------------------------------------------------


def test_the_version_and_input_files_come_from_the_summary(positive_summary) -> None:
    """`cli_report` used to pass both in; the generator reads them itself, which
    is why passing them was both wrong and unnecessary."""
    html = render(positive_summary)
    assert "9.9.9" in html
    assert "sample.bam" in html


# ---------------------------------------------------------------------------
# The effective BWA reference selection reaches the report (#163)
# ---------------------------------------------------------------------------


def _labeled_value(html: str, label: str) -> str:
    """The text rendered after ``<strong>{label}:</strong>`` in the report."""
    match = re.search(rf"<strong>{re.escape(label)}:</strong>\s*([^<]*)", html)
    assert match, f"label {label!r} not found in report"
    return match.group(1).strip()


def test_the_effective_bwa_reference_reaches_the_report_when_a_fallback_was_taken(tmp_path) -> None:
    """`select_bwa_reference` (Task 7, #163) can substitute a UCSC-family reference for
    a requested NCBI/ENSEMBL one; the report is the operator's only later evidence that
    happened, so all four fields `pipeline.py` records must actually be rendered."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        reference_assembly_requested="hg38_ensembl",
        reference_key_used="bwa_reference_hg38",
        reference_path="/refs/hg38.fa",
        reference_source_effective="ucsc",
    )

    html = render(tmp_path)

    assert _labeled_value(html, "Reference assembly requested") == "hg38_ensembl"
    assert _labeled_value(html, "Reference key used") == "bwa_reference_hg38"
    assert _labeled_value(html, "Reference path") == "/refs/hg38.fa"
    assert _labeled_value(html, "Effective reference source") == "ucsc"


def test_the_requested_and_effective_reference_sources_are_shown_as_different(tmp_path) -> None:
    """A fallback must read as a substitution, not as a match.

    `reference_assembly_requested` names the label actually asked for (`hg38_ensembl`,
    whose own source is `ensembl`); when the run fell back, `reference_source_effective`
    is `ucsc` instead. Both must be visible and distinct so an operator cannot read the
    requested label as proof of what was actually used.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        reference_assembly_requested="hg38_ensembl",
        reference_key_used="bwa_reference_hg38",
        reference_path="/refs/hg38.fa",
        reference_source_effective="ucsc",
    )

    html = render(tmp_path)

    requested = _labeled_value(html, "Reference assembly requested")
    effective = _labeled_value(html, "Effective reference source")
    assert requested != effective
    assert requested == "hg38_ensembl"
    assert effective == "ucsc"


def test_an_older_summary_without_reference_selection_fields_still_renders(positive_summary) -> None:
    """A summary written before this change (or one with no BWA reference at all, e.g.
    a BAM-only run) simply omits the four fields; the section must not appear rather
    than rendering empty labels."""
    html = render(positive_summary)

    assert "Reference assembly requested" not in html
    assert "Reference key used" not in html
    assert "Reference path" not in html
    assert "Effective reference source" not in html


def test_an_unrecognised_pipeline_version_leaves_the_coverage_figures_alone() -> None:
    """The legacy correction requires positive identification, and defaults to inaction.

    A missing `coverage_qc` column alone is not proof of a pre-2.0.8 summary - a
    hand-built or third-party one may simply omit it - and silently rescaling its mean
    would be worse than the problem being fixed.
    """
    from vntyper.scripts.generate_report import _predates_region_wide_mean

    assert _predates_region_wide_mean("2.0.7") is True
    assert _predates_region_wide_mean("v2.0.6") is True
    assert _predates_region_wide_mean("2.0.8") is False
    assert _predates_region_wide_mean("2.1.0") is False
    for unknown in (None, "", "garbage", "not.a.version"):
        assert _predates_region_wide_mean(unknown) is False, f"{unknown!r} must not trigger a rescale"


def test_a_pre_2_0_8_summary_is_judged_on_its_corrected_mean(tmp_path) -> None:
    """The whole legacy path, end to end through the rendered report.

    A summary recorded by 2.0.7 carries a mean over *covered* positions. Stored 150.0 at
    40% uncovered stands for a region-wide 90.0, which fails the 100x threshold - judging
    the stored figure would pass it and print a QC verdict the data does not support.
    """
    legacy_row = {**COVERAGE_ROW, "mean": 150.0, "percent_uncovered": 40.0}
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [legacy_row]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        version="2.0.7",
    )

    html = render(tmp_path)

    assert _coverage_qc_cell(html) == "FAIL", "a pre-2.0.8 mean must be corrected before it is judged"


def test_the_same_summary_recorded_by_this_version_is_taken_at_face_value(tmp_path) -> None:
    """The correction must fire only for versions that actually predate the change."""
    current_row = {**COVERAGE_ROW, "mean": 150.0, "percent_uncovered": 40.0}
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [current_row]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        version="2.0.8",
    )

    html = render(tmp_path)

    assert _coverage_qc_cell(html) == "PASS", "a current summary's mean is already region-wide"
