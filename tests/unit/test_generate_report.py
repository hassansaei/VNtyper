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
from datetime import datetime, timezone
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
    assert "<h1>MUC1 VNTR report — sample</h1>" in html


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
    assert "<h1>MUC1 VNTR report — unnamed sample</h1>" in html
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


def test_a_truncated_placeholder_row_is_still_not_a_result(tmp_path) -> None:
    """A placeholder need not carry all ten columns to be a non-result.

    This case used to pin that the negative placeholder kept its `Motif` column
    (contract C3): the row carries `Motif` and no `Motifs` at all, which is the
    shape that made the missing-column defect invisible to a positive run. The
    column contract is now pinned where it is declared
    (`test_the_kestrel_display_columns_key_on_the_annotated_motif`) and on a
    positive run (`test_the_motif_column_reaches_the_report`), because a negative
    run no longer renders a row at all - so what is left to check here is that a
    partial placeholder is recognised as one rather than tabulated.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [{"Motif": "None", "Confidence": "Negative"}]),
    )

    html = render(tmp_path)

    assert 'id="kestrel_table"' not in html, "a placeholder row was tabulated as a result"
    assert "No variant detected by Kestrel" in visible_text(html)


def test_the_kestrel_display_columns_key_on_the_annotated_motif() -> None:
    """Contract C3, pinned at the declaration: the fix is on the report side.
    `motif_processing.py` keeps both names and neither may be renamed -- the
    Kestrel flagging rules `eval` against `Motif` and the duplicate ordering
    against `Motifs`, and a missing name evaluates to False rather than raising."""
    from vntyper.scripts.report_formatting import KESTREL_DISPLAY_COLUMNS

    assert KESTREL_DISPLAY_COLUMNS["Motif"] == "Motif"
    assert "Motifs" not in KESTREL_DISPLAY_COLUMNS


def test_the_confidence_column_is_colour_coded(positive_summary) -> None:
    """The literal carries an underline since #242's accessibility pass.

    Red and orange are the same grey in print, in greyscale and to a reader with a
    red-green deficiency, and both confidence values were bold - so the hue was the
    whole of the difference between a high-precision and a low-precision call. The
    underline style is the non-colour half of that distinction; the hues themselves
    are unchanged. See `tests/unit/test_report_presentation.py`.
    """
    html = render(positive_summary)
    assert '<span style="color:red;font-weight:bold;text-decoration:underline solid;">High_Precision</span>' in html


# ---------------------------------------------------------------------------
# The Kestrel empty states - issue #242
# ---------------------------------------------------------------------------

#: What `kestrel_genotyping.output_empty_result` writes when Kestrel ran and called
#: nothing: the literal string "None" in all nine value columns and `Negative` in
#: `Confidence`. `summary.parse_tsv` splits the file on tabs and coerces nothing
#: (AGENTS.md trap 5: all `parsed_result` values are strings), so this is exactly the
#: row that reaches the report.
NEGATIVE_KESTREL_ROW = {
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

#: Script and style bodies, which are not text the reader sees.
_SCRIPT_OR_STYLE = re.compile(r"<(script|style)\b[^>]*>.*?</\1>", re.DOTALL | re.IGNORECASE)

#: Any tag.
_TAG = re.compile(r"<[^>]+>")


def visible_text(html: str) -> str:
    """Return what the report says, with its markup removed.

    The defect this file's empty-state cases describe is a *reading*: the cells
    ``None None None None None None None None NaN Negative`` under a heading
    claiming to list identified variants. Asserting that against the raw HTML is
    vacuous - every cell is separated by tags - so the tags come out first.

    Args:
        html: The rendered report.

    Returns:
        str: The document's text, whitespace collapsed.
    """
    return " ".join(_TAG.sub(" ", _SCRIPT_OR_STYLE.sub(" ", html)).split())


def negative_summary(output_dir: Path) -> Path:
    """Write a summary whose Kestrel step carries only the empty-result placeholder."""
    write_summary(
        output_dir,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [NEGATIVE_KESTREL_ROW]),
        input_files={"bam": "NEG1.bam"},
    )
    return output_dir


def test_a_negative_kestrel_run_renders_a_sentence_not_a_none_row(tmp_path) -> None:
    """The commonest report in any cohort, and it read as a crashed pipeline.

    ``output_empty_result`` writes the literal string ``"None"`` into all nine value
    columns, and ``build_kestrel_frames`` coerced ``Depth_Score`` through
    ``pd.to_numeric(errors="coerce")`` - so a negative sample's report tabulated
    ``None None None None None None None None NaN Negative`` under the heading
    "Kestrel Identified Variants".
    """
    text = visible_text(render(negative_summary(tmp_path)))

    assert "None None" not in text, "the empty-result placeholder is still tabulated as a variant"
    assert "NaN" not in text, "a coerced placeholder depth score is still displayed"
    assert "No variant detected by Kestrel" in text


def test_a_negative_kestrel_run_makes_no_row_count_claim(tmp_path) -> None:
    """The count line says nothing was withheld; counting a non-result defeats it.

    Rendered over the placeholder it read "Showing 1 of 1 Kestrel row; none flagged."
    above a table containing no variant. The adVNTR side has always suppressed its
    count line when it has no table, and this mirrors it.
    """
    html = render(negative_summary(tmp_path))

    assert "Showing 1 of 1 Kestrel row" not in html
    assert "Kestrel row" not in html, "a report with no Kestrel table still makes a Kestrel row-count claim"


def test_an_absent_kestrel_step_does_not_emit_a_zero_column_table(tmp_path) -> None:
    """``vntyper report`` renders a supplied summary, which need not have the step.

    ``to_html`` on a frame with no columns produces a headerless, bodyless table -
    a stray empty box beneath a heading. ``escaped_table_html`` returns "" for an
    empty frame instead, which is the hook the authored sentence hangs on.
    """
    write_summary(tmp_path, tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]))

    html = render(tmp_path)

    assert 'id="kestrel_table"' not in html, "a report with no Kestrel step still emits an empty table"
    assert '<tr style="text-align: right;">' not in html
    assert "Kestrel genotyping was not performed" in visible_text(html)


#: What ``record_step`` writes when a stage's result file is absent (#212). Built by
#: calling the real recorder rather than by hand: the flag is the only structural
#: evidence that "the step produced nothing" is not "the step found nothing", and a
#: hand-written imitation of it would let the writer change shape without this noticing.
_STEP_START = datetime(2026, 1, 1, 12, 0, 0, tzinfo=timezone.utc).replace(tzinfo=None)
_STEP_END = datetime(2026, 1, 1, 12, 0, 30, tzinfo=timezone.utc).replace(tzinfo=None)


def unreadable_kestrel_step(missing_result_file: Path) -> dict:
    """Record a Kestrel step whose result file does not exist, the way a run would.

    Args:
        missing_result_file: A path that is not there.

    Returns:
        dict: The step mapping ``record_step`` produced.
    """
    from vntyper.scripts.summary import record_step

    summary: dict = {"steps": []}
    record_step(
        summary,
        summary_steps.STEP_KESTREL,
        str(missing_result_file),
        "tsv",
        "java -jar kestrel.jar",
        _STEP_START,
        _STEP_END,
    )
    return summary["steps"][0]


def test_a_kestrel_result_file_that_is_missing_is_not_reported_as_a_negative(tmp_path) -> None:
    """#212's other half, closed on the report side.

    ``record_step`` flags the step ``result_file_missing`` because ``md5sum`` swallows
    the ``FileNotFoundError`` and ``parse_tsv`` turns it into a comment and an empty
    ``data`` list - which is exactly what a run that legitimately found nothing
    produces. Nothing read that flag, so the report rendered the two identically; and
    once the empty state was authored, the failure rendered as the *sentence* "No
    variant detected by Kestrel in this sample."

    A report that states a negative the run never established is the defect this whole
    issue exists to remove, so the third state is rendered as its own.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        unreadable_kestrel_step(tmp_path / "kestrel" / "kestrel_result.tsv"),
        input_files={"bam": "NEG1.bam"},
    )

    text = visible_text(render(tmp_path))

    assert "No variant detected by Kestrel" not in text, "a step that produced nothing is reported as a negative"
    assert "Kestrel genotyping was not performed" not in text, "the step did run; saying otherwise is a second claim"
    assert "result file is missing or could not be read" in text
    assert "this is not a negative" in text


@pytest.mark.parametrize(
    ("parsed_result", "case"),
    [
        (None, "record_step's initial value, left in place when parsing never ran"),
        ({"error": "Error parsing file: boom"}, "record_step's parse-failure shape"),
        ({"error": "Unsupported file type for result parsing: bed"}, "record_step's unsupported-type shape"),
    ],
    ids=["null", "parse-error", "unsupported-type"],
)
def test_a_kestrel_step_that_could_not_be_read_is_not_reported_as_a_negative(
    tmp_path, parsed_result, case: str
) -> None:
    """Every shape ``record_step`` can leave behind when it did not get a result.

    The missing-file flag is the common one; these are the others, and each is
    recognised structurally rather than by the wording of a message.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        {"step": summary_steps.STEP_KESTREL, "parsed_result": parsed_result},
    )

    text = visible_text(render(tmp_path))

    assert "No variant detected by Kestrel" not in text, f"{case} is reported as a negative"
    assert "result file is missing or could not be read" in text


def test_the_three_kestrel_states_do_not_read_alike(tmp_path) -> None:
    """Ran and called nothing, did not run, could not be read: three facts, three states.

    This is the same distinction #223 drew for an unreadable derived VCF and the one
    ``screening_summary.NOT_PERFORMED`` draws for adVNTR: a stage that was never asked
    to run has said nothing, a stage whose result could not be read has said nothing
    either, and a report that renders either as a negative is asserting something the
    run never established.

    Each sentence is checked against *both* of the others, so collapsing any two of the
    three back into one fails here rather than passing quietly.
    """
    states = {}
    for name in ("ran", "absent", "unreadable"):
        directory = tmp_path / name
        directory.mkdir()
        states[name] = directory
    negative_summary(states["ran"])
    write_summary(states["absent"], tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]))
    write_summary(
        states["unreadable"],
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        unreadable_kestrel_step(tmp_path / "gone.tsv"),
    )

    sentences = {
        "ran": "No variant detected by Kestrel",
        "absent": "Kestrel genotyping was not performed",
        "unreadable": "result file is missing or could not be read",
    }
    rendered = {name: visible_text(render(directory)) for name, directory in states.items()}

    for name, text in rendered.items():
        assert sentences[name] in text, f"the {name} state does not say so"
        for other, sentence in sentences.items():
            if other != name:
                assert sentence not in text, f"the {name} state also reads as the {other} state"


def test_suppressing_the_placeholder_leaves_the_screening_state_alone(tmp_path) -> None:
    """The state is computed from the rows, and a non-result was never one.

    With the placeholder in the frame every configured Kestrel rule broke on its
    ``Confidence`` condition and the block's ``default`` - "negative" - was returned;
    an empty frame returns the same default by the shortest path in
    ``compute_algorithm_result``. This pins that the two agree, so removing the row
    from the table cannot move the screening verdict.
    """
    html = render(negative_summary(tmp_path))

    assert 'data-state="no-finding"' in html
    assert "Kestrel: negative" in visible_text(html)


def test_a_row_that_names_a_variant_is_never_suppressed(tmp_path) -> None:
    """Adversarial: a real call whose ``Confidence`` happens to be the placeholder token.

    Suppression keys on the whole row being empty, not on one cell. A rule reading
    ``Confidence == "Negative"`` alone would delete this variant from the report.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, "Confidence": "Negative"}]),
    )

    html = render(tmp_path)

    assert 'id="kestrel_table"' in html, "a row naming a position, a REF and an ALT was suppressed"
    assert ">67</td>" in html
    assert "Showing 1 of 1 Kestrel row; none flagged." in html


def test_a_negative_run_still_names_its_sample_and_its_coverage(tmp_path) -> None:
    """The empty state replaces the table, not the report.

    A negative run is the commonest one in a cohort, so everything that makes the
    file a record - who it is about, what was measured - has to survive the branch.
    """
    html = render(negative_summary(tmp_path))

    assert "<h1>MUC1 VNTR report — NEG1</h1>" in html
    assert _coverage_qc_cell(html) == "PASS"


def test_kestrel_conversion_failure_preserves_both_frames(monkeypatch, caplog) -> None:
    """A formatting conversion failure keeps matching evidence and display data.

    The display frame carries the sample's own string **unescaped** here, and that is
    deliberate: since the table is rendered through ``escaped_table_html`` the escaping
    happens once, at render time, over every column not named in ``html_columns``.
    Escaping here as well would double-escape it, so a motif sequence containing a ``<``
    would reach the reader as the literal text ``&lt;``. What must stay true is that
    nothing sample-derived reaches the HTML unescaped, which is asserted end to end by
    ``test_every_kestrel_cell_but_the_two_we_build_is_escaped``.
    """
    monkeypatch.setattr(generate_report.pd, "to_numeric", Mock(side_effect=ValueError("bad depth")))
    caplog.set_level(logging.WARNING, logger=generate_report.logger.name)
    caplog.clear()

    display_frame, matching_frame = generate_report.build_kestrel_frames(
        [{**KESTREL_ROW, "Motif_sequence": "<untrusted>"}]
    )

    assert len(display_frame) == len(matching_frame) == 1
    assert matching_frame.loc[0, "Confidence"] == "High_Precision"
    assert matching_frame.loc[0, "Motif Sequence"] == "<untrusted>"
    assert display_frame.loc[0, "Confidence"] == (
        '<span style="color:red;font-weight:bold;text-decoration:underline solid;">High_Precision</span>'
    )
    assert display_frame.loc[0, "Motif Sequence"] == "<untrusted>"
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


def test_the_screening_state_reaches_the_report(positive_summary) -> None:
    """I2: the computed screening state must reach the template, not just `is_positive`.

    `positive_summary` is Kestrel-High_Precision, adVNTR-absent, well-covered -- so
    this also pins the exact provenance line's wording for the common case.
    """
    html = render(positive_summary)
    assert 'data-state="finding"' in html
    assert "Kestrel: High_Precision" in html
    assert "adVNTR: not performed" in html
    assert "Coverage QC: PASS" in html


def test_a_negative_screening_state_carries_the_no_finding_state(tmp_path) -> None:
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, []),
        tabular_step(summary_steps.STEP_ADVNTR, []),
    )
    html = render(tmp_path)
    assert 'data-state="no-finding"' in html
    assert "Kestrel: negative" in html
    assert "adVNTR: negative" in html


def test_the_provenance_line_never_prints_the_not_performed_token_raw(positive_summary) -> None:
    """`advntr_result` is the literal `"none"` when the stage never ran -- distinct
    from `"negative"`, which means it ran and found nothing. The provenance line must
    render this as words a reader understands, not the raw internal token."""
    html = render(positive_summary)
    assert "adVNTR: none" not in html


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


def test_the_advntr_not_performed_state_is_worded_once(positive_summary) -> None:
    """One state, one sentence - and it is the template's.

    ``generate_summary_report`` also built ``<p>adVNTR genotyping was not performed.</p>``
    for this state, which the template's ``{% if advntr_available and advntr_highlight %}``
    can never reach: when adVNTR is unavailable the guard is false and the ``{% else %}``
    below it prints its own, differently worded line. Two sentences for one state, one of
    them unreachable. This counts what the reader actually gets.
    """
    text = visible_text(render(positive_summary))

    assert "adVNTR genotyping was not performed or no adVNTR results are available." in text
    assert text.count("adVNTR genotyping was not performed") == 2, (
        "the report should say this exactly twice - once in the configured screening message and once "
        f"in the adVNTR section - and says it {text.count('adVNTR genotyping was not performed')} times"
    )


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


#: The adVNTR row the escaping checks below plant a payload into. Kept beside them
#: rather than reused from the numeric section, so a change to the precision
#: specimen cannot quietly change what the trust-boundary checks are testing.
ADVNTR_ESCAPING_ROW = {
    "VID": 25561,
    "Variant": "I22_G_LEN1",
    "NumberOfSupportingReads": 14,
    "MeanCoverage": 98.5,
    "Pvalue": 1e-09,
    "RU": "CGCGG",
    "POS": 67,
    "REF": "G",
    "ALT": "GG",
    "Flag": "Not flagged",
}

#: Every adVNTR display column except the one exemption, derived from the display
#: table so a column added later is covered without editing this list. These are the
#: cells a widened ``html_columns`` would expose.
ADVNTR_ESCAPED_COLUMNS = tuple(column for column in generate_report.ADVNTR_DISPLAY_COLUMNS if column != "Flag")


def advntr_summary(tmp_path: Path, **overrides) -> Path:
    """Write a summary whose adVNTR step carries one row with ``overrides`` applied."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        tabular_step(summary_steps.STEP_ADVNTR, [{**ADVNTR_ESCAPING_ROW, **overrides}]),
    )
    return tmp_path


def test_a_malicious_flag_in_an_advntr_row_is_server_escaped(tmp_path) -> None:
    """The adVNTR table's one escaping exemption, planted against.

    #242 moved this table from a blanket ``to_html(escape=True)`` to the same
    per-column model the Kestrel table uses: ``escape=False`` for the whole table,
    with ``escape_frame_cells`` escaping every cell except ``Flag``, whose markup
    ``flag_html`` builds *and escapes* itself. The two halves have to stay together.

    **The concrete state this catches**: drop or reorder the ``flag_html`` call in
    ``generate_summary_report`` while ``html_columns=("Flag",)`` stays. Nothing
    raises, the table still renders, every other test still passes - and a
    sample-derived ``Flag`` reaches the HTML unescaped. That is this codebase's
    signature failure mode: a silently wrong call, not a crash.
    """
    html = render(advntr_summary(tmp_path, Flag=PAYLOAD))

    assert PAYLOAD not in html
    assert ESCAPED in html


@pytest.mark.parametrize("column", ADVNTR_ESCAPED_COLUMNS)
def test_every_other_advntr_cell_is_escaped(tmp_path, column: str) -> None:
    """The columns the exemption does *not* cover, including the numeric ones.

    A payload in a numeric column is not absurd: the value comes out of a supplied
    ``pipeline_summary.json``, the formatters pass a non-numeric value through
    untouched by design, and the whole table is rendered with ``escape=False``. So
    each of these is escaped only because it is *not* named in ``html_columns`` -
    which is exactly what widening that tuple would undo.
    """
    html = render(advntr_summary(tmp_path, **{column: PAYLOAD}))

    assert PAYLOAD not in html, f"an adVNTR {column} value reached the HTML unescaped"
    assert ESCAPED in html


def test_an_advntr_flagged_row_states_its_reason_in_the_table(tmp_path) -> None:
    """The adVNTR flag cell is rendered server-side too, not just Kestrel's.

    Both tables lost their client-side flag renderer in #242, so both have to gain
    the server-side one; the escaping check above passes vacuously if this table
    stopped carrying a flag cell at all.
    """
    from vntyper.scripts.report_formatting import FLAG_WARNING_GLYPH

    html = render(advntr_summary(tmp_path, Flag="Low_Depth"))

    assert FLAG_WARNING_GLYPH in html
    assert "Low_Depth" in html


#: Every Kestrel display column except the two whose markup VNtyper builds itself and
#: the one whose value never survives to be escaped. Derived from the display table so
#: a column added later is covered without editing this list.
#:
#: ``Depth_Score`` is excluded because ``build_kestrel_frames`` runs it through
#: ``pd.to_numeric(errors="coerce")``, so a non-numeric value is NaN before it reaches
#: any escaping at all - covered by
#: ``test_a_payload_in_the_depth_score_is_coerced_rather_than_escaped``.
KESTREL_ESCAPED_COLUMNS = tuple(
    column for column in generate_report.KESTREL_DISPLAY_COLUMNS if column not in ("Confidence", "Flag", "Depth_Score")
)


@pytest.mark.parametrize("column", KESTREL_ESCAPED_COLUMNS)
def test_every_kestrel_cell_but_the_two_we_build_is_escaped(tmp_path, column: str) -> None:
    """The Kestrel table's escaping, asserted where it is observable.

    The table is rendered with ``escape=False`` - the ``Confidence`` span and the
    ``Flag`` cell are markup VNtyper builds - so every other cell is escaped only
    because ``escaped_table_html`` escapes it. #242 routed this table through that
    helper instead of calling ``to_html`` directly; before that the frame was escaped a
    step earlier, in ``build_kestrel_frames``. Either arrangement is correct and only
    one of them may be in force, because doing both renders ``&lt;`` to the reader as
    text. This is the assertion that survives the choice.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, column: PAYLOAD}]),
    )

    html = render(tmp_path)

    assert PAYLOAD not in html, f"a Kestrel {column} value reached the HTML unescaped"
    assert ESCAPED in html
    assert "&amp;lt;" not in html, "the Kestrel table is escaped twice, so the reader sees the escape sequence"


def test_a_payload_in_the_depth_score_is_coerced_rather_than_escaped(tmp_path) -> None:
    """The one Kestrel column whose value cannot reach the escaping at all.

    ``build_kestrel_frames`` sorts on ``Depth Score`` and runs the column through
    ``pd.to_numeric(errors="coerce")`` first, so anything that is not a number is NaN
    before the table is built. Asserted rather than assumed: it is the reason that
    column is left out of the parametrisation above, and if the coercion moved, the
    exclusion would silently stop covering anything.
    """
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, "Depth_Score": PAYLOAD}]),
    )

    html = render(tmp_path)

    assert PAYLOAD not in html
    assert ESCAPED not in html, "the payload survived the coercion, so this column does need escaping"


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
    autoescaping on without marking them would print the span markup as text.

    The literal gained `role="img"` and a name in #242's accessibility pass: a bare
    `&#10004;` is announced as its code point or skipped, so the `Status` column
    read as empty to a screen reader while looking complete on screen.
    """
    html = render(positive_summary)
    assert '<span style="color:green;font-weight:bold;" role="img" aria-label="No warning">&#10004;</span>' in html
    assert "&lt;span style=" not in html


def test_escaping_does_not_neuter_the_screening_message(positive_summary) -> None:
    """Configured messages carry `<br>` line breaks."""
    html = render(positive_summary)
    assert "<br>" in html


def test_escaping_does_not_neuter_the_results_tables(positive_summary) -> None:
    html = render(positive_summary)
    assert 'id="kestrel_table"' in html
    assert '<span style="color:red;font-weight:bold;text-decoration:underline solid;">High_Precision</span>' in html


def test_escaping_does_not_neuter_the_igv_script_block(positive_summary) -> None:
    html = render(positive_summary)
    assert 'const tableJson = {"headers": [], "rows": []};' in html


#: Every context value the template is allowed to interpolate unescaped, and why.
#:
#: This is the audit record of the report's trust model, so each entry has to say
#: what is actually true of the value today. **Neither results table is escaped by
#: pandas.** Both go out through ``to_html(escape=False)``, because both carry
#: markup VNtyper built, and both therefore rely on ``escape_frame_cells`` having
#: escaped every cell first *except* the columns named in ``html_columns``. That
#: exemption is per column and each exempted column escapes its own value; widening
#: ``html_columns`` is what would expose a sample's string, not editing this dict.
SAFE_BY_DESIGN = {
    "kestrel_highlight": (
        "pandas table rendered with escape=False; escape_frame_cells escapes every cell except "
        "`Confidence` and `Flag`, whose markup confidence_html and flag_html build and escape themselves"
    ),
    "advntr_highlight": (
        "pandas table rendered with escape=False through escaped_table_html; escape_frame_cells escapes "
        "every cell except `Flag`, whose markup flag_html builds and escapes itself - or one of two fixed <p> "
        "sentences when adVNTR did not run or found nothing"
    ),
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
    a BAM-only run) simply omits the four fields.

    Three of them stay conditional: they name a *file* the run opened, and a
    missing label is the honest rendering of a run that opened none. The fourth
    moved into the provenance block, where an absent value is stated rather than
    left blank - the requested assembly is what selects the region even for a run
    that reads no reference, so its absence is itself a fact about the run (#242).
    """
    html = render(positive_summary)

    assert _labeled_value(html, "Reference assembly requested") == "not recorded by this run"
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


# ---------------------------------------------------------------------------
# No per-sample result row is ever hidden - issue #242
# ---------------------------------------------------------------------------

PER_SAMPLE_TEMPLATE = TEMPLATE_DIR / "report_template.html"

#: Constructs that take a row out of the reader's view. Deliberately shape-based
#: rather than name-based: banning the literal `toggleFlagged` would be satisfied by
#: renaming it. ``.remove(`` is matched only with empty parentheses, so
#: ``classList.remove("selected")`` - which removes a class, not an element - does not
#: register.
_ROW_HIDING_VERB = re.compile(
    r"""
      removeChild\s*\(
    | \.remove\s*\(\s*\)
    | \.detach\s*\(
    | \.hide\s*\(
    | \.filter\s*\(
    | style\.display\s*=\s*['"]\s*none
    | classList\.add\s*\(\s*['"](?:d-none|hidden|invisible)
    | setAttribute\s*\(\s*['"]hidden
    | \.hidden\s*=\s*true
    """,
    re.VERBOSE,
)

#: Words that mean the surrounding statement is talking about a table row. A verb
#: from the list above is only an offence when it is applied to one of these.
_ROW_SUBJECT = re.compile(r"\b(?:tr|nTr|aoData|row|rows|tbody|dataIndex)\b", re.IGNORECASE)

#: A CSS declaration block that makes its subject unreadable.
_CSS_HIDES = re.compile(r"display\s*:\s*none|visibility\s*:\s*hidden")

#: ``tr`` as an element selector - not as part of ``.tr-thing`` or ``#tr``.
_CSS_ROW_SELECTOR = re.compile(r"(?<![\w.#-])tr(?![\w-])")


def _style_blocks(source: str) -> list[str]:
    """Return the contents of every ``<style>`` element."""
    return re.findall(r"<style[^>]*>(.*?)</style>", source, re.DOTALL)


def test_no_per_sample_result_row_enters_a_hiding_path() -> None:
    """No construct in the per-sample template can take a results row off the page.

    The defect issue #242 is named after is not that a flag is styled badly; it is
    that the row the screening summary narrates is *removed from the DOM* by a
    client-side DataTables predicate before the reader sees the table. This is the
    invariant that stops it coming back under another name.

    **What this can see.** Four shapes, over the template source text: DataTables'
    ``ext.search`` row-visibility hook; its ``paging`` and ``searching`` options,
    which hide rows just as effectively; a JavaScript statement pairing a removal or
    hiding verb with a word meaning "row"; and a CSS rule that hides ``tr``.

    **What this cannot see.** It is a tripwire, not a behavioural test - the unit
    tier has no JavaScript engine. It cannot evaluate an expression, so a hiding
    call assembled from string fragments, spread across lines, or reached through an
    alias (``var f = el.remove; f.call(row)``) escapes it, as does any construct
    whose shape is not one of the four. ``tests/browser/test_flagged_rows.py``
    measures the visible row count in a real browser and is what actually proves the
    rows are there; this exists so that a *renamed* filter fails in the tier
    everybody runs.
    """
    source = PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")

    assert "$.fn.dataTable" in source or "DataTable(" in source, (
        "the template no longer initialises DataTables at all; these assertions would be vacuous"
    )
    assert "ext.search" not in source, (
        "the per-sample report registers a DataTables row-visibility filter, which removes rows from the DOM"
    )
    assert '"paging": false' in source, "DataTables paging hides every row past the first page"
    assert '"searching": false' in source, "DataTables searching hides every row that does not match"

    offenders = [
        (number, line.strip())
        for number, line in enumerate(source.splitlines(), start=1)
        if _ROW_HIDING_VERB.search(line) and _ROW_SUBJECT.search(line)
    ]
    assert offenders == [], f"the per-sample template applies a hiding construct to a row: {offenders}"

    hidden_rows = [
        (selector.strip(), body.strip())
        for block in _style_blocks(source)
        for selector, body in re.findall(r"([^{}]+)\{([^{}]*)\}", block)
        if _CSS_ROW_SELECTOR.search(selector) and _CSS_HIDES.search(body)
    ]
    assert hidden_rows == [], f"the per-sample template hides table rows with CSS: {hidden_rows}"

    inline = re.findall(r"<tr[^>]*style=\"[^\"]*(?:display\s*:\s*none|visibility\s*:\s*hidden)", source)
    assert inline == [], f"the per-sample template hides a row with an inline style: {inline}"


def test_the_flag_switch_highlights_rather_than_filters() -> None:
    """P4: hiding flagged rows is defensible for cohort triage and indefensible for
    a single-patient read, so the per-sample switch changes emphasis only."""
    source = PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")

    assert 'id="highlightFlagged"' in source
    assert "Highlight flagged values" in source
    assert "toggleFlagged" not in source, "the per-sample report no longer has a show/hide flagged switch"
    assert "Show flagged values" not in source


def test_the_highlight_switch_survives_a_script_that_never_loaded() -> None:
    """The handler must not share a ``<script>`` element with the jQuery code.

    A ``$`` that never resolved throws at the top of its block and takes every
    statement after it down with it, so a handler appended to that block works only
    for readers whose browser reached the CDN - which is the shape of defect this
    whole change is removing.
    """
    source = PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")
    blocks = re.findall(r"<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>", source, re.DOTALL)

    owning = [block for block in blocks if "highlightFlagged" in block]
    assert len(owning) == 1, f"expected exactly one script block to handle the switch, found {len(owning)}"
    assert "$(" not in owning[0], "the highlight handler shares its block with jQuery code that can throw first"


def test_the_cohort_report_keeps_its_own_filter() -> None:
    """Scope boundary, pinned: this change is per-sample only (precondition P4)."""
    cohort = (TEMPLATE_DIR / "cohort_summary_template.html").read_text(encoding="utf-8")

    assert "ext.search" in cohort, "the cohort filter was removed; that is a separate, reviewed decision"
    assert "toggleFlagged" in cohort


# ---------------------------------------------------------------------------
# Every displayed number is computed on the server - issue #242
# ---------------------------------------------------------------------------

ADVNTR_ROW = {
    "VID": 25561,
    "Variant": "I22_G_LEN1",
    "NumberOfSupportingReads": 14,
    "MeanCoverage": 98.5,
    "Pvalue": 1e-09,
    "RU": "CGCGG",
    "POS": 67,
    "REF": "G",
    "ALT": "GG",
    "Flag": "Not flagged",
}

#: A Kestrel row whose every numeric column discriminates between the rounding rules.
PRECISE_KESTREL_ROW = {**KESTREL_ROW, "Depth_Score": 0.010012}


def test_no_displayed_number_is_computed_in_the_browser() -> None:
    """``applyRounding`` rewrote every numeric cell of every initialised table."""
    source = PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")

    assert "applyRounding" not in source
    assert "roundValue" not in source
    assert "toFixed" not in source


def test_the_reader_is_not_shown_a_row_count_the_browser_computed() -> None:
    """DataTables' "Showing 1 to 3 of 3 entries" footer is a second, contradictory
    count that only exists when the CDNs resolve."""
    assert '"info": false' in PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")


@pytest.fixture
def both_tables(tmp_path):
    """A run with both algorithms reporting, at discriminating precision."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [PRECISE_KESTREL_ROW]),
        tabular_step(summary_steps.STEP_ADVNTR, [ADVNTR_ROW]),
    )
    return tmp_path


@pytest.mark.parametrize(
    ("column", "displayed"),
    [
        ("POS", ">67<"),
        ("Estimated_Depth_AlternateVariant", ">120<"),
        ("Estimated_Depth_Variant_ActiveRegion", ">12000<"),
        ("Depth_Score", ">0.010012<"),
        ("VID", ">25561<"),
        ("NumberOfSupportingReads", ">14<"),
        ("MeanCoverage", ">98.50<"),
        ("Pvalue", ">1e-09<"),
    ],
)
def test_each_numeric_column_reaches_the_html_already_formatted(both_tables, column, displayed) -> None:
    """The string the reader sees is in the file, not assembled by a script."""
    assert displayed in render(both_tables), f"{column} is not rendered as {displayed!r}"


def test_a_p_value_is_not_destroyed_by_rounding(both_tables) -> None:
    """``parseFloat((1e-9).toFixed(4)).toString()`` is ``"0"``: the online report
    displayed a highly significant adVNTR p-value as zero."""
    html = render(both_tables)

    assert ">0<" not in html.split("<h2>adVNTR Identified Variants</h2>")[1].split("</table>")[0]


# ---------------------------------------------------------------------------
# The flag cell says why - issue #242
# ---------------------------------------------------------------------------


def _kestrel_table(html: str) -> str:
    """Return the Kestrel table's markup."""
    start = html.index('id="kestrel_table"')
    return html[start : html.index("</table>", start)]


def test_a_flagged_row_states_its_reason_in_the_table(tmp_path) -> None:
    """The reason used to live only in a Bootstrap ``title`` on a glyph, so it was
    invisible in print, invisible to a screen reader and absent when jQuery failed."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [{**KESTREL_ROW, "Flag": "Low_Depth"}]),
    )

    from vntyper.scripts.report_formatting import FLAG_WARNING_GLYPH

    table = _kestrel_table(render(tmp_path))

    assert "Low_Depth" in table
    assert FLAG_WARNING_GLYPH in table


def test_an_unflagged_row_is_marked_as_clean(positive_summary) -> None:
    from vntyper.scripts.report_formatting import FLAG_OK_GLYPH

    table = _kestrel_table(render(positive_summary))

    assert "Not flagged" in table
    assert FLAG_OK_GLYPH in table


def test_the_kestrel_table_states_how_many_rows_are_shown(tmp_path) -> None:
    """Rendered in Python from the frame, so it cannot contradict the table."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(
            summary_steps.STEP_KESTREL,
            [KESTREL_ROW, {**KESTREL_ROW, "Motif": "6", "Flag": "Low_Depth", "Depth_Score": 0.008}],
        ),
    )

    assert "Showing 2 of 2 Kestrel rows; 1 flagged." in render(tmp_path)


def test_the_advntr_table_states_how_many_rows_are_shown(both_tables) -> None:
    assert "Showing 1 of 1 adVNTR row; none flagged." in render(both_tables)


def test_a_run_with_no_advntr_table_makes_no_advntr_count_claim(positive_summary) -> None:
    assert "adVNTR row" not in render(positive_summary)


# ---------------------------------------------------------------------------
# The report identifies itself (#242)
# ---------------------------------------------------------------------------

#: The four shapes `resolve_pipeline_input` produces (`pipeline_inputs.py:162-170`).
#: The template branched on two of them, so a CRAM run and a single-end FASTQ run
#: rendered an empty `Input Files:` line.
INPUT_SHAPES = [
    {"cram": "S1.cram"},
    {"fastq1": "S1.fq.gz"},
    {"bam": "S1.bam"},
    {"fastq1": "a.fq", "fastq2": "b.fq"},
]


def _provenance_block(html: str) -> str:
    """The provenance section of the rendered report, and nothing else.

    The extraction stops at the first ``</div>``, which is only the end of the
    block while the block contains no nested ``<div>``. If one is ever added, this
    would silently return a prefix and every negative assertion made against it -
    ``"hg19" not in ...`` above all - would pass by looking at less and less of the
    report. That is a guard quietly becoming a no-op, so it is checked rather than
    assumed.

    Args:
        html: The rendered report.

    Returns:
        str: The markup between the block's opening tag and its closing ``</div>``.
    """
    start = html.index('id="provenance"')
    end = html.index("</div>", start)
    block = html[start:end]
    assert "<div" not in block, (
        "the provenance block now contains a nested <div>, so this extraction stops at the "
        "inner closing tag; every negative assertion made against the block is now vacuous. "
        "Parse the block properly rather than deleting this check."
    )
    return block


@pytest.mark.parametrize("inputs", INPUT_SHAPES)
def test_every_input_shape_names_its_files(tmp_path, inputs) -> None:
    """A report that does not say what it was run on is not a record of anything."""
    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        input_files=inputs,
    )

    listed = _labeled_value(render(tmp_path), "Input Files")

    for name in inputs.values():
        assert name in listed, f"{name!r} is missing from the Input Files line"


def test_the_title_carries_the_sample_name(tmp_path) -> None:
    write_summary(tmp_path, input_files={"cram": "S1.cram"})

    assert "<title>MUC1 VNTR report — S1" in render(tmp_path)


def test_the_heading_carries_the_sample_name(tmp_path) -> None:
    """Two tabs and a printed page are all indistinguishable without this."""
    write_summary(tmp_path, input_files={"cram": "S1.cram"})

    assert "<h1>MUC1 VNTR report — S1</h1>" in render(tmp_path)


def test_the_header_names_the_assay_and_the_version(tmp_path) -> None:
    write_summary(tmp_path, input_files={"bam": "S1.bam"}, version="2.0.18")

    html = render(tmp_path)

    assert _labeled_value(html, "Assay") == "MUC1 coding VNTR genotyping"
    assert _labeled_value(html, "Sample") == "S1"
    assert _labeled_value(html, "VNtyper Version") == "2.0.18"


def test_the_printed_header_line_names_the_run_and_its_use(tmp_path) -> None:
    """The line at the head of the printed record, rendered end to end.

    It is ``display: none`` on screen and the only identity a filed sheet carries, so
    every field it states has to survive the render - and, like everything else in the
    report, it is escaped rather than interpolated into the stylesheet (see
    ``tests/unit/test_report_presentation.py::test_no_value_is_interpolated_into_a_stylesheet``).
    """
    from vntyper.scripts.report_identity import RESEARCH_USE_STATEMENT

    write_summary(
        tmp_path,
        tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]),
        input_files={"bam": "S1.bam"},
        version="2.0.18",
        reference_assembly_requested="hg19",
    )

    header = re.search(r'<div class="print-header">(.*?)</div>', render(tmp_path), re.DOTALL)

    assert header, "the report has no printed header line"
    line = " ".join(header.group(1).split())
    for value in ("S1", "MUC1 coding VNTR genotyping", "hg19", "VNtyper 2.0.18", RESEARCH_USE_STATEMENT):
        assert value in line, f"the printed header line does not state {value!r}: {line!r}"


def test_an_explicit_sample_name_wins_over_the_input_files(tmp_path) -> None:
    write_summary(tmp_path, input_files={"cram": "S1.cram"})

    html = render(tmp_path, sample_name="PATIENT_042")

    assert "<title>MUC1 VNTR report — PATIENT_042" in html
    assert _labeled_value(html, "Sample") == "PATIENT_042"


def test_the_report_uses_the_name_the_run_recorded_for_itself(tmp_path) -> None:
    """`vntyper pipeline -s PATIENT_042 --bam foo.bam` labels Kestrel's outputs
    and its VCF header `PATIENT_042`. Until the run recorded that name, the report
    could not see it and titled itself `foo` -- one run, two identities, in the
    artefact that gets forwarded."""
    write_summary(tmp_path, input_files={"bam": "foo.bam"}, sample_name="PATIENT_042")

    html = render(tmp_path)

    assert "<title>MUC1 VNTR report — PATIENT_042</title>" in html
    assert _labeled_value(html, "Sample") == "PATIENT_042"
    # The file it ran on is still stated; the two are different facts.
    assert _labeled_value(html, "Input Files") == "foo.bam"


def test_an_explicit_sample_name_beats_the_one_the_run_recorded(tmp_path) -> None:
    write_summary(tmp_path, input_files={"bam": "foo.bam"}, sample_name="PATIENT_042")

    html = render(tmp_path, sample_name="RENAMED")

    assert _labeled_value(html, "Sample") == "RENAMED"


def test_the_recorded_placeholder_name_does_not_displace_the_derivation(tmp_path) -> None:
    """`cli_handlers` records the literal `"sample"` when it resolved nothing.
    A report titled `sample` beside an input file that names the sample perfectly
    well is a placeholder winning over a fact."""
    write_summary(tmp_path, input_files={"fastq1": "S1_R1.fastq.gz"}, sample_name="sample")

    assert _labeled_value(render(tmp_path), "Sample") == "S1"


def test_a_legacy_summary_with_no_recorded_name_derives_one_as_before(tmp_path) -> None:
    """Every archived run predates the field. Adding a level above the derivation
    must not change what those reports are called."""
    write_summary(tmp_path, input_files={"fastq1": "S1_R1.fastq.gz", "fastq2": "S1_R2.fastq.gz"})

    assert _labeled_value(render(tmp_path), "Sample") == "S1"


def test_a_summary_with_no_input_files_still_names_the_report(tmp_path) -> None:
    """A report with nothing to derive a name from must still be a report."""
    write_summary(tmp_path)

    html = render(tmp_path)

    assert "<title>MUC1 VNTR report — unnamed sample" in html
    assert _labeled_value(html, "Input Files") == "not recorded by this run"


def test_a_malformed_input_files_value_costs_the_line_and_not_the_report(tmp_path) -> None:
    """`vntyper report` renders a *supplied* summary (#207), so a wrong-typed
    `input_files` must not replace the report with a traceback."""
    write_summary(tmp_path, tabular_step(summary_steps.STEP_KESTREL, [KESTREL_ROW]), input_files=["sample.bam"])

    html = render(tmp_path)

    assert _labeled_value(html, "Input Files") == "not recorded by this run"
    assert "<title>MUC1 VNTR report — unnamed sample" in html
    assert _kestrel_table(html), "the rest of the report must still render"


#: An injection payload carrying no ``/``. ``PAYLOAD`` does, and the sample name
#: is taken from a *basename*, so ``</script>`` alone is stripped by the path
#: split before any escaping happens -- which would make an escaping assertion on
#: the title pass for the wrong reason.
TITLE_PAYLOAD = "<img src=x onerror=alert(1)>"
TITLE_ESCAPED = "&lt;img src=x onerror=alert(1)&gt;"


def test_a_sample_name_derived_from_an_input_file_is_escaped(tmp_path) -> None:
    """The name is derived from a sample-supplied basename, so it reaches
    ``<title>`` and ``<h1>`` as attacker-influenced text."""
    write_summary(tmp_path, input_files={"bam": f"{TITLE_PAYLOAD}.bam"})

    html = render(tmp_path)

    assert TITLE_PAYLOAD not in html
    assert f"<title>MUC1 VNTR report — {TITLE_ESCAPED}</title>" in html
    assert f"<h1>MUC1 VNTR report — {TITLE_ESCAPED}</h1>" in html


def test_a_sample_name_in_the_printed_header_line_is_escaped(tmp_path) -> None:
    """The printed header line states the name too, so it is the third place to check.

    It is also the reason that line is a block in the document: the alternative was a
    running header in the page margin, which would have put this value inside a
    ``<style>`` element where HTML escaping means nothing.
    """
    write_summary(tmp_path, input_files={"bam": f"{TITLE_PAYLOAD}.bam"})

    header = re.search(r'<div class="print-header">(.*?)</div>', render(tmp_path), re.DOTALL)

    assert header, "the report has no printed header line"
    assert TITLE_PAYLOAD not in header.group(1)
    assert TITLE_ESCAPED in header.group(1)


# ---------------------------------------------------------------------------
# Provenance: recorded, or said to be absent -- never guessed (#242, P5)
# ---------------------------------------------------------------------------


def test_a_legacy_summary_renders_not_recorded(tmp_path) -> None:
    """A summary written before this change carries none of the provenance
    fields, and the report must say so rather than reading the config default.

    ``config["default_values"]["reference_assembly"]`` is ``hg19``. Printing it
    would mislabel every ``--reference-assembly`` override, and it cannot
    reconstruct ``--custom-regions`` at all.
    """
    write_summary(tmp_path, version="2.0.11")

    html = render(tmp_path)

    assert "not recorded by this run" in html
    assert "hg19" not in _provenance_block(html)


def test_a_current_summary_prints_the_resolved_region(tmp_path) -> None:
    write_summary(tmp_path, schema_version=1, region_resolved="chr1:155184000-155194000")

    assert "chr1:155,184,000-155,194,000" in render(tmp_path)


def test_the_provenance_block_reads_the_assembly_fields_already_recorded(tmp_path) -> None:
    """Two of the four provenance rows are not new keys.

    ``assembly_declared`` is ``reference_assembly_requested``, written by
    ``start_summary``; ``assembly_detected`` is the ``BAM Header Parsing`` step's
    ``assembly_text``. Recording either again under a second name would be the
    divergent-source problem ``cli_report.py``'s docstring warns about.
    """
    write_summary(
        tmp_path,
        {
            "step": summary_steps.STEP_BAM_HEADER,
            "parsed_result": {"assembly_text": "GRCh38", "assembly_contig": "chr1"},
        },
        reference_assembly_requested="hg38_ensembl",
        schema_version=1,
        region_resolved="chr1:155184000-155194000",
    )

    block = _provenance_block(render(tmp_path))

    assert "hg38_ensembl" in block
    assert "GRCh38" in block
    assert "chr1:155,184,000-155,194,000" in block
    assert "not recorded by this run" not in block


def test_the_summary_schema_version_is_shown(tmp_path) -> None:
    write_summary(tmp_path, schema_version=1)

    assert _labeled_value(render(tmp_path), "Summary schema version") == "1"


# ---------------------------------------------------------------------------
# Run time is not render time
# ---------------------------------------------------------------------------


def test_the_run_time_and_the_render_time_are_both_shown(tmp_path) -> None:
    write_summary(tmp_path, pipeline_start="2020-01-02T03:04:05.678901")

    html = render(tmp_path)

    assert _labeled_value(html, "Pipeline run started") == "2020-01-02 03:04:05 UTC"
    # Both carry a zone, so a reader can subtract one from the other.
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \S+", _labeled_value(html, "This report rendered"))


def test_re_rendering_an_archived_run_does_not_restamp_the_run_time(tmp_path) -> None:
    """`vntyper report -o <finished run>` produced an artefact claiming to be
    newer than the analysis, because the only date on it was `datetime.now()`."""
    write_summary(tmp_path, pipeline_start="2020-01-02T03:04:05.678901")

    first = _labeled_value(render(tmp_path), "Pipeline run started")
    second = _labeled_value(render(tmp_path), "Pipeline run started")

    assert first == second == "2020-01-02 03:04:05 UTC"


def test_a_summary_with_no_start_time_says_so(tmp_path) -> None:
    write_summary(tmp_path)

    assert _labeled_value(render(tmp_path), "Pipeline run started") == "not recorded by this run"


# ---------------------------------------------------------------------------
# What the pipeline actually wrote, read back by the real report
# ---------------------------------------------------------------------------


def test_the_pipeline_puts_the_resolved_region_on_disk_before_the_report_reads_it(tmp_path) -> None:
    """The write-ordering trap, pinned end to end.

    ``pipeline_summary.json`` is written incrementally by ``record_step``, and
    ``generate_summary_report`` reads it back **from disk** -- while the final
    ``write_summary`` runs after the report. The resolved region does not exist
    until the coverage stage, so a key set after the last ``record_step`` would
    never reach the report at all.

    This renders from the bytes that were on disk at the instant the pipeline
    called the report generator, not from a hand-built summary, so it fails if
    the ordering regresses even when the finished file is correct.
    """
    from tests.support.pipeline_harness import run_pipeline_under_harness

    run_dir = tmp_path / "run"
    captured: dict[str, str] = {}

    def _capture(*args, **kwargs):
        """Read the summary at the instant the pipeline asked for the report."""
        captured["summary"] = (run_dir / "pipeline_summary.json").read_text(encoding="utf-8")

    harness = run_pipeline_under_harness(run_dir, stage_side_effects={"generate_summary_report": _capture})
    assert harness.error is None

    on_disk = json.loads(captured["summary"])
    assert on_disk["schema_version"] == 1
    assert on_disk["region_resolved"] == "chr1:155158000-155163000"

    report_dir = tmp_path / "rendered"
    report_dir.mkdir()
    (report_dir / "pipeline_summary.json").write_text(captured["summary"], encoding="utf-8")

    html = render(report_dir)

    assert "chr1:155,158,000-155,163,000" in _provenance_block(html)
    # The harness passes the literal `"sample"` placeholder, so this is also the
    # end-to-end proof that it falls through to the input basename (`in.bam`).
    assert "<title>MUC1 VNTR report — in</title>" in html


def test_the_operators_sample_name_survives_the_run_into_the_report(tmp_path) -> None:
    """`-s` on `vntyper pipeline`, end to end, through what the pipeline wrote.

    ``start_summary`` runs before any step, so the name is on disk from the first
    ``record_step`` onwards -- the same write-ordering the resolved region needed,
    with no new plumbing.
    """
    from tests.support.pipeline_harness import run_pipeline_under_harness

    run_dir = tmp_path / "run"
    captured: dict[str, str] = {}

    def _capture(*args, **kwargs):
        captured["summary"] = (run_dir / "pipeline_summary.json").read_text(encoding="utf-8")

    harness = run_pipeline_under_harness(
        run_dir,
        sample_name="PATIENT_042",
        stage_side_effects={"generate_summary_report": _capture},
    )
    assert harness.error is None

    assert json.loads(captured["summary"])["sample_name"] == "PATIENT_042"
    # The same string the run handed Kestrel, so the report and the VCF agree.
    assert harness.kwargs("run_kestrel")["sample_name"] == "PATIENT_042"

    report_dir = tmp_path / "rendered"
    report_dir.mkdir()
    (report_dir / "pipeline_summary.json").write_text(captured["summary"], encoding="utf-8")

    assert "<title>MUC1 VNTR report — PATIENT_042</title>" in render(report_dir)
