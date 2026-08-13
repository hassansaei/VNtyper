"""The archived PDF, printed by a real engine and read back as text.

A report is forwarded and filed as a PDF, and that PDF outlives the HTML it came
from. Until issue #242 neither template carried a single author media rule -
measured in the browser, ``document.styleSheets`` reported ``mediaRules: []`` -
so everything the printed page did was whatever the screen stylesheet happened to
imply. Three of those were defects, and none of them is visible to a test that
reads the stylesheet:

* ``.table td`` clamps to ``max-width: 150px`` with ``overflow: hidden`` and
  reveals the rest on ``:hover``. **The motif sequence is 121 bp in real data**
  and paper has no hover, so the archived record carried an ellipsis where the
  evidence was. The 13-character fixture in ``tests/unit/test_generate_report.py``
  fits inside the clamp and proves nothing about this;
* a collapsed ``<details>`` printed as its summary and nothing else;
* nothing after page 1 said whose report it was, and no page was numbered.

The unit tier can assert that a rule exists. It cannot assert that the rule
reached the paper: a stylesheet can satisfy every source-text assertion while the
PDF still omits closed content, page numbers or long cells. So this prints the
report and reads the artefact.

Supported print renderer
------------------------
**Chromium 151.0.7922.34**, driven through Playwright's ``page.pdf()`` with
``prefer_css_page_size=True`` so the document's own ``@page`` rule decides the
sheet, and text extracted with **poppler ``pdftotext`` 26.01.0** (a system binary,
not a Python dependency - the test skips where it is absent). This matters
because page counters and running headers differ materially between engines:
``@page`` margin boxes with ``counter(page)``/``counter(pages)`` and ``string()``
are honoured here and are ignored outright by Firefox and WebKit, neither of which
exposes a PDF API to Playwright at all. Measured, not assumed:

* ``page.pdf()`` **does** dispatch ``beforeprint``, so the handler that reopens a
  reader's collapsed sections runs on this path exactly as it does behind Ctrl+P;
* the ``@media print`` rule ``details > summary ~ *`` does **not** open a closed
  ``<details>`` in Chromium - the UA hides closed-details content in a way no
  author rule reaches. It is the no-JS fallback for engines that honour it, and
  the handler is what does the work here. The JavaScript-disabled case below
  therefore asserts what is actually true rather than what would be convenient.

What a reader with scripting off can still lose
-----------------------------------------------
Measured in this file, both states with one gesture between them: a ``<details>``
served ``open`` prints its contents with no script at all, and the same disclosure
after a click on its ``<summary>`` prints its heading and nothing else. Nothing in
CSS can undo that click - ``open`` is an attribute, not a style - so the mitigation
is that the report is *served* with every print-relevant disclosure open, and
collapsing one is a deliberate act rather than the state the file arrives in. The
residual limitation is stated in ``docs/pipeline/reports.md`` in the same words.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
from playwright.sync_api import Browser

from tests.browser.conftest import COVERAGE_ROW, TEMPLATE_DIR
from vntyper.cli import load_config
from vntyper.scripts import summary_steps
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.report_identity import ASSAY_NAME, RESEARCH_USE_STATEMENT

pytestmark = pytest.mark.browser

#: How long the ``Motif_sequence`` cell is in real data: **121 bp**, against the
#: 13-character stub the unit fixtures use. The length is the whole point - anything
#: that fits inside the 150px screen clamp cannot exhibit the defect this file exists
#: for - so the bases are the fixtures' own motif repeated to that length rather than
#: a sequence invented to look real.
MOTIF_SEQUENCE_LENGTH = 121
MOTIF_SEQUENCE = ("GGCCACCACCCTG" * 10)[:MOTIF_SEQUENCE_LENGTH]

#: The sample this report is about. Distinctive enough that finding it in the
#: extracted text cannot be an accident.
SAMPLE_NAME = "PRINTCHECK_042"

#: The rest of the identity the running header carries, and what the run recorded for
#: each. They are constants because every one of them is asserted **per page**.
PIPELINE_VERSION = "9.9.9"
ASSEMBLY = "hg19"
RUN_STARTED = "2020-01-02T03:04:05"

#: How the run time reads once ``format_run_timestamp`` has rendered it, with the label
#: that keeps it from being mistaken for the time the report was made. Those are two
#: different facts and the report is at pains to keep them apart.
RUN_TIME_LABEL = "Run 2020-01-02 03:04:05 UTC"

#: A log line that must never reach the paper: the log prints as a one-line pointer
#: to the HTML original, not as a wall of DEBUG output.
DEBUG_LOG_LINE = "DEBUG kestrel_genotyping: writing kestrel_result.tsv"

#: The sentence that replaces it.
LOG_POINTER = "The full pipeline log is in the HTML original of this report."

#: Three Kestrel rows. Each depth score is unique and printed to six decimal places
#: by the server, so counting them counts rows that actually reached the page.
KESTREL_ROWS: tuple[dict[str, Any], ...] = (
    {
        "Motif": "5",
        "Variant": "Insertion",
        "POS": 67,
        "REF": "G",
        "ALT": "GG",
        "Motif_sequence": MOTIF_SEQUENCE,
        "Estimated_Depth_AlternateVariant": 120,
        "Estimated_Depth_Variant_ActiveRegion": 12000,
        "Depth_Score": 0.010012,
        "Confidence": "High_Precision",
        "Flag": "Not flagged",
    },
    {
        "Motif": "6",
        "Variant": "Insertion",
        "POS": 68,
        "REF": "C",
        "ALT": "CC",
        "Motif_sequence": MOTIF_SEQUENCE,
        "Estimated_Depth_AlternateVariant": 41,
        "Estimated_Depth_Variant_ActiveRegion": 5100,
        "Depth_Score": 0.008034,
        "Confidence": "Low_Precision",
        "Flag": "Low_Depth",
    },
    {
        "Motif": "7",
        "Variant": "Deletion",
        "POS": 69,
        "REF": "GG",
        "ALT": "G",
        "Motif_sequence": MOTIF_SEQUENCE,
        "Estimated_Depth_AlternateVariant": 18,
        "Estimated_Depth_Variant_ActiveRegion": 3000,
        "Depth_Score": 0.006001,
        "Confidence": "Low_Precision",
        "Flag": "Depth_Score_Below_Threshold",
    },
)

#: What each row is recognised by in the extracted text, as the server renders it.
ROW_MARKERS = ("0.010012", "0.008034", "0.006001")


def _render_report(output_dir: Path, *, sample_name: str = SAMPLE_NAME) -> Path:
    """Render a report whose evidence is the size real evidence is.

    Args:
        output_dir: Directory to render into.
        sample_name: What the run recorded for itself, explicitly - so the report
            prints it verbatim and a hostile name is not derived away before it can
            be tested.

    Returns:
        Path: The rendered ``summary_report.html``.
    """
    log_file = output_dir / "pipeline.log"
    log_file.write_text(f"{DEBUG_LOG_LINE}\n", encoding="utf-8")
    payload = {
        "version": PIPELINE_VERSION,
        "sample_name": sample_name,
        "sample_name_is_explicit": True,
        "input_files": {"bam": f"{sample_name}.bam"},
        "reference_assembly_requested": ASSEMBLY,
        "pipeline_start": RUN_STARTED,
        "steps": [
            {
                "step": summary_steps.STEP_COVERAGE,
                "parsed_result": {"comments": [], "data": [COVERAGE_ROW]},
            },
            {
                "step": summary_steps.STEP_KESTREL,
                "parsed_result": {"comments": [], "data": list(KESTREL_ROWS)},
            },
        ],
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(payload), encoding="utf-8")

    generate_summary_report(
        output_dir=str(output_dir),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=str(log_file),
        config=load_config(None),
    )
    return output_dir / "summary_report.html"


@pytest.fixture
def printable_report(tmp_path: Path) -> Path:
    """The specimen every case below prints.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``summary_report.html``.
    """
    return _render_report(tmp_path)


def _printed_pages(
    browser: Browser,
    report: Path,
    *,
    javascript: bool = True,
    before: str = "",
    click: str = "",
) -> list[str]:
    """Print a report to PDF and return the text of each page.

    Args:
        browser: Playwright browser, supplied by ``pytest-playwright``.
        report: The rendered report to open.
        javascript: Whether the page may run scripts. False is the no-JS case: a
            reader with scripting off, and the only way to see what the print
            stylesheet does on its own.
        before: JavaScript evaluated after load and before printing, used to put the
            page into the state a reader might print it from. Needs ``javascript``.
        click: A selector clicked after load and before printing. This is the *reader's*
            way into the same states, and it is the only one available with scripting
            off: toggling a ``<details>`` from its ``<summary>`` is user-agent behaviour,
            not script.

    Returns:
        list[str]: One string per printed page, whitespace as ``pdftotext`` gives it.
    """
    if browser.browser_type.name != "chromium":
        pytest.skip("printing to PDF is a Chromium capability; see this module's docstring")
    if shutil.which("pdftotext") is None:
        pytest.skip("poppler's pdftotext is not installed, so the printed artefact cannot be read back")

    context = browser.new_context(java_script_enabled=javascript)
    try:
        page = context.new_page()
        page.goto(report.as_uri(), wait_until="load")
        if before:
            page.evaluate(before)
        if click:
            page.click(click)
        pdf = report.with_suffix(".pdf")
        page.pdf(path=str(pdf), prefer_css_page_size=True, print_background=True)
    finally:
        context.close()

    # `-raw` keeps each drawn run together instead of sorting the page by physical
    # position. That matters for exactly one assertion and it is the important one:
    # a wrapped 121 bp cell is drawn as a dozen short lines, and in the default mode
    # poppler interleaves the neighbouring columns' single lines between them by
    # y-coordinate - so the sequence is on the page and is still not a substring of
    # the extraction. Measured both ways before choosing this one.
    extracted = subprocess.run(
        ["pdftotext", "-raw", str(pdf), "-"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return extracted.split("\f")[:-1] if extracted.endswith("\f") else extracted.split("\f")


def _squashed(text: str) -> str:
    """Return text with every whitespace character removed.

    A cell wrapped across two printed lines arrives with a newline inside it, so a
    121-character sequence is never one substring of the extraction. Removing the
    whitespace is what makes "the whole sequence is on the page" answerable.

    Args:
        text: Extracted PDF text.

    Returns:
        str: The same text with no whitespace at all.
    """
    return "".join(text.split())


def test_the_specimen_carries_the_sequence_length_real_data_carries() -> None:
    """The measurement the print budget is derived from, pinned where it is stated.

    Every assertion below about the sequence surviving the page is worth exactly as
    much as this constant is long. A fixture quietly shortened to 13 characters would
    pass every one of them against an unmodified template.
    """
    assert len(MOTIF_SEQUENCE) == MOTIF_SEQUENCE_LENGTH == 121
    assert " " not in MOTIF_SEQUENCE, "a sequence with a space in it has a break opportunity real data does not"


def test_the_printed_pdf_is_a_complete_record(printable_report: Path, browser: Browser) -> None:
    """Everything the reader is owed reaches the paper, and the log does not.

    Renderer: Chromium 151.0.7922.34 via ``page.pdf(prefer_css_page_size=True)``;
    text read back with poppler ``pdftotext`` 26.01.0. See the module docstring on
    why the engine is named.
    """
    pages = _printed_pages(browser, printable_report)
    text = "\n".join(pages)
    squashed = _squashed(text)

    assert SAMPLE_NAME in squashed, "the printed record does not say whose report it is"
    assert MOTIF_SEQUENCE in squashed, "the 121 bp motif sequence is still truncated in print"
    for marker in ROW_MARKERS:
        assert squashed.count(marker) == 1, f"the variant row identified by {marker} did not print exactly once"
    assert "Researchuseonly" in squashed, "the printed record makes no statement about its use"
    assert "Page1of" in squashed, "the printed pages are not numbered"
    assert _squashed(LOG_POINTER) in squashed, "the printed log points nowhere"
    assert _squashed(DEBUG_LOG_LINE) not in squashed, "the pipeline log printed as a wall of DEBUG output"


def test_every_printed_page_is_numbered_and_says_whose_report_it_is(printable_report: Path, browser: Browser) -> None:
    """A sheet on its own must still say what it is a sheet of.

    Page 2 of a separated or rescanned printout used to carry no identity at all: the
    header was ordinary document content, so it appeared once, and the ``@page``
    margin boxes carried only ``Page N of M``. The identity is now in the margin
    boxes too, which is why the escaper exists - Chromium 151 drops a margin-box
    ``content`` list containing ``string()``, so the document cannot hand the margin a
    value and the only route is writing it into the stylesheet.

    Every element of the running header is asserted on **page 2 specifically**, not on
    the concatenation: a page-1-only header passes the concatenated form, which is how
    the defect survived the earlier version of this case. The specimen is asserted to
    be longer than one page for the same reason.
    """
    pages = _printed_pages(browser, printable_report)

    assert len(pages) >= 2, f"the specimen printed on {len(pages)} page(s), so a running header proves nothing"
    for number, page_text in enumerate(pages, start=1):
        squashed = _squashed(page_text)
        assert f"Page{number}of{len(pages)}" in squashed, f"page {number} is not numbered"
        assert SAMPLE_NAME in squashed, f"page {number} does not say whose report it is"
        assert _squashed(ASSAY_NAME) in squashed, f"page {number} does not say what was assayed"
        assert "hg19" in squashed, f"page {number} does not say which assembly the run asked for"
        assert "VNtyper9.9.9" in squashed, f"page {number} does not say which version produced it"
        assert _squashed(RUN_TIME_LABEL) in squashed, f"page {number} does not say when the run started"
        assert _squashed(RESEARCH_USE_STATEMENT) in squashed, f"page {number} makes no statement about its use"


def test_a_hostile_sample_name_cannot_escape_the_running_header(tmp_path: Path, browser: Browser) -> None:
    """The identity is written into a stylesheet, so the escaper is load-bearing here.

    A sample name is a file name, and a file name is not a controlled vocabulary. This
    renders a real report for a sample whose name tries to close the CSS string, the
    declaration, the rule and the ``<style>`` element, then prints it and reads the
    paper back: the name must arrive in the margin as text, and the page must still be
    a page - numbered, with its tables on it - rather than a document that fell apart
    at the point the name was interpolated.
    """
    hostile = 'S1" ; } </style><script>x</script> {'
    report = _render_report(tmp_path, sample_name=hostile)
    pages = _printed_pages(browser, report, javascript=False)
    squashed = _squashed("\n".join(pages))

    assert _squashed(hostile) in squashed, "the hostile name did not reach the margin as text"
    assert "Page1of" in squashed, "the page counter is gone, so the stylesheet was broken by the name"
    assert MOTIF_SEQUENCE in squashed, "the document did not survive the interpolated name"
    assert "<script>" not in report.read_text(encoding="utf-8").split("</head>")[0], (
        "the name reopened a script element inside the document head"
    )


def test_a_section_the_reader_collapsed_still_prints(printable_report: Path, browser: Browser) -> None:
    """Printing is not a reading of what happens to be open.

    A reader collapses a section to get it out of the way and then prints; the
    archived record must not be missing it. ``page.pdf()`` dispatches ``beforeprint``,
    so this exercises the same handler Ctrl+P does. The log is the deliberate
    exception and is checked in the same pass, because "open everything" and "open
    everything except the log" differ only there.

    The collapse is a **click on the summary**, not an assignment to ``.open``: that is
    the reader's own gesture, it is user-agent behaviour rather than script, and it is
    the one thing this case and
    :func:`test_a_disclosure_the_reader_collapsed_does_not_print_without_scripting`
    have in common - so the two differ in exactly one variable.
    """
    collapsed = _printed_pages(browser, printable_report, click="#variantsToggle")
    squashed = _squashed("\n".join(collapsed))

    assert "IGVvarianttable" in squashed or "novariant" in squashed.lower() or "Variants" in squashed
    assert _squashed("This run produced no IGV variant table.") in squashed, (
        "the collapsed disclosure's body did not print, so a reader's collapse is a hole in the record"
    )
    assert _squashed(DEBUG_LOG_LINE) not in squashed, "the log was opened along with everything else"


def test_a_log_the_reader_expanded_still_prints_as_a_pointer(printable_report: Path, browser: Browser) -> None:
    """The one case in which the log rule is the only thing standing between a reader
    and a wall of DEBUG output in the archived PDF.

    Every other specimen here leaves ``#logDisclosure`` closed, and a closed disclosure
    is hidden by the UA whether or not the print block says anything - so the rule
    ``.log-section > .details-body { display: none !important; }`` could be deleted with
    all six of these cases and all 34 unit cases still passing. A reader who expanded
    the log to look at it and then printed is what it is for, and it is why this case
    exists rather than a stronger assertion on one of the others.
    """
    pages = _printed_pages(
        browser,
        printable_report,
        before='() => { document.getElementById("logDisclosure").open = true; }',
    )
    squashed = _squashed("\n".join(pages))

    assert _squashed(DEBUG_LOG_LINE) not in squashed, (
        "a log the reader had expanded printed its body, so the archived PDF carries the DEBUG output"
    )
    assert _squashed(LOG_POINTER) in squashed, "the expanded log printed neither its body nor the pointer"


#: Text planted inside a disclosure body by the server, so that "did this disclosure's
#: contents reach the paper" is answerable at all. The shipped report's one open
#: disclosure is filled in by ``initTable()``, so with scripting off it is empty
#: whatever its state - which would make the two cases below indistinguishable and both
#: of them vacuous.
SERVER_RENDERED_EVIDENCE = "SERVERRENDEREDEVIDENCE42"

#: The empty container ``initTable()`` writes into, and where the marker goes instead.
VARIANT_TABLE_CONTAINER = '<div id="tableSelectorDiv"></div>'


def _with_evidence_inside_the_disclosure(report: Path) -> Path:
    """Write a copy of a rendered report with server-rendered text in a disclosure.

    This is specimen preparation, not a fixture of the product: the marker stands for
    any evidence a future change might put inside a ``<details>``. The disclosure, its
    ``open`` attribute and every print rule are the shipped ones.

    Args:
        report: The rendered ``summary_report.html``.

    Returns:
        Path: A sibling file whose Variants disclosure has a body the server wrote.
    """
    html = report.read_text(encoding="utf-8")
    assert VARIANT_TABLE_CONTAINER in html, "the variant-table container moved; this specimen plants nothing"
    # A click on a `<summary>` toggles, so "collapsed" is only what it means if the
    # served state was open. That is the mitigation itself, pinned in the source text
    # by `tests/unit/test_report_presentation.py::
    # test_every_print_relevant_disclosure_is_served_open`, and restated here so the
    # case below cannot quietly become "the reader expanded it".
    assert re.search(r"<details[^>]*\bid=\"variantsDisclosure\"[^>]*\bopen\b", html), (
        "the Variants disclosure is not served open, so clicking its summary expands rather than collapses"
    )
    marked = report.with_name("marked_report.html")
    marked.write_text(
        html.replace(VARIANT_TABLE_CONTAINER, f'<div id="tableSelectorDiv"><p>{SERVER_RENDERED_EVIDENCE}</p></div>'),
        encoding="utf-8",
    )
    return marked


def test_a_disclosure_is_served_open_so_a_reader_with_no_script_gets_its_contents(
    printable_report: Path, browser: Browser
) -> None:
    """The mitigation, measured: printing works with no script because nothing is closed.

    Every print-relevant disclosure carries ``open`` in the HTML the server writes, so
    the reachable no-script case is the one that works. Collapsing is then something a
    reader has to choose, rather than the state the file arrives in.
    """
    pages = _printed_pages(browser, _with_evidence_inside_the_disclosure(printable_report), javascript=False)

    assert SERVER_RENDERED_EVIDENCE in _squashed("\n".join(pages)), (
        "a disclosure the server left open did not print its contents with scripting off"
    )


def test_a_disclosure_the_reader_collapsed_does_not_print_without_scripting(
    printable_report: Path, browser: Browser
) -> None:
    """The residual limitation, stated rather than papered over.

    A reader with scripting off clicks a ``<summary>`` - which is user-agent behaviour
    and needs no script - and prints. **Chromium 151 does not put that section's
    contents on the paper.** The ``@media print`` rule ``details > summary ~ *
    { display: block !important }`` does not reach them: the UA hides a closed
    disclosure's content in a way no author rule overrides, and ``open`` is an attribute
    rather than a style, so no stylesheet can set it. The ``beforeprint`` handler that
    reopens collapsed sections is the fix, and it is exactly what this reader does not
    have.

    So the archived PDF is short a section, and only the ``<summary>`` line says one was
    there. Both halves are asserted, because "the reader can see something is missing"
    is the whole of what survives. The mitigation is the sibling case above: the file is
    served with nothing collapsed, so reaching this state takes a deliberate click.
    ``docs/pipeline/reports.md`` says the same in the operator's words.
    """
    pages = _printed_pages(
        browser,
        _with_evidence_inside_the_disclosure(printable_report),
        javascript=False,
        click="#variantsToggle",
    )
    squashed = _squashed("\n".join(pages))

    assert SERVER_RENDERED_EVIDENCE not in squashed, (
        "Chromium now prints a reader-collapsed disclosure without script - the limitation "
        "documented here and in docs/pipeline/reports.md is gone, so update both"
    )
    assert "Variants" in squashed, "not even the collapsed section's heading printed, so nothing marks the gap"
    assert "Page1of" in squashed, "the rest of the print stylesheet stopped working in this state"


def test_the_record_survives_a_reader_with_scripting_off(printable_report: Path, browser: Browser) -> None:
    """What the print stylesheet does on its own, asserted honestly.

    Everything structural - the page numbers, the running header, the un-clamped
    cells, the whole coverage table and the log pointer - is CSS and needs no script.
    Reopening a collapsed ``<details>`` is the one thing that does, in this engine,
    and this pins that so the fallback cannot be quietly believed to do more than it
    does.
    """
    pages = _printed_pages(browser, printable_report, javascript=False)
    squashed = _squashed("\n".join(pages))

    assert MOTIF_SEQUENCE in squashed
    assert "Page1of" in squashed
    assert SAMPLE_NAME in squashed
    assert _squashed(LOG_POINTER) in squashed
    assert _squashed(DEBUG_LOG_LINE) not in squashed


def test_the_printed_record_carries_the_coverage_verdict(printable_report: Path, browser: Browser) -> None:
    """The QC verdict lives only in the view a switch reveals, and paper has no switch.

    ``#detailedCoverageView`` is ``display: none`` until the reader ticks "Show
    detailed coverage", so the printed report used to carry the two headline coverage
    figures and neither the eight statistics nor the verdict.
    """
    squashed = _squashed("\n".join(_printed_pages(browser, printable_report)))

    assert "CoverageQCPASS" in squashed, "the printed record omits the coverage QC verdict"
    assert "UncoveredBases" in squashed, "the printed record omits the detailed coverage statistics"
