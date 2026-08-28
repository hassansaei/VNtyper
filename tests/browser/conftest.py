"""Fixtures for the browser tier: a real report, opened by a real browser engine.

Why this tier exists
--------------------
``tests/unit/test_generate_report.py`` renders the shipped template and asserts on
the HTML that comes out, which is everything the report *is* on the server side.
It cannot see what the report *does*, because the unit tier has no JavaScript
engine (``tests/unit/test_template_escaping.py`` says so in as many words). Three
of the report's behaviours existed only once a browser had run its scripts, and
all three were defects (#242):

* a DataTables custom search predicate hid every row whose ``Flag`` was anything
  other than ``Not flagged``/``Not applicable``/empty, unless the reader ticked a
  toggle they had no reason to know about;
* ``applyRounding`` rewrote every numeric cell to four decimal places *in the
  DOM*, so the figure on screen was not the figure in the file;
* both lived behind jQuery, DataTables and Bootstrap loaded from three CDNs, so a
  reader with no network got neither - and got no warning either.

Together that meant the same file read differently depending on the reader's
network. This tier renders one report and opens it twice, once with the network
reachable and once with everything but ``file://`` blocked, and holds the report
to reading the same either way.

How "was this pass really online?" is answered now
--------------------------------------------------
It used to be answered by asking the page whether jQuery and DataTables had
initialised, and the file it was asked about is the reason: an online pass that
silently became an offline one would have made ``online == offline`` pass against
a report nobody had checked.

The CDN tags are gone (#242), so there is no jQuery left to ask - and the question
has been replaced rather than dropped. :func:`external_requests` reads the
**request log**, and the log is populated from Playwright's ``request`` event,
which fires *before* routing: an attempted fetch is recorded whether the network
would have answered it, the route handler aborted it, or a CDN 404'd it. So "this
document asked for nothing off-machine" is a claim with the same force in both
passes, and it no longer depends on the network having worked.
``test_offline_document.py`` proves the recorder itself is not vacuous, by
pointing a control document at a host and watching the recorder see it.

Purity
------
Rendering goes through the real ``generate_summary_report`` against the real
shipped template and the real ``report_config.json`` - nothing here mocks the
report - but it stays filesystem-only: ``bed_file`` is left unset, so
``create_report`` (igv-reports) is never invoked. The *browser*, on the other
hand, deliberately does reach the network in the online pass. That is the
behaviour under test, and it is why this tier is not the unit tier and is not in
``make check-all``.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Callable, Iterator
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
from playwright.sync_api import Browser, BrowserContext, Page, Request, Response, Route

import vntyper
from vntyper.cli import load_config
from vntyper.scripts import summary_steps
from vntyper.scripts.generate_report import generate_summary_report

logger = logging.getLogger(__name__)

#: The shipped Jinja2 templates, rendered as-is. Same source as the unit tier.
TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"

#: The rows of the Kestrel table, as CSS sees them. The id is ``kestrel_table``
#: (``generate_report.py``: ``to_html(table_id="kestrel_table", ...)``), *not*
#: ``kestrel``: a selector that matches nothing makes both snapshots empty and
#: turns the determinism test green against an unmodified template.
KESTREL_ROW_SELECTOR = "#kestrel_table tbody tr"

#: How many Kestrel rows the fixture puts in the report.
EXPECTED_KESTREL_ROWS = 3


def kestrel_column_text(heading: str) -> str:
    """Return JS that reads one visible Kestrel column, keyed on its heading.

    Positional selectors are what this replaces, and they were wrong the moment
    ``KESTREL_DISPLAY_COLUMNS`` moved ``Motif_sequence`` to the end (#242): the
    ``Flag`` cell stopped being ``td:last-child`` and ``Depth Score`` stopped being
    ``td:nth-child(9)``, and both checks went on running against whatever cell had
    taken the position. Reading the heading row makes the guard say what it means and
    survive the next reorder, which is also what the migration note in
    ``docs/user-guide/output-files.md`` tells anyone scraping the report to do.

    Args:
        heading: The column's heading text, as the report renders it.

    Returns:
        str: A JavaScript expression for ``eval_on_selector_all`` over
        :data:`KESTREL_ROW_SELECTOR`, evaluating to one whitespace-normalised string
        per visible row. A row with no such column yields ``""``.
    """
    return f"""
els => {{
    const table = els.length ? els[0].closest('table') : null;
    const headings = table && table.tHead
        ? Array.from(table.tHead.rows[0].cells).map(cell => cell.textContent.trim())
        : [];
    const index = headings.indexOf({heading!r});
    return els.filter(row => row.offsetParent !== null)
              .map(row => {{
                  const cells = row.querySelectorAll('td');
                  return index >= 0 && cells.length > index
                      ? cells[index].textContent.replace(/\\s+/g, ' ').trim()
                      : '';
              }});
}}
"""


#: A ``Flag`` value the removed client-side predicate treated as clean, so the row
#: carrying it was the one that survived online.
#:
#: These two constants describe the specimen, and a test that asserts on them
#: alone asserts nothing - it would still pass if they were changed to values the
#: predicate never hid. Every guard that depends on which side of the predicate a
#: value falls therefore spells the allowlist out as a literal of its own and
#: cross-checks these against it.
CLEAN_FLAG = "Not flagged"

#: ``Flag`` values the removed client-side predicate hid, so the rows carrying them
#: disappeared online while staying in the file. No shipped configuration produces
#: one - ``report_config.json`` declares no ``flag_rules`` key at all - so the
#: fixture states them directly rather than hoping a rule fires.
HIDDEN_FLAGS = ("Low_Depth", "Depth_Score_Below_Threshold")

#: Three Kestrel rows: one clean, two flagged. Every ``Depth_Score`` carries more
#: than four decimal places, which is what made the removed ``applyRounding``
#: observable (``0.010012`` in the file, ``0.01`` on screen).
KESTREL_ROWS: tuple[dict[str, Any], ...] = (
    {
        "Motifs": "X-5",
        "Motif": "5",
        "Variant": "Insertion",
        "POS": 67,
        "REF": "G",
        "ALT": "GG",
        "Motif_sequence": "GGCCACCACCCTG",
        "Estimated_Depth_AlternateVariant": 120,
        "Estimated_Depth_Variant_ActiveRegion": 12000,
        "Depth_Score": 0.010012,
        "Confidence": "High_Precision",
        "Flag": CLEAN_FLAG,
    },
    {
        "Motifs": "X-6",
        "Motif": "6",
        "Variant": "Insertion",
        "POS": 68,
        "REF": "C",
        "ALT": "CC",
        "Motif_sequence": "GGCCACCACCCTG",
        "Estimated_Depth_AlternateVariant": 41,
        "Estimated_Depth_Variant_ActiveRegion": 5100,
        "Depth_Score": 0.008034,
        "Confidence": "Low_Precision",
        "Flag": HIDDEN_FLAGS[0],
    },
    {
        "Motifs": "X-7",
        "Motif": "7",
        "Variant": "Deletion",
        "POS": 69,
        "REF": "GG",
        "ALT": "G",
        "Motif_sequence": "GGCCACCACCCTG",
        "Estimated_Depth_AlternateVariant": 18,
        "Estimated_Depth_Variant_ActiveRegion": 3000,
        "Depth_Score": 0.006001,
        "Confidence": "Low_Precision",
        "Flag": HIDDEN_FLAGS[1],
    },
)

#: A well-covered sample, so nothing in the coverage section competes for the
#: reader's attention with the table this tier is about.
COVERAGE_ROW: dict[str, Any] = {
    "mean": 250.0,
    "median": 248.0,
    "stdev": 12.5,
    "min": 100,
    "max": 400,
    "region_length": 1000,
    "uncovered_bases": 5,
    "percent_uncovered": 0.5,
}

#: Collect the visible rows of the Kestrel table and normalise their whitespace.
#: The removed DataTables predicate took filtered rows out of the DOM rather than
#: hiding them, so the ``offsetParent`` check is belt and braces against a future
#: filter that hides them with CSS instead.
_VISIBLE_ROW_TEXT = """
els => els.filter(e => e.offsetParent !== null)
          .map(e => e.textContent.replace(/\\s+/g, ' ').trim())
"""

#: The two libraries the report used to load from a CDN. Nothing may put them back:
#: this tier's whole subject is that the document stands on its own, and a page that
#: has jQuery in it again has a CDN tag in it again.
_REMOVED_LIBRARIES = """() => ({
    jquery: (window.jQuery && window.jQuery.fn) ? window.jQuery.fn.jquery : null,
    dataTable: !!(window.jQuery && window.jQuery.fn && window.jQuery.fn.dataTable),
    bootstrap: !!(window.bootstrap || (window.jQuery && window.jQuery.fn && window.jQuery.fn.tooltip)),
})"""


def _tabular_step(name: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Shape a tsv/csv-derived step the way ``summary.py`` records it.

    Args:
        name: The step name, matched verbatim by ``generate_report.py``.
        rows: The step's parsed rows.

    Returns:
        dict: The step mapping as it appears in ``pipeline_summary.json``.
    """
    return {"step": name, "parsed_result": {"comments": [], "data": rows}}


def _block_everything_but_local_files(route: Route) -> None:
    """Abort any request that is not a ``file://`` read.

    Offline is the whole point of this tier, so it is enforced here rather than
    trusted: an unplugged network is not reproducible and a warm HTTP cache would
    quietly make the offline pass online.

    Args:
        route: The intercepted route.
    """
    if route.request.url.startswith("file://"):
        route.continue_()
    else:
        route.abort()


@pytest.fixture
def rendered_report(tmp_path: Path) -> Path:
    """Render a real report from a specimen ``pipeline_summary.json``.

    The specimen is chosen to be able to exhibit the defect: three Kestrel rows,
    two of them flagged with values the shipped client-side filter hides, and
    depth scores with more than four decimal places.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``summary_report.html``.
    """
    payload = {
        "version": "9.9.9",
        "input_files": {"bam": "sample.bam"},
        "steps": [
            _tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
            _tabular_step(summary_steps.STEP_KESTREL, list(KESTREL_ROWS)),
        ],
    }
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(payload), encoding="utf-8")

    generate_summary_report(
        output_dir=str(tmp_path),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=None,
        config=load_config(None),
    )
    return tmp_path / "summary_report.html"


@pytest.fixture
def rendered_report_with_custom_fastp_cutoffs(tmp_path: Path) -> Path:
    """Render a report whose fastp values sit exactly on custom cutoffs.

    The values are intentionally unlike the shipped defaults. A correct report
    must show these labels and judge every value as passing, so the test using
    this fixture detects a label and icon decision taking their cutoffs from
    different places.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``summary_report.html``.
    """
    payload = {
        "version": "9.9.9",
        "input_files": {"bam": "sample.bam"},
        "steps": [
            _tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
            _tabular_step(summary_steps.STEP_KESTREL, list(KESTREL_ROWS)),
        ],
    }
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(payload), encoding="utf-8")
    fastp_dir = tmp_path / "fastq_bam_processing"
    fastp_dir.mkdir()
    (fastp_dir / "output.json").write_text(
        json.dumps(
            {
                "summary": {
                    "before_filtering": {"total_reads": 10000},
                    "after_filtering": {"q20_rate": 0.7555, "q30_rate": 0.6543},
                },
                "duplication": {"rate": 0.1234},
                "filtering_result": {"passed_filter_reads": 7765},
            }
        ),
        encoding="utf-8",
    )
    config = deepcopy(load_config(None))
    config["thresholds"].update(
        {
            "duplication_rate": 0.1234,
            "q20_rate": 0.7555,
            "q30_rate": 0.6543,
            "passed_filter_reads_rate": 0.7765,
        }
    )

    generate_summary_report(
        output_dir=str(tmp_path),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=None,
        config=config,
    )
    return tmp_path / "summary_report.html"


@pytest.fixture
def rendered_report_with_fastp_half_ties(tmp_path: Path) -> Path:
    """Render every fastp metric at a Decimal half-tie configured cutoff.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``summary_report.html``.
    """
    payload = {
        "version": "9.9.9",
        "input_files": {"bam": "sample.bam"},
        "steps": [
            _tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
            _tabular_step(summary_steps.STEP_KESTREL, list(KESTREL_ROWS)),
        ],
    }
    (tmp_path / "pipeline_summary.json").write_text(json.dumps(payload), encoding="utf-8")
    fastp_dir = tmp_path / "fastq_bam_processing"
    fastp_dir.mkdir()
    (fastp_dir / "output.json").write_text(
        json.dumps(
            {
                "summary": {
                    "before_filtering": {"total_reads": 100000},
                    "after_filtering": {"q20_rate": 0.77645, "q30_rate": 0.70045},
                },
                "duplication": {"rate": 0.05045},
                "filtering_result": {"passed_filter_reads": 80045},
            }
        ),
        encoding="utf-8",
    )
    config = deepcopy(load_config(None))
    config["thresholds"].update(
        {
            "duplication_rate": 0.05045,
            "q20_rate": 0.77645,
            "q30_rate": 0.70045,
            "passed_filter_reads_rate": 0.80045,
        }
    )

    generate_summary_report(
        output_dir=str(tmp_path),
        template_dir=str(TEMPLATE_DIR),
        report_file="summary_report.html",
        log_file=None,
        config=config,
    )
    return tmp_path / "summary_report.html"


#: Every URL each page asked for, keyed by page. Module-level rather than a fixture so
#: :func:`external_requests` can be called as a plain function from a test that already
#: has the page - the same shape ``removed_libraries`` has, and one less fixture to
#: thread through every signature.
_REQUEST_LOG: dict[Page, list[str]] = {}


#: A page in the shape ``create_report`` writes: the container marker, the two
#: JavaScript literals, and a body end. ``run_igv_report`` is stubbed out with this
#: because igv-reports needs a real BAM and a real reference, and this tier has
#: neither - but everything downstream of the stub is the shipped code path, including
#: the fragment extraction and the "is there a session" decision that gates the
#: half-megabyte payload.
IGV_REPORTS_PAGE = (
    "<html><body>\n"
    '<div id="container">\n<div id="igvDiv"></div>\n</div>\n<script>\n'
    'const tableJson = {"headers": ["unique_id", "CHROM", "POS", "REF", "ALT"], '
    '"rows": [["0", "chr1", "155188205", "G", "GG"], ["1", "chr1", "155188305", "GG", "G"]]}\n'
    'const sessionDictionary = {"0": "session0.json", "1": "session1.json"}\n'
    "</script>\n</body></html>\n"
)

#: The two spliced variants, as the report renders them. Named here so the tests that
#: read the variant table do not each restate the fixture.
IGV_VARIANT_IDS = ("0", "1")


@pytest.fixture
def rendered_report_with_alignments(tmp_path: Path, monkeypatch) -> Callable[..., Path]:
    """Return a callable that renders a report carrying a real alignment session.

    This is the shape that costs 497 KB and the only one in which the embedded library
    is written at all, so every test of the alignment panel needs it. It replaces the
    older approach of splicing the two JavaScript literals into an already-rendered
    document: a splice cannot produce the *markup* the template authors from
    ``igv_session_available``, so it would have tested the script against a page the
    generator never writes.

    Args:
        tmp_path: Pytest's per-test temporary directory.
        monkeypatch: Pytest's patcher, used to stand in for igv-reports.

    Yields:
        Callable: ``render(mode="embedded") -> Path`` to the rendered report.
    """
    from vntyper.scripts import generate_report

    def _stub(bed_file, bam_file, fasta_file, output_html, **kwargs) -> None:
        Path(output_html).write_text(IGV_REPORTS_PAGE, encoding="utf-8")

    monkeypatch.setattr(generate_report, "run_igv_report", _stub)

    def _render(mode: str = "embedded") -> Path:
        run = tmp_path / mode
        run.mkdir(exist_ok=True)
        payload = {
            "version": "9.9.9",
            "input_files": {"bam": "sample.bam"},
            "steps": [
                _tabular_step(summary_steps.STEP_COVERAGE, [COVERAGE_ROW]),
                _tabular_step(summary_steps.STEP_KESTREL, list(KESTREL_ROWS)),
            ],
        }
        (run / "pipeline_summary.json").write_text(json.dumps(payload), encoding="utf-8")
        bed = run / "regions.bed"
        bed.write_text("chr1\t155188100\t155188300\n", encoding="utf-8")

        generate_summary_report(
            output_dir=str(run),
            template_dir=str(TEMPLATE_DIR),
            report_file="summary_report.html",
            log_file=None,
            config=load_config(None),
            bed_file=str(bed),
            report_igv=mode,
        )
        report = run / "summary_report.html"
        logger.info("Rendered %s (--report-igv %s): %d bytes", report, mode, report.stat().st_size)
        return report

    return _render


@pytest.fixture
def failed_requests() -> dict[Page, list[str]]:
    """Per-page record of the requests that failed at the **transport** layer.

    This is a necessary but nowhere near sufficient signal that a pass was online:
    ``requestfailed`` fires for a dropped connection or a blocked route, and not
    for an HTTP 404 or 503. Pair it with :data:`error_responses` and
    :func:`external_requests`, which sees the cases they cannot.

    Returns:
        dict: Page -> the URLs that failed, filled in by :func:`open_report`.
    """
    return {}


@pytest.fixture
def error_responses() -> dict[Page, list[str]]:
    """Per-page record of the responses that came back with an error status.

    A CDN that answers 404 or 503 completes the request as far as the transport
    layer is concerned, so :data:`failed_requests` stays empty while the page is
    left unenhanced - which makes an online pass indistinguishable from an offline
    one. This is the half of the evidence that sees those.

    Returns:
        dict: Page -> ``"<status> <url>"`` for every response at 400 or above.
    """
    return {}


@pytest.fixture
def page_errors() -> dict[Page, list[str]]:
    """Per-page JavaScript exceptions raised while the document executes."""
    return {}


@pytest.fixture
def open_report(
    browser: Browser,
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
    page_errors: dict[Page, list[str]],
) -> Iterator[Callable[..., Page]]:
    """Return a callable that opens a report in a browser, online or offline.

    Each call gets its **own browser context**. Two independent page loads are
    the whole point - if the offline pass left its route handler, or its empty
    HTTP cache, on a page the online pass then reused, the online pass would be
    offline too and the comparison would compare nothing.

    Args:
        browser: Playwright browser, supplied by ``pytest-playwright``.
        failed_requests: Registry to record each page's transport failures in.
        error_responses: Registry to record each page's error statuses in.
        page_errors: Registry to record JavaScript exceptions in each page.

    Yields:
        Callable: ``open_report(path, *, offline=False) -> Page``.
    """
    contexts: list[BrowserContext] = []

    def _open(path: Path, *, offline: bool = False, before_load: str | None = None) -> Page:
        context = browser.new_context()
        contexts.append(context)
        page = context.new_page()
        failures: list[str] = []
        errors: list[str] = []
        requested: list[str] = []
        failed_requests[page] = failures
        error_responses[page] = errors
        page_errors[page] = []
        _REQUEST_LOG[page] = requested
        # Registered before `goto`, and before the route handler: `request` fires for
        # every attempt, aborted or not, which is what makes the log evidence about the
        # *document* rather than about the network.
        page.on("request", lambda request: requested.append(request.url))
        page.on("requestfailed", lambda request: failures.append(_describe_failure(request)))
        page.on("response", lambda response: _record_error_response(response, errors))
        page.on("pageerror", lambda error: page_errors[page].append(str(error)))
        if offline:
            page.route("**/*", _block_everything_but_local_files)
        if before_load is not None:
            # Runs in a fresh document before any of the page's own script does, which
            # is the only way to take a global away from the page under test.
            page.add_init_script(before_load)
        page.goto(path.as_uri(), wait_until="load")
        logger.info(
            "Opened %s (offline=%s); %d transport failure(s) %s; %d error response(s) %s",
            path,
            offline,
            len(failures),
            failures,
            len(errors),
            errors,
        )
        return page

    yield _open

    for context in contexts:
        context.close()
    _REQUEST_LOG.clear()


def _describe_failure(request: Request) -> str:
    """Render one failed request as ``<reason> <url>``.

    Args:
        request: The request the browser could not complete.

    Returns:
        str: A one-line description for an assertion message.
    """
    return f"{request.failure} {request.url}"


def _record_error_response(response: Response, errors: list[str]) -> None:
    """Record a response whose status says the browser did not get what it asked for.

    Args:
        response: The response the page received.
        errors: The page's registry, appended to in place.
    """
    if response.status >= 400:
        errors.append(f"{response.status} {response.url}")


def removed_libraries(page: Page) -> dict[str, Any]:
    """Report whether any of the three removed third-party libraries is back.

    Args:
        page: An already-loaded page.

    Returns:
        dict: ``jquery`` (the version string, or ``None``), ``dataTable`` and
        ``bootstrap`` (whether either is registered).
    """
    state: dict[str, Any] = page.evaluate(_REMOVED_LIBRARIES)
    logger.info("Removed-library state: %s", state)
    return state


def external_requests(page: Page) -> list[str]:
    """Every request this page issued that was not a local file read.

    The registry behind this is filled from Playwright's ``request`` event, which
    fires before the route handler runs - so an attempted fetch is recorded whether
    it succeeded, was aborted by the offline route, or came back 404. That is what
    makes "the document asked for nothing off-machine" a claim worth making in
    either pass rather than a restatement of "the network was unplugged".

    Args:
        page: A page opened through :func:`open_report`.

    Returns:
        list[str]: The non-``file://`` URLs, in the order they were requested.
    """
    requested = _REQUEST_LOG.get(page, [])
    external = [url for url in requested if not url.startswith("file://")]
    logger.info("%d request(s), %d of them off-machine: %s", len(requested), len(external), external)
    return external


@pytest.fixture
def visible_kestrel_rows() -> Callable[[Page], list[str]]:
    """Return a callable reading the Kestrel rows a reader can actually see.

    Returns:
        Callable: ``visible_kestrel_rows(page) -> list[str]``, one whitespace-
        normalised string per visible row.
    """

    def _rows(page: Page) -> list[str]:
        rows: list[str] = page.eval_on_selector_all(KESTREL_ROW_SELECTOR, _VISIBLE_ROW_TEXT)
        logger.info("Visible Kestrel rows (%d): %s", len(rows), rows)
        return rows

    return _rows
