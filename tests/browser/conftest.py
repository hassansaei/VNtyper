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
to reading the same either way. The first two behaviours are gone; the CDN tags
are still there, which is why the online pass still has to prove it was online.

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

#: Whether the third-party scripts the report still loads actually initialised.
#: ``requestfailed`` cannot answer this: a 404, a 503 and a 200 carrying unusable
#: JavaScript all leave the transport layer happy and the page unenhanced.
_ENHANCEMENT_STATE = """() => ({
    jquery: (window.jQuery && window.jQuery.fn) ? window.jQuery.fn.jquery : null,
    dataTable: !!(window.jQuery && window.jQuery.fn && window.jQuery.fn.dataTable),
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
def failed_requests() -> dict[Page, list[str]]:
    """Per-page record of the requests that failed at the **transport** layer.

    This is a necessary but nowhere near sufficient signal that a pass was online:
    ``requestfailed`` fires for a dropped connection or a blocked route, and not
    for an HTTP 404 or 503. Pair it with :data:`error_responses` and
    :func:`enhancement_state`, which see the cases it cannot.

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
def open_report(
    browser: Browser,
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
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

    Yields:
        Callable: ``open_report(path, *, offline=False) -> Page``.
    """
    contexts: list[BrowserContext] = []

    def _open(path: Path, *, offline: bool = False) -> Page:
        context = browser.new_context()
        contexts.append(context)
        page = context.new_page()
        failures: list[str] = []
        errors: list[str] = []
        failed_requests[page] = failures
        error_responses[page] = errors
        page.on("requestfailed", lambda request: failures.append(_describe_failure(request)))
        page.on("response", lambda response: _record_error_response(response, errors))
        if offline:
            page.route("**/*", _block_everything_but_local_files)
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


def enhancement_state(page: Page) -> dict[str, Any]:
    """Report whether the page's third-party scripts actually initialised.

    The report still loads jQuery and DataTables from CDNs, and every one of the
    behaviours this tier compares depends on them having run. Asking the page
    itself is the only way to know: a script served as a 200 of unusable
    JavaScript leaves no trace in either request registry.

    Args:
        page: An already-loaded page.

    Returns:
        dict: ``jquery`` (the version string, or ``None``) and ``dataTable``
        (whether the DataTables plugin is registered on jQuery).
    """
    state: dict[str, Any] = page.evaluate(_ENHANCEMENT_STATE)
    logger.info("Third-party enhancement state: %s", state)
    return state


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
