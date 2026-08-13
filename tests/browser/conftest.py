"""Fixtures for the browser tier: a real report, opened by a real browser engine.

Why this tier exists
--------------------
``tests/unit/test_generate_report.py`` renders the shipped template and asserts on
the HTML that comes out, which is everything the report *is* on the server side.
It cannot see what the report *does*, because the unit tier has no JavaScript
engine (``tests/unit/test_template_escaping.py`` says so in as many words). Three
of the report's behaviours only exist once a browser has run its scripts:

* the DataTables custom search filter hides every row whose ``Flag`` is anything
  other than ``Not flagged``/``Not applicable``/empty, unless the reader ticks a
  toggle they have no reason to know about;
* ``applyRounding`` rewrites every numeric cell to four decimal places *in the
  DOM*, so the figure on screen is not the figure in the file;
* both of those live behind jQuery, DataTables and Bootstrap loaded from three
  CDNs, so a reader with no network gets neither - and gets no warning either.

Together that means the same file reads differently depending on the reader's
network, which is issue #242. This tier renders one report and opens it twice,
once with the network reachable and once with everything but ``file://`` blocked.

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
from playwright.sync_api import Browser, BrowserContext, Page, Request, Route

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

#: A ``Flag`` value the shipped client-side filter treats as clean, so the row
#: carrying it survives online.
CLEAN_FLAG = "Not flagged"

#: ``Flag`` values the shipped client-side filter hides, so the rows carrying them
#: disappear online while staying in the file. No shipped configuration produces
#: one - ``report_config.json`` declares no ``flag_rules`` key at all - so the
#: fixture states them directly rather than hoping a rule fires.
HIDDEN_FLAGS = ("Low_Depth", "Depth_Score_Below_Threshold")

#: Three Kestrel rows: one clean, two flagged. Every ``Depth_Score`` carries more
#: than four decimal places, which is what makes the browser's ``applyRounding``
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
#: DataTables removes filtered rows from the DOM rather than hiding them, so the
#: ``offsetParent`` check is belt and braces against a future filter that hides
#: them with CSS instead.
_VISIBLE_ROW_TEXT = """
els => els.filter(e => e.offsetParent !== null)
          .map(e => e.textContent.replace(/\\s+/g, ' ').trim())
"""


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
    """Per-page record of the requests the browser could not complete.

    An online pass run on a machine with no network is indistinguishable from an
    offline pass, and would make the determinism check pass while proving
    nothing. This is how a test tells the difference.

    Returns:
        dict: Page -> the URLs that failed, filled in by :func:`open_report`.
    """
    return {}


@pytest.fixture
def open_report(browser: Browser, failed_requests: dict[Page, list[str]]) -> Iterator[Callable[..., Page]]:
    """Return a callable that opens a report in a browser, online or offline.

    Each call gets its **own browser context**. Two independent page loads are
    the whole point - if the offline pass left its route handler, or its empty
    HTTP cache, on a page the online pass then reused, the online pass would be
    offline too and the comparison would compare nothing.

    Args:
        browser: Playwright browser, supplied by ``pytest-playwright``.
        failed_requests: Registry to record each page's failed requests in.

    Yields:
        Callable: ``open_report(path, *, offline=False) -> Page``.
    """
    contexts: list[BrowserContext] = []

    def _open(path: Path, *, offline: bool = False) -> Page:
        context = browser.new_context()
        contexts.append(context)
        page = context.new_page()
        failures: list[str] = []
        failed_requests[page] = failures
        page.on("requestfailed", lambda request: failures.append(_describe_failure(request)))
        if offline:
            page.route("**/*", _block_everything_but_local_files)
        page.goto(path.as_uri(), wait_until="load")
        logger.info("Opened %s (offline=%s); %d request(s) failed: %s", path, offline, len(failures), failures)
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
