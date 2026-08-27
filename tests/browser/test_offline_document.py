"""The rendered report resolves nothing off-machine, measured in a browser.

Why a browser and not a grep
----------------------------
A ``src``/``href`` scan over the source is not sufficient and never was. It misses
CSS ``url()``, ``@import``, ``srcset``, anything a script fetches, and the session
URLs igv.js loads - and the cohort report makes the point on its own: Plotly's inlined
bundle carries dozens of ``https://`` string literals (basemap tile templates, licence
links, XML namespaces) that it never requests for a donut chart, so a substring scan of
that file reports 67 "remote references" and every one of them is a false positive.

The browser is where an unapproved request is *observable*. Playwright's ``request``
event fires for every attempt before routing, so this sees a fetch whether it would
have succeeded, was aborted, or came back 404.

Eleven CDN tags existed across the two templates before this change, over four origins,
with no ``integrity`` attribute and no fallback of any kind.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
from playwright.sync_api import Page, expect

from tests.browser.conftest import external_requests, removed_libraries
from vntyper.cli import load_config
from vntyper.scripts.cohort_summary import generate_cohort_summary_report

pytestmark = pytest.mark.browser

#: The four hosts the two reports used to reach for. Named as literals so a failure
#: message says which one came back rather than only that something did.
FORMER_CDN_HOSTS = (
    "stackpath.bootstrapcdn.com",
    "cdn.datatables.net",
    "code.jquery.com",
    "cdn.jsdelivr.net",
)

#: A host nothing is listening on, in the loopback range, used only to prove the
#: recorder below can see a request at all. Port 9 is `discard`, and the request is
#: made from an offline context, so nothing leaves the machine either way.
UNREACHABLE = "http://127.0.0.1:9/probe.js"


@pytest.fixture
def rendered_cohort_report(tmp_path: Path) -> Path:
    """Render a real cohort report through the real renderer.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``cohort_summary.html``.
    """
    generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=pd.DataFrame(
            [
                {"Sample": "sample_one", "Motif": "5", "Confidence": "High_Precision", "Flag": "Not flagged"},
                {"Sample": "sample_two", "Motif": "6", "Confidence": "Low_Precision", "Flag": "Low_Coverage"},
            ]
        ),
        advntr_df=pd.DataFrame([{"Sample": "sample_one", "VID": "25561", "Flag": "Not flagged"}]),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )
    return tmp_path / "cohort_summary.html"


@pytest.mark.parametrize("mode", ["embedded", "sidecar", "off"])
def test_no_unapproved_request_leaves_the_per_sample_page(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
    mode: str,
) -> None:
    """Every ``--report-igv`` mode, opened with the network available.

    The mode matters: ``embedded`` is the one where a library gets loaded at all, and
    it is the one where a mistake would be a fetch rather than an omission. The other
    two must not reach out either - a report that says the browser is elsewhere has no
    business asking a CDN for one.
    """
    page = open_report(rendered_report_with_alignments(mode), offline=False)
    page.wait_for_timeout(500)

    external = external_requests(page)

    assert external == [], f"--report-igv {mode} fetched {len(external)} resource(s) off-machine: {external}"
    for host in FORMER_CDN_HOSTS:
        assert not any(host in url for url in external), f"{host} is back"


def test_no_unapproved_request_leaves_the_cohort_page(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The cohort report too, Plotly's two inlined bundles included.

    Plotly is what makes this worth measuring rather than reasoning about: its bundle
    is 4.8 MB of JavaScript full of URLs, and whether a donut chart touches any of them
    is a question about Plotly's code, not about ours.
    """
    page = open_report(rendered_cohort_report, offline=False)
    page.wait_for_timeout(500)

    external = external_requests(page)

    assert external == [], f"the cohort report fetched {len(external)} resource(s) off-machine: {external}"


def test_the_shared_plotly_bundle_renders_both_cohort_charts_without_page_errors(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
    page_errors: dict[Page, list[str]],
) -> None:
    """Both chart fragments execute against the one embedded Plotly library."""
    page = open_report(rendered_cohort_report, offline=True)
    graphs = page.locator(".plotly-graph-div")

    expect(graphs).to_have_count(2)
    for index in range(2):
        expect(graphs.nth(index).locator("svg.main-svg").first).to_be_visible()
    assert page_errors[page] == []


def test_the_page_error_recorder_sees_an_exception(
    tmp_path: Path,
    open_report: Callable[..., Page],
    page_errors: dict[Page, list[str]],
) -> None:
    """The chart assertion cannot pass vacuously because page errors were unobserved."""
    control = tmp_path / "page-error.html"
    control.write_text("<script>throw new Error('page-error-control')</script>", encoding="utf-8")

    page = open_report(control, offline=True)

    assert any("page-error-control" in error for error in page_errors[page])


@pytest.mark.parametrize("mode", ["embedded", "sidecar", "off"])
def test_none_of_the_three_removed_libraries_is_back(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
    mode: str,
) -> None:
    """jQuery, DataTables and Bootstrap are gone from the page, not merely from a tag.

    The request check above would also catch a re-added ``<script src>``; this catches
    the other route back in, which is inlining one of them into the document. That is
    not hypothetical - the alignment library is inlined, and "inline the others too"
    is the obvious next idea.
    """
    page = open_report(rendered_report_with_alignments(mode), offline=False)

    state = removed_libraries(page)

    assert state == {"jquery": None, "dataTable": False, "bootstrap": False}, (
        f"a removed library is back in the page: {state}"
    )


def test_the_request_recorder_sees_a_request_that_is_made(
    tmp_path: Path,
    open_report: Callable[..., Page],
) -> None:
    """The control that stops every assertion above from passing vacuously.

    If the recorder saw nothing at all - a listener registered too late, a page that
    never loaded, a Playwright change - "no external request" would be the same green
    for a report full of CDN tags. So: a document that *does* reference an external
    host, opened through the same fixture, must be recorded.

    Offline, so the request is aborted at the route rather than attempted against the
    world; ``request`` fires before routing, which is exactly the property being
    proved. The host is loopback and the port is one nothing serves, so no packet
    leaves the machine even if the abort were to fail.
    """
    control = tmp_path / "control.html"
    control.write_text(f'<!DOCTYPE html><html><body><script src="{UNREACHABLE}"></script></body></html>', "utf-8")

    page = open_report(control, offline=True)

    assert external_requests(page) == [UNREACHABLE], (
        "the request recorder did not see a request the document plainly makes, so every "
        "'no external request' assertion in this file would pass against any report"
    )
