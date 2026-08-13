"""A control that focuses and does nothing is not operable, and a fallback nobody
sees is not a fallback.

``test_keyboard_operability.py`` proves every control can be *reached*. That is half
of a control. This file is the other half: what happens when it is used, and what
the reader is given when the thing behind it never arrived.

Both behaviours here were fixed in the same commit that made the controls reachable,
and both were fixed for the same reason - a ``$`` that never loaded throws at the
top of whichever ``<script>`` element it is in and takes every statement after it in
that block down with it:

* **"Show detailed coverage"** lived in the jQuery block, so offline it focused,
  ticked, and did nothing. That matters beyond hygiene: the detailed view is the
  **only** route to the Coverage QC verdict, the minimum and maximum, and the
  uncovered-base count. An inert switch takes that evidence away from every reader
  whose browser cannot reach three CDNs.
* **the alignment panel** was built inside ``igv.createBrowser(...).then(...)``, in
  the same block, so offline neither the variant table nor any message about the
  missing viewer was ever written. The reader got blank space where an explanation
  belonged.

Neither fix had a test. Moving the handler back into the jQuery block, or deleting
the ``typeof igv === "undefined"`` guard, left every suite in the repository green -
which is the precise regression each fix exists to prevent. AGENTS.md's "touch a
file, add tests for it" is not satisfied by a fix whose removal is silent.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import enhancement_state
from vntyper.scripts.report_formatting import EMPTY_SESSION_DICTIONARY, EMPTY_TABLE_JSON

pytestmark = pytest.mark.browser

#: A distinctive fragment of the message the report writes when igv.js is not there.
#: A literal, not an import: a test that asks the template what it says agrees with
#: whatever the template says, including nothing.
VIEWER_UNAVAILABLE = "The alignment viewer could not be loaded"

#: What the report says when the run produced no IGV report at all. A different
#: fact, and the two must not be confused - hence asserting each is absent from the
#: other's case.
NO_IGV_REPORT = "No IGV alignment visualization is available for this sample"

#: Evidence that exists only in the detailed coverage view. If the switch does not
#: swap the views, this is what the reader never sees.
DETAILED_ONLY_EVIDENCE = ("Coverage QC", "Minimum Coverage", "Uncovered Bases")

#: A variant table and a session per row, spliced into the rendered report. This
#: tier never runs ``create_report`` (igv-reports needs a BED and a BAM), so a
#: report carrying an alignment session has to be constructed. The specimen is local
#: to this file on purpose: ``test_keyboard_operability.py`` carries its own, and a
#: shared one would move both files' subjects at once when it changed.
SPLICED_HEADERS: list[str] = ["unique_id", "CHROM", "POS", "REF", "ALT"]
SPLICED_VARIANTS: list[list[str]] = [["0", "chr1", "155188205", "G", "GG"]]
SPLICED_TABLE_JSON: dict[str, list] = {"headers": SPLICED_HEADERS, "rows": SPLICED_VARIANTS}
SPLICED_SESSION_DICTIONARY: dict[str, str] = {"0": "session0.json"}

#: Whether the alignment library is actually present in the page. This is a better
#: proof that a pass was online than asking jQuery its version, because it is the
#: exact global the guard under test branches on.
_IGV_PRESENT = "() => typeof igv !== 'undefined'"


@pytest.fixture
def report_with_alignment_session(rendered_report: Path) -> Path:
    """Return a copy of the rendered report that claims an alignment session.

    ``initIGV`` is only reached when ``sessionDictionary`` is non-empty, so without
    this the ``typeof igv === "undefined"`` branch is unreachable and a test of it
    would pass against a report that never ran it.

    Both replacements are asserted. A splice that missed would leave the original
    document behind, the "no IGV report was generated" branch would run instead, and
    the tests below would quietly assert something else entirely.

    Args:
        rendered_report: The report rendered by the shared fixture.

    Returns:
        Path: A sibling file carrying a one-row variant table and its session.
    """
    source = rendered_report.read_text(encoding="utf-8")
    empty_table = f"const tableJson = {EMPTY_TABLE_JSON};"
    empty_sessions = f"const sessionDictionary = {EMPTY_SESSION_DICTIONARY};"
    assert empty_table in source, "the rendered report does not carry the empty tableJson literal to replace"
    assert empty_sessions in source, "the rendered report does not carry the empty sessionDictionary literal"

    spliced = source.replace(empty_table, f"const tableJson = {json.dumps(SPLICED_TABLE_JSON)};").replace(
        empty_sessions, f"const sessionDictionary = {json.dumps(SPLICED_SESSION_DICTIONARY)};"
    )

    target = rendered_report.with_name("summary_report_with_session.html")
    target.write_text(spliced, encoding="utf-8")
    return target


# ---------------------------------------------------------------------------
# The detailed coverage switch actually swaps the views
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("offline", [False, True], ids=["online", "offline"])
def test_the_detailed_coverage_switch_swaps_the_two_views(
    rendered_report: Path,
    open_report: Callable[..., Page],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
    offline: bool,
) -> None:
    """Ticking the switch hides the basic view and shows the detailed one, and back.

    Parametrised over the network because the network is what broke it. The
    ``offline`` case is the one that fails if the handler goes back into the jQuery
    block; the ``online`` case is what says the move did not break the working
    path, and it proves it was online rather than assuming it.
    """
    page = open_report(rendered_report, offline=offline)
    if not offline:
        assert not failed_requests[page], f"the online pass could not complete every request: {failed_requests[page]}"
        assert not error_responses[page], f"the online pass got error statuses back: {error_responses[page]}"
        state = enhancement_state(page)
        assert state["jquery"] and state["dataTable"], (
            f"the online pass loaded no working jQuery/DataTables, so it is a second offline pass: {state}"
        )

    basic = page.locator("#basicCoverageView")
    detailed = page.locator("#detailedCoverageView")

    assert basic.is_visible(), "the basic coverage view is not shown by default"
    assert not detailed.is_visible(), "the detailed coverage view is shown before the switch is ticked"

    page.check("#toggleDetailedCoverage")

    assert detailed.is_visible(), (
        f"ticking 'Show detailed coverage' did not show the detailed view (offline={offline}). "
        "The switch focuses and does nothing, so the Coverage QC verdict, the minimum and maximum "
        "and the uncovered-base count are unreachable."
    )
    assert not basic.is_visible(), "both coverage views are shown at once, so the mean is stated twice"

    page.uncheck("#toggleDetailedCoverage")

    assert basic.is_visible(), "unticking the switch did not restore the basic view"
    assert not detailed.is_visible(), "unticking the switch left the detailed view showing"


@pytest.mark.parametrize("offline", [False, True], ids=["online", "offline"])
def test_the_detailed_view_is_the_only_route_to_the_coverage_verdict(
    rendered_report: Path,
    open_report: Callable[..., Page],
    offline: bool,
) -> None:
    """What the switch is *for*, asserted as evidence rather than as a style change.

    A test that only watched ``display`` swap would still pass if the detailed view
    were emptied of the figures that justify its existence.
    """
    page = open_report(rendered_report, offline=offline)
    page.check("#toggleDetailedCoverage")

    shown = page.locator("#detailedCoverageView").inner_text()
    missing = [label for label in DETAILED_ONLY_EVIDENCE if label not in shown]
    assert missing == [], f"the detailed coverage view does not carry {missing} (offline={offline}): {shown!r}"

    basic = page.locator("#basicCoverageView").inner_text()
    assert "Coverage QC" not in basic, (
        "the basic view now states the QC verdict too, so this test no longer describes why the switch exists"
    )


# ---------------------------------------------------------------------------
# The alignment panel explains itself when the viewer is not there
# ---------------------------------------------------------------------------


def test_the_panel_says_the_viewer_could_not_load_when_it_did_not(
    report_with_alignment_session: Path,
    open_report: Callable[..., Page],
) -> None:
    """igv.js is absent offline, and the reader is told so instead of shown blank space.

    Without the ``typeof igv === "undefined"`` guard this throws a ``ReferenceError``
    and ``#igvContainer`` keeps an empty ``<div id="igvDiv">`` - a panel-shaped hole
    with no explanation in it.
    """
    page = open_report(report_with_alignment_session, offline=True)

    assert page.evaluate("() => typeof igv === 'undefined'"), (
        "igv.js is present, so this pass never reaches the branch under test"
    )

    shown = page.locator("#igvContainer").inner_text()
    assert VIEWER_UNAVAILABLE in shown, f"the reader is given no explanation for the missing viewer: {shown!r}"
    assert NO_IGV_REPORT not in shown, (
        "the report claims no IGV data was produced, which is a different fact from a viewer that failed to load"
    )


def test_the_panel_carries_no_failure_message_when_the_viewer_loaded(
    report_with_alignment_session: Path,
    open_report: Callable[..., Page],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
) -> None:
    """The other half: the message is conditional, not decoration printed always.

    A guard that wrote its fallback unconditionally would satisfy the check above
    and tell every reader the viewer had failed. The precondition is asserted
    positively - igv.js really is in the page - so a run on a machine with no
    network fails here rather than passing for the wrong reason.
    """
    page = open_report(report_with_alignment_session, offline=False)

    assert not failed_requests[page], f"the online pass could not complete every request: {failed_requests[page]}"
    assert not error_responses[page], f"the online pass got error statuses back: {error_responses[page]}"
    assert page.evaluate(_IGV_PRESENT), "igv.js did not load, so this is a second offline pass and proves nothing"

    shown = page.locator("#igvContainer").inner_text()
    assert VIEWER_UNAVAILABLE not in shown, f"the report says the viewer failed while the viewer is loaded: {shown!r}"


def test_a_run_with_no_alignment_data_says_that_instead(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Two different absences, two different sentences.

    The shipped fixture produces no IGV report at all, which is not the same as a
    viewer that failed to load, and a reader who is told the wrong one draws the
    wrong conclusion about their run.
    """
    page = open_report(rendered_report, offline=True)

    shown = page.locator("#igvContainer").inner_text()
    assert NO_IGV_REPORT in shown, f"a run with no alignment data says nothing about it: {shown!r}"
    assert VIEWER_UNAVAILABLE not in shown, "a run that produced no IGV report blames the viewer for it"
