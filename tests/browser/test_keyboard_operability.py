"""Every control in the per-sample report must be reachable with the Tab key.

Why this is a browser test and cannot be anything else
-----------------------------------------------------
Focusability is *computed*, not declared. Nothing in the template says an element
is out of the focus order; ``input[type='checkbox'] { display: none; }`` said it,
from a stylesheet, about four inputs at once - and Bootstrap kept drawing the
switch from ``::before``/``::after`` on the sibling label, so the page still
*looked* right. Measured in Chromium against the unmodified template, online and
offline alike:

    #highlightFlagged        width 0, height 0, display none, focusable false
    #toggleDetailedCoverage  width 0, height 0, display none, focusable false
    #collapsible             width 0, height 0, display none, focusable false
    #logToggle               width 0, height 0, display none, focusable false

A source assertion cannot see any of that, and the unit tier has no browser
(``tests/unit/test_template_escaping.py`` says so). So this presses Tab in a real
engine and records what the engine focused.

What the assertion has to be
----------------------------
The **complete** expected set, not "some control was reachable". A traversal that
only checks elements which already happened to be focusable proves nothing: run
against the unmodified template it would find whatever DataTables added to the
column headers and pass. Every control the report offers is named below, and the
set of Tab stops must contain all of them.

Offline, deliberately
---------------------
The report has to stand on its own, and the offline pass is the one where it must.
There is nothing left to load (#242), so every Tab stop measured here is one the
document itself put there - which used to be true only offline, because online
DataTables added header stops of its own and the measurement became about DataTables.
The online/offline *equivalence* of the report is the subject of
``test_report_determinism.py``; this is about the document's own keyboard path.

The two IGV cases
-----------------
The shipped fixture renders the no-alignment case: no BED, so no ``create_report``,
empty ``tableJson`` and ``sessionDictionary``, and an authored sentence where the
alignment view would be. That case must not take the rest of the page's keyboard path
down with it.

The IGV-present case renders a real report with a real session through
``rendered_report_with_alignments``, in ``--report-igv sidecar`` mode. Sidecar rather
than embedded on purpose: it produces exactly the same variant table and row controls,
and it leaves no half-megabyte library to expand asynchronously and add controls of its
own to the tab order halfway through a traversal. What igv.js does once it *is*
expanded is ``test_alignment_states.py``'s subject, not this file's.

This replaces splicing two JavaScript literals into an already-rendered file. The
splice could not produce the markup the template authors from ``igv_session_available``,
so it tested the script against a page the generator never writes.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

pytestmark = pytest.mark.browser

#: The controls the per-sample report offers when no IGV table was produced. Every
#: one of these was unreachable by keyboard before this change.
#:
#: The last two are the sortable column headings, which are new interactive controls
#: (#242): they replace DataTables' own header controls, they are built by the report's
#: own script, and a control that a script creates is exactly the kind that gets built
#: as a ``<th onclick>`` and reaches nobody's keyboard. First and last of the Kestrel
#: heading row are named here so the traversal below covers the whole row's span;
#: :func:`test_every_sortable_heading_is_an_operable_control` then holds *all* of them
#: to being focusable, named, keyboard-operable and reflected in ``aria-sort``.
BASE_CONTROLS: tuple[str, ...] = (
    "#highlightFlagged",
    "#toggleDetailedCoverage",
    "#kestrel_table thead th:first-child button",
    "#kestrel_table thead th:last-child button",
    "#variantsToggle",
    "#logToggle",
)

#: Every sortable column heading of the Kestrel table.
SORT_HEADERS = "#kestrel_table thead th button"

#: The per-row controls of the IGV variant table, one per row of the shared fixture.
VARIANT_CONTROLS: tuple[str, ...] = ("#variantSelect_0", "#variantSelect_1")

#: Tag every expected control with the selector it was named by, so a Tab stop can be
#: identified as "the thing the test asked for" rather than by whatever text it happens
#: to contain. The sortable headings have no ids and their labels are column names that
#: move, so a descriptor built from the element's own text would have to be rewritten
#: every time a column is renamed - and would silently stop matching instead of failing
#: on the thing it was written about.
_TAG_EXPECTED = """
(selectors) => selectors.map(selector => {
    const el = document.querySelector(selector);
    if (el) {
        el.dataset.tabProbe = selector;
    }
    return !!el;
})
"""

#: Describe whatever the browser has focused, or ``null`` once focus has left the
#: document. The tag above wins where it is present, then the element's id; anything
#: else is described well enough to read in a failure message.
_ACTIVE_DESCRIPTOR = """
() => {
    const el = document.activeElement;
    if (!el || el === document.body || el === document.documentElement) {
        return null;
    }
    if (el.dataset && el.dataset.tabProbe) {
        return el.dataset.tabProbe;
    }
    if (el.id) {
        return "#" + el.id;
    }
    const text = (el.textContent || "").replace(/\\s+/g, " ").trim().slice(0, 40);
    return el.tagName.toLowerCase() + "[" + text + "]";
}
"""

#: Whether one control is present, has a box a pointer could hit, and takes focus
#: when asked. A control that is focusable but has no box is the same defect in a
#: different costume, so both are measured.
_PROBE = """
(selector) => {
    const el = document.querySelector(selector);
    if (!el) {
        return {present: false};
    }
    const rect = el.getBoundingClientRect();
    const style = getComputedStyle(el);
    el.focus();
    return {
        present: true,
        width: rect.width,
        height: rect.height,
        display: style.display,
        focusable: document.activeElement === el,
    };
}
"""


@pytest.fixture
def report_with_variant_table(rendered_report_with_alignments: Callable[..., Path]) -> Path:
    """A real report carrying a real two-row variant table and a session per row.

    ``--report-igv sidecar``: the same variant table and the same row controls as
    ``embedded``, without a half-megabyte library expanding asynchronously and adding
    igv.js's own controls to the tab order partway through a traversal.

    Args:
        rendered_report_with_alignments: The shared renderer.

    Returns:
        Path: The rendered report.
    """
    return rendered_report_with_alignments("sidecar")


def _tab_stops(page: Page, limit: int = 40) -> list[str]:
    """Press Tab until focus leaves the document, recording each stop.

    Args:
        page: An already-loaded page.
        limit: How many presses to make before giving up, so a document that traps
            focus fails with a readable list rather than hanging.

    Returns:
        list[str]: One descriptor per Tab stop, in order, stopping at the first
        repeat (headless Chromium wraps around rather than leaving for browser
        chrome).
    """
    stops: list[str] = []
    for _ in range(limit):
        page.keyboard.press("Tab")
        stop: str | None = page.evaluate(_ACTIVE_DESCRIPTOR)
        if stop is None or stop in stops:
            break
        stops.append(stop)
    return stops


def _assert_reachable(page: Page, expected: tuple[str, ...]) -> None:
    """Fail unless every expected control is a Tab stop, and say which were not.

    Args:
        page: An already-loaded page.
        expected: The complete set of control selectors that must be reachable.
    """
    found = page.evaluate(_TAG_EXPECTED, list(expected))
    absent = [selector for selector, present in zip(expected, found, strict=True) if not present]
    assert absent == [], f"these controls are not in the document at all, so no traversal can reach them: {absent}"

    stops = _tab_stops(page)
    missing = [selector for selector in expected if selector not in stops]
    assert not missing, f"these controls cannot be reached with the Tab key: {missing}; the Tab order was {stops}"


def test_every_control_is_reachable_by_tab_when_no_igv_table_was_produced(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The IGV-failure case: no variant table, and the rest of the page still works.

    ``bed_file`` is unset for this tier, so ``create_report`` never runs and the
    report falls back to its empty literals. That is a real shape - a run with no
    alignment input produces exactly this file - and the reader of one still needs
    the two switches and both disclosures.
    """
    page = open_report(rendered_report, offline=True)

    assert page.locator("#variant_table").count() == 0, "this fixture is supposed to be the IGV-failure case"
    _assert_reachable(page, BASE_CONTROLS)


def test_every_control_is_reachable_by_tab_when_an_igv_table_is_present(
    report_with_variant_table: Path,
    open_report: Callable[..., Page],
) -> None:
    """The IGV-present case: one keyboard-operable control per variant row.

    The rows were ``<tr onclick>`` with no keyboard path at all, so this is the
    half of "every control" that no amount of fixing the stylesheet would have
    reached.
    """
    page = open_report(report_with_variant_table, offline=True)

    assert page.locator("#variant_table tbody tr").count() == len(VARIANT_CONTROLS), (
        "the spliced variant table did not render, so this is not the IGV-present case"
    )
    _assert_reachable(page, BASE_CONTROLS + VARIANT_CONTROLS)


def test_every_control_has_a_box_and_takes_focus(
    report_with_variant_table: Path,
    open_report: Callable[..., Page],
) -> None:
    """The measurement the defect was found with, run as an assertion.

    ``display: none`` took the four inputs out of the focus order *and* out of the
    accessibility tree while the label kept drawing the switch, which is why visual
    QA never saw it. Pinning the box as well as the focus is what stops the fix
    being "give it a tabindex and leave it invisible".
    """
    page = open_report(report_with_variant_table, offline=True)

    measured = {selector: page.evaluate(_PROBE, selector) for selector in BASE_CONTROLS + VARIANT_CONTROLS}

    for selector, state in measured.items():
        assert state["present"], f"{selector} is not in the document at all"
        assert state["display"] != "none", f"{selector} is display:none, so it is out of the accessibility tree"
        assert state["width"] > 0 and state["height"] > 0, f"{selector} has no box a reader could see or hit: {state}"
        assert state["focusable"], f"{selector} does not take focus: {state}"


def test_each_control_is_announced_by_the_name_its_label_gives_it(
    report_with_variant_table: Path,
    open_report: Callable[..., Page],
) -> None:
    """A focus stop with no accessible name is a control nobody can identify.

    Every lookup here goes through Playwright's role/name resolution, which is the
    accessibility tree rather than the DOM - so it fails if the name is carried
    only by an adjacent ``<span>`` the browser never associates with the control.
    """
    page = open_report(report_with_variant_table, offline=True)

    assert page.get_by_role("checkbox", name="Highlight flagged values").count() == 1
    assert page.get_by_role("checkbox", name="Show detailed coverage").count() == 1
    assert page.get_by_role("group", name="Variants").count() == 1
    assert page.get_by_role("group", name="Pipeline log").count() == 1
    # The row control names itself with the variant it shows, so a screen-reader
    # user moving between rows hears which one they are on rather than "Show".
    assert page.get_by_role("button", name="Showing in the alignment view: CHROM chr1, POS 155188205").count() == 1
    assert page.get_by_role("button", name="Show in the alignment view: CHROM chr1, POS 155188305").count() == 1


def test_enter_and_space_move_the_variant_selection(
    report_with_variant_table: Path,
    open_report: Callable[..., Page],
) -> None:
    """Keyboard activation, and the current row saying so in the accessibility tree.

    The first row is selected on render. Focusing the second row's control and
    pressing Enter must move ``aria-current``; Space must do the same on the way back.
    Neither key did anything at all against ``<tr onclick>``.

    ``aria-current`` rather than ``aria-selected`` inside a ``role="grid"`` (#242): the
    grid role promises Left/Right/Home/End movement this table does not implement, and
    a role that over-promises is worse than no role. ``aria-current`` states the same
    fact with no obligation attached, and it is *removed* rather than set to "false" on
    the rows that are not current, which is how the attribute is defined.
    """
    page = open_report(report_with_variant_table, offline=True)

    assert page.locator("#row_0").get_attribute("aria-current") == "true"

    page.focus("#variantSelect_1")
    page.keyboard.press("Enter")
    assert page.locator("#row_1").get_attribute("aria-current") == "true"
    assert page.locator("#row_0").get_attribute("aria-current") is None

    page.focus("#variantSelect_0")
    page.keyboard.press(" ")
    assert page.locator("#row_0").get_attribute("aria-current") == "true"
    assert page.locator("#row_1").get_attribute("aria-current") is None


def test_the_selected_variant_row_says_so_in_words(
    report_with_variant_table: Path,
    open_report: Callable[..., Page],
) -> None:
    """Which row is showing must not be carried by the highlight colour alone.

    The selected row is tinted and outlined, and that was the whole of the signal.
    The control's own text now reads ``Showing`` on the selected row and ``Show``
    on the others, so the selection survives greyscale, print and a screen reader.
    """
    page = open_report(report_with_variant_table, offline=True)

    assert page.locator("#variantSelect_0").inner_text().strip() == "Showing"
    assert page.locator("#variantSelect_1").inner_text().strip() == "Show"

    page.focus("#variantSelect_1")
    page.keyboard.press("Enter")

    assert page.locator("#variantSelect_0").inner_text().strip() == "Show"
    assert page.locator("#variantSelect_1").inner_text().strip() == "Showing"


# ---------------------------------------------------------------------------
# The sortable column headings, which are new controls (#242)
# ---------------------------------------------------------------------------


def test_every_sortable_heading_is_an_operable_control(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Every heading of a sortable table is focusable, named, and has a box.

    These replace DataTables' own header controls. The report's script builds them,
    which is precisely the shape that arrives as a ``<th onclick>`` with no keyboard
    path, no accessible name and no state - the same defect the variant rows had.

    The accessible name is the column's own heading text, so it needs no label and
    cannot drift from what the reader sees.
    """
    page = open_report(rendered_report, offline=True)

    headings = page.locator(SORT_HEADERS)
    count = headings.count()
    assert count >= 8, f"the Kestrel table has {count} sortable headings; the fixture renders eleven columns"

    for index in range(count):
        heading = headings.nth(index)
        label = heading.inner_text().strip()
        assert label, f"sortable heading {index} has no text, so it has no accessible name"
        box = heading.bounding_box()
        assert box and box["width"] > 0 and box["height"] > 0, f"heading {label!r} has no box: {box}"
        heading.focus()
        assert page.evaluate("() => document.activeElement.tagName.toLowerCase()") == "button", (
            f"heading {label!r} does not take focus"
        )


def test_every_sortable_heading_is_reachable_by_tab(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """All eleven of them, not just the two named in ``BASE_CONTROLS``.

    A traversal that stopped at the first and last would pass against a middle column
    whose heading had been left as plain text.
    """
    page = open_report(rendered_report, offline=True)

    count = page.locator(SORT_HEADERS).count()
    selectors = tuple(f"#kestrel_table thead th:nth-child({index + 1}) button" for index in range(count))

    _assert_reachable(page, selectors)


def test_enter_and_space_sort_the_table_and_say_which_way(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Keyboard activation, and ``aria-sort`` reflecting the result.

    A ``<button>`` answers both keys by itself, which is why it is a button - but
    "the control responds" and "the accessibility tree knows what it did" are two
    claims, and only the second survives a change to how the arrow is drawn.

    **The rows are counted before and after.** Sorting reorders; it must never remove.
    That is the invariant that makes a sort control acceptable in a single-patient
    report at all, and it is the one a future "sort and filter" would break.
    """
    page = open_report(rendered_report, offline=True)

    heading = page.locator(SORT_HEADERS).first
    column = page.locator("#kestrel_table thead th").first
    before = page.locator("#kestrel_table tbody tr").count()

    assert column.get_attribute("aria-sort") == "none"

    heading.focus()
    page.keyboard.press("Enter")
    assert column.get_attribute("aria-sort") == "ascending"

    page.keyboard.press(" ")
    assert column.get_attribute("aria-sort") == "descending"

    assert page.locator("#kestrel_table tbody tr").count() == before, "sorting removed a row"
    assert page.locator("#kestrel_table tbody tr:visible").count() == before, "sorting hid a row"


def test_sorting_reorders_the_rows_it_was_asked_to_reorder(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The control does the thing, and the thing is a permutation.

    Asserting ``aria-sort`` alone would pass against a header that announces a
    direction and moves nothing. The set of rows is compared before and after, so a
    sort that dropped or duplicated one fails here rather than looking like a
    reordering.
    """
    page = open_report(rendered_report, offline=True)

    read = "els => els.map(e => e.textContent.replace(/\\s+/g, ' ').trim())"
    heading = page.locator("#kestrel_table thead th button").nth(2)
    before: list[str] = page.eval_on_selector_all("#kestrel_table tbody tr", read)

    heading.click()
    ascending: list[str] = page.eval_on_selector_all("#kestrel_table tbody tr", read)
    heading.click()
    descending: list[str] = page.eval_on_selector_all("#kestrel_table tbody tr", read)

    # Ascending against descending, rather than either against the render order: the
    # fixture's rows may already be in a column's ascending order, and a test that
    # demanded movement from a *sorted* starting point would be asserting that the sort
    # is wrong. Three distinct values cannot be in both orders at once.
    assert sorted(ascending) == sorted(before), "sorting changed which rows are in the table, not only their order"
    assert sorted(descending) == sorted(before), "sorting changed which rows are in the table, not only their order"
    assert ascending != descending, "the sort control announced two directions and produced one order"
    assert ascending == list(reversed(descending)), f"{ascending} is not the reverse of {descending}"
