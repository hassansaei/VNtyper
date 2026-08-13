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
The report has to stand on its own, and the offline pass is the one where it must:
no jQuery, no DataTables, no Bootstrap, so every Tab stop measured here is one the
template itself put in the document. Online the same stops are present plus
DataTables' own header stops, which would make the measurement about DataTables.
The online/offline *equivalence* of the report is the subject of
``test_report_determinism.py``; this is about the document's own keyboard path.

The two IGV cases
-----------------
``create_report`` is not run by this tier, so the shipped fixture renders the
IGV-failure case: empty ``tableJson``, empty ``sessionDictionary``, an authored
message where the alignment view would be. That case must not take the rest of the
page's keyboard path down with it. The IGV-present case is produced by splicing a
real variant table into the same rendered file, which is the only way to get one
without igv-reports and a BAM; the fixture fails loudly if the splice does not
land, because an IGV-present case that silently degraded to the failure case would
assert nothing about the row controls.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from vntyper.scripts.report_formatting import EMPTY_SESSION_DICTIONARY, EMPTY_TABLE_JSON

pytestmark = pytest.mark.browser

#: The controls the per-sample report offers when no IGV table was produced. Every
#: one of these was unreachable by keyboard before this change.
BASE_CONTROLS: tuple[str, ...] = (
    "#highlightFlagged",
    "#toggleDetailedCoverage",
    "#variantsToggle",
    "#logToggle",
)

#: The per-row controls of the spliced-in IGV variant table, one per row.
VARIANT_CONTROLS: tuple[str, ...] = ("#variantSelect_0", "#variantSelect_1")

#: The variant table spliced into the rendered report for the IGV-present case.
#: Shaped exactly like igv-reports' own ``tableJson``: column 0 is the unique id
#: the session dictionary is keyed by and is not displayed.
SPLICED_HEADERS: list[str] = ["unique_id", "CHROM", "POS", "REF", "ALT"]
SPLICED_VARIANTS: list[list[str]] = [
    ["0", "chr1", "155188205", "G", "GG"],
    ["1", "chr1", "155188305", "GG", "G"],
]
SPLICED_TABLE_JSON: dict[str, list] = {"headers": SPLICED_HEADERS, "rows": SPLICED_VARIANTS}

#: A session per row, so ``initIGV`` takes its "a session exists" branch and the
#: page exercises the same code path a real IGV report does. igv.js itself is
#: absent offline; the report is expected to say so and keep the table.
SPLICED_SESSION_DICTIONARY = {"0": "session0.json", "1": "session1.json"}

#: Describe whatever the browser has focused, or ``null`` once focus has left the
#: document. The key is the element's id where it has one, because that is what the
#: expected set is written in; anything else is described well enough to read in a
#: failure message.
_ACTIVE_DESCRIPTOR = """
() => {
    const el = document.activeElement;
    if (!el || el === document.body || el === document.documentElement) {
        return null;
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
def report_with_variant_table(rendered_report: Path) -> Path:
    """Return a copy of the rendered report carrying a real IGV variant table.

    ``create_report`` needs igv-reports and a BAM, neither of which this tier has,
    so the two JavaScript literals the template interpolates are replaced in the
    rendered file. Both replacements are asserted: a splice that missed would leave
    the IGV-failure document behind and the IGV-present test would quietly assert
    the failure case twice.

    Args:
        rendered_report: The report rendered by the shared fixture.

    Returns:
        Path: A sibling file with a two-row variant table and a session per row.
    """
    source = rendered_report.read_text(encoding="utf-8")
    empty_table = f"const tableJson = {EMPTY_TABLE_JSON};"
    empty_sessions = f"const sessionDictionary = {EMPTY_SESSION_DICTIONARY};"
    assert empty_table in source, "the rendered report does not carry the empty tableJson literal to replace"
    assert empty_sessions in source, "the rendered report does not carry the empty sessionDictionary literal"

    spliced = source.replace(empty_table, f"const tableJson = {json.dumps(SPLICED_TABLE_JSON)};").replace(
        empty_sessions, f"const sessionDictionary = {json.dumps(SPLICED_SESSION_DICTIONARY)};"
    )

    target = rendered_report.with_name("summary_report_with_igv.html")
    target.write_text(spliced, encoding="utf-8")
    return target


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

    assert page.locator("#variant_table tbody tr").count() == len(SPLICED_VARIANTS), (
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
    """Keyboard activation, and the selected row saying so in the accessibility tree.

    The first row is selected on render. Focusing the second row's control and
    pressing Enter must move ``aria-selected``; Space must do the same on the way
    back. Neither key did anything at all against ``<tr onclick>``.
    """
    page = open_report(report_with_variant_table, offline=True)

    assert page.locator("#row_0").get_attribute("aria-selected") == "true"

    page.focus("#variantSelect_1")
    page.keyboard.press("Enter")
    assert page.locator("#row_1").get_attribute("aria-selected") == "true"
    assert page.locator("#row_0").get_attribute("aria-selected") == "false"

    page.focus("#variantSelect_0")
    page.keyboard.press(" ")
    assert page.locator("#row_0").get_attribute("aria-selected") == "true"
    assert page.locator("#row_1").get_attribute("aria-selected") == "false"


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
