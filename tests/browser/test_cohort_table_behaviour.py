"""The cohort report's three affordances, reimplemented in vanilla JS and measured.

Removing DataTables removed the *implementation* of the cohort's flag filter, its
search and its paging in one go. All three are reviewed behaviour for triage across
samples (#242, precondition P4), so all three were rebuilt rather than quietly dropped,
and this file is what makes "two of three" fail rather than look finished.

Why the cohort may hide a row and the per-sample report may not
--------------------------------------------------------------
It is a difference in what the document is for, not in what the code can do. Cohort
triage across hundreds of samples is a search for the ones worth opening, and hiding
rows is what that task consists of. A per-sample report is a record of one patient's
run: a control that removes evidence from it is the defect this whole issue is named
after. ``tests/unit/test_generate_report.py::test_no_per_sample_result_row_enters_a_hiding_path``
holds the other side of the line, and it does so by naming the constructs *this* file
exercises - so the two cannot both be satisfied by one implementation shared across
both reports.

What withholding a row obliges the report to say
------------------------------------------------
That it is doing it. DataTables' "Showing 1 to 10 of 57 entries" was removed from the
per-sample report because it counted the rows left *after* the filter, which is a
second count contradicting the server's. Here it is the opposite: the reader has
controls that remove rows on purpose, so a count of what is left - and of what is being
withheld - is the only thing standing between a filtered table and a wrong conclusion.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
from playwright.sync_api import Page

from vntyper.cli import load_config
from vntyper.scripts.cohort_summary import generate_cohort_summary_report

pytestmark = pytest.mark.browser

#: The default page size the script uses. A literal here rather than an import: a test
#: that asks the template what its page size is agrees with whatever it says.
PAGE_SIZE = 25

#: Enough samples to page. 30 clean rows and 4 flagged ones, so the first page is full,
#: the second is short, and the flag filter has something to withhold on both.
CLEAN_SAMPLES = 30
FLAGGED_SAMPLES = 4

#: Read the rows a reader can actually see, from one table by index. The vanilla filter
#: sets ``hidden`` and leaves the row in the document, so presence is not visibility -
#: and the index matters because the cohort report has three tables and a bare
#: ``table tbody tr`` would pool the Kestrel table's rows with the adVNTR table's.
_VISIBLE_ROWS = """
(index) => Array.from(document.querySelectorAll("table")[index].tBodies[0].rows)
    .filter(e => e.offsetParent !== null)
    .map(e => e.textContent.replace(/\\s+/g, ' ').trim())
"""

#: The Kestrel table's toolbar and pager, addressed by the table index its controls
#: carry. Built by the script, so a page with no script has none of them - which is the
#: correct fallback, and the last test here.
SEARCH = "#tableSearch_0"
PAGE_SIZE_SELECT = "#tablePageSize_0"
STATUS = "#tableStatus_0"
NEXT = "#tableNext_0"
PREVIOUS = "#tablePrevious_0"


@pytest.fixture
def rendered_cohort_report(tmp_path: Path) -> Path:
    """A cohort big enough to page, with flagged rows to filter.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``cohort_summary.html``.
    """
    rows = [
        {
            "Sample": f"clean_{index:03d}",
            "Motif": str(index),
            "Confidence": "High_Precision",
            "Flag": "Not flagged",
        }
        for index in range(CLEAN_SAMPLES)
    ] + [
        {
            "Sample": f"flagged_{index:03d}",
            "Motif": str(index),
            "Confidence": "Low_Precision",
            "Flag": "Low_Coverage",
        }
        for index in range(FLAGGED_SAMPLES)
    ]
    generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=pd.DataFrame(rows),
        advntr_df=pd.DataFrame([{"Sample": "clean_000", "VID": "25561", "Flag": "Not flagged"}]),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )
    return tmp_path / "cohort_summary.html"


def _visible(page: Page) -> list[str]:
    """The Kestrel table's visible rows.

    Args:
        page: An already-loaded page.

    Returns:
        list[str]: One whitespace-normalised string per visible row.
    """
    rows: list[str] = page.evaluate(_VISIBLE_ROWS, 0)
    return rows


def test_the_report_pages_its_rows(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Paging: a full first page, a shorter second, and Previous back to the first.

    The fixture has 30 clean rows and a default page size of 25, so this measures a
    boundary rather than a round number.
    """
    page = open_report(rendered_cohort_report, offline=True)

    first = _visible(page)
    assert len(first) == PAGE_SIZE, f"the first page shows {len(first)} rows, not {PAGE_SIZE}"

    page.click(NEXT)
    second = _visible(page)
    assert len(second) == CLEAN_SAMPLES - PAGE_SIZE, f"the second page shows {len(second)} rows"
    assert set(second).isdisjoint(first), "the second page repeats rows from the first"

    page.click(PREVIOUS)
    assert _visible(page) == first, "Previous did not return to the first page"


def test_the_page_size_is_the_readers_choice(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Changing the page size changes the page, and returns to the first one."""
    page = open_report(rendered_cohort_report, offline=True)

    page.select_option(PAGE_SIZE_SELECT, "10")

    assert len(_visible(page)) == 10, f"the page size control did not take: {len(_visible(page))} rows"


def test_the_search_box_narrows_the_table(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Search: a substring over the row, including the flag reason it no longer shows.

    The flag cell is rebuilt as a glyph, so its reason is no longer in the row's text.
    Searching it anyway is what keeps "find the samples flagged for coverage" working
    after that rebuild - the reason is read from the ``data-original`` the script keeps.
    """
    page = open_report(rendered_cohort_report, offline=True)
    page.check("#toggleFlagged")

    page.fill(SEARCH, "flagged_00")
    narrowed = _visible(page)

    assert len(narrowed) == FLAGGED_SAMPLES, f"search matched {len(narrowed)} rows: {narrowed}"
    assert all("flagged_00" in row for row in narrowed)

    page.fill(SEARCH, "Low_Coverage")
    by_reason = _visible(page)
    assert len(by_reason) == FLAGGED_SAMPLES, (
        f"searching the flag reason matched {len(by_reason)} rows; the reason is drawn as a glyph, so it has "
        "to be searched from the value the script remembered"
    )

    page.fill(SEARCH, "")
    assert len(_visible(page)) == PAGE_SIZE, "clearing the search did not restore the table"


def test_the_flag_filter_hides_flagged_rows_until_the_reader_asks(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The reviewed cohort behaviour, unchanged in effect and rewritten in implementation.

    Default is unticked, and flagged rows are withheld. This is the one place in either
    report where that is allowed, and it is allowed because triage is a search.
    """
    page = open_report(rendered_cohort_report, offline=True)

    page.fill(SEARCH, "flagged_")
    assert _visible(page) == [], "a flagged row is visible before the reader asked for one"

    page.check("#toggleFlagged")
    assert len(_visible(page)) == FLAGGED_SAMPLES, f"ticking the switch showed {len(_visible(page))} flagged rows"

    page.uncheck("#toggleFlagged")
    assert _visible(page) == [], "unticking the switch did not hide the flagged rows again"


def test_the_switch_says_what_it_will_do_next(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The label is rewritten, and the input it wraps survives the rewrite.

    Writing to the ``<label>`` rather than to its ``<span>`` would delete the control
    the label names, which is why the span exists at all.
    """
    page = open_report(rendered_cohort_report, offline=True)

    assert page.inner_text("#toggleFlaggedLabel").strip() == "Show flagged values"

    page.check("#toggleFlagged")
    assert page.inner_text("#toggleFlaggedLabel").strip() == "Hide flagged values"
    assert page.locator("#toggleFlagged").count() == 1, "rewriting the label deleted the control"


def test_the_reader_is_told_how_many_rows_are_being_withheld(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """A filtered table that does not say it is filtered is a wrong answer.

    The count is the whole justification for allowing the controls: a reader who sees
    25 of 34 rows and is told so can ask for the rest, and a reader who is not told has
    no reason to suspect there are any.
    """
    page = open_report(rendered_cohort_report, offline=True)

    status = page.inner_text(STATUS)
    assert "Showing 1 to 25 of 30" in status, f"the count line reads {status!r}"
    assert "hidden by the controls above" in status, f"the count line does not say rows are withheld: {status!r}"

    page.check("#toggleFlagged")
    assert "of 34" in page.inner_text(STATUS), f"the count did not follow the filter: {page.inner_text(STATUS)!r}"


def test_a_search_that_matches_nothing_says_how_many_rows_are_in_the_file(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """An empty table is the state most likely to be misread as "there is nothing".

    So it says the opposite in words: nothing matched, and this many rows are here.
    """
    page = open_report(rendered_cohort_report, offline=True)

    page.fill(SEARCH, "no-such-sample")

    assert _visible(page) == []
    status = page.inner_text(STATUS)
    assert "No row matches" in status and "34 are in this file" in status, f"the empty state reads {status!r}"


def test_no_control_exists_at_all_when_no_script_ran(
    rendered_cohort_report: Path,
    browser,
) -> None:
    """With scripting off, every row of every table is visible and nothing is offered.

    This is the correct fallback and it is the reason the toolbar is built by the script
    rather than written into the markup: a search box that cannot search and a pager
    that cannot page are controls that lie about what the reader can do, and a printed
    copy would carry them.
    """
    context = browser.new_context(java_script_enabled=False)
    try:
        page = context.new_page()
        page.goto(rendered_cohort_report.as_uri(), wait_until="load")

        assert page.locator(SEARCH).count() == 0, "a search box exists that cannot search"
        assert page.locator(".table-pager").count() == 0, "a pager exists that cannot page"
        assert len(_visible(page)) == CLEAN_SAMPLES + FLAGGED_SAMPLES, (
            "some rows are withheld from a reader who has no way to ask for them"
        )
    finally:
        context.close()
