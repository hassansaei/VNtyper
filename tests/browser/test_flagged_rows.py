"""A flagged variant row must reach every reader of the file, and say why it is flagged.

Issue #242 is named after this defect: a ``High_Precision_flagged`` report narrated
a flagged pathogenic variant in its screening summary and then rendered a table the
reader could not find it in. The row was in the file -- ``to_html`` wrote it -- but
a client-side DataTables search predicate removed it from the DOM before the reader
saw anything, unless they ticked a "Show flagged values" switch they had no reason
to know existed.

That made the defect invisible to the unit tier by construction: asserting the flag
reason is in the rendered HTML passed against the unmodified code, because the
hiding happened after the HTML was parsed. It was only observable as a **visible
row count after initialisation**, which is why these checks live here. The
predicate is gone and they are green; they stay so the next change cannot bring it
back under another name.

The specimen comes from ``conftest.py``: three Kestrel rows, one clean and two
carrying ``Flag`` values that predicate hid. No shipped configuration produces a
flagged row -- ``report_config.json`` declares no ``flag_rules`` key at all -- so
the fixture states them directly rather than hoping a rule fires.

Measured before the fix: **1 of 3 rows online, 3 of 3 offline**.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import (
    CLEAN_FLAG,
    EXPECTED_KESTREL_ROWS,
    HIDDEN_FLAGS,
    KESTREL_ROW_SELECTOR,
    enhancement_state,
)

pytestmark = pytest.mark.browser

#: Read the Flag cell -- the last cell of each Kestrel row -- as the reader sees it.
_FLAG_CELL_TEXT = """
els => els.filter(e => e.offsetParent !== null)
          .map(e => {
              const cells = e.querySelectorAll('td');
              return cells.length ? cells[cells.length - 1].textContent.replace(/\\s+/g, ' ').trim() : '';
          })
"""

#: Read the Depth Score cell (index 8 of the eleven display columns).
_DEPTH_SCORE_TEXT = """
els => els.filter(e => e.offsetParent !== null)
          .map(e => {
              const cells = e.querySelectorAll('td');
              return cells.length > 8 ? cells[8].textContent.trim() : '';
          })
"""


def _require_online(
    page: Page,
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
    *,
    offline: bool,
) -> None:
    """Fail an online pass that did not actually get its third-party scripts.

    Proving a pass was online takes a POSITIVE check. The absence of a transport
    failure is not one: ``requestfailed`` fires for a dropped connection and a
    blocked route, and not for a 404, a 503, or a 200 carrying JavaScript that
    never defines ``jQuery``. In any of those the page is left unenhanced, the
    online pass silently becomes a second offline pass, and the online
    parametrisations below pass while proving nothing about the online path.

    So: nothing failed at the transport layer, nothing came back with an error
    status, and the enhancement is live in the page itself. The same three-part
    guard ``test_report_determinism`` uses, reading the same registries - see the
    SUCCESSOR note there for what replaces it once the CDN tags are gone and there
    is no jQuery left to ask.

    Args:
        page: The page the pass opened.
        failed_requests: Transport-layer failure registry, keyed by page.
        error_responses: Error-status registry, keyed by page.
        offline: Whether the pass was deliberately offline. Nothing is asserted
            then - an offline pass is *expected* to fail every non-``file://``
            request, and that is what the fixture's route handler makes it do.
    """
    if offline:
        return

    assert not failed_requests[page], (
        f"the online pass could not complete every request, so it was not online: {failed_requests[page]}"
    )
    assert not error_responses[page], (
        f"the online pass got error statuses back, so its scripts did not run: {error_responses[page]}"
    )
    state = enhancement_state(page)
    assert state["jquery"] and state["dataTable"], (
        "the online pass loaded no working jQuery/DataTables, so it rendered the same unenhanced page the "
        f"offline pass does and this check proves nothing about the online path: {state}"
    )


@pytest.mark.parametrize("offline", [False, True], ids=["online", "offline"])
def test_every_flagged_row_is_visible_after_initialisation(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
    offline: bool,
) -> None:
    """Every Kestrel row the pipeline wrote is visible, whatever the reader's network.

    Failed before the fix in the ``online`` case only, where the removed
    ``$.fn.dataTable.ext.search`` predicate left 1 of 3 rows standing.
    """
    page = open_report(rendered_report, offline=offline)
    _require_online(page, failed_requests, error_responses, offline=offline)

    rows = visible_kestrel_rows(page)

    assert len(rows) == EXPECTED_KESTREL_ROWS, (
        f"the reader sees {len(rows)} of {EXPECTED_KESTREL_ROWS} Kestrel rows (offline={offline}): {rows}"
    )
    for flag in HIDDEN_FLAGS:
        assert any(flag in row for row in rows), f"the row flagged {flag!r} is not visible (offline={offline}): {rows}"
    assert any(CLEAN_FLAG in row for row in rows), f"the unflagged row is not visible (offline={offline}): {rows}"


@pytest.mark.parametrize("offline", [False, True], ids=["online", "offline"])
def test_the_flag_reason_is_readable_text_rather_than_a_tooltip(
    rendered_report: Path,
    open_report: Callable[..., Page],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
    offline: bool,
) -> None:
    """The reason a row is flagged is text in the cell, not a hover-only attribute.

    ``updateFlagColumn`` used to replace the cell with a bare glyph carrying the
    reason in a Bootstrap ``title``, which is invisible in print, invisible to a
    screen reader once Bootstrap moves ``title`` into ``data-original-title``, and
    absent entirely when the script does not run.
    """
    page = open_report(rendered_report, offline=offline)
    _require_online(page, failed_requests, error_responses, offline=offline)

    flags: list[str] = page.eval_on_selector_all(KESTREL_ROW_SELECTOR, _FLAG_CELL_TEXT)

    assert len(flags) == EXPECTED_KESTREL_ROWS, f"expected every row's Flag cell, got {flags} (offline={offline})"
    for reason in (CLEAN_FLAG, *HIDDEN_FLAGS):
        assert any(reason in cell for cell in flags), (
            f"no Flag cell reads {reason!r} as text (offline={offline}): {flags}"
        )


@pytest.mark.parametrize("offline", [False, True], ids=["online", "offline"])
def test_the_flag_switch_emphasises_without_removing_anything(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
    offline: bool,
) -> None:
    """Ticking "Highlight flagged values" changes styling and nothing else.

    The switch is operated through its ``<label>``, which is how a pointer reaches
    it. It used to be the *only* way: the template hid every bare checkbox so that
    its CSS-only accordions could work, which left the input itself unreachable by
    keyboard and absent from the accessibility tree. That is fixed, and the
    keyboard path is measured in ``test_keyboard_operability.py``; this check is
    about what the switch does, not how it is reached, so it keeps clicking the
    label.
    """
    page = open_report(rendered_report, offline=offline)
    _require_online(page, failed_requests, error_responses, offline=offline)
    before = visible_kestrel_rows(page)

    page.click('label[for="highlightFlagged"]')

    assert visible_kestrel_rows(page) == before, "the switch changed which rows are visible"
    assert (
        page.eval_on_selector(KESTREL_ROW_SELECTOR, "el => el.closest('table').className").find("highlight-flagged")
        >= 0
    ), "the switch did not put the highlight class on the results table"


def test_the_depth_score_reads_the_same_online_and_offline(
    rendered_report: Path,
    open_report: Callable[..., Page],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
) -> None:
    """A displayed number is a property of the run, not of the reader's network.

    ``applyRounding`` rewrote every numeric cell of every initialised table to four
    decimal places in the DOM, so the same file showed ``0.01`` online and
    ``0.010012`` offline.
    """
    online_page = open_report(rendered_report, offline=False)
    _require_online(online_page, failed_requests, error_responses, offline=False)

    online: list[str] = online_page.eval_on_selector_all(KESTREL_ROW_SELECTOR, _DEPTH_SCORE_TEXT)
    offline: list[str] = open_report(rendered_report, offline=True).eval_on_selector_all(
        KESTREL_ROW_SELECTOR, _DEPTH_SCORE_TEXT
    )

    assert online, "no Depth Score cells were found online - the selector or the fixture is wrong"
    assert online == offline, f"the depth score depends on the network: online {online}, offline {offline}"
