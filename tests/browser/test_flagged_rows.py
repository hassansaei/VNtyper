"""A flagged variant row must reach every reader of the file, and say why it is flagged.

Issue #242 is named after this defect: a ``High_Precision_flagged`` report narrates
a flagged pathogenic variant in its screening summary and then renders a table the
reader cannot find it in. The row is in the file -- ``to_html`` wrote it -- but the
client-side DataTables search predicate removes it from the DOM before the reader
sees anything, unless they tick a "Show flagged values" switch they have no reason
to know exists.

That makes this defect invisible to the unit tier by construction: asserting the
flag reason is in the rendered HTML passes against the unmodified code, because the
hiding happens after the HTML is parsed. It is only observable as a **visible row
count after initialisation**, which is why these checks live here.

The specimen comes from ``conftest.py``: three Kestrel rows, one clean and two
carrying ``Flag`` values the shipped predicate hides. No shipped configuration
produces a flagged row -- ``report_config.json`` declares no ``flag_rules`` key at
all -- so the fixture states them directly rather than hoping a rule fires.

Measured before the fix: **1 of 3 rows online, 3 of 3 offline**.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import CLEAN_FLAG, EXPECTED_KESTREL_ROWS, HIDDEN_FLAGS, KESTREL_ROW_SELECTOR

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


def _require_online(failed: list[str], offline: bool) -> None:
    """Fail an online pass that could not reach the network.

    A machine with no network turns the online pass into a second offline pass, so
    a check that only ever ran offline would report green against the unmodified
    template.

    Args:
        failed: The URLs this page could not fetch.
        offline: Whether the pass was deliberately offline.
    """
    if not offline:
        assert not failed, f"the online pass could not complete every request, so it proves nothing: {failed}"


@pytest.mark.parametrize("offline", [False, True], ids=["online", "offline"])
def test_every_flagged_row_is_visible_after_initialisation(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
    failed_requests: dict[Page, list[str]],
    offline: bool,
) -> None:
    """Every Kestrel row the pipeline wrote is visible, whatever the reader's network.

    Fails before the fix in the ``online`` case only: the shipped
    ``$.fn.dataTable.ext.search`` predicate leaves 1 of 3 rows standing.
    """
    page = open_report(rendered_report, offline=offline)
    _require_online(failed_requests[page], offline)

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
    offline: bool,
) -> None:
    """The reason a row is flagged is text in the cell, not a hover-only attribute.

    ``updateFlagColumn`` used to replace the cell with a bare glyph carrying the
    reason in a Bootstrap ``title``, which is invisible in print, invisible to a
    screen reader once Bootstrap moves ``title`` into ``data-original-title``, and
    absent entirely when the script does not run.
    """
    page = open_report(rendered_report, offline=offline)
    _require_online(failed_requests[page], offline)

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
    offline: bool,
) -> None:
    """Ticking "Highlight flagged values" changes styling and nothing else.

    The switch is operated through its ``<label>`` rather than the input: the
    template hides every bare checkbox so that its CSS-only accordions work, and a
    label click is how the reader operates it in both passes.
    """
    page = open_report(rendered_report, offline=offline)
    _require_online(failed_requests[page], offline)
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
) -> None:
    """A displayed number is a property of the run, not of the reader's network.

    ``applyRounding`` rewrote every numeric cell of every initialised table to four
    decimal places in the DOM, so the same file showed ``0.01`` online and
    ``0.010012`` offline.
    """
    online_page = open_report(rendered_report, offline=False)
    _require_online(failed_requests[online_page], offline=False)

    online: list[str] = online_page.eval_on_selector_all(KESTREL_ROW_SELECTOR, _DEPTH_SCORE_TEXT)
    offline: list[str] = open_report(rendered_report, offline=True).eval_on_selector_all(
        KESTREL_ROW_SELECTOR, _DEPTH_SCORE_TEXT
    )

    assert online, "no Depth Score cells were found online - the selector or the fixture is wrong"
    assert online == offline, f"the depth score depends on the network: online {online}, offline {offline}"
