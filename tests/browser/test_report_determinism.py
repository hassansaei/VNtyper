"""The report must say the same thing to every reader who opens the same file.

This is the standing acceptance check for issue #242. A VNtyper report is a single
HTML file that gets forwarded, archived and reopened months later; two people
looking at the same file must see the same rows and the same numbers. They did
not, and the difference was not in the file - it was in whether the reader's
browser could reach three CDNs:

* **online**, jQuery and DataTables loaded, a custom search predicate ran, and
  every row whose ``Flag`` was not ``Not flagged``/``Not applicable``/empty was
  removed from the DOM before the reader ever saw it, while ``applyRounding``
  rewrote every numeric cell to four decimal places;
* **offline**, none of that happened, so every row stayed, at the precision the
  pipeline wrote.

Measured on the unmodified template: **1 of 3 rows online, 3 of 3 offline**, with
the surviving row's depth score reading ``0.01`` against ``0.010012`` in the file.
Both client-side passes are now gone and this check is **green**. It is not a
record of a past defect - it is what stops the next change reintroducing one, so
it must keep failing loudly for the same reason it once did.

Why it is not marked ``xfail``
-----------------------------
It passes, so a ``strict=True`` xfail would XPASS and fail the tier. Even while it
was red, an unqualified ``xfail`` would have been the wrong tool: it reports an
*expected failure* for any outcome that is not a pass, so the narrower the reason
is pinned the better, and a bare marker pins nothing. Nothing gates on this tier
either way - ``make check-all`` runs ``format-check lint type-check-all test-unit
check-integration-compatibility`` and CI runs ``pytest -m unit tests/unit``, so
neither reaches ``tests/browser``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import CLEAN_FLAG, EXPECTED_KESTREL_ROWS, HIDDEN_FLAGS, enhancement_state

pytestmark = pytest.mark.browser

#: The allowlist the removed client-side predicate used, and which the server-side
#: flag glyph still uses to decide clean from flagged. Spelled out here as a
#: literal on purpose: a guard that asks ``conftest`` which of its own values are
#: flagged is self-referential and passes whatever they contain - setting
#: ``HIDDEN_FLAGS = ("Not applicable", "")`` would satisfy every membership check
#: while hiding no row at all.
_CLEAN_FLAG_ALLOWLIST = frozenset({"", "Not flagged", "Not applicable"})

#: The specimen's ``Flag`` values, in the order the report renders its rows, as
#: literals. Cross-checked against ``conftest``'s constants below, so a change
#: there fails this guard rather than quietly redefining what it demands.
_EXPECTED_FLAG_REASONS = ("Not flagged", "Low_Depth", "Depth_Score_Below_Threshold")

#: The specimen's ``Motif`` column, in render order.
_EXPECTED_MOTIFS = ("5", "6", "7")

#: The specimen's ``Depth Score`` column, in render order and at the precision the
#: pipeline recorded. ``applyRounding`` used to show ``0.01`` for the first of
#: these to any reader whose browser reached the CDNs.
_EXPECTED_DEPTH_SCORES = ("0.010012", "0.008034", "0.006001")

#: The flag reason as words. ``report_formatting.flag_html`` renders a glyph span
#: and a reason span; this reads the reason, so the guard does not depend on which
#: glyph is chosen.
_FLAG_REASON_SELECTOR = "#kestrel_table tbody tr .flag-reason"

#: The first cell of each Kestrel row (``Motif``).
_MOTIF_SELECTOR = "#kestrel_table tbody tr td:first-child"

#: The ninth cell of each Kestrel row (``Depth Score``), 1-indexed by CSS.
_DEPTH_SCORE_SELECTOR = "#kestrel_table tbody tr td:nth-child(9)"

_CELL_TEXT = "els => els.map(e => e.textContent.replace(/\\s+/g, ' ').trim())"


def _texts(page: Page, selector: str) -> list[str]:
    """Read the normalised text of everything matching ``selector``.

    Args:
        page: An already-loaded page.
        selector: A CSS selector.

    Returns:
        list[str]: One whitespace-normalised string per match, in document order.
    """
    matched: list[str] = page.eval_on_selector_all(selector, _CELL_TEXT)
    return matched


def test_the_fixture_is_a_specimen_that_can_show_the_divergence(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
) -> None:
    """The specimen really is one the removed predicate would have cut down.

    This guards the check below from passing vacuously, and it does so without
    asking the specimen to grade itself: the allowlist is a literal here, the
    expected rows are literals here, and ``conftest``'s constants are checked
    against them rather than trusted.
    """
    # 1. The specimen's constants are the ones this guard was written against.
    assert (CLEAN_FLAG, *HIDDEN_FLAGS) == _EXPECTED_FLAG_REASONS, (
        f"conftest now describes the specimen as {(CLEAN_FLAG, *HIDDEN_FLAGS)}, which is not what this "
        f"guard checks ({_EXPECTED_FLAG_REASONS}). Update both together, deliberately."
    )
    assert len(_EXPECTED_FLAG_REASONS) == EXPECTED_KESTREL_ROWS

    # 2. Those values really do fall on the sides of the predicate the specimen
    #    needs them to. Judged against the literal allowlist, not against itself.
    assert _EXPECTED_FLAG_REASONS[0] in _CLEAN_FLAG_ALLOWLIST, (
        f"{_EXPECTED_FLAG_REASONS[0]!r} is not a clean value, so no row survives the predicate and the "
        "online/offline comparison has nothing to be different about"
    )
    for reason in _EXPECTED_FLAG_REASONS[1:]:
        assert reason not in _CLEAN_FLAG_ALLOWLIST, (
            f"{reason!r} is in the clean allowlist {sorted(_CLEAN_FLAG_ALLOWLIST)}, so the predicate would "
            "have hidden nothing and this specimen cannot exhibit the defect"
        )

    # 3. The report renders exactly that specimen, with no JavaScript involved.
    page = open_report(rendered_report, offline=True)
    rows = visible_kestrel_rows(page)

    assert len(rows) == EXPECTED_KESTREL_ROWS, (
        f"the specimen must render {EXPECTED_KESTREL_ROWS} Kestrel rows with no JavaScript, got {rows}"
    )
    assert _texts(page, _FLAG_REASON_SELECTOR) == list(_EXPECTED_FLAG_REASONS)
    assert _texts(page, _MOTIF_SELECTOR) == list(_EXPECTED_MOTIFS)
    assert _texts(page, _DEPTH_SCORE_SELECTOR) == list(_EXPECTED_DEPTH_SCORES)


def test_the_report_reads_identically_online_and_offline(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
    failed_requests: dict[Page, list[str]],
    error_responses: dict[Page, list[str]],
) -> None:
    """One file, two readers, one answer - issue #242.

    Green since the client-side row predicate and the client-side rounding pass
    were removed. Before that the online reader was shown one of three Kestrel
    rows, with its depth score rounded in the DOM, and was told nothing about it.
    """
    online_page = open_report(rendered_report, offline=False)
    online = visible_kestrel_rows(online_page)
    offline = visible_kestrel_rows(open_report(rendered_report, offline=True))

    # Proving the online pass was really online takes a POSITIVE check. The
    # absence of a transport failure is not one: `requestfailed` fires for a
    # dropped connection and a blocked route, and not for a 404, a 503, or a 200
    # carrying unusable JavaScript. In any of those the third-party scripts never
    # initialise, both passes render the same raw table, and `online == offline`
    # passes against a report nobody has checked. So: nothing failed, nothing came
    # back with an error status, and the enhancement is live in the page itself.
    #
    # SUCCESSOR, for the task that removes the CDN tags: once the report fetches
    # nothing, "the online pass was really online" is not a meaningful thing to
    # assert, and jQuery will not be there to assert it with. The guard that
    # replaces this one is the opposite claim - that the online pass issued ZERO
    # non-`file://` requests, which is the same evidence read the other way round.
    # It would fail today, which is why it is a note here and not code.
    assert not failed_requests[online_page], (
        f"the online pass could not complete every request, so it was not online: {failed_requests[online_page]}"
    )
    assert not error_responses[online_page], (
        f"the online pass got error statuses back, so its scripts did not run: {error_responses[online_page]}"
    )
    state = enhancement_state(online_page)
    assert state["jquery"] and state["dataTable"], (
        "the online pass loaded no working jQuery/DataTables, so it rendered the same unenhanced page the "
        f"offline pass does and this comparison proves nothing: {state}"
    )

    # An empty selector would make `[] == []` pass against an unmodified report,
    # which is how this check silently stops testing anything.
    assert online, "the online pass found no Kestrel rows at all - the selector or the fixture is wrong"
    assert offline, "the offline pass found no Kestrel rows at all - the selector or the fixture is wrong"

    assert online == offline, (
        "the same report file reads differently depending on the reader's network:\n"
        f"  online  ({len(online)} rows): {online}\n"
        f"  offline ({len(offline)} rows): {offline}"
    )
