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

from tests.browser.conftest import (
    CLEAN_FLAG,
    EXPECTED_KESTREL_ROWS,
    HIDDEN_FLAGS,
    KESTREL_ROW_SELECTOR,
    external_requests,
    kestrel_column_text,
)

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

#: The first cell of each Kestrel row (``Motif``). Still first after the column-order
#: change, and still the column the table is read from left to right by.
_MOTIF_SELECTOR = "#kestrel_table tbody tr td:first-child"

#: The ``Depth Score`` column, keyed on its heading. It used to be
#: ``td:nth-child(9)``, which is ``Confidence`` now that ``Motif_sequence`` is last -
#: a positional selector goes on matching a cell after the column it named has moved,
#: so it fails on the value rather than on the position and reads as a data defect.
_DEPTH_SCORE_TEXT = kestrel_column_text("Depth Score")

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
    depth_scores: list[str] = page.eval_on_selector_all(KESTREL_ROW_SELECTOR, _DEPTH_SCORE_TEXT)
    assert depth_scores == list(_EXPECTED_DEPTH_SCORES)


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

    # **This is the successor guard the previous version of this comment described.**
    #
    # Proving the online pass was really online used to take a positive check that
    # jQuery and DataTables had initialised: the absence of a transport failure is not
    # one, because `requestfailed` fires for a dropped connection and a blocked route
    # and not for a 404, a 503, or a 200 carrying unusable JavaScript. In any of those
    # the third-party scripts never initialised, both passes rendered the same raw
    # table, and `online == offline` passed against a report nobody had checked.
    #
    # The CDN tags are gone (#242) and there is no jQuery left to ask, so the guard is
    # the opposite claim over the same evidence: the online pass issued **zero**
    # non-`file://` requests. That is a stronger statement than the one it replaces -
    # it is about the document rather than about the network, and it holds in the
    # offline pass too, because Playwright records a request before the route handler
    # can abort it. `tests/browser/test_offline_document.py` proves the recorder is
    # not vacuous by pointing a control document at a host and watching it be seen.
    assert external_requests(online_page) == [], (
        f"the online pass fetched something off-machine: {external_requests(online_page)}"
    )
    assert not failed_requests[online_page], (
        f"the online pass could not complete every request: {failed_requests[online_page]}"
    )
    assert not error_responses[online_page], f"the online pass got error statuses back: {error_responses[online_page]}"

    # An empty selector would make `[] == []` pass against an unmodified report,
    # which is how this check silently stops testing anything.
    assert online, "the online pass found no Kestrel rows at all - the selector or the fixture is wrong"
    assert offline, "the offline pass found no Kestrel rows at all - the selector or the fixture is wrong"

    assert online == offline, (
        "the same report file reads differently depending on the reader's network:\n"
        f"  online  ({len(online)} rows): {online}\n"
        f"  offline ({len(offline)} rows): {offline}"
    )


@pytest.mark.parametrize("offline", (False, True), ids=("online", "offline"))
def test_fastp_labels_and_statuses_share_configured_cutoffs_in_a_browser(
    rendered_report_with_custom_fastp_cutoffs: Path,
    open_report: Callable[..., Page],
    offline: bool,
) -> None:
    """A reader sees each custom cutoff paired with its matching decision.

    The fixture puts every displayed rate exactly on its own non-default
    threshold. A stale label or a status icon judged with a different cutoff
    therefore makes this browser-level report contract fail.
    """
    page = open_report(rendered_report_with_custom_fastp_cutoffs, offline=offline)

    for label, cutoff in (
        ("Duplication Rate", "12.34%"),
        ("Q20 Rate", "75.55%"),
        ("Q30 Rate", "65.43%"),
        ("Passed Filter Rate", "77.65%"),
    ):
        row = page.locator("tr").filter(has_text=f"{label} (Cutoff: {cutoff})")
        assert row.count() == 1, f"the report has no {label!r} row with its configured {cutoff} cutoff"
        assert row.locator('[aria-label="No warning"]').count() == 1, (
            f"the {label!r} status does not pass at its configured {cutoff} cutoff"
        )


@pytest.mark.parametrize("offline", (False, True), ids=("online", "offline"))
def test_fastp_half_tie_values_and_icons_agree_in_a_browser(
    rendered_report_with_fastp_half_ties: Path,
    open_report: Callable[..., Page],
    offline: bool,
) -> None:
    """The visible half-tie text, cutoff, and decision remain one browser contract."""
    page = open_report(rendered_report_with_fastp_half_ties, offline=offline)

    for label, displayed in (
        ("Duplication Rate", "5.05%"),
        ("Q20 Rate", "77.65%"),
        ("Q30 Rate", "70.05%"),
        ("Passed Filter Rate", "80.05%"),
    ):
        row = page.locator("tr").filter(has_text=f"{label} (Cutoff: {displayed})")
        assert row.count() == 1, f"the report has no {label!r} half-tie row at {displayed}"
        assert row.locator("td").nth(1).inner_text().strip() == displayed
        assert row.locator('[aria-label="No warning"]').count() == 1


@pytest.mark.parametrize("offline", (False, True), ids=("online", "offline"))
def test_exact_fastp_json_boundaries_keep_visible_values_and_icons_together(
    rendered_report_with_exact_fastp_boundaries: Path,
    open_report: Callable[..., Page],
    offline: bool,
) -> None:
    """Exact JSON operands remain below their cutoffs in either browser mode."""
    page = open_report(rendered_report_with_exact_fastp_boundaries, offline=offline)

    for label, displayed, cutoff in (
        ("Q20 Rate", "60.04%", "60.05%"),
        ("Passed Filter Rate", "77.64%", "77.65%"),
    ):
        row = page.locator("tr").filter(has_text=f"{label} (Cutoff: {cutoff})")
        assert row.count() == 1, f"the report has no {label!r} row with its exact {cutoff} cutoff"
        assert row.locator("td").nth(1).inner_text().strip() == displayed
        assert row.locator('[aria-label="Warning"]').count() == 1
