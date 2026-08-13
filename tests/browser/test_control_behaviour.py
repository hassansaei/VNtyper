"""A control that focuses and does nothing is not operable.

``test_keyboard_operability.py`` proves every control can be *reached*. That is half
of a control. This file is the other half: what happens when it is used.

The behaviour here was fixed for a reason that has since gone away, and it is worth
keeping the reason on the record because the fix outlives it. "Show detailed coverage"
used to live in the jQuery ``<script>`` block, and a ``$`` that never loaded throws at
the top of whichever block it is in and takes every statement after it down with it -
so offline the switch focused, ticked, and did nothing. That matters beyond hygiene:
the detailed view is the **only** route to the Coverage QC verdict, the minimum and
maximum, and the uncovered-base count.

There is no jQuery left to throw (#242), so the failure mode this guards against is now
a different one: a handler that shares a block with something that throws for any other
reason. The check is unchanged, and it is still not satisfied by a fix whose removal is
silent - moving the handler back in with the sort code would leave every other suite in
this repository green.

The alignment panel's states used to be tested here too. They are now
``tests/browser/test_alignment_states.py``, which covers all five of them rather than
the two that existed when this file was written.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import external_requests

pytestmark = pytest.mark.browser

#: Evidence that exists only in the detailed coverage view. If the switch does not
#: swap the views, this is what the reader never sees.
DETAILED_ONLY_EVIDENCE = ("Coverage QC", "Minimum Coverage", "Uncovered Bases")

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
        # The online pass used to prove it was online by asking whether jQuery had
        # initialised. There is nothing left to load, so the honest claim is the
        # opposite one, read off the same evidence: the document asked for nothing.
        assert external_requests(page) == [], "the report fetched something off-machine"

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
