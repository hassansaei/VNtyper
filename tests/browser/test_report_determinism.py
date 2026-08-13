"""The report must say the same thing to every reader who opens the same file.

This is the acceptance check for issue #242. A VNtyper report is a single HTML
file that gets forwarded, archived and reopened months later; two people looking
at the same file must see the same rows and the same numbers. Today they do not,
and the difference is not in the file - it is in whether the reader's browser
could reach three CDNs:

* **online**, jQuery and DataTables load, the custom search filter runs, and
  every row whose ``Flag`` is not ``Not flagged``/``Not applicable``/empty is
  removed from the DOM before the reader ever sees it;
* **offline**, none of that happens, so every row stays - and every number keeps
  the precision the pipeline wrote rather than the four decimal places
  ``applyRounding`` imposes.

``test_the_report_reads_identically_online_and_offline`` therefore **fails
today**, on purpose. It is the evidence the issue is real, and the check that
says when it is fixed. It is left failing rather than marked ``xfail``: a strict
xfail also swallows an *error*, so a browser tier that could not launch a browser
at all would report green, which is exactly the dishonest RED this tier exists to
prevent. Nothing gates on it - ``make check-all`` runs ``format-check lint
type-check-all test-unit check-integration-compatibility`` and CI runs
``pytest -m unit tests/unit``, so neither reaches ``tests/browser``.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import CLEAN_FLAG, EXPECTED_KESTREL_ROWS, HIDDEN_FLAGS

pytestmark = pytest.mark.browser


def test_the_fixture_is_a_specimen_that_can_show_the_divergence(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
) -> None:
    """Without JavaScript, the report shows every Kestrel row the pipeline wrote.

    This guards the check below from passing vacuously. If the fixture ever stops
    producing rows, or stops producing a row the shipped filter hides, then
    ``online == offline`` becomes true for a reason that has nothing to do with
    the report being fixed.
    """
    rows = visible_kestrel_rows(open_report(rendered_report, offline=True))

    assert len(rows) == EXPECTED_KESTREL_ROWS, (
        f"the specimen must render {EXPECTED_KESTREL_ROWS} Kestrel rows with no JavaScript, got {rows}"
    )
    assert any(CLEAN_FLAG in row for row in rows), f"no unflagged row to survive the online filter: {rows}"
    for flag in HIDDEN_FLAGS:
        assert any(flag in row for row in rows), f"the specimen no longer carries the hidden flag {flag!r}: {rows}"


def test_the_report_reads_identically_online_and_offline(
    rendered_report: Path,
    open_report: Callable[..., Page],
    visible_kestrel_rows: Callable[[Page], list[str]],
    failed_requests: dict[Page, list[str]],
) -> None:
    """One file, two readers, two different answers - issue #242.

    Fails today. The online reader is shown one of three Kestrel rows, with its
    depth score rounded in the DOM; the offline reader is shown all three, at the
    precision the pipeline recorded. Neither is told anything is missing.
    """
    online_page = open_report(rendered_report, offline=False)
    online = visible_kestrel_rows(online_page)
    offline = visible_kestrel_rows(open_report(rendered_report, offline=True))

    # A machine with no network turns the online pass into a second offline pass,
    # and the comparison below then compares a report with itself.
    assert not failed_requests[online_page], (
        "the online pass could not complete every request, so it was not online and this check "
        f"proves nothing: {failed_requests[online_page]}"
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
