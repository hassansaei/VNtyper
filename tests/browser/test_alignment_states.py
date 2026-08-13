"""The alignment panel is never blank space, in any state it can be in.

Until #242 a failed load threw before any fallback was reached, so the reader got a
panel-shaped hole with no explanation in it, and the console entry saying why is not
somewhere a reader of an archived report will ever look. There are five states now and
every one of them is authored:

1. **no alignment session** - the run produced none;
2. **``--report-igv off``** - the operator asked for none;
3. **``--report-igv sidecar``** - the browser is the file beside this one;
4. **``DecompressionStream`` missing** - the library is in this file and this engine
   cannot expand it (anything older than Chrome 80 / Safari 16.4 / Firefox 113);
5. **the expansion or the browser failed** - said with the error, in the page.

The first three are written into the *markup*, from facts Python has before the file is
written, which is why they survive a reader with scripting off. The last two can only
be known in the browser.

The success path is here too, and it is the one worth stating plainly: the embedded
library expands and defines ``window.igv`` **with the network blocked**. That is the
whole point of the change, and nothing else in this repository measures it.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest
from playwright.sync_api import Page

from tests.browser.conftest import external_requests

pytestmark = pytest.mark.browser

#: The authored heading of each state, as a literal. Not imported from anywhere: a test
#: that asks the template what it says agrees with whatever the template says, including
#: with nothing at all.
NO_SESSION = "No alignment visualisation is available for this sample."
SWITCHED_OFF = "Alignment visualisation was switched off for this run."
SIDECAR = "The alignment browser is a separate file."
CANNOT_EXPAND = "This browser cannot expand the embedded alignment browser."
COULD_NOT_START = "The alignment browser could not be started."
EXPANDING = "igv.js 3.0.2 is embedded in this file, compressed"

#: Every authored heading, so each test can assert that the state it expects is the one
#: the reader got *and* that no other state is showing beside it.
ALL_STATES = (NO_SESSION, SWITCHED_OFF, SIDECAR, CANNOT_EXPAND, COULD_NOT_START, EXPANDING)

#: Take the global away before the document's own script runs. This is the only way to
#: reach state 4 in a modern engine, and it is exactly what an engine from before 2023
#: presents to the page.
DELETE_DECOMPRESSION_STREAM = "delete window.DecompressionStream"

#: Whether the embedded library actually expanded into a usable global.
IGV_PRESENT = "() => typeof window.igv !== 'undefined' && typeof window.igv.createBrowser === 'function'"


def _panel(page: Page) -> str:
    """The alignment panel's visible text.

    Args:
        page: An already-loaded page.

    Returns:
        str: Whatever the reader can read in the panel, whitespace normalised.
    """
    return " ".join(page.locator("#igvContainer").inner_text().split())


def _assert_only_state(page: Page, expected: str) -> None:
    """Fail unless exactly the expected state is showing.

    Two halves, and the second is the one that catches a fallback printed
    unconditionally: a panel that always says "the viewer could not be loaded" would
    satisfy any single-state check and tell every reader the viewer had failed.

    Args:
        page: An already-loaded page.
        expected: The authored heading the reader should see.
    """
    shown = _panel(page)
    assert expected in shown, f"the reader is not told {expected!r}; the panel says {shown!r}"
    others = [state for state in ALL_STATES if state != expected and state in shown]
    assert others == [], f"the panel also claims {others}, which are different facts about the same run"


def test_a_run_with_no_alignment_session_says_so(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The shipped fixture: no BED, so no session. The commonest report there is."""
    page = open_report(rendered_report, offline=True)

    _assert_only_state(page, NO_SESSION)


def test_a_run_with_no_alignment_session_says_so_with_no_script_at_all(
    rendered_report: Path,
    browser,
) -> None:
    """And it says it in the markup, so scripting off is not blank space either.

    This is what "authored in the template rather than written by the script" buys, and
    it is measured rather than asserted about the source: a context with JavaScript
    disabled is the strongest available statement of "the file says this on its own".
    """
    context = browser.new_context(java_script_enabled=False)
    try:
        page = context.new_page()
        page.goto(rendered_report.as_uri(), wait_until="load")

        assert NO_SESSION in _panel(page)
    finally:
        context.close()


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("off", SWITCHED_OFF), ("sidecar", SIDECAR)],
)
def test_each_mode_that_carries_no_library_says_which_one_it_is(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
    mode: str,
    expected: str,
) -> None:
    """A run *with* alignments, in the two modes that deliberately embed nothing.

    Both are states the operator chose, and telling a reader "no alignment
    visualisation is available" would report a choice as an absence.
    """
    page = open_report(rendered_report_with_alignments(mode), offline=True)

    _assert_only_state(page, expected)


def test_an_engine_without_decompression_stream_is_told_why(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
) -> None:
    """State 4, the one that used to be blank space.

    Deleting the global before load is what an engine older than Chrome 80 / Safari
    16.4 / Firefox 113 presents to this page. Before #242 the load threw a
    ``ReferenceError`` and ``#igvContainer`` kept an empty ``<div id="igvDiv">``: a
    panel-shaped hole, and the reader had no way to know whether the alignment view was
    missing because the sample had none or because their browser was too old.
    """
    page = open_report(
        rendered_report_with_alignments("embedded"),
        offline=True,
        before_load=DELETE_DECOMPRESSION_STREAM,
    )

    assert page.evaluate("() => typeof DecompressionStream === 'undefined'"), (
        "DecompressionStream is still defined, so this pass never reaches the branch under test"
    )
    _assert_only_state(page, CANNOT_EXPAND)

    shown = _panel(page)
    assert "Chrome 80" in shown and "Safari 16.4" in shown and "Firefox 113" in shown, (
        f"the reader is not told which engines can expand it: {shown!r}"
    )
    assert "Nothing has been hidden" in shown, (
        f"the reader is not told the tables above are complete without it: {shown!r}"
    )


def test_the_variant_table_survives_an_engine_that_cannot_expand_the_library(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
) -> None:
    """The evidence beside the frame is not collateral damage.

    The variant list is built from a literal already in the document, so it has no
    business depending on the library expanding - and it used to, because the whole
    block threw at the first line.
    """
    page = open_report(
        rendered_report_with_alignments("embedded"),
        offline=True,
        before_load=DELETE_DECOMPRESSION_STREAM,
    )

    assert page.locator("#variant_table tbody tr").count() == 2, "the variant table was lost with the viewer"
    assert page.locator("#kestrel_table tbody tr").count() == 3, "the Kestrel table was lost with the viewer"


def test_the_embedded_library_expands_with_the_network_blocked(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
) -> None:
    """The success path, offline. This is the change, measured.

    ``window.igv`` exists because 497 KB of base64 in this file was decoded, gunzipped
    and evaluated - with every request that is not a local file read aborted at the
    route. Nothing else in this repository says that the vendored payload is a working
    library rather than merely 497 KB of the right length.
    """
    page = open_report(rendered_report_with_alignments("embedded"), offline=True)
    page.wait_for_function(IGV_PRESENT, timeout=15_000)

    assert page.evaluate(IGV_PRESENT)
    assert external_requests(page) == [], "expanding the embedded library reached off-machine"


def test_a_browser_that_cannot_start_says_so_rather_than_leaving_the_placeholder(
    rendered_report_with_alignments: Callable[..., Path],
    open_report: Callable[..., Page],
) -> None:
    """State 5. The session file this tier produces does not exist, so this is real.

    ``create_report`` is not run here, so ``session0.json`` is named by the extracted
    session dictionary and is not on disk. ``igv.createBrowser`` therefore rejects, and
    what the reader must not be left with is the "expanding now" placeholder, forever,
    with the reason in a console they will never open.
    """
    page = open_report(rendered_report_with_alignments("embedded"), offline=True)
    page.wait_for_function(
        f"() => document.getElementById('igvState') === null || "
        f"document.getElementById('igvState').innerText.includes({COULD_NOT_START!r})",
        timeout=15_000,
    )

    _assert_only_state(page, COULD_NOT_START)


def test_the_placeholder_is_what_the_reader_sees_first(
    rendered_report_with_alignments: Callable[..., Path],
    browser,
) -> None:
    """Expansion is asynchronous, so the first paint is a state and needs authoring.

    Measured with scripting off, which is the same document the reader sees for the
    interval before the script has finished - and a permanent state for a reader who
    has scripting off for good. It is deliberately not a spinner: a spinner promises
    that something is coming, and this sentence says what is here and what it costs.
    """
    report = rendered_report_with_alignments("embedded")
    context = browser.new_context(java_script_enabled=False)
    try:
        page = context.new_page()
        page.goto(report.as_uri(), wait_until="load")
        shown = _panel(page)

        assert EXPANDING in shown, f"the alignment frame's first paint is not authored: {shown!r}"
        assert "makes no network request" in shown
    finally:
        context.close()
