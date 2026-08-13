"""The cohort report's flag mark must say, in words, what its glyph and colour say.

``updateFlagColumn`` replaces every ``Flag`` cell of the cohort report with a green
``✓`` or a red ``✖``. That is two visual channels and no textual one: to a screen
reader the cell announced as a bare character, and in greyscale the two glyph shapes
were the only difference left. The mark now carries ``role="img"`` and an
``aria-label`` holding the reason.

Why this is a browser test
--------------------------
The cell is rebuilt **in the browser**, so nothing about the rendered file says what
the reader ends up with. ``tests/unit/test_cohort_summary_oracle.py`` cannot see it
either - ``_skeleton()`` replaces every script body with ``<SCRIPT-BODY>`` before
hashing - so reverting the accessible name was silent in every tier this repository
had. This file is what makes it loud.

The cohort report keeps its filter, deliberately (#242, precondition P4): hiding
flagged rows is defensible for triage across samples and indefensible for a
single-patient read. So the flagged specimen is reached by ticking the switch, which
is also the only way a reader reaches it.

**What changed when DataTables left.** The mark used to be built on DataTables'
``preDrawCallback``, which meant it existed only for a reader whose browser had
reached three CDNs - so "offline" was the no-script case and the last test here used
it as one. The vanilla script runs for every reader, online or offline, so the
no-script case is now stated as what it always meant: a context with JavaScript
disabled. That is a stricter test than the one it replaces, because it no longer
depends on a network condition to stand in for a scripting condition.
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

#: The name given to a ``Flag`` cell that carries no value at all. It has to be a
#: literal rather than the empty string: ``role="img"`` with an empty ``aria-label``
#: is an image with **no accessible name**, which a screen reader skips entirely -
#: strictly worse than the bare ``✓`` it replaced, because that at least announced
#: as a character.
UNRECORDED_FLAG_NAME = "No flag recorded"

#: Three ``Flag`` states, because they are announced three different ways. The empty
#: one is the case the first version of this change got wrong.
CLEAN_FLAG = "Not flagged"
FLAGGED_REASON = "Low_Coverage"
UNRECORDED_FLAG = ""

#: A two-column cohort is enough: ``updateFlagColumn`` keys off the *last* column
#: being headed ``Flag`` and does nothing else with the frame's shape.
KESTREL_ROWS = pd.DataFrame(
    [
        {"Sample": "sample_one", "Motif": "5", "Confidence": "High_Precision", "Flag": CLEAN_FLAG},
        {"Sample": "sample_two", "Motif": "6", "Confidence": "Low_Precision", "Flag": FLAGGED_REASON},
        {"Sample": "sample_three", "Motif": "7", "Confidence": "High_Precision", "Flag": UNRECORDED_FLAG},
    ]
)

#: adVNTR gets one clean row, so the second table proves the treatment is not
#: specific to the Kestrel table.
ADVNTR_ROWS = pd.DataFrame([{"Sample": "sample_one", "VID": "25561", "Flag": CLEAN_FLAG}])

#: Every rebuilt flag cell holds exactly one ``<span>``, and the selector below
#: identifies it by **position** rather than by any of the attributes under test.
#: That is deliberate: selecting on ``[role='img']`` would make a mark that lost its
#: role invisible to the query, and the resulting "no marks were built" would blame
#: DataTables for a defect in the mark. Restricting to the last cell keeps the
#: colour-coded ``Confidence`` span out of the result.
MARK_SELECTOR = "table tbody td:last-child span"

#: Read every flag mark the reader can actually see, with the row it sits in for the
#: failure message.
#:
#: **Visible, not merely present.** DataTables took a filtered row out of the DOM, so
#: counting elements and counting what a reader sees were the same number; the vanilla
#: filter sets ``hidden`` on the row and leaves it in the document, which is a better
#: implementation and a worse thing to count. The ``offsetParent`` check is what keeps
#: this measuring the reader's experience rather than the DOM's contents.
_MARKS = """
els => els.filter(e => e.offsetParent !== null).map(e => ({
    row: e.closest('tr').textContent.replace(/\\s+/g, ' ').trim(),
    role: e.getAttribute('role') || '',
    label: e.getAttribute('aria-label'),
    glyph: e.textContent.trim(),
}))
"""

#: The ``Flag`` cell as the server wrote it, before any script touched it.
_LAST_CELL_TEXT = """
els => els.map(e => {
    const cells = e.querySelectorAll('td');
    return cells.length ? cells[cells.length - 1].textContent.trim() : '';
})
"""


@pytest.fixture
def rendered_cohort_report(tmp_path: Path) -> Path:
    """Render a real cohort report through the real renderer.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``cohort_summary.html``.
    """
    generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=KESTREL_ROWS.copy(),
        advntr_df=ADVNTR_ROWS.copy(),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )
    return tmp_path / "cohort_summary.html"


def _marks(page: Page) -> list[dict[str, str]]:
    """Return every flag mark in the page.

    Args:
        page: An already-loaded page whose scripts have run.

    Returns:
        list[dict[str, str]]: One mapping per mark, with its row text, its
        ``aria-label`` and its glyph.
    """
    marks: list[dict[str, str]] = page.eval_on_selector_all(MARK_SELECTOR, _MARKS)
    return marks


def _require_the_marks_were_built(page: Page, expected: int) -> list[dict[str, str]]:
    """Fail unless the browser actually rebuilt the flag column.

    A mark existing at all means the exact code path under test executed, so this is
    what stops every assertion below from holding vacuously over an empty list.

    Args:
        page: An already-loaded page.
        expected: How many marks the specimen should have produced.

    Returns:
        list[dict[str, str]]: The marks.
    """
    marks = _marks(page)
    assert len(marks) == expected, (
        f"expected {expected} flag marks, found {len(marks)}: {marks}. The script did not rebuild the flag "
        "column, so nothing below tests the accessible name."
    )
    return marks


def test_every_flag_mark_has_an_accessible_name(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """No mark may be an image with no name, including the empty-``Flag`` case.

    Two rows are visible by default - the clean one and the one with no ``Flag``
    value at all - plus adVNTR's clean row. ``isClean`` treats an empty string the
    same as ``Not flagged``, so the empty cell got a ``✓`` and, before this was
    fixed, ``aria-label=""``: an image a screen reader announces as nothing.
    """
    page = open_report(rendered_cohort_report, offline=False)
    marks = _require_the_marks_were_built(page, expected=3)

    unexposed = [mark for mark in marks if mark["role"] != "img"]
    assert unexposed == [], (
        f"these flag marks are not exposed as images, so their aria-label is ignored entirely: {unexposed}"
    )

    unnamed = [mark for mark in marks if not (mark["label"] or "").strip()]
    assert unnamed == [], f"these flag marks are images with no accessible name: {unnamed}"


def test_the_mark_is_announced_as_the_reason_it_stands_for(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The name is the ``Flag`` value, so colour and glyph carry nothing extra.

    The flagged row is behind the switch, because the cohort report filters it out
    by default and that is the reviewed behaviour for triage across samples. It is
    ticked here for the same reason a reader would tick it.
    """
    page = open_report(rendered_cohort_report, offline=False)
    _require_the_marks_were_built(page, expected=3)

    page.check("#toggleFlagged")
    marks = _require_the_marks_were_built(page, expected=4)

    by_label = {mark["label"]: mark for mark in marks}
    assert CLEAN_FLAG in by_label, f"no mark is announced as {CLEAN_FLAG!r}: {marks}"
    assert FLAGGED_REASON in by_label, f"no mark is announced as {FLAGGED_REASON!r}: {marks}"
    assert UNRECORDED_FLAG_NAME in by_label, f"no mark names the unrecorded flag: {marks}"

    # The glyph and the name have to agree, or the colour says one thing and the
    # screen reader another.
    assert by_label[CLEAN_FLAG]["glyph"] == "✓"
    assert by_label[FLAGGED_REASON]["glyph"] == "✖"


def test_the_reason_is_still_readable_when_no_script_ran(
    rendered_cohort_report: Path,
    browser,
) -> None:
    """The mark is an enhancement over text, not the only carrier of the reason.

    With scripting off nothing rebuilds the cell, so it keeps what the server wrote
    into it. If that were not the reason in words, the accessible name above would be
    papering over a cell that says nothing at all to a reader whose browser did not
    run the script - and it would also mean every row of the file was unreadable in
    a text extractor, which is how an archived cohort table gets read years later.

    This used to be measured *offline*, because offline was when jQuery failed to
    arrive and the rebuild never happened. The rebuild no longer depends on the
    network, so the condition is now stated directly.
    """
    context = browser.new_context(java_script_enabled=False)
    try:
        page = context.new_page()
        page.goto(rendered_cohort_report.as_uri(), wait_until="load")

        assert _marks(page) == [], "a mark was built with scripting off, so this is not the no-script case"

        cells: list[str] = page.eval_on_selector_all("table tbody tr", _LAST_CELL_TEXT)
        assert CLEAN_FLAG in cells, f"the clean reason is not readable without a script: {cells}"
        assert FLAGGED_REASON in cells, f"the flag reason is not readable without a script: {cells}"
    finally:
        context.close()
