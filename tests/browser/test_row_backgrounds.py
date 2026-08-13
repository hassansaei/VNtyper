"""Row backgrounds, measured. Two findings whose fixes had no test at all.

**No zebra in the per-sample report.** It was removed by not emitting Bootstrap's
``table-striped`` class, and a comment asserted it - but ``git grep table-striped``
across ``tests/`` returned nothing and no test measured a row background anywhere, so
re-adding the class passed the entire suite. The class is now meaningless without
Bootstrap, which makes the claim harder to break by accident and *impossible* to check
by reading the class list: whether rows alternate is now a property of this
repository's own stylesheet. So it is measured, in the only place a computed background
exists.

The rule this protects is not decorative. The per-sample tables have exactly one row
treatment that carries meaning - the flagged-value highlight - and a stripe underneath
it is a second, meaningless one competing for the same signal. It was also the
``#f2f2f2`` that put the old inline ``orange`` at 1.76:1.

**The cohort keeps its stripes**, and that is asserted here too, in the same terms.
The two reports disagree on purpose, and a test that only pinned one of them would let
the other drift into matching it.

**The flag glyph has a size of its own.** A recent change dropped ``font-size: 16px``
from the mark and left ``.flag-glyph`` setting only ``font-weight``, so the tick and
the cross rendered at the table's 13px. Nothing measured it, because a font size is
exactly the kind of thing a source scan reads as present when it is inherited from
somewhere else entirely.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest
from playwright.sync_api import Page

from vntyper.cli import load_config
from vntyper.scripts.cohort_summary import generate_cohort_summary_report

pytestmark = pytest.mark.browser

#: Read the computed background of every row of one table.
_ROW_BACKGROUNDS = "els => els.map(e => getComputedStyle(e).backgroundColor)"

#: Read the computed font size of every flag mark.
_GLYPH_SIZES = "els => els.map(e => getComputedStyle(e).fontSize)"

#: The token layer's ``--step-0``. 16px, the report's body size, and what the mark was
#: set to inline before the size was lost.
EXPECTED_GLYPH_SIZE = "16px"


@pytest.fixture
def rendered_cohort_report(tmp_path: Path) -> Path:
    """Render a real cohort report with four Kestrel rows.

    Four rather than two, so "alternating" is distinguishable from "the second row is
    different": a two-row table cannot tell a zebra from an off-by-one.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        Path: The rendered ``cohort_summary.html``.
    """
    generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=pd.DataFrame(
            [
                {"Sample": f"sample_{index}", "Motif": str(index), "Confidence": "High_Precision", "Flag": flag}
                for index, flag in enumerate(("Not flagged", "Not flagged", "Not flagged", "Not flagged"))
            ]
        ),
        advntr_df=pd.DataFrame([{"Sample": "sample_0", "VID": "25561", "Flag": "Not flagged"}]),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )
    return tmp_path / "cohort_summary.html"


def test_the_per_sample_rows_all_have_the_same_background(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """Three Kestrel rows, one background. This is the finding, pinned.

    Measured rather than read off the class list, because the class list stopped being
    the answer when Bootstrap left: ``table-striped`` could be re-added tomorrow and
    would draw nothing until a rule for it appeared in the per-sample report's own
    stylesheet. What a reader sees is the computed value, so that is what is compared.
    """
    page = open_report(rendered_report, offline=True)

    backgrounds: list[str] = page.eval_on_selector_all("#kestrel_table tbody tr", _ROW_BACKGROUNDS)

    assert len(backgrounds) == 3, f"expected the fixture's three Kestrel rows, measured {backgrounds}"
    assert len(set(backgrounds)) == 1, (
        f"the per-sample rows alternate: {backgrounds}. A stripe competes with the flagged-value "
        "highlight, which is the one row treatment in these tables that carries meaning."
    )


def test_the_per_sample_table_carries_no_striping_class(
    rendered_report: Path,
) -> None:
    """The other half, in the source, so the failure names the cause rather than a colour.

    Kept alongside the measurement rather than instead of it: this one says *why* the
    backgrounds are equal and fails the moment somebody re-adds the class, before any
    rule exists to make the measurement notice.

    The scan is over ``<table>`` tags, not over the whole document: the shared token
    layer explains at length why the per-sample report has no stripe and why the cohort
    does, and the word appears in that prose. Scanning it would make an accurate comment
    fail the build and reward deleting it.
    """
    tags = re.findall(r"<table\b[^>]*>", rendered_report.read_text(encoding="utf-8"))

    assert tags, "the report rendered no tables at all, so this guard would pass vacuously"
    striped = [tag for tag in tags if "table-striped" in tag]
    assert striped == [], f"a per-sample table carries the striping class: {striped}"


def test_the_cohort_rows_alternate(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """And the cohort's do, deliberately.

    Reading along a wide row across many samples is what a stripe is for, and the
    cohort's flagged rows are filtered rather than tinted, so nothing competes with it.
    The rule is in the cohort template's own stylesheet now - Bootstrap used to draw it.
    """
    page = open_report(rendered_cohort_report, offline=True)

    backgrounds: list[str] = page.eval_on_selector_all("table tbody tr", _ROW_BACKGROUNDS)

    assert len(backgrounds) >= 4, f"the cohort fixture rendered {len(backgrounds)} rows; four are needed"
    assert backgrounds[0] != backgrounds[1], f"the cohort report lost its stripes: {backgrounds}"
    assert backgrounds[0] == backgrounds[2], f"the cohort stripe is not alternating: {backgrounds}"
    assert backgrounds[1] == backgrounds[3], f"the cohort stripe is not alternating: {backgrounds}"


def test_the_cohort_flag_glyph_has_an_intentional_size(
    rendered_cohort_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """The glyph is 16px, not the 13px it inherits from the table.

    The size used to be an inline ``font-size: 16px`` on the element the script built.
    It was lost when the glyph's colour moved out of the script and into the token
    layer's classes: the weight came across and the size did not, so the mark quietly
    shrank to the table's own step. It is in the token layer now, in the scale, for both
    reports.
    """
    page = open_report(rendered_cohort_report, offline=True)

    sizes: list[str] = page.eval_on_selector_all(".flag-glyph", _GLYPH_SIZES)

    assert sizes, "no flag mark was built, so this measures nothing"
    assert set(sizes) == {EXPECTED_GLYPH_SIZE}, (
        f"the flag glyph renders at {sorted(set(sizes))} rather than {EXPECTED_GLYPH_SIZE}; it is inheriting "
        "the table's font size instead of carrying one of its own"
    )


def test_the_per_sample_flag_glyph_has_the_same_size(
    rendered_report: Path,
    open_report: Callable[..., Page],
) -> None:
    """One rule, both reports. The per-sample mark is built by Python and the cohort's
    by a script, so the same declaration has to reach two different element shapes."""
    page = open_report(rendered_report, offline=True)

    sizes: list[str] = page.eval_on_selector_all("#kestrel_table .flag-glyph", _GLYPH_SIZES)

    assert sizes, "no flag mark was rendered, so this measures nothing"
    assert set(sizes) == {EXPECTED_GLYPH_SIZE}, f"the per-sample flag glyph renders at {sorted(set(sizes))}"
