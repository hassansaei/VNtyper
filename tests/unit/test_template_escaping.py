"""Neither report template may build markup by concatenating a sample-derived value.

SPECIFICATION (#207, #216): the server escapes every ``Flag`` cell and every value
reaching the two IGV table cells, and until this milestone the client-side
JavaScript undid that escaping in two different ways:

* ``updateFlagColumn()`` -- present, verbatim, in both templates -- read the cell
  back with ``.text()`` (which *decodes* the server's escaping), concatenated the
  decoded value into a ``title="..."`` attribute inside an HTML string, and handed
  that string to ``.html()``, which reparses it as markup. A ``Flag`` of
  ``"></span><img src=x onerror=alert(1)><span title="`` renders safely in the
  initial document and executes on the first DataTables redraw.
* the IGV variant table (``report_template.html`` only) wrote each igv-reports
  cell with ``cell.innerHTML = ...``, which parses its argument as markup instead
  of writing it as a text node.

This is a **tripwire over the template source text**, not a behavioural test. The
unit tier has no JavaScript engine -- no network, no browser, no node (AGENTS.md)
-- so nothing here can *run* either template and prove the page is safe. What this
module proves is that the exact sinks which made the page unsafe are gone and have
not come back in one of a few adjacent shapes (a different quote style, a template
literal, ``insertAdjacentHTML``, a bare ``innerHTML``/``outerHTML`` assignment).
It **cannot** prove there is no other, differently-shaped sink somewhere else in
either file, and it cannot prove anything about runtime behaviour at all -- this
repo has no browser test tier, and AGENTS.md forbids adding one at the unit tier.
The residual risk is a future sink whose shape these patterns do not match.
"""

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

#: Every template file this package ships, the shared token layer included -- it is
#: template source and no linter reaches it either, so the checks below apply to it as
#: much as to the two reports that ``{% include %}`` it.
TEMPLATES = (
    Path("vntyper/templates/_report_base.html"),
    Path("vntyper/templates/cohort_summary_template.html"),
    Path("vntyper/templates/report_template.html"),
)

#: A value spliced into a markup string that is then handed to ``.html()`` or
#: ``insertAdjacentHTML()`` for (re)parsing -- the #207 sink, plus the adjacent
#: shapes the original single-quote-only regex missed (Finding 3): double-quoted
#: concatenation, a template literal argument, and ``insertAdjacentHTML``.
_UNSAFE_HTML_SINK = re.compile(
    r"""
    \.html\(\s*(['"])(?:(?!\1).)*?\1\s*\+   # .html('...' + or .html("..." +, quotes matched consistently
    | \.html\(\s*`                          # .html(`template literal`)
    | \.insertAdjacentHTML\(                # insertAdjacentHTML(...)
    """,
    re.VERBOSE,
)


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_template_concatenates_a_variable_into_parsed_markup(template: Path) -> None:
    source = template.read_text(encoding="utf-8")
    offenders = _UNSAFE_HTML_SINK.findall(source)
    assert offenders == [], f"{template} splices a variable into markup that gets reparsed: {offenders}"


#: Templates that still build the flag cell in the browser. The per-sample report no
#: longer does: #242 moved that rendering to ``report_formatting.flag_html``, so the
#: reason is text in the file rather than a tooltip assembled by a script, and the
#: two checks below have nothing left to look at there. They are not weakened for the
#: cohort report, which keeps ``updateFlagColumn`` and therefore keeps the sink.
CLIENT_SIDE_FLAG_TEMPLATES = (Path("vntyper/templates/cohort_summary_template.html"),)


@pytest.mark.parametrize("template", CLIENT_SIDE_FLAG_TEMPLATES, ids=lambda p: p.name)
def test_the_flag_tooltip_title_is_set_through_attr(template: Path) -> None:
    source = template.read_text(encoding="utf-8")
    assert ".attr('title', originalText)" in source, (
        f"{template} must set the tooltip title with .attr(), which takes a value rather than a fragment"
    )


@pytest.mark.parametrize("template", CLIENT_SIDE_FLAG_TEMPLATES, ids=lambda p: p.name)
def test_the_flag_cell_is_emptied_before_the_mark_is_appended(template: Path) -> None:
    """``.empty().append()`` replaces the cell's content without an HTML parse."""
    source = template.read_text(encoding="utf-8")
    assert "$flagCell.empty().append(" in source, f"{template} must rebuild the flag cell with DOM APIs"


def test_the_per_sample_flag_cell_is_no_longer_built_in_the_browser() -> None:
    """The other half of narrowing the two checks above: the sink is gone, not moved.

    ``report_formatting.flag_html`` escapes the value and the table names ``Flag`` in
    ``html_columns``, so ``escape_frame_cells`` escapes every other cell -- the same
    per-column exemption the ``Confidence`` span already had.
    """
    source = Path("vntyper/templates/report_template.html").read_text(encoding="utf-8")

    assert "updateFlagColumn" not in source
    assert "data-original" not in source
    assert "tooltip" not in source


def test_the_igv_variant_table_writes_text_not_markup() -> None:
    """#216: both cells carry igv-reports values derived from the sample's reads.

    ``innerHTML`` parses; ``textContent`` does not. Neither cell needs markup.
    Verified this is not a display regression: VNtyper's BED is three columns
    (``kestrel_genotyping.py:682-686``), so igv-reports' ``has_name`` is false and
    the one column it pre-escapes (``bedtable.py:32``) is never rendered here.
    """
    source = Path("vntyper/templates/report_template.html").read_text(encoding="utf-8")
    assert "cell.innerHTML" not in source
    assert "cell.textContent = headers[j];" in source
    assert "cell.textContent = rowData[j];" in source


#: ``innerHTML``/``outerHTML`` assigned from anything other than a plain string
#: literal (Finding 3) -- the general shape of the sink A3 removes. A plain
#: literal is deliberately left alone: ``report_template.html`` sets
#: ``#igvContainer``'s ``innerHTML`` to a fixed, VNtyper-authored ``"<p>...`` string
#: with no interpolation, which is not a sink because nothing sample-derived ever
#: reaches it.
_UNSAFE_INNER_OUTER_HTML_ASSIGNMENT = re.compile(r"\b(?:innerHTML|outerHTML)\s*=(?!\s*['\"])")


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_inner_or_outer_html_is_assigned_from_a_non_literal(template: Path) -> None:
    source = template.read_text(encoding="utf-8")
    offenders = _UNSAFE_INNER_OUTER_HTML_ASSIGNMENT.findall(source)
    assert offenders == [], f"{template} assigns innerHTML/outerHTML from something that is not a string literal"


#: AGENTS.md sets line length 120 for this repository. Ruff enforces it for `.py` only -
#: these two files are HTML, so nothing gated them and one line had reached 127 characters.
TEMPLATE_LINE_LENGTH = 120


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_template_line_exceeds_the_repository_line_length(template: Path) -> None:
    """The templates carry the report's JavaScript, and no linter reaches it.

    Trailing whitespace is stripped before measuring, matching how ruff measures a
    Python line: a line that is only long because of trailing blanks is a whitespace
    problem, not a length one, and reporting it here would be noise.
    """
    over = [
        (number, len(line.rstrip()))
        for number, line in enumerate(template.read_text(encoding="utf-8").splitlines(), start=1)
        if len(line.rstrip()) > TEMPLATE_LINE_LENGTH
    ]
    assert over == [], f"{template} has lines over {TEMPLATE_LINE_LENGTH} characters: {over}"
