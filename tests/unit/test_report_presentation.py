"""The report's presentation layer, asserted over the shipped templates.

This is the tier that can read the template source but cannot run it. That split
is load-bearing here, because the defect this file was opened for was invisible on
both sides of it on its own:

* ``input[type='checkbox'] { display: none; }`` was **one unscoped rule** that took
  all four of the report's interactive inputs out of the focus order and out of the
  accessibility tree at once, while Bootstrap kept drawing the switch from
  ``::before``/``::after`` on the sibling ``<label>`` - so the page looked and
  clicked exactly right and only a keyboard or a screen reader found it;
* nothing in the markup *said* any of that. Focusability is computed, so the
  measurement lives in ``tests/browser/test_keyboard_operability.py``.

What is asserted here is the other half: the constructs that make the computed
result possible. A rule that removes an interactive element from the page, a
disclosure built out of a checkbox rather than ``<details>``, a missing viewport,
a missing landmark, an animation with no ``prefers-reduced-motion`` escape, and a
status glyph with no name. Each of those is a property of the source text, and
each one is cheap to reintroduce by accident.

Later tasks in this plan add contrast, print-block and token assertions to this
file; it is the presentation tier's home, not this change's scratchpad.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import vntyper
from vntyper.scripts import report_formatting as rf

pytestmark = pytest.mark.unit

TEMPLATE_DIR = Path(vntyper.__file__).resolve().parent / "templates"
PER_SAMPLE_TEMPLATE = TEMPLATE_DIR / "report_template.html"
COHORT_TEMPLATE = TEMPLATE_DIR / "cohort_summary_template.html"
TEMPLATES = (PER_SAMPLE_TEMPLATE, COHORT_TEMPLATE)

#: Element names that are focusable and operable in their own right. A stylesheet
#: rule that hides one of these hides a control, whatever the author intended.
INTERACTIVE_ELEMENTS = ("input", "button", "summary", "select", "textarea", "details")

#: The unscoped rule this change removed, in the shape the brief names it.
_HIDDEN_CHECKBOX = re.compile(r"input\[type=.checkbox.\][^{]*\{[^}]*display:\s*none")

#: A declaration that takes its subject off the page entirely.
_HIDES = re.compile(r"display\s*:\s*none|visibility\s*:\s*hidden")

#: Layout properties whose animation the design detector reports as
#: ``layout-transition``, and which the ``<details>`` rewrite exists to remove.
_LAYOUT_TRANSITION = re.compile(r"transition[^;{}]*:[^;{}]*\b(?:max-height|height|width|padding|margin)\b")


#: HTML, Jinja2 and C-style block comments. Everything below asks what the
#: template *does*, and both templates document at length what they used to do -
#: including, verbatim, the rule and the handler these tests ban. Scanning the
#: prose would make an accurate comment fail the build and reward deleting it.
_COMMENT = re.compile(r"<!--.*?-->|\{#.*?#\}|/\*.*?\*/", re.DOTALL)


def _markup(template: Path) -> str:
    """Read a template with its comments removed.

    Args:
        template: The template file.

    Returns:
        str: The source text with every HTML, Jinja2 and block comment stripped.
    """
    return _COMMENT.sub("", template.read_text(encoding="utf-8"))


def _style_blocks(source: str) -> list[str]:
    """Return the contents of every ``<style>`` element.

    Args:
        source: A template's source text.

    Returns:
        list[str]: One string per ``<style>`` element.
    """
    return re.findall(r"<style[^>]*>(.*?)</style>", source, re.DOTALL)


def _rules(source: str) -> list[tuple[str, str]]:
    """Return every ``(selector, declarations)`` pair in the template's stylesheets.

    Nested at-rules are flattened by the same crude split, which is fine for what
    is asked of it: the questions below are all "does any rule anywhere say X".

    Args:
        source: A template's source text.

    Returns:
        list[tuple[str, str]]: Selector text and declaration text, both stripped.
    """
    return [
        (selector.strip(), body.strip())
        for block in _style_blocks(source)
        for selector, body in re.findall(r"([^{}]+)\{([^{}]*)\}", block)
    ]


def _declarations(style: str) -> list[str]:
    """Split an inline style string into its declarations.

    Args:
        style: An inline ``style`` attribute value.

    Returns:
        list[str]: Whitespace-stripped, non-empty declarations.
    """
    return [part.strip() for part in style.split(";") if part.strip()]


# ---------------------------------------------------------------------------
# Nothing removes an interactive element from the accessibility tree
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_interactive_input_is_removed_from_the_accessibility_tree(template: Path) -> None:
    """The exact rule the defect was: unscoped, and four controls wide.

    ``display: none`` is not a styling choice for a control. It removes the element
    from the focus order *and* from the accessibility tree, so the switch could be
    operated only by clicking the label a mouse user happened to find.
    """
    assert not _HIDDEN_CHECKBOX.search(_markup(template)), (
        f"{template.name} hides every checkbox on the page, which takes each one out of the focus order"
    )


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_stylesheet_rule_hides_an_interactive_element(template: Path) -> None:
    """The general shape, so the fix cannot come back under a different selector.

    Rules whose subject is a **pseudo-element** are exempt: hiding
    ``summary::-webkit-details-marker`` replaces the disclosure triangle, it does
    not hide the disclosure. That is why the exemption is written as "the selector
    contains ``::``" rather than as a list of blessed selectors.
    """
    offenders = [
        (selector, body)
        for selector, body in _rules(_markup(template))
        if "::" not in selector
        and _HIDES.search(body)
        and any(re.search(rf"(?<![\w.#-]){element}(?![\w-])", selector) for element in INTERACTIVE_ELEMENTS)
    ]
    assert offenders == [], f"{template.name} hides an interactive element with CSS: {offenders}"


# ---------------------------------------------------------------------------
# Disclosures are native
# ---------------------------------------------------------------------------


def test_every_disclosure_is_a_details_element() -> None:
    """A checkbox plus a label is not a disclosure; it only looks like one.

    ``<details>`` is focusable, is announced as a disclosure with its state, is
    opened by find-in-page and is opened by printing. The checkbox accordion had
    none of those, and its ``max-height`` transition was the one finding the design
    detector reported on this template.
    """
    source = _markup(PER_SAMPLE_TEMPLATE)

    assert source.count("<details") == source.count("<summary")
    assert source.count("<details") >= 2, "both the variant table and the pipeline log are disclosures"
    assert 'class="toggle"' not in source
    assert "lbl-toggle" not in source
    assert "collapsible-content" not in source


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_rule_animates_a_layout_property(template: Path) -> None:
    """``transition: max-height`` was how the accordion opened, and it is gone.

    Animating a layout property makes the browser re-lay-out on every frame. It is
    also the shape that survives a rewrite unnoticed, because the result still
    looks fine on a fast machine.
    """
    offenders = [(selector, body) for selector, body in _rules(_markup(template)) if _LAYOUT_TRANSITION.search(body)]
    assert offenders == [], f"{template.name} animates a layout property: {offenders}"


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_motion_is_switched_off_for_a_reader_who_asked_for_that(template: Path) -> None:
    """Every remaining transition needs an escape, and the escape is a media query.

    A reader with a vestibular disorder sets this at the operating system. The
    report animated its disclosures unconditionally.
    """
    source = _markup(template)
    if "transition" not in source and "animation" not in source:
        pytest.skip(f"{template.name} declares no motion at all")
    assert "prefers-reduced-motion" in source, f"{template.name} animates without honouring prefers-reduced-motion"


# ---------------------------------------------------------------------------
# Viewport, landmarks and announcements
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_the_template_declares_a_viewport(template: Path) -> None:
    """Without this a phone lays the page out at 980px and scales it down.

    Measured: a 390px device rendered the 12.8px body text at an effective ~5.1px.
    Both templates were missing it, which is why both are checked here.
    """
    source = _markup(template)
    assert re.search(r'<meta\s+name="viewport"\s+content="[^"]*width=device-width', source), (
        f"{template.name} declares no viewport, so a phone lays it out at 980px"
    )


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_the_report_body_is_a_main_landmark(template: Path) -> None:
    """A screen reader user skips to the content by landmark; there was none."""
    source = _markup(template)
    assert "<main" in source and "</main>" in source, f"{template.name} has no main landmark"


def test_the_header_warning_is_announced() -> None:
    """The BAM-header warning is the report's one genuinely urgent message.

    It was a red paragraph, which is colour plus a word. ``role="alert"`` is what
    puts it in front of a reader who is not looking at it.
    """
    source = _markup(PER_SAMPLE_TEMPLATE)
    block = re.search(r'<p[^>]*role="alert"[^>]*>.*?</p>', source, re.DOTALL)
    assert block, "the header warning does not announce itself"
    assert "header_warning" in block.group(0), "role=alert is on something other than the header warning"


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_a_keyboard_reader_can_see_where_the_focus_is(template: Path) -> None:
    """Focus that cannot be seen is focus that cannot be used.

    Neither template styled ``:focus`` or ``:focus-visible`` at all, so the visible
    ring was whatever the engine's default happened to be - and on the switch there
    was no box to draw it around.
    """
    rules = _rules(_markup(template))
    focus_rules = [(selector, body) for selector, body in rules if ":focus-visible" in selector]
    assert focus_rules, f"{template.name} styles no focus ring"
    assert any("outline" in body for _, body in focus_rules), (
        f"{template.name} has a :focus-visible rule that draws no outline: {focus_rules}"
    )


# ---------------------------------------------------------------------------
# Every meaning carried by colour or a glyph is also carried by text
# ---------------------------------------------------------------------------


def test_both_status_glyphs_carry_an_accessible_name() -> None:
    """The coverage and fastp columns are a bare ⚠ or ✔ with a colour.

    A glyph with no name is announced as its code point or skipped entirely, so the
    Status column of the fastp table read as empty. The two names must also differ,
    or the column says the same thing whatever the measurement was.
    """
    names = []
    for icon in (rf.WARNING_ICON, rf.OK_ICON):
        assert 'role="img"' in icon, f"{icon} is not exposed as an image, so its aria-label is ignored"
        match = re.search(r'aria-label="([^"]+)"', icon)
        assert match, f"{icon} has no accessible name"
        names.append(match.group(1))

    assert names[0] != names[1], f"both status glyphs are announced as {names[0]!r}"


def test_the_confidence_styling_does_not_rest_on_colour_alone() -> None:
    """Two confidence values distinguished only by hue are one value in greyscale.

    ``High_Precision`` was red and ``Low_Precision`` orange, both bold - so the
    only difference between them, beyond the label itself, was the hue. Every pair
    the colours separate has to be separated by something else as well.
    """
    styles = rf.CONFIDENCE_STYLES
    colours = {next(d for d in _declarations(style) if d.startswith("color:")) for style in styles.values()}
    other = {frozenset(d for d in _declarations(style) if not d.startswith("color:")) for style in styles.values()}

    assert len(other) >= len(colours), (
        f"the confidence styles distinguish {len(colours)} colours but only {len(other)} non-colour treatments, "
        f"so some of the distinction exists in hue alone: {styles}"
    )


@pytest.mark.parametrize("value", sorted(rf.CONFIDENCE_STYLES))
def test_the_confidence_label_is_always_present_as_text(value: str) -> None:
    """Whatever the styling does, the value itself is words in the cell."""
    assert f">{value}<" in rf.confidence_html(value)


# ---------------------------------------------------------------------------
# The IGV variant rows are real controls
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Print is a build target - issue #242
# ---------------------------------------------------------------------------


def _at_rule_body(source: str, prelude: str) -> str:
    """Return the body of one at-rule, brace-matched so nested at-rules survive.

    ``@media print`` contains ``@page``-style nesting and the crude
    ``([^{}]+)\\{([^{}]*)\\}`` split used elsewhere in this file cannot see it. This
    walks the braces instead.

    Args:
        source: A template's source text.
        prelude: The at-rule's prelude, e.g. ``"@media print"``.

    Returns:
        str: Everything between the rule's braces, or "" when it is absent.
    """
    start = source.find(prelude)
    if start == -1:
        return ""
    opening = source.find("{", start)
    depth = 0
    for index in range(opening, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[opening + 1 : index]
    return ""


def _print_block() -> str:
    """The per-sample report's ``@media print`` body, comments stripped."""
    return _at_rule_body(_markup(PER_SAMPLE_TEMPLATE), "@media print")


def test_the_report_has_a_print_stylesheet() -> None:
    """The archived PDF is the artefact that outlives the HTML.

    Measured in a browser rather than grepped, before this change:
    ``document.styleSheets`` reported ``mediaRules: []`` for both templates - not one
    author media rule, so the printed page inherited every screen affordance.
    """
    assert "@media print" in PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")
    assert _print_block(), "the @media print block is empty"


def test_print_unsets_the_truncation_clamp() -> None:
    """``.table td`` clamps to 150px with ``overflow: hidden`` and reveals on hover.

    The motif sequence is 121 bp in real data and paper has no hover, so the archived
    record carried a truncated sequence with an ellipsis where the evidence was.
    """
    block = _print_block()

    assert "max-width: none" in block
    assert "white-space: normal" in block
    assert "overflow: visible" in block, "the cell is still clipped by overflow: hidden"


def test_print_gives_an_unbroken_sequence_somewhere_to_break() -> None:
    """Unsetting the clamp is not enough on its own.

    A 121-character motif sequence contains no space, so with the clamp gone it still
    overflows its cell and is clipped by the table's own width. Measured in Chromium:
    ``max-width: none`` plus ``white-space: normal`` printed 100 of the 121 bases;
    adding a break opportunity printed all 121.
    """
    block = _print_block()

    assert "overflow-wrap: anywhere" in block or "word-break: break-all" in block, (
        "print unsets the clamp but gives a long unbroken sequence no break opportunity, so it is clipped anyway"
    )


def test_print_hides_the_controls_that_do_nothing_on_paper() -> None:
    """The two switches print as ticked or unticked boxes that nobody can operate."""
    assert ".controls" in _print_block()


def test_the_printed_page_is_numbered() -> None:
    """A dropped sheet is unrecoverable without this, and a truncated print is invisible.

    ``@page`` margin boxes, which the supported renderer (Chromium) honours;
    ``tests/browser/test_printed_record.py`` is what proves they reach the PDF.
    """
    page_rule = _at_rule_body(_markup(PER_SAMPLE_TEMPLATE), "@page")

    assert "counter(page)" in page_rule
    assert "counter(pages)" in page_rule, "the page number does not say how many pages there are"


def test_the_printed_record_states_its_identity_in_the_document() -> None:
    """Who, what, which assembly, which version - at the head of the printed sheet.

    It is a block in the document rather than an ``@page`` margin box because the
    supported engine cannot put a document value in one: measured, Chromium 151 drops
    a ``content`` list containing ``string()`` whole. See the note in the template and
    ``test_no_value_is_interpolated_into_a_stylesheet`` for the alternative that was
    refused.
    """
    source = _markup(PER_SAMPLE_TEMPLATE)
    block = re.search(r'<div class="print-header">(.*?)</div>', source, re.DOTALL)

    assert block, "the printed record has no header line"
    for value in ("sample_name", "assay_name", "assembly_declared", "pipeline_version", "research_use_statement"):
        assert value in block.group(1), f"the printed header line does not state {value}"


def test_no_value_is_interpolated_into_a_stylesheet() -> None:
    """The printed identity is a document value, and it stays out of the CSS.

    A ``<style>`` element is a raw text element and CSS is not an HTML context:
    autoescaping turns ``&`` into ``&amp;`` there, which nothing decodes, and does
    nothing about the characters that matter in a CSS string - a sample name reaching
    ``content:`` could close the element outright. The sample name is derived from a
    sample-supplied basename and the report is a file people forward, so a running
    header in the page margin is worth less than this guarantee.
    """
    for style in _style_blocks(PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")):
        assert "{{" not in style, f"a value is interpolated into a stylesheet: {style[:200]}"


def test_the_log_prints_as_a_pointer_rather_than_a_wall_of_debug() -> None:
    """A DEBUG log is pages of it, and it is in the HTML original either way.

    This is the one disclosure that does *not* print open: it is the exception the
    print block states explicitly, so the rule that opens the others cannot be read as
    an accident.
    """
    block = _print_block()

    assert "log-section" in block, "the print block does not treat the log differently from any other disclosure"
    assert "HTML original" in block, "the printed log does not point anywhere"
    assert 'class="log-section"' in _markup(PER_SAMPLE_TEMPLATE) or "log-section" in _markup(PER_SAMPLE_TEMPLATE)


def test_every_other_disclosure_prints_open() -> None:
    """The no-JS half of the belt and braces.

    CSS alone cannot open a closed ``<details>`` in Chromium - measured: with this rule
    and scripting off, the body of a closed disclosure did not print - so the handler
    below is what does the work there. The rule is still what an engine honouring the
    author's cascade over the UA's closed-details shadow content needs, and it costs
    nothing.
    """
    block = _print_block()

    assert "details > summary ~ *" in block
    assert "display: block !important" in block


def test_a_section_the_reader_collapsed_is_restored_after_printing() -> None:
    """Printing must not be a mutation. Both halves of the handler, and its own block.

    A ``$`` that never resolved throws at the top of whichever block it is in and takes
    every statement after it down, so the print handler shares its block with no jQuery
    - the same reason the switches were moved out of it.
    """
    source = PER_SAMPLE_TEMPLATE.read_text(encoding="utf-8")
    blocks = re.findall(r"<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>", source, re.DOTALL)

    owning = [block for block in blocks if "beforeprint" in block]
    assert len(owning) == 1, f"expected exactly one script block to handle printing, found {len(owning)}"
    assert "afterprint" in owning[0], "the handler opens sections for printing and never closes them again"
    assert "log-section" in owning[0], "the handler would open the log as well, which prints as a pointer"
    assert "$(" not in owning[0], "the print handler shares its block with jQuery code that can throw first"


def test_the_printed_record_carries_the_whole_coverage_table() -> None:
    """The detailed view is hidden behind a switch, and paper has no switches.

    ``#detailedCoverageView`` carries eight statistics and the QC verdict that the
    basic view does not, and it is ``display: none`` until a reader ticks a control the
    print block removes. Printing the basic view alone archives less than the screen
    could show.
    """
    block = _print_block()

    assert "#detailedCoverageView" in block
    assert "#basicCoverageView" in block, "both coverage views print, so the mean is stated twice"


def test_the_report_states_that_it_is_research_use_only() -> None:
    """The one thing a forwarded artefact has to say about itself.

    Quoted from ``README.md`` and declared as a constant beside the assay name rather
    than written into the template, for the same reason every other fixed string in
    ``report_identity`` is: interpretive text is not invented at a call site.
    """
    from vntyper.scripts.report_identity import RESEARCH_USE_STATEMENT

    assert "research use only" in RESEARCH_USE_STATEMENT.lower()
    assert "research_use_statement" in _markup(PER_SAMPLE_TEMPLATE)


def test_the_igv_variant_rows_are_operated_by_a_real_control() -> None:
    """``<tr onclick>`` has no keyboard path, no name and no state.

    A ``<button>`` brings Enter and Space, a focus ring and a name for free. The
    row carries ``aria-selected`` so which variant is showing is a fact in the
    accessibility tree rather than a background colour.
    """
    source = _markup(PER_SAMPLE_TEMPLATE)

    assert "row.onclick" not in source, "the variant row is still activated by a click handler on the row"
    assert 'createElement("button")' in source, "no per-row control is created"
    assert 'setAttribute("aria-selected"' in source, "the selected variant row is not announced as selected"
    assert 'setAttribute("role", "grid")' in source, (
        "aria-selected is only meaningful on a row inside a grid, so the table has to declare one"
    )
