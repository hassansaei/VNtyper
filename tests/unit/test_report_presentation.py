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
