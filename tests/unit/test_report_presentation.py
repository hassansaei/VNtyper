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

#: The token layer both reports include. It is a partial rather than a base both
#: templates ``extend``: the two reports share a stylesheet and nothing else - no block
#: structure, no shared body - and ``{% include %}`` says exactly that. Both renderers
#: build a ``FileSystemLoader`` over the configured template directory
#: (``generate_report.py`` and ``cohort_summary.py``), so a partial beside the two
#: templates is reachable from both, including from a third-party directory configured
#: through ``config.json``'s ``paths.template_dir``.
SHARED_PARTIAL = TEMPLATE_DIR / "_report_base.html"

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


#: ``{% include "name.html" %}``, with or without whitespace control. Both report
#: templates pull the token layer in this way, so "what does this template say" stopped
#: being a question about one file: a rule that used to sit in ``report_template.html``
#: now lives in the partial and applies to both reports. Everything below therefore reads
#: a template's *effective* source rather than its own bytes.
_INCLUDE = re.compile(r"\{%-?\s*include\s+['\"]([^'\"]+)['\"].*?-?%\}", re.DOTALL)


def _source(template: Path, _seen: tuple[str, ...] = ()) -> str:
    """Read a template with every ``{% include %}`` resolved, recursively.

    Args:
        template: The template file.
        _seen: Templates already being expanded on this branch, so a cycle terminates
            instead of recursing until the interpreter gives up.

    Returns:
        str: The source text with each included partial spliced in where it was
        included.
    """
    if template.name in _seen:
        return ""
    seen = (*_seen, template.name)
    return _INCLUDE.sub(
        lambda match: _source(TEMPLATE_DIR / match.group(1), seen),
        template.read_text(encoding="utf-8"),
    )


def _markup(template: Path) -> str:
    """Read a template with its includes resolved and its comments removed.

    Args:
        template: The template file.

    Returns:
        str: The effective source text with every HTML, Jinja2 and block comment
        stripped.
    """
    return _COMMENT.sub("", _source(template))


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


# ---------------------------------------------------------------------------
# Contrast, computed from the stylesheet the report ships
#
# The point of this section is that nothing in it transcribes a colour. Every
# value is parsed out of `_report_base.html`, so a palette edit is measured rather
# than believed - a table of hardcoded pairs would go on passing while the token
# it claims to describe drifted underneath it.
#
# What it replaces, all measured with the same function below:
#
#   1.72:1  white on `lightskyblue`   - every table heading in the per-sample report
#   1.97:1  `orange` on white         - `Low_Precision`, the Kestrel confidence column
#   1.76:1  `orange` on a striped row - the same value in the cohort report
#   4.00:1  `red` on white            - `High_Precision`, and the fastp warning glyph
#   3.57:1  `red` on a striped row    - the same value in the cohort report
#   3.80:1  `red` on `#f9f9f9`        - the BAM-header warning, in `.header-info`
#
# Every one of those carried meaning a reader has to be able to read. The only
# colour pair in the whole report that met AA was the green tick that means
# nothing is wrong.
# ---------------------------------------------------------------------------

#: WCAG 2.2 SC 1.4.3 for normal-size text. Nothing here is large-scale text: the
#: report's headings are the only type over 18pt and they carry no colour of their own.
WCAG_AA_NORMAL_TEXT = 4.5

#: Tokens that colour text, and the surfaces text is drawn on. The product of the two
#: is what is asserted, because a token layer does not say which pairs occur - any
#: foreground may end up on any of the three surfaces, and `--sunken` is a table head
#: in one report and a disclosure summary in the other.
TEXT_TOKENS = (
    "--ink",
    "--ink-muted",
    "--accent",
    "--state-finding",
    "--state-caution",
    "--state-ok",
    "--state-none",
)
SURFACE_TOKENS = ("--paper", "--surface", "--sunken")

#: The CSS colour keywords this repository's stylesheets have actually used, so a
#: literal one can be recognised where it appears. It is deliberately not the full CSS
#: named-colour list: the check below is a tripwire over the shapes that shipped here,
#: not a CSS parser, and it says so rather than implying completeness.
_NAMED_COLOURS = {
    "white": "#ffffff",
    "black": "#000000",
    "red": "#ff0000",
    "green": "#008000",
    "orange": "#ffa500",
    "lightskyblue": "#87cefa",
    "lightblue": "#add8e6",
    "darkblue": "#00008b",
    "gray": "#808080",
    "grey": "#808080",
}

#: A custom property declaration.
_CUSTOM_PROPERTY = re.compile(r"(--[\w-]+)\s*:\s*([^;{}]+)")

#: A ``:root`` rule and its body.
_ROOT_RULE = re.compile(r":root\s*\{([^{}]*)\}")

#: A reference to a custom property.
_VAR_REFERENCE = re.compile(r"var\(\s*--[\w-]+\s*(?:,[^()]*)?\)")

#: A literal colour: hex, a functional notation, or one of the keywords above.
_LITERAL_COLOUR = re.compile(r"#[0-9a-fA-F]{3,8}\b|\b(?:rgb|rgba|hsl|hsla)\(|\b(?:" + "|".join(_NAMED_COLOURS) + r")\b")

#: A declaration that paints text or a surface. ``border-color`` and friends are
#: excluded by the lookbehind: a hairline is not text and is held to 3:1, not 4.5:1.
_PAINT_DECLARATION = re.compile(r"(?<![\w-])(color|background|background-color)\s*:\s*([^;{}]+)")


def _channels(colour: str) -> tuple[int, int, int]:
    """Return the 8-bit sRGB channels of one colour.

    Args:
        colour: A hex colour (``#rgb`` or ``#rrggbb``) or one of :data:`_NAMED_COLOURS`.

    Returns:
        tuple[int, int, int]: The red, green and blue channels.
    """
    text = colour.strip().lower()
    text = _NAMED_COLOURS.get(text, text).lstrip("#")
    if len(text) == 3:
        text = "".join(channel * 2 for channel in text)
    return (int(text[0:2], 16), int(text[2:4], 16), int(text[4:6], 16))


def _relative_luminance(colour: str) -> float:
    """Return the WCAG relative luminance of one colour.

    Args:
        colour: Any colour :func:`_channels` understands.

    Returns:
        float: The relative luminance, 0 for black and 1 for white.
    """

    def linear(value: int) -> float:
        srgb = value / 255
        return srgb / 12.92 if srgb <= 0.04045 else ((srgb + 0.055) / 1.055) ** 2.4

    red, green, blue = (linear(channel) for channel in _channels(colour))
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def contrast_ratio(foreground: str, background: str) -> float:
    """Return the WCAG contrast ratio between two colours.

    Args:
        foreground: The text colour.
        background: The colour behind it.

    Returns:
        float: The ratio, between 1 and 21. The order of the arguments does not matter.
    """
    lighter, darker = sorted(
        (_relative_luminance(foreground), _relative_luminance(background)),
        reverse=True,
    )
    return (lighter + 0.05) / (darker + 0.05)


def _root_tokens(css: str) -> dict[str, str]:
    """Return the custom properties every ``:root`` rule in ``css`` declares.

    Args:
        css: A stylesheet, or the body of one at-rule.

    Returns:
        dict[str, str]: Property name -> value, later declarations winning.
    """
    return {name: value.strip() for body in _ROOT_RULE.findall(css) for name, value in _CUSTOM_PROPERTY.findall(body)}


def _palettes() -> dict[str, dict[str, str]]:
    """Return every palette the shared token layer declares, keyed by its name.

    The dark and print palettes redefine a subset of the light one, so each is the
    light palette updated with its own overrides - which is also the only way to notice
    a redefinition that *removes* a token the light palette had.

    Returns:
        dict[str, dict[str, str]]: ``light``, ``dark`` and ``print`` -> token -> colour.
    """
    source = SHARED_PARTIAL.read_text(encoding="utf-8") if SHARED_PARTIAL.is_file() else ""
    dark = _at_rule_body(source, "@media (prefers-color-scheme: dark)")
    printed = _at_rule_body(source, "@media print")

    base = source
    for body in (dark, printed):
        if body:
            base = base.replace(body, "")

    light = _root_tokens(base)
    return {
        "light": light,
        "dark": {**light, **_root_tokens(dark)},
        "print": {**light, **_root_tokens(printed)},
    }


def _paint_declarations() -> list[tuple[str, str, str]]:
    """Every colour a rendered report paints text or a surface with.

    Two sources, because the report's colour comes from two places: the stylesheets in
    the templates (the shared partial included), and the inline ``style`` attributes
    ``report_formatting`` and ``cohort_tables`` build into table cells.

    ``@page`` bodies are excluded. Their margin boxes sit outside the document tree and
    whether a ``:root`` custom property reaches them was not measured in the supported
    engine, so the page number keeps a literal grey rather than a token that might
    resolve to nothing.

    Returns:
        list[tuple[str, str, str]]: ``(where, property, value)`` for each declaration.
    """
    found: list[tuple[str, str, str]] = []

    for template in (*TEMPLATES, SHARED_PARTIAL):
        if not template.is_file():
            continue
        for block in _style_blocks(_COMMENT.sub("", template.read_text(encoding="utf-8"))):
            page_rule = _at_rule_body(block, "@page")
            scanned = block.replace(page_rule, "") if page_rule else block
            found += [(template.name, name, value) for name, value in _PAINT_DECLARATION.findall(scanned)]

    for where, markup in _rendered_inline_styles():
        for style in re.findall(r'style="([^"]*)"', markup):
            found += [(where, name, value) for name, value in _PAINT_DECLARATION.findall(style)]

    return found


#: The page background both reports had before the token layer existed: the user-agent
#: default, which Bootstrap restates as ``background-color: #fff``. A literal that has
#: escaped the token layer is measured against it, because that is what a reader saw.
_ASSUMED_PAGE_BACKGROUND = "#ffffff"


def _literal_colours(value: str) -> list[str]:
    """Return the literal colours in one declaration value.

    ``var()`` references are removed first, so ``color-mix(in oklab, var(--x) 8%,
    transparent)`` names no colour of its own and neither does ``inherit``.

    Args:
        value: The declaration's value.

    Returns:
        list[str]: The literal colours, empty when the value resolves through tokens.
    """
    return _LITERAL_COLOUR.findall(_VAR_REFERENCE.sub("", value))


def _describe_literal(where: str, name: str, value: str, surface: str) -> str:
    """Describe one escaped literal, with its measured contrast where that is meaningful.

    Args:
        where: The file or constant the declaration came from.
        name: The property name.
        value: The declaration's value.
        surface: The colour to measure a foreground literal against.

    Returns:
        str: A one-line description for the failure message.
    """
    literal = _literal_colours(value)[0]
    declaration = f"{where}: {name}: {value.strip()}"
    if name != "color" or literal.endswith("("):
        return declaration
    return f"{declaration} ({contrast_ratio(literal, surface):.2f}:1 on {surface})"


def _rendered_inline_styles() -> list[tuple[str, str]]:
    """The fragments this codebase builds with an inline ``style`` attribute.

    Returns:
        list[tuple[str, str]]: ``(where, markup)`` for each fragment.
    """
    from vntyper.scripts.cohort_tables import confidence_span

    return [
        ("report_formatting.WARNING_ICON", rf.WARNING_ICON),
        ("report_formatting.OK_ICON", rf.OK_ICON),
        ("report_formatting.flag_html('Low_Depth')", rf.flag_html("Low_Depth")),
        ("report_formatting.flag_html('Not flagged')", rf.flag_html("Not flagged")),
        ("report_formatting.confidence_html('High_Precision')", rf.confidence_html("High_Precision")),
        ("report_formatting.confidence_html('Low_Precision')", rf.confidence_html("Low_Precision")),
        ("cohort_tables.confidence_span('High_Precision')", str(confidence_span("High_Precision"))),
        ("cohort_tables.confidence_span('Low_Precision')", str(confidence_span("Low_Precision"))),
    ]


def test_the_contrast_function_agrees_with_the_measurements_this_change_records() -> None:
    """Calibration for everything below, against the pairs the token layer replaces.

    Without this, a broken ratio function would make the assertions in this section
    vacuous in the direction that matters - too *generous* - and nothing would say so.
    The six values are the ones recorded in the commit message and in the section
    comment above, so the two cannot drift apart either.
    """
    measured = {
        "white on lightskyblue": contrast_ratio("white", "lightskyblue"),
        "orange on white": contrast_ratio("orange", "#ffffff"),
        "orange on a striped row": contrast_ratio("orange", "#f2f2f2"),
        "red on white": contrast_ratio("red", "#ffffff"),
        "red on a striped row": contrast_ratio("red", "#f2f2f2"),
        "red on the header-info panel": contrast_ratio("red", "#f9f9f9"),
    }

    assert {name: round(ratio, 2) for name, ratio in measured.items()} == {
        "white on lightskyblue": 1.72,
        "orange on white": 1.97,
        "orange on a striped row": 1.76,
        "red on white": 4.00,
        "red on a striped row": 3.57,
        "red on the header-info panel": 3.80,
    }
    assert contrast_ratio("#ffffff", "#000000") == pytest.approx(21.0)
    assert contrast_ratio("#777777", "#777777") == pytest.approx(1.0)


def test_every_token_pair_meets_wcag_aa() -> None:
    """Every foreground token on every surface token, in all three palettes.

    Parsed from the stylesheet rather than transcribed, so editing a token is what this
    measures. Asserting the *product* rather than the pairs that happen to occur today
    is deliberate: which surface a state colour lands on is a markup decision, and a
    later task moves that markup.
    """
    palettes = _palettes()
    assert palettes["light"], f"{SHARED_PARTIAL.name} declares no colour tokens at all"

    missing = sorted(
        f"{name}: {token}"
        for name, palette in palettes.items()
        for token in (*TEXT_TOKENS, *SURFACE_TOKENS)
        if token not in palette
    )
    assert missing == [], f"the token layer does not declare: {missing}"

    failures = [
        f"{name}: {foreground} ({palette[foreground]}) on {background} ({palette[background]}) "
        f"is {contrast_ratio(palette[foreground], palette[background]):.2f}:1"
        for name, palette in palettes.items()
        for foreground in TEXT_TOKENS
        for background in SURFACE_TOKENS
        if contrast_ratio(palette[foreground], palette[background]) < WCAG_AA_NORMAL_TEXT
    ]
    assert failures == [], f"token pairs below {WCAG_AA_NORMAL_TEXT}:1 for normal text: {failures}"


def test_no_report_colour_is_declared_outside_the_token_layer() -> None:
    """The other half: a token layer nothing uses proves nothing.

    Every ``color`` and ``background`` the two reports paint with has to resolve through
    ``:root``, or the test above measures a palette the page does not use. This is the
    check that fails when a literal comes back - which is how ``lightskyblue`` and
    ``orange`` survived every previous accessibility pass.
    """
    surface = _palettes()["light"].get("--surface", _ASSUMED_PAGE_BACKGROUND)
    offenders = [
        _describe_literal(where, name, value, surface)
        for where, name, value in _paint_declarations()
        if _literal_colours(value)
    ]
    assert offenders == [], f"colours declared outside the token layer: {offenders}"


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


@pytest.mark.parametrize(
    "renderer",
    ["report_formatting.confidence_html", "cohort_tables.confidence_span"],
)
@pytest.mark.parametrize("value", sorted(rf.CONFIDENCE_CLASSES))
def test_no_confidence_value_is_distinguished_by_colour(renderer: str, value: str) -> None:
    """The distinction between two confidence values is the word, in both reports.

    ``High_Precision`` was ``color:red`` and ``Low_Precision`` ``color:orange``, both
    bold, so beyond the label the hue was the whole of the difference - and red on the
    most trustworthy call sat one column from a red ``Flag`` glyph meaning the opposite.
    A transitional underline was added so the two were not separated by hue alone; it
    goes with the hue, because the label was always the honest channel and the label is
    text in the cell either way.

    Both renderers are asserted. They used to be near-copies of one another and the
    cohort one kept its own colour table, which is how the same value came to measure
    1.76:1 there and 1.97:1 in the per-sample report.
    """
    from vntyper.scripts.cohort_tables import confidence_span

    render = rf.confidence_html if renderer.startswith("report_formatting") else confidence_span
    markup = str(render(value))

    assert "style=" not in markup, f"{renderer} still carries an inline style: {markup}"
    assert "color" not in markup.replace(rf.CONFIDENCE_CLASSES[value], ""), f"{renderer} still names a colour: {markup}"
    assert "text-decoration" not in markup, f"{renderer} still distinguishes by underline: {markup}"
    assert rf.CONFIDENCE_CLASSES[value] in markup, f"{renderer} drops the class the stylesheet hangs on: {markup}"


@pytest.mark.parametrize("value", sorted(rf.CONFIDENCE_CLASSES))
def test_the_confidence_label_is_always_present_as_text(value: str) -> None:
    """Whatever the styling does, the value itself is words in the cell."""
    assert f">{value}<" in rf.confidence_html(value)


def test_the_shared_stylesheet_gives_the_confidence_class_no_hue() -> None:
    """The other end of the same fact: the class exists and it paints nothing.

    A class the stylesheet colours would put the hue straight back, one indirection
    further away from the test above.
    """
    coloured = [
        (selector, body)
        for selector, body in _rules(_COMMENT.sub("", SHARED_PARTIAL.read_text(encoding="utf-8")))
        if "confidence" in selector and re.search(r"(?<![\w-])color\s*:", body)
    ]
    assert coloured == [], f"the confidence class is colour-coded after all: {coloured}"


# ---------------------------------------------------------------------------
# The IGV variant rows are real controls
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Print is a build target - issue #242
# ---------------------------------------------------------------------------


def _print_block() -> str:
    """The per-sample report's ``@media print`` body, comments stripped."""
    return _at_rule_body(_markup(PER_SAMPLE_TEMPLATE), "@media print")


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_the_report_has_a_print_stylesheet(template: Path) -> None:
    """The archived PDF is the artefact that outlives the HTML.

    Measured in a browser rather than grepped, before this change:
    ``document.styleSheets`` reported ``mediaRules: []`` for both templates - not one
    author media rule, so the printed page inherited every screen affordance.

    Both templates are checked because there is now one print block and both include
    it. It used to exist in the per-sample report only, so a printed cohort report
    inherited the screen's clipped cells and had no page numbers.
    """
    assert "@media print" in _source(template)
    assert _at_rule_body(_markup(template), "@media print"), f"{template.name}'s @media print block is empty"


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_each_template_reaches_the_token_layer_by_include(template: Path) -> None:
    """The composition mechanism, pinned where a later change would notice.

    ``{% include %}``, not ``{% extends %}``: the two reports share a stylesheet and
    have no common body. A resolver that walks the shipped templates - to check that
    every context key is referenced, say - has to follow includes recursively to see
    the partial at all, and the partial has to be in the same loader directory for
    either renderer to find it.
    """
    assert SHARED_PARTIAL.is_file(), f"{SHARED_PARTIAL.name} is missing"
    assert f'include "{SHARED_PARTIAL.name}"' in template.read_text(encoding="utf-8"), (
        f"{template.name} does not include the shared token layer"
    )
    assert "{% extends" not in template.read_text(encoding="utf-8"), (
        f"{template.name} extends a base; the resolver in a later task is written for includes"
    )


def test_the_rules_both_reports_need_are_declared_once() -> None:
    """The duplication this partial removes, asserted so it cannot creep back.

    The truncation block was byte-identical in the two templates and the switch was a
    near-copy; ``.table-container`` was used by the per-sample report and defined only
    in the cohort one, so that wrapper styled nothing at all. Each of these now has
    exactly one declaration, in the file both reports include.
    """
    # Comments are stripped: both templates document at length what they used to do,
    # and scanning the prose would make an accurate comment fail the build.
    files = {t.name: _COMMENT.sub("", t.read_text(encoding="utf-8")) for t in (*TEMPLATES, SHARED_PARTIAL)}

    declared_once = (
        "--state-finding:",
        ".switch input",
        ":focus-visible",
        "@media print",
        "prefers-reduced-motion",
        'meta name="viewport"',
    )
    for name in declared_once:
        owners = [where for where, text in files.items() if name in text]
        assert owners == [SHARED_PARTIAL.name], f"{name!r} is declared in {owners}, not only in the token layer"

    # The clamp itself is gone rather than moved. `max-width: 150px` with
    # `text-overflow: ellipsis` put the 121 bp motif sequence behind a `:hover`, and
    # the exception that let the `Flag` reason escape it was `td:last-child` - which
    # moving `Motif_sequence` last would have handed to the motif and taken off the
    # flag. Cells wrap instead, so no column needs an exception.
    for gone in ("text-overflow", "max-width: 150px"):
        owners = [where for where, text in files.items() if gone in text]
        assert owners == [], f"{gone!r} came back in {owners}"

    # `.table-container` is the one hook both reports already used and only one
    # declared. The cohort keeps its own width and centring on top; what it no longer
    # owns is whether the wrapper does anything at all.
    assert ".table-container {" in SHARED_PARTIAL.read_text(encoding="utf-8")
    for template in TEMPLATES:
        assert 'class="table-container"' in template.read_text(encoding="utf-8"), (
            f"{template.name} no longer uses the wrapper the token layer styles"
        )


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


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_value_is_interpolated_into_a_stylesheet(template: Path) -> None:
    """The printed identity is a document value, and it stays out of the CSS.

    A ``<style>`` element is a raw text element and CSS is not an HTML context:
    autoescaping turns ``&`` into ``&amp;`` there, which nothing decodes, and does
    nothing about the characters that matter in a CSS string - a sample name reaching
    ``content:`` could close the element outright. The sample name is derived from a
    sample-supplied basename and the report is a file people forward, so a running
    header in the page margin is worth less than this guarantee.

    Read through the include, so the shared token layer is held to it too - it is the
    file that now carries the ``@page`` margin boxes the temptation belongs to.
    """
    for style in _style_blocks(_source(template)):
        assert "{{" not in style, f"a value is interpolated into {template.name}'s stylesheet: {style[:200]}"


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
