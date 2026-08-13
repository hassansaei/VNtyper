"""The escaper that lets a sample name into a stylesheet at all.

Interpolating a value into CSS was refused outright for as long as there was no
escaper, and the refusal was correct: a ``<style>`` is a raw text element, Jinja2's
autoescaping does nothing useful inside one, and the sample name is derived from a
sample-supplied basename in a file people forward. Chromium 151 drops a ``@page``
margin box whose ``content`` contains ``string()``, so a running header carrying the
sample name on every printed sheet has no other route.

:func:`~vntyper.scripts.css_escaping.css_string_literal` is therefore the whole of
the defence, and these cases are written against the attacks rather than against the
implementation: the parsers that read this text are the CSS tokenizer *and* the HTML
tokenizer, and only one of them cares about quotes.

The properties asserted, each one an escape route closed:

* the literal cannot be closed early (``"``), and the escape mechanism cannot be
  turned back on itself (``\\``);
* no raw control character survives, so a newline cannot terminate the string and a
  NUL cannot confuse the parse;
* no ``<`` survives, so ``</style>`` cannot end the raw text element - the one
  attack no amount of CSS-correct quoting stops;
* nothing outside printable ASCII survives, so the ``<style>`` element's encoding
  cannot decide what the header says;
* and the escapes are terminated in a way that neither eats the next character nor
  merges with it, which is where a *correct-looking* escaper goes wrong quietly.
"""

from __future__ import annotations

import re

import pytest

from vntyper.scripts.css_escaping import css_string_literal

pytestmark = pytest.mark.unit

#: Every character that must never appear in the produced literal's body at all. A
#: single membership test catches an escaper that handles the case a test was written
#: for and not the class it belongs to.
FORBIDDEN_IN_BODY = ('"', "<", ">", "&")

#: Printable ASCII minus ``"``, ``&``, ``<``, ``>`` and ``\\`` - written out as ranges
#: rather than as a negated class, so that "what may appear literally" is stated
#: positively and a character added to the escaper's denylist has to be removed here
#: too.
_LITERAL_CHARACTER = r" !#-%'-;=?-\[\]-~"

#: The grammar the body must match: those characters, and backslashes that introduce
#: exactly six hex digits and a space. Escaping ``\\`` and ``"`` as hex rather than as
#: ``\\\\`` and ``\\"`` is what makes this one pattern the whole contract instead of a
#: list of cases.
BODY_GRAMMAR = re.compile(rf"^(?:[{_LITERAL_CHARACTER}]|\\[0-9a-f]{{6}} )*$")

#: The attacks. Each one is a sample name somebody could actually be handed - a file
#: name is not a controlled vocabulary - and each one exits a different parser.
ADVERSARIAL_NAMES = (
    'S1" ; }',  # closes the literal, ends the declaration, ends the rule
    "S1\\",  # a trailing backslash escapes whatever follows it
    "S1\nS2",  # a raw newline is a CSS string terminator, with an error
    "S1\x00S2",  # NUL
    "S1\x7fS2",  # DEL, a control character on the printable side of the table
    "</style><script>alert(1)</script>",  # leaves CSS entirely, via the HTML tokenizer
    "</StYlE>",  # the tag match is case-insensitive; the escaper must not rely on case
    "patiënt-042",  # ordinary non-ASCII, which is not an attack and must still be safe
    "‮sample",  # a bidi override, which rewrites what a reader sees
)


@pytest.mark.parametrize("name", ADVERSARIAL_NAMES)
def test_no_adversarial_name_escapes_the_literal(name: str) -> None:
    """The whole contract, over every attack at once.

    Args:
        name: A hostile or merely awkward sample name.
    """
    literal = css_string_literal(name)
    body = literal[1:-1]

    assert literal.startswith('"') and literal.endswith('"')
    for character in FORBIDDEN_IN_BODY:
        assert character not in body, f"{character!r} survived into the literal body: {literal!r}"
    assert all(0x20 <= ord(character) <= 0x7E for character in body), f"a raw non-printable survived: {literal!r}"
    assert BODY_GRAMMAR.match(body), f"the body is not made only of literal characters and complete escapes: {body!r}"


def test_a_backslash_cannot_be_left_to_escape_the_closing_quote() -> None:
    """``S1\\`` immediately before the closing quote is the classic way out.

    If the backslash were passed through, the literal would end ``\\"`` and the quote
    would be escaped rather than closing - so the string would run on into the rest
    of the stylesheet, and everything after it in the file would be inside a string.
    """
    assert css_string_literal("S1\\") == '"S1\\00005c "'


def test_a_quote_cannot_close_the_literal() -> None:
    assert css_string_literal('a"b') == '"a\\000022 b"'


@pytest.mark.parametrize(
    "character,expected",
    [
        ("\n", "\\00000a "),
        ("\x00", "\\000000 "),
        ("\x7f", "\\00007f "),
        ("<", "\\00003c "),
        (">", "\\00003e "),
        ("&", "\\000026 "),
        ("é", "\\0000e9 "),
        ("·", "\\0000b7 "),
        ("😀", "\\01f600 "),
    ],
)
def test_each_escaped_class_is_written_as_six_hex_digits_and_a_space(character: str, expected: str) -> None:
    """Six digits is the maximum a CSS escape may have, so the parser stops there.

    A shorter escape runs into a following hex digit and decodes as something else
    entirely: ``\\b7`` followed by the letter ``a`` is U+B7A, not ``·a``.

    Args:
        character: One character of the class under test.
        expected: The escape it must produce.
    """
    assert css_string_literal(character) == f'"{expected}"'


def test_an_escape_does_not_eat_the_space_after_it() -> None:
    """The bug a *correct-looking* escaper ships with.

    A CSS escape consumes one whitespace character as its terminator. The escaper
    emits that terminator itself, so a space in the value survives as a space; an
    escaper that relied on the value's own space would silently render
    ``sample · assay`` as ``sample ·assay``.
    """
    literal = css_string_literal("A · B")

    assert literal == '"A \\0000b7  B"'
    assert literal.count("  ") == 1, "the terminator and the value's own space must both be there"


def test_a_hex_escape_cannot_merge_with_the_character_after_it() -> None:
    """``é`` followed by a hex digit, which is the ambiguity six digits removes."""
    assert css_string_literal("éa") == '"\\0000e9 a"'


@pytest.mark.parametrize("name", ADVERSARIAL_NAMES)
def test_the_literal_never_terminates_the_style_element(name: str) -> None:
    """The HTML tokenizer's rule, checked the way the HTML tokenizer applies it.

    A raw text element ends at ``</`` followed by its tag name, case-insensitively,
    and CSS quoting is not consulted. This asserts on a whole ``<style>`` element
    built the way the template builds one.

    Args:
        name: A hostile or merely awkward sample name.
    """
    element = f"<style>@page {{ @top-left {{ content: {css_string_literal(name)}; }} }}</style>"

    assert len(re.findall(r"</\s*style", element, re.IGNORECASE)) == 1, f"the element is closed twice: {element!r}"


def test_ordinary_text_is_left_alone() -> None:
    """An escaper that mangles a normal name is one nobody will keep."""
    assert css_string_literal("PATIENT_042 (hg19)") == '"PATIENT_042 (hg19)"'


def test_an_empty_value_is_still_a_literal() -> None:
    """A ``content`` of nothing is valid; a ``content`` of *no value* is a parse error
    that takes the whole declaration with it."""
    assert css_string_literal("") == '""'
