"""
css_escaping.py

Module Purpose:
---------------
One function, and it is a security boundary: turn an arbitrary string into a CSS
string literal that cannot escape the literal, the declaration, the rule, or the
``<style>`` element the rule is written into.

Why this exists at all. The printed report needs a running header in the page
margin, and Chromium 151 drops a ``@page`` margin box whose ``content`` list
contains ``string()`` -- so the only way to put the sample name on every sheet is to
write it into the stylesheet. That was refused for as long as there was no escaper,
and the refusal was right: a ``<style>`` is a **raw text element**, so Jinja2's
autoescaping does nothing useful there (it turns ``&`` into ``&amp;``, which CSS
then renders literally, and leaves every character that actually matters alone),
and the sample name is derived from a sample-supplied basename in a file people
forward. This module is that escaper. It is what makes the interpolation
defensible, and it is the whole of the defence, so it is tested against the attacks
before anything calls it.

What is escaped, and why each one:

* ``\\``, which starts a CSS escape sequence, and ``"``, which closes the literal;
* control characters -- a raw newline terminates a CSS string with an error, and a
  NUL is a parse hazard in its own right;
* everything non-ASCII, because the ``<style>`` element's encoding is not
  guaranteed to be the document's and a mis-decoded byte is a character nobody
  chose;
* ``<``, ``>`` and ``&``, which are not CSS-special at all. They are escaped
  because CSS is not the only parser reading this text: the HTML tokenizer ends a
  raw text element at ``</style``, and no amount of CSS-correct quoting stops it.
  A sample named ``</style><script>...`` is the whole reason.

**All of them take the same ``\\XXXXXX `` hex form**, including the two that could
have been written ``\\\\`` and ``\\"``. One form rather than two is what makes the
result checkable in a sentence -- *the body is printable ASCII, contains no
``"``/``<``/``>``/``&``, and every backslash in it introduces exactly six hex digits
and a space* -- and a property a test can assert as a grammar is worth more here
than two characters of brevity, because this function is the entire defence.

Functions:
    css_string_literal: One arbitrary string to one quoted CSS string literal
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

#: The quote the produced literal is delimited with. Single quotes are equally legal
#: CSS; one is chosen so that exactly one quote character needs escaping inside.
QUOTE = '"'

#: Characters escaped despite being ordinary printable ASCII. ``\\`` and ``"`` are
#: CSS-special. ``<``, ``>`` and ``&`` are not special to CSS at all: they are here
#: because the ``<style>`` element is raw text and the HTML tokenizer stops at
#: ``</style`` whatever CSS thinks of the surrounding quotes. ``>`` and ``&`` follow
#: ``<`` rather than being separately necessary - the cost is nothing and "no markup
#: character survives" is a rule a reader can check.
ALWAYS_ESCAPED = '\\"<>&'

#: The highest code point written verbatim. Above this is escaped; ``\x7f`` (DEL) is
#: below it only by accident of the table and is handled as a control character.
LAST_PRINTABLE_ASCII = 0x7E

#: The lowest code point written verbatim.
FIRST_PRINTABLE_ASCII = 0x20


def css_string_literal(value: str) -> str:
    """Render a string as a CSS string literal, delimiters included.

    Every escape is written as six hex digits followed by a space. Six digits is the
    maximum a CSS escape may have, so the parser stops there of its own accord, and
    the trailing space is then consumed as the escape's terminator -- which is what
    makes a space in ``value`` *immediately after* an escaped character survive
    instead of being eaten. Both halves matter: without the padding, ``\\b7`` runs
    into a following hex digit and decodes as a different character; without the
    space, ``\\0000b7`` swallows the space after it.

    The result is safe to write into a ``<style>`` element as-is: it contains no
    ``<``, ``>``, ``&``, no ``"`` and no raw control character between its
    delimiters, and nothing outside printable ASCII.

    Args:
        value: The text to render. Any string, including one a sample named.

    Returns:
        str: The quoted literal, e.g. ``'"S1"'``. Never an empty string: an empty
        input renders as ``'""'``.

    Example:
        >>> css_string_literal('a"b')
        '"a\\\\000022 b"'
    """
    pieces: list[str] = [QUOTE]
    for character in value:
        code_point = ord(character)
        if character in ALWAYS_ESCAPED or code_point < FIRST_PRINTABLE_ASCII or code_point > LAST_PRINTABLE_ASCII:
            pieces.append(f"\\{code_point:06x} ")
            continue
        pieces.append(character)
    pieces.append(QUOTE)
    return "".join(pieces)
