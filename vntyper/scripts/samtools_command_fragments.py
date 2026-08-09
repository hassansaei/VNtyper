"""Shared shell fragments for samtools command construction."""

from __future__ import annotations

import shlex
from pathlib import Path


def quote_path(value: str | Path | int) -> str:
    """Quote a user-supplied command value.

    Args:
        value: Path, region, count, or other non-tool command operand.

    Returns:
        The shell-quoted value.
    """
    return shlex.quote(str(value))


def thread_flag(threads: int) -> str:
    """Return the optional samtools thread flag.

    Args:
        threads: Requested thread count.

    Returns:
        A trailing-space-terminated ``-@`` fragment, or an empty string.
    """
    if threads <= 1:
        return ""
    return f"-@ {quote_path(threads)} "


def reference_flag(reference_path: str | Path | None) -> str:
    """Return the optional samtools reference flag.

    Args:
        reference_path: Reference FASTA path, if one was resolved.

    Returns:
        A trailing-space-terminated ``-T`` fragment, or an empty string.
    """
    if reference_path is None:
        return ""
    return f"-T {quote_path(reference_path)} "


def customized_index_input(alignment_path: str | Path, index_path: str | Path | None) -> str:
    """Return alignment operands with an optional explicit custom index.

    Args:
        alignment_path: BAM or CRAM input.
        index_path: Exact BAI or CRAI selected with samtools ``-X``.

    Returns:
        The quoted alignment alone, or ``-X alignment index``.
    """
    if index_path is None:
        return quote_path(alignment_path)
    return f"-X {quote_path(alignment_path)} {quote_path(index_path)}"
