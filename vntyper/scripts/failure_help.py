"""Operator-facing help for pipeline failures: what went wrong, and where to read more (#338).

A failed run used to end in a Python traceback. That is right for a bug report and wrong
as the last thing an operator reads, especially for setup errors that are the same on every
sample. This module holds the pure decisions: the one-line summary a failed run ends with,
and the fix hints for the Kestrel errors whose meaning is known.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

#: The published troubleshooting page. Every new operator-facing failure message links it.
TROUBLESHOOTING_URL = "https://hassansaei.github.io/VNtyper/user-guide/troubleshooting/"

#: Known Kestrel ERROR texts and what to do about each, checked in order. Each marker is
#: a substring of the message the pinned Kestrel 1.0.1 logs, measured against the jar.
KESTREL_ERROR_HINTS: tuple[tuple[str, str], ...] = (
    (
        "Error reading reference sequence",
        "The MUC1 motif reference is missing, empty or corrupt. Reinstall it from the repository "
        "root with `vntyper install-references -d reference` and run VNtyper from that directory.",
    ),
    (
        "indexed k-mer count (IKC)",
        "Kestrel could not read the k-mer counts KAnalyze wrote. Check free disk space and write "
        "permission in the output directory, then re-run.",
    ),
    (
        "Error setting k-mer counts",
        "Kestrel could not load the k-mer counts KAnalyze wrote. Check free disk space and write "
        "permission in the output directory, then re-run.",
    ),
    (
        "File not found while variant writer",
        "Kestrel could not write its VCF. Check that the output directory still exists and is writable.",
    ),
)


def kestrel_error_hint(errors: list[str]) -> str | None:
    """Return the fix for the first Kestrel error whose meaning is known.

    Args:
        errors: Kestrel ERROR lines, in log order.

    Returns:
        str | None: The hint for the earliest line that matches a known marker, or None
        when none does.
    """
    for line in errors:
        for marker, hint in KESTREL_ERROR_HINTS:
            if marker in line:
                return hint
    return None


def summarize_failure(exc: BaseException) -> str:
    """Render the last line a failed run prints: the cause, and where to get help.

    Args:
        exc: The exception that ended the run.

    Returns:
        str: ``VNtyper failed: <first line of the message>.`` plus the troubleshooting link.
        The first line is used because several messages continue with quoted detail that
        is already in the log above.
    """
    lines = str(exc).strip().splitlines()
    # A first line that introduces detail ends in ":"; normalise it to one sentence.
    headline = (lines[0] if lines else type(exc).__name__).rstrip(" :.")
    return (
        f"VNtyper failed: {headline}. "
        f"The full message and traceback are above and in pipeline.log. Help: {TROUBLESHOOTING_URL}"
    )
