"""Read the Kestrel call log for the failures its exit status does not report (#338).

The vendored Kestrel 1.0.1 exits ``0`` after a fatal error. Measured against the pinned
``kestrel.jar``:

| Condition                     | Exit | VCF          | Log                                              |
| ----------------------------- | ---- | ------------ | ------------------------------------------------ |
| motif reference missing/empty | 0    | none         | ``ERROR ... Error reading reference sequence(s)`` |
| IKC missing                   | 0    | header-only  | ``ERROR ... Cannot open indexed k-mer count``     |
| IKC empty (zero k-mers)       | 0    | header-only  | ``ERROR ... Error setting k-mer counts``          |

The first row reached users as "Kestrel produced no usable VCF" on every sample, with the
cause visible only in ``kestrel_kmer_<k>.log``. The last two are worse: a header-only VCF
is, correctly, a usable empty result to ``describe_unusable_vcf``, so the sample would be
reported as a confident negative. The exit status cannot tell these apart from success;
the log can, so the log is the contract.

KAnalyze does not need this: it exits ``2`` on the same kinds of input failure.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from pathlib import Path

from vntyper.scripts.failure_help import TROUBLESHOOTING_URL, kestrel_error_hint

logger = logging.getLogger(__name__)

#: One Kestrel logback line at ERROR level: ``HH:MM:SS [thread] ERROR logger - message``.
#: Anchored on the timestamp and the thread bracket so that a message which merely
#: *contains* the word ERROR (a sample named ``ERROR``, a path) cannot match.
ERROR_LINE = re.compile(r"^\d{2}:\d{2}:\d{2}(?:\.\d+)? \[[^\]]*\]\s+ERROR\s")

#: How many ERROR lines a raised message quotes. Kestrel repeats one cause through
#: several layers (``CountMap`` then ``KestrelRunner``), so the first few carry it.
MAX_QUOTED_ERRORS = 3


def find_kestrel_errors(lines: Iterable[str]) -> list[str]:
    """Return the ERROR-level lines of a Kestrel log, in order.

    Args:
        lines: The log's lines, with or without trailing newlines.

    Returns:
        list[str]: Each matching line with its trailing whitespace removed; empty when
        the log records no error.
    """
    return [line.rstrip() for line in lines if ERROR_LINE.match(line)]


def read_kestrel_errors(log_file: str | Path) -> list[str]:
    """Return the ERROR-level lines of the Kestrel log at ``log_file``.

    A missing log yields no errors. ``run_command`` opens the log before it starts the
    process, so a real run always leaves one; its absence means nothing ran that could
    have logged. Any other read failure raises, because a log that exists and cannot be
    read cannot vouch for a run that exited ``0``.

    Args:
        log_file: The Kestrel call log.

    Returns:
        list[str]: The ERROR lines, possibly empty.

    Raises:
        RuntimeError: If the log exists but cannot be read.
    """
    try:
        with open(log_file, encoding="utf-8", errors="replace") as handle:
            return find_kestrel_errors(handle)
    except FileNotFoundError:
        return []
    except OSError as exc:
        msg = f"Could not read the Kestrel log {log_file} to check it for errors ({exc})."
        logger.error(msg)
        raise RuntimeError(msg) from exc


def describe_kestrel_errors(errors: list[str], *, kmer_size: int, log_file: str | Path) -> str:
    """Render the message raised when Kestrel logged errors but exited ``0``.

    Args:
        errors: The ERROR lines from :func:`find_kestrel_errors`; must be non-empty.
        kmer_size: The attempt's k-mer size.
        log_file: The Kestrel call log, named so the operator can read the rest.

    Returns:
        str: One message quoting up to :data:`MAX_QUOTED_ERRORS` lines.

    Raises:
        ValueError: If ``errors`` is empty; there is nothing to describe.
    """
    if not errors:
        raise ValueError("describe_kestrel_errors needs at least one error line.")
    quoted = "\n  ".join(errors[:MAX_QUOTED_ERRORS])
    more = len(errors) - MAX_QUOTED_ERRORS
    tail = f"\n  ... and {more} more" if more > 0 else ""
    hint = kestrel_error_hint(errors)
    fix = f"Fix: {hint}\n" if hint else ""
    return (
        f"Kestrel reported {len(errors)} error(s) for k-mer size {kmer_size} but exited 0, "
        f"so this attempt is treated as failed rather than as a result:\n  {quoted}{tail}\n"
        f"{fix}Full log: {log_file}\nHelp: {TROUBLESHOOTING_URL}"
    )
