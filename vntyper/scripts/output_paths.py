"""
Keep a caller-named report file inside the output directory it was meant for.

``vntyper report --report-file`` and ``vntyper cohort --summary-file`` are documented as
*names* -- "Name of the output report file" -- but both are plain strings joined onto
``--output-dir`` with :func:`os.path.join`, which honours a parent-relative or absolute
value and quietly writes wherever it points. ``--report-file ../../summary.html``
overwrites a sibling run's report; ``--report-file /etc/cron.d/vntyper`` ignores
``--output-dir`` altogether.

This module states the rule those two options were always documented as having: a single
path segment, resolving inside the output directory. It is pure -- a name, a directory
and an answer -- so it can be tested without writing anything, and it lives outside
``generate_report.py`` and ``cohort_summary.py`` because both are over the size limit and
both need the same check (AGENTS.md rule 3).

**Rejected, never repaired.** ``os.path.basename("../../report.html")`` is
``"report.html"``, so silently taking the basename would write a file the caller did not
name and say nothing about it. A caller who passed a path meant something by it, and the
useful answer is the one that says so.
"""

import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


def contained_output_path(output_dir: str | Path, name: str, option: str) -> Path:
    """
    Resolve a caller-supplied report file name inside the output directory.

    Args:
        output_dir (str | Path): The directory the run was told to write into.
        name (str): The value of ``--report-file`` / ``--summary-file``. Untrusted.
        option (str): The option the value came from, named in any error so the message
            says which of the two was wrong.

    Returns:
        Path: ``output_dir / name``, once it is known to be a single name that lands
            inside ``output_dir``.

    Raises:
        ValueError: If the value is empty, is not a single path segment, or resolves
            outside the output directory.
    """
    candidate = (name or "").strip()
    if not candidate or candidate in {".", ".."} or os.path.basename(candidate) != candidate:
        msg = (
            f"{option} must be a single file name inside the output directory, not a path: {name!r}. "
            "Use --output-dir to choose where the report is written."
        )
        logger.error(msg)
        raise ValueError(msg)

    directory = Path(output_dir)
    destination = directory / candidate
    # Resolved, because a single segment can still be a symlink pointing elsewhere and
    # writing through it would overwrite the target rather than a file in this directory.
    if os.path.dirname(os.path.realpath(destination)) != os.path.realpath(directory):
        msg = f"{option} does not resolve inside the output directory: {name!r}"
        logger.error(msg)
        raise ValueError(msg)

    return destination
