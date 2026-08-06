"""
Working-tree guards for the advisory mutation harness.

``mutation_test.py`` rewrites production source **in place**: it reads a target file,
writes one mutant after another over it, and finally writes back the text it read at
the start. That contract is only safe on a clean tree, because the text it read at the
start is whatever happened to be on disk - not what is committed. Sweeping a tree with
an uncommitted edit in a target file therefore goes wrong three ways at once:

* the published score and the committed evidence page describe **uncommitted** code
  while presenting it as the committed baseline;
* the end-of-sweep "is everything back?" check sees the maintainer's own edit and
  reports it as a failed restoration; and
* whatever that error suggests doing about it is aimed at code the harness never
  touched.

This module holds the two things that follow from that: asking git which targets are
dirty, and wording the two failures. Both are separated from the harness because they
are the parts a test can drive without a sweep - the parsing and the messages are pure,
and the one subprocess call has a single well-defined command.

.. warning::

   Nothing here may ever recommend ``git checkout --``, ``git restore``,
   ``git reset --hard`` or ``git clean``. The situation in which a user reads these
   messages is precisely the one where uncommitted work is at stake, and a recovery
   instruction that discards it is worse than no instruction at all.
   ``tests/unit/test_mutation_test.py`` asserts the absence of every one of them.
"""

from __future__ import annotations

import subprocess
from collections.abc import Iterable, Sequence
from pathlib import Path

#: ``git status --porcelain`` (v1) prefixes every entry with a two-character status
#: code and one space - ``" M path"``, ``"?? path"``, ``"R  old -> new"``. The format is
#: explicitly documented as stable for scripts, which is why it is parsed rather than
#: the human-readable output.
_PORCELAIN_PREFIX = 3


def parse_porcelain(output: str) -> list[str]:
    """
    Extract the paths named by ``git status --porcelain`` output.

    Args:
        output (str): Raw stdout of ``git status --porcelain``.

    Returns:
        list[str]: One repo-relative path per entry, in the order reported. Renames and
            copies contribute their destination - the path that is actually on disk.
    """
    paths: list[str] = []
    for line in output.splitlines():
        entry = line[_PORCELAIN_PREFIX:].strip()
        if not entry:
            continue
        # `R  old -> new` and `C  old -> new`: the file on disk is the right-hand side.
        if " -> " in entry:
            entry = entry.split(" -> ", 1)[1]
        paths.append(entry.strip('"'))
    return paths


def dirty_paths(repo_root: Path, paths: Iterable[str]) -> list[str]:
    """
    Report which of the given paths have uncommitted changes.

    Args:
        repo_root (Path): Repository root; the command runs here.
        paths (Iterable[str]): Repo-relative paths to ask about. An empty iterable asks
            about the whole tree, so callers must not pass one by accident.

    Returns:
        list[str]: The subset of ``paths`` git reports as modified, staged or untracked.
            Empty when the tree is clean, and also empty when git is unavailable or this
            is not a repository - the guard exists to protect a maintainer's checkout,
            so it declines to block a working sweep over a missing tool.
    """
    try:
        completed = subprocess.run(
            ["git", "status", "--porcelain", "--", *paths],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return []
    if completed.returncode != 0:
        return []
    return parse_porcelain(completed.stdout)


def format_dirty_tree_refusal(dirty: Sequence[str]) -> str:
    """
    Word the refusal to start a sweep over uncommitted work.

    Args:
        dirty (Sequence[str]): Target paths with uncommitted changes.

    Returns:
        str: A message naming every dirty file and saying what to do, without naming a
            command that would throw the changes away.
    """
    listed = "\n".join(f"  {path}" for path in dirty)
    return (
        "ERROR: refusing to start - these target files have uncommitted changes:\n"
        f"{listed}\n"
        "\n"
        "The sweep rewrites each target in place and puts back the text it read at the\n"
        "start, so running it now would measure and publish uncommitted code as the\n"
        "committed baseline, and would then report your own edit as a failed restore.\n"
        "Commit or stash these changes and re-run. Nothing has been modified."
    )


def format_unrestored_warning(dirty: Sequence[str]) -> str:
    """
    Word the end-of-sweep failure where a target was left differing from the tree.

    The tree was verified clean before the sweep, so a difference here is a mutant the
    harness failed to put back. The message still stops short of prescribing a command:
    it cannot prove the diff is only the mutant, and the cost of being wrong is somebody
    else's work.

    Args:
        dirty (Sequence[str]): Target paths that still differ.

    Returns:
        str: A message naming every file and how to inspect it.
    """
    listed = "\n".join(f"  {path}" for path in dirty)
    return (
        "ERROR: the sweep may not have restored these files:\n"
        f"{listed}\n"
        "\n"
        "They were clean when the sweep started, so the difference is most likely a\n"
        "mutant left on disk. Review it with `git diff -- <path>`, and do not commit,\n"
        "build, install or package from this tree until you have."
    )
