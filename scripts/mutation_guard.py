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

The generated files are in the same position. ``write_outputs()`` rewrites
``docs/development/mutation-testing.md`` and the raw results wholesale, on every run
including ``--render-only``, so an uncommitted edit to either is simply overwritten. They
are guarded alongside the sources; see :func:`writable_paths`.

This module holds what follows from that: listing what a run would write, asking git
which of those are dirty, and wording the three failures. All of it is separated from the
harness because it is the part a test can drive without a sweep - the parsing and the
messages are pure, and the one subprocess call has a single well-defined command.

The guard **fails closed**. When git cannot answer, :func:`dirty_paths` raises rather than
assuming the tree is clean: unknown is not clean, and the asymmetry is total - refusing
wrongly costs a message, proceeding wrongly costs uncommitted work.

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


def writable_paths(repo_root: Path, targets: Iterable[str], outputs: Iterable[Path | None]) -> list[str]:
    """
    List every repo-relative path a run would write, sources and generated files alike.

    The sweep rewrites its mutation targets, and ``write_outputs()`` overwrites the
    report and the raw results unconditionally - including on ``--render-only``, which
    rewrites the docs page without mutating anything. Guarding only the sources therefore
    left an uncommitted edit to ``docs/development/mutation-testing.md`` to be destroyed
    by the next ``make mutation``.

    Args:
        repo_root (Path): Repository root. Outputs are expressed relative to it.
        targets (Iterable[str]): Repo-relative source paths the sweep will mutate.
        outputs (Iterable[Path | None]): Output destinations, ``None`` for the ones this
            run was not asked to write.

    Returns:
        list[str]: Sorted, de-duplicated repo-relative paths. Outputs that fall outside
            the repository are omitted - git has nothing to say about them, and a
            ``--output /tmp/report.md`` is a legitimate way to run the harness.
    """
    paths = set(targets)
    root = repo_root.resolve()
    for output in outputs:
        if output is None:
            continue
        try:
            paths.add(str(Path(output).resolve().relative_to(root)))
        except ValueError:
            continue
    return sorted(paths)


def format_indeterminate_refusal(reason: str) -> str:
    """
    Word the refusal to proceed when cleanliness could not be determined at all.

    Args:
        reason (str): What git said, or why it could not be run.

    Returns:
        str: A message giving the reason and saying that nothing was touched.
    """
    return (
        "ERROR: refusing to start - cannot determine whether the working tree is clean:\n"
        f"  {reason}\n"
        "\n"
        "The sweep rewrites production source in place and puts back the text it read at\n"
        "the start, and it overwrites the generated report, so it is only safe on a clean\n"
        "tree. Without an answer from git that precondition cannot be established, and\n"
        "assuming the tree is clean is the assumption whose cost is somebody's uncommitted\n"
        "work. Run inside a git working tree, with git on PATH.\n"
        "Nothing has been modified."
    )


def dirty_paths(repo_root: Path, paths: Iterable[str]) -> list[str]:
    """
    Report which of the given paths have uncommitted changes.

    Args:
        repo_root (Path): Repository root; the command runs here.
        paths (Iterable[str]): Repo-relative paths to ask about. An empty iterable asks
            about the whole tree, so callers must not pass one by accident.

    Returns:
        list[str]: The subset of ``paths`` git reports as modified, staged or untracked.
            Empty when, and only when, git reports the tree clean.

    Raises:
        RuntimeError: When git cannot answer - it is missing, not executable, or exits
            non-zero. **Unknown is not clean.** An earlier version returned ``[]`` here
            on the grounds that the guard should not block a sweep over a missing tool;
            that made the guard fail open, which is the pattern it was written to remove.
            A corrupted index, a ``safe.directory`` refusal and a genuine non-repository
            are indistinguishable from in here, and only one of them is harmless.
    """
    try:
        completed = subprocess.run(
            ["git", "status", "--porcelain", "--", *paths],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise RuntimeError(format_indeterminate_refusal(f"git could not be run: {exc}")) from exc
    if completed.returncode != 0:
        detail = completed.stderr.strip() or f"git status exited {completed.returncode}"
        raise RuntimeError(format_indeterminate_refusal(detail))
    return parse_porcelain(completed.stdout)


def format_dirty_tree_refusal(dirty: Sequence[str]) -> str:
    """
    Word the refusal to start a sweep over uncommitted work.

    Args:
        dirty (Sequence[str]): Paths the run would write - mutation targets or generated
            outputs - that have uncommitted changes.

    Returns:
        str: A message naming every dirty file and saying what to do, without naming a
            command that would throw the changes away.
    """
    listed = "\n".join(f"  {path}" for path in dirty)
    return (
        "ERROR: refusing to start - these files the run would overwrite have uncommitted\n"
        "changes:\n"
        f"{listed}\n"
        "\n"
        "The sweep rewrites each target in place and puts back the text it read at the\n"
        "start, so running it now would measure and publish uncommitted code as the\n"
        "committed baseline, and would then report your own edit as a failed restore.\n"
        "The generated report and results are rewritten wholesale, so an uncommitted edit\n"
        "to either is simply lost.\n"
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
