"""The checks that decide whether a gate run is admissible as evidence.

Comparing two sides is only meaningful once three separate things are true, and none of
them follows from the comparison itself:

1. **Each side did what its case declared it would do.** Every case in the matrix carries
   an ``expect_exit``, and a case expected to exit zero also carries the artefacts it must
   have produced. Without this, two runs that both die at exit 1 before writing a single
   genotype artefact compare *equal*: :func:`golden_cohort.compare._presence` classifies
   ``(None, None)`` as ``absent_both``, :func:`golden_cohort.compare.diff_case` excludes
   ``absent_both`` from the delta list, and the verdict comes back ``IDENTICAL``. The gate
   then exits 0 having measured nothing at all. That is the worst failure a gate has, and
   for the first four recorded runs nothing in the harness prevented it - ``expect_exit``
   was written seven times in :mod:`golden_cohort.matrix` and read nowhere.

2. **The two sides really are two different, oppositely-labelled trees.** The launcher
   proves "this process ran the tree the caller named", which is sound and is the hard
   part. It does not prove the caller named two *different* trees. Passing the same run
   root twice, or two run roots whose sides both recorded ``before``, or two sides that
   both expected the marker present, all produce a confident ``IDENTICAL`` over a tree
   compared against itself - the exact failure the whole wrapper exists to prevent,
   reached by misconfiguration rather than by import resolution.

3. **The revision each side ran is recorded, not asserted.** See :mod:`golden_cohort`'s
   entry point: before this module the harness never ran ``git rev-parse``, never looked
   at the working tree's cleanliness and never accepted an expected SHA, so ``side.json``
   recorded a *path*. A path is not a revision: the same path is a different commit ten
   minutes later. Every "attests commit X and nothing after it" sentence on the gate page
   was therefore an operator's assertion that the instrument could neither check nor
   contradict.

The exit and artefact rules are pure functions of the case and what is on disk, so they
are unit-testable without a cohort. The revision recording shells out to ``git`` and takes
its runner as an argument for the same reason.
"""

from __future__ import annotations

import logging
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

#: The only exit codes ``vntyper`` itself produces. AGENTS.md: "Only exit codes 0 and 1
#: are used." Anything else out of a launched case is a harness-level event - 97 is the
#: launcher's refusal to start, 98 an unhandled exception inside the CLI, 99 the runner's
#: timeout - and never satisfies a case's expectation, whichever way that expectation
#: points. A probe that is *supposed* to fail must fail the way the product fails.
CLI_EXIT_CODES: frozenset[int] = frozenset({0, 1})

#: Paths whose dirtiness can change a genotype, and therefore the paths whose working-tree
#: state is recorded per side. Mirrors the gate page's own
#: ``git diff --stat <candidate>..HEAD -- vntyper/ docker/ scripts/``.
REVISION_PATHS: tuple[str, ...] = ("vntyper", "docker", "scripts")

#: Artefacts every zero-exit ``vntyper pipeline`` case must have written, relative to its
#: ``--output-dir``. Measured against run 4 rather than assumed: all 59 zero-expected
#: pipeline cases produced all five on both sides, and the two nonzero probes produced
#: none of them. ``advntr/output_adVNTR_result.tsv`` is deliberately **not** here - adVNTR
#: is an optional module and its absence is legitimate on the 56 cases that do not run it,
#: so requiring it would make the rule case-dependent for no gain. The consequence is
#: named: an adVNTR case that silently produces no adVNTR output is not caught by this
#: rule, only by the row-set comparison finding the artefact missing on one side.
PIPELINE_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "pipeline_summary.json",
    "kestrel/kestrel_result.tsv",
    "kestrel/kestrel_pre_result.tsv",
    "coverage/coverage_summary.tsv",
    "summary_report.html",
)

#: Artefacts every ``vntyper cohort`` case that aggregates at least one sample must have
#: written. ``cohort_empty`` aggregates none, writes only its log and still exits 0 by
#: design, so it declares an empty requirement instead of this one.
COHORT_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cohort_summary.html",
    "cohort_kestrel.csv",
    "cohort_kestrel.tsv",
    "cohort_kestrel.json",
    "cohort_advntr.csv",
    "cohort_advntr.tsv",
    "cohort_advntr.json",
)


def check_exit(expect_exit: str, exit_code: int | None, *, aborted: bool, timed_out: bool) -> str:
    """Judge one case's exit code against what its matrix entry declared.

    Args:
        expect_exit: ``"zero"`` or ``"nonzero"``, from the case.
        exit_code: What the launched process returned, or None if it never ran.
        aborted: Whether the launcher refused to start this run.
        timed_out: Whether the runner killed it.

    Returns:
        str: An empty string when the expectation is met, otherwise a one-line reason.

    Raises:
        ValueError: If ``expect_exit`` is neither ``"zero"`` nor ``"nonzero"``. A typo
            there would otherwise silently exempt the case from the whole check.
    """
    if expect_exit not in ("zero", "nonzero"):
        msg = f"Unknown expect_exit {expect_exit!r}; the matrix declares only 'zero' or 'nonzero'."
        logger.error(msg)
        raise ValueError(msg)
    if aborted:
        return "the launcher refused to start this run, so its exit code is not a pipeline result"
    if timed_out:
        return "the run was killed on the harness timeout, so its exit code is not a pipeline result"
    if exit_code is None:
        return "no exit code was recorded, so this case never ran"
    if exit_code not in CLI_EXIT_CODES:
        return f"exit {exit_code} is not one of vntyper's own exit codes {sorted(CLI_EXIT_CODES)}"
    if expect_exit == "zero" and exit_code != 0:
        return f"expected exit 0, got {exit_code}"
    if expect_exit == "nonzero" and exit_code == 0:
        return "expected a nonzero exit, got 0"
    return ""


def missing_artifacts(output_dir: Path, required: Sequence[str]) -> list[str]:
    """List the declared artefacts this case did not write.

    Args:
        output_dir: The case's ``--output-dir``.
        required: Paths relative to it, from the case's ``required_artifacts``.

    Returns:
        list[str]: The relative paths that are not files, in declaration order.
    """
    return [name for name in required if not (output_dir / name).is_file()]


def check_case(case: dict[str, Any], record: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Judge one completed case against its declared exit code and artefact set.

    Args:
        case: The matrix entry, carrying ``expect_exit`` and ``required_artifacts``.
        record: The runner's record for this case, carrying ``exit_code``, ``aborted``
            and ``timed_out``.
        output_dir: Where the case wrote.

    Returns:
        dict[str, Any]: ``expectation_met``, ``exit_problem``, ``missing_artifacts`` and
        ``expectation_problems`` - the last being every reason as a list, empty when the
        case is admissible.
    """
    expect_exit = case.get("expect_exit", "zero")
    exit_problem = check_exit(
        expect_exit,
        record.get("exit_code"),
        aborted=bool(record.get("aborted")),
        timed_out=bool(record.get("timed_out")),
    )
    required = list(case.get("required_artifacts", ()))
    missing = missing_artifacts(output_dir, required) if not exit_problem else []

    problems = [problem for problem in (exit_problem,) if problem]
    if missing:
        problems.append(
            f"exited as expected but did not write {len(missing)} of {len(required)} required artefact(s): "
            f"{', '.join(missing)}"
        )
    return {
        "expect_exit": expect_exit,
        "required_artifacts": required,
        "missing_artifacts": missing,
        "exit_problem": exit_problem,
        "expectation_problems": problems,
        "expectation_met": not problems,
    }


def describe_tree(
    tree: Path,
    *,
    paths: Sequence[str] = REVISION_PATHS,
    run: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> dict[str, Any]:
    """Record which revision a source tree is at, and whether it is clean.

    Args:
        tree: The side's source tree.
        paths: The paths whose working-tree state is recorded, relative to ``tree``.
        run: The subprocess runner, injectable so this is testable without a git repo.

    Returns:
        dict[str, Any]: ``head``, ``branch``, ``dirty`` (the whole tree), ``dirty_paths``
        (the porcelain lines under ``paths``), ``dirty_relevant`` and ``error``. Every
        field is None or empty when git could not answer; this function never raises,
        because a tree that is not a git checkout is a fact to record about the run, not a
        reason to lose the run.
    """
    info: dict[str, Any] = {
        "head": None,
        "branch": None,
        "dirty": None,
        "dirty_paths": [],
        "dirty_relevant": None,
        "revision_paths": list(paths),
        "error": None,
    }

    def _git(*args: str) -> str | None:
        try:
            completed = run(
                ["git", "-C", str(tree), *args],
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            info["error"] = f"{type(exc).__name__}: {exc}"
            return None
        if completed.returncode != 0:
            stderr = (completed.stderr or "").strip()[:300]
            info["error"] = stderr or f"git {' '.join(args)} exited {completed.returncode}"
            return None
        # Deliberately not stripped: `git status --porcelain` encodes the index and
        # worktree states in the first two characters, so " M path" and "M  path" are
        # different facts and a leading space is data.
        return completed.stdout

    head = _git("rev-parse", "HEAD")
    if head is None:
        return info
    info["head"] = head.strip()
    branch = _git("rev-parse", "--abbrev-ref", "HEAD")
    info["branch"] = branch.strip() if branch is not None else None

    whole = _git("status", "--porcelain")
    if whole is not None:
        info["dirty"] = bool(whole.strip())
    scoped = _git("status", "--porcelain", "--", *paths)
    if scoped is not None:
        lines = [line for line in scoped.splitlines() if line.strip()]
        info["dirty_paths"] = lines
        info["dirty_relevant"] = bool(lines)
    return info


def verify_revision(revision: dict[str, Any] | None, expected_sha: str | None, *, side: str) -> list[str]:
    """Check one side's recorded revision against what the caller expected.

    Args:
        revision: The mapping from :func:`describe_tree`, or None if the side predates
            revision recording.
        expected_sha: The SHA the caller asserts this side is at, full or abbreviated, or
            None to record without checking.
        side: ``before`` or ``after``, for the messages.

    Returns:
        list[str]: The problems, empty when there are none.
    """
    if expected_sha is None:
        return []
    if not revision or not revision.get("head"):
        return [
            f"{side}: --expect-{side}-sha {expected_sha} was given but this side recorded no revision "
            f"(error: {(revision or {}).get('error')}). A side whose revision is unknown cannot attest a commit."
        ]
    head = str(revision["head"])
    if not head.startswith(expected_sha) and not expected_sha.startswith(head):
        return [f"{side}: expected {expected_sha}, but the side ran at {head}"]
    return []


def check_sides_opposed(
    before_root: Path,
    after_root: Path,
    before_side: dict[str, Any],
    after_side: dict[str, Any],
) -> list[str]:
    """Check that two side records describe two genuinely opposed runs.

    The launcher proves each process ran the tree its caller named. Nothing proved the
    caller named two different trees, or labelled them the right way round, or asked for
    opposite marker states - so "tree compared against itself" stayed reachable by
    misconfiguration and would report a serene ``IDENTICAL``.

    Args:
        before_root: The baseline run root, as resolved.
        after_root: The candidate run root, as resolved.
        before_side: The baseline ``side.json``.
        after_side: The candidate ``side.json``.

    Returns:
        list[str]: The problems, empty when the two sides are properly opposed.
    """
    problems: list[str] = []

    if before_root == after_root:
        problems.append(f"--before-root and --after-root are the same directory ({before_root})")
    if before_side.get("side") != "before":
        problems.append(f"--before-root holds a side labelled {before_side.get('side')!r}, not 'before'")
    if after_side.get("side") != "after":
        problems.append(f"--after-root holds a side labelled {after_side.get('side')!r}, not 'after'")

    before_tree, after_tree = before_side.get("tree"), after_side.get("tree")
    if before_tree and before_tree == after_tree:
        problems.append(f"both sides ran the same source tree ({before_tree}); there is nothing to compare")

    before_marker, after_marker = before_side.get("marker"), after_side.get("marker")
    if before_marker != after_marker:
        problems.append(
            f"the two sides used different marker modules ({before_marker!r} vs {after_marker!r}), so neither "
            "side's resolution proof says anything about the other"
        )
    before_expect, after_expect = before_side.get("expect_marker"), after_side.get("expect_marker")
    if before_expect == after_expect:
        problems.append(
            f"both sides expected the marker {before_expect!r}; the marker distinguishes the two trees only when "
            "the expectations are opposed"
        )

    before_head = (before_side.get("revision") or {}).get("head")
    after_head = (after_side.get("revision") or {}).get("head")
    if before_head and after_head and before_head == after_head:
        problems.append(f"both sides ran the same commit ({before_head}); there is nothing to compare")

    # A side that ran nothing compares clean against another side that ran nothing, for the
    # same reason two failed cases do: there is no delta because there is no measurement.
    # `build_matrix` refuses a zero-case matrix and `run_side` reports a zero-case side as
    # unverified; this is the last of the three locks, and the only one that applies to a
    # run root the harness did not write.
    for label, side in (("before", before_side), ("after", after_side)):
        launched = len(side.get("pipeline_results") or {}) + len(side.get("cohort_results") or {})
        if launched == 0:
            problems.append(f"{label}: this side has no case results at all, so it measured nothing")

    return problems
