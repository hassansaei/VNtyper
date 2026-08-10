#!/usr/bin/env python3
"""
Advisory mutation scoring for VNtyper's pure-decision modules.

Line coverage answers "did a test execute this line?". It does not answer "would a
test have noticed if this line were wrong?" - and for this codebase the second
question is the one that matters, because the characteristic failure here is not a
crash but a silently wrong genotype call. ``confidence_assignment.py`` was the proof:
100% line coverage and a 21% mutation score, i.e. four out of five deliberate defects
in it went undetected by a fully green build.

This harness makes that measurement reproducible. It:

0. Refuses to start unless each selected target and requested real output is clean. The
   selected committed target bytes define the measurement, while real outputs are
   replaced only after the measurement is complete. See ``mutation_guard``.
1. Captures HEAD, creates a disposable detached worktree, and overlays the current
   non-ignored working state except selected targets and requested output paths.
2. Proves imports resolve inside that pinned workspace, then requires both a green
   unmutated baseline and a known-killed canary before measuring ordinary mutants.
3. Tokenizes each selected module and discards mutants that do not compile (for example,
   ``*args`` -> ``/args``).
4. Writes one mutant at a time only inside the disposable workspace, runs the tests,
   and restores the exact post-overlay baseline after every attempt.
5. Counts a mutant as **killed** if the tests fail and **survived** if they pass, installs
   complete requested outputs atomically, removes the worktree, and verifies that real
   production source retained its startup digests.

A survivor is a defect the suite cannot see. The score is ``killed / total``.

The green-baseline preflight
----------------------------
Step 4 is the whole measurement and it reads a single bit: pytest's return code. It
cannot tell "this mutant broke a test" from "this repository was already broken" - an
unrelated failing test, a collection error, a missing dependency, an import-time crash
or a timeout all return non-zero, and every one of them would be recorded as a kill.
A sufficiently broken environment therefore scores **100%**, prints a triumphant report
and overwrites the committed evidence page with it.

So :func:`baseline_refusal` runs, before a single byte is mutated, exactly the pytest
invocations the sweep will use to judge mutants:

* each selected target's scoped test files, serially, as the kill check does; and
* the whole unit tier in parallel, as the survivor escalation does - a broken full-tier
  run is the *more* dangerous of the two, because it converts every genuine survivor
  into a phantom "killed (full tier)".

If any of them fails the run aborts with the captured pytest output, and nothing is
written. This is not hypothetical: the ``.pyc`` warning below records two sweeps on this
branch whose results were entirely fictional.

Scoping and accuracy
--------------------
Running all 1559 unit tests for every mutant would take hours, so each target declares
the test files that exercise it and those run first with ``-x``. That alone would bias
the score *downward* - a mutant killed only by some other file's test would be
miscounted as a survivor - so any mutant that survives its scoped tests is **escalated
to the full unit tier** before being recorded as a survivor. Killed mutants (the common
case) stay fast; survivors are confirmed against everything.

This is ADVISORY. It is not gated in CI and it must not be: equivalent mutants - those
that change the source without changing behaviour, such as a bound that no input can
reach - have not been hand-classified, so the true score is somewhat higher than the
number printed here. Use it to find untested decisions, not as a pass/fail threshold.

.. warning::

   **The ``.pyc`` trap - read this before changing how mutants are written.**

   Two mutation sweeps on this branch produced fictional results before this was found.

   **The invariant: no child process may load bytecode generated for a different
   revision of a target module.** Everything below exists to hold that inside the
   disposable workspace, and nothing below is worth defending on its own terms.

   Why it is easy to break: CPython validates a cached ``.pyc`` against the source's
   **(mtime, size)** pair, and mtime has one-second granularity. A mutant that is
   written within the same second as the file it replaces *and* is the same byte length
   therefore has an identical (mtime, size), the interpreter considers the stale
   ``.pyc`` valid, loads the **unmutated** bytecode, and the tests pass. Every such
   mutant "survives" and the score is garbage.

   How much of the mutant population that covers, measured rather than asserted: over
   the six modules in :data:`TARGETS` at the time of writing, **46 of 121** generated
   mutants (38%) were byte-length preserving. Not most - but far more than enough, and
   they are the threshold flips (``==`` -> ``!=``, ``-`` -> ``+``, ``1`` -> ``2``) that
   matter most. The exact split moves with the source; the hazard does not.

   Two examples that were previously given here are wrong and are corrected so nobody
   re-derives the rule from them: :data:`OPERATOR_MUTATIONS` maps ``<`` -> ``>=``, so
   the harness never performs ``<`` -> ``>``; and ``and`` -> ``or`` is three characters
   to two, which *changes* the size and so is caught by the cache validator anyway.

   Two defences, both required and both applied below:
     * ``python -B`` **and** ``PYTHONDONTWRITEBYTECODE=1`` in the child, so no ``.pyc``
       is written during the run at all. Either flag alone already disables bytecode
       writing; :func:`run_pytest` sets both, unconditionally, so the defence does not
       depend on how the harness was invoked. The ``PYTHONDONTWRITEBYTECODE=1`` on the
       ``make mutation`` recipe is defence in depth for the *parent*, which never
       imports a target module - it is not what holds the invariant.
     * every ``__pycache__`` under ``vntyper/`` is deleted before the sweep starts *and*
       again after each mutant is written, immediately before its tests run - so no
       ``.pyc`` from an earlier run or an earlier mutant can be picked up either.

   Neither alone is sufficient: the flags stop new caches, the deletion stops old ones.
   Do not remove either, and do not "optimise" the sweep by leaving caches in place.

Usage:
    python scripts/mutation_test.py                       # all targets
    python scripts/mutation_test.py --module scoring.py   # one target
    python scripts/mutation_test.py --output report.txt   # also write to a file
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import tokenize
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import FrameType
from typing import NoReturn

from mutation_guard import dirty_paths, format_dirty_tree_refusal, writable_paths
from mutation_output import atomic_write_text, validate_disjoint_paths
from mutation_source import WorkspaceRoot, atomic_replace_source, read_source_bytes
from mutation_workspace import detached_head_workspace, verify_import_provenance
from mutation_workspace_fs import git_capability_path, using_root_capability

REAL_REPO_ROOT = Path(__file__).resolve().parents[1]

#: Target module -> the unit test files that exercise it.
#:
#: These five are the pure-decision modules: they turn variants and motifs into the
#: confidence, score and flag values that end up in the report. They are the modules
#: where a wrong branch is invisible, which is exactly where mutation testing pays.
#: The scoped files are a fast path only - see "Scoping and accuracy" above.
#: Keep these lists WIDE. Every file that so much as imports the module belongs here,
#: found with `grep -rln <module> tests/unit/`. A mutant that survives the scoped run is
#: escalated to the full tier, which is correct but costs ~20s, so a missing entry does
#: not corrupt the score - it just makes the sweep crawl.
TARGETS: dict[str, tuple[str, ...]] = {
    "vntyper/scripts/confidence_assignment.py": (
        "tests/unit/test_confidence_assignment.py",
        "tests/unit/test_confidence_boundaries.py",
        "tests/unit/test_calibration_consistency.py",
        "tests/unit/test_scoring.py",
        "tests/unit/test_builders.py",
    ),
    "vntyper/scripts/scoring.py": (
        "tests/unit/test_scoring.py",
        "tests/unit/test_calibration_consistency.py",
        "tests/unit/test_confidence_assignment.py",
    ),
    "vntyper/scripts/motif_processing.py": (
        "tests/unit/test_motif_filtering_issue_136.py",
        "tests/unit/test_motif_preprocessing.py",
        "tests/unit/test_motif_decisions.py",
        "tests/unit/test_motif_config_guard.py",
        "tests/unit/test_kestrel_filtering.py",
        "tests/unit/test_generate_report.py",
    ),
    # The pure decision layer extracted out of motif_correction_and_annotation.
    # Registered in the same commit as the extraction, deliberately: a module absent
    # from TARGETS is not mutated at all, so extracting the hard decisions without
    # adding them here would have made motif_processing.py's score rise for the worst
    # possible reason - the decisions leaving the measurement rather than becoming
    # tested. The two scores are only comparable to the 30.9% baseline together.
    "vntyper/scripts/motif_decisions.py": (
        "tests/unit/test_motif_decisions.py",
        "tests/unit/test_motif_filtering_issue_136.py",
        "tests/unit/test_motif_config_guard.py",
        "tests/unit/test_kestrel_filtering.py",
        "tests/unit/test_generate_report.py",
    ),
    "vntyper/scripts/variant_parsing.py": ("tests/unit/test_variant_parsing.py",),
    "vntyper/scripts/flagging.py": (
        "tests/unit/test_flagging.py",
        "tests/unit/test_advntr_flagging_rules.py",
        "tests/unit/test_advntr_command.py",
        "tests/unit/test_screening_summary.py",
        "tests/unit/test_generate_report.py",
    ),
}

#: Known high-signal mutant that the scoring unit tests must kill. The harness refuses
#: a measurement unless this exact identity produces its exact assertion failure and no
#: infrastructure/plugin error while pytest exits exactly 1.
CANARY_KEY = ("vntyper/scripts/scoring.py", 74, "/", "*")
CANARY_EXPECTED_NODE = "tests/unit/test_scoring.py::test_split_depth_and_calculate_frame_score_no_frameshift"
CANARY_EXPECTED_ASSERTION = "E       AssertionError: Frame_Score should be 1.0"
CANARY_INFRASTRUCTURE_MARKERS = (
    "INTERNALERROR",
    "ERROR collecting",
    "ERROR at setup",
    "ERROR at teardown",
    "ERROR:",
    "ERROR ",
    "ImportError",
    "ModuleNotFoundError",
    "PluginValidationError",
    "Traceback (most recent call last)",
)

#: Operator token substitutions. Each flips a decision without changing what the
#: surrounding statement *is*, which is what makes an undetected one a real defect
#: rather than a crash: `>=` -> `<` moves a threshold, it does not raise.
OPERATOR_MUTATIONS: dict[str, str] = {
    "<": ">=",
    "<=": ">",
    ">": "<=",
    ">=": "<",
    "==": "!=",
    "!=": "==",
    "+": "-",
    "-": "+",
    "*": "/",
    "/": "*",
    "+=": "-=",
    "-=": "+=",
}

#: Keyword substitutions. `and`/`or` flips rule composition; `True`/`False` flips a
#: literal outcome. Both are the shape of bug this codebase actually ships.
#:
#: `not` -> `` (deletion) is the highest-signal operator of the set for this codebase:
#: it inverts a guard while leaving the statement's shape intact, which is precisely
#: the "stages mark, they do not filter" logic in flagging.py and motif_processing.py.
#: The empty replacement leaves a harmless double space (`if not x` -> `if  x`); the
#: compile check below rejects the cases where it would not parse.
KEYWORD_MUTATIONS: dict[str, str] = {
    "and": "or",
    "or": "and",
    "True": "False",
    "False": "True",
    "not": "",
    "in": "not in",
    "break": "continue",
    "continue": "break",
}


def mutate_number(literal: str) -> str | None:
    """
    Return a numeric literal one greater than the given one, or None if unparseable.

    Off-by-one on a numeric literal is the classic threshold defect, and these five
    modules are threshold logic. Underscored, hex and imaginary literals are skipped
    rather than guessed at.

    Args:
        literal (str): The source text of a NUMBER token.

    Returns:
        str | None: The incremented literal, or None if it is not a plain int/float.
    """
    try:
        if any(c in literal.lower() for c in "xob_j"):
            return None
        if "." in literal or "e" in literal.lower():
            return repr(float(literal) + 1.0)
        return str(int(literal) + 1)
    except ValueError:
        return None


#: Survivors hand-classified as EQUIVALENT: the mutation changes the source without
#: changing any behaviour a test could observe, so no test could ever kill it and it is
#: not a gap in the suite.
#:
#: Keyed by ``(repo-relative path, line, original token, replacement)`` and valued with
#: the reason. Keeping this as data rather than prose means a re-run re-applies the
#: classification automatically instead of silently reverting to a raw score.
#:
#: Only add an entry you have reasoned through and can justify in one line. An
#: unjustified entry here is worse than a low score: it is a real gap papered over.
#: Line numbers are checked against the live sweep, so a stale entry is reported rather
#: than silently ignored.
EQUIVALENT_MUTANTS: dict[tuple[str, int, str, str], str] = {
    # --- confidence_assignment.py -------------------------------------------------
    # Six entries used to live here, one per `<subdict>.get(<key>, <default>)` read of a
    # calibration constant. They are gone because the code is: the six defaults were
    # DELETED (#184 follow-up), not reclassified. The constants are now read as direct
    # subscripts, so a missing key raises KeyError instead of silently substituting a
    # second, wrong calibration - and there is no default operand left to mutate.
    # --- variant_parsing.py -------------------------------------------------------
    ("vntyper/scripts/variant_parsing.py", 114, "0.0", "1.0"): (
        "`.get()` default for `alt_filtering.gg_depth_score_threshold`; the shipped config supplies 0.00469"
    ),
    # --- motif_processing.py ------------------------------------------------------
    # Line moved 315 -> 342 when 11e2300 extracted the decision layer. The sweep's stale-
    # entry check caught the drift; the mutant is the same one and is still equivalent.
    ("vntyper/scripts/motif_processing.py", 342, "60", "61"): (
        "`.get()` default for `motif_filtering.position_threshold`; the shipped config supplies 60"
    ),
    # --- flagging.py --------------------------------------------------------------
    # NOTE the deliberately narrow scope here. The shipped config sets
    # `duplicate_flagging.enabled = false`, so the whole duplicate-marking block is dead
    # by default - but `enabled` is a SUPPORTED TOGGLE, so that code is reachable by
    # configuration and its mutants are real gaps, not equivalents. Only the `.get()`
    # default operand itself is classified here.
    ("vntyper/scripts/flagging.py", 150, "False", "True"): (
        "`.get()` default for `duplicate_flagging.enabled`; the shipped config supplies it explicitly"
    ),
    # `itertuples(index=True)` only prepends an `Index` field to the namedtuple; the loop
    # body reads `row_tuple.Flag` by name and never touches position, so both forms yield
    # the same value. Checked against the awkward cases too - a column literally named
    # `Index`, duplicate `Flag` columns, a string index and a MultiIndex all resolve
    # `.Flag` to the same column either way, because namedtuple's `rename=True` renames
    # the colliding field and not `Flag`.
    # Line moved 238 -> 243 with the duplicate-flagging changes in 6e7cda2.
    ("vntyper/scripts/flagging.py", 243, "False", "True"): (
        "`itertuples(index=)` only adds an `Index` field; the loop reads `row_tuple.Flag` by name"
    ),
}


@dataclass
class Mutant:
    """One candidate mutation of a single source file."""

    path: str
    line: int
    original: str
    replacement: str
    source: str
    #: Set once the mutant has been run.
    killed: bool | None = None
    killed_by: str = ""

    def describe(self) -> str:
        """Return a one-line human-readable identifier for this mutant."""
        return f"{self.path}:{self.line}  {self.original!r} -> {self.replacement!r}"

    @property
    def key(self) -> tuple[str, int, str, str]:
        """Return this mutant's key into :data:`EQUIVALENT_MUTANTS`."""
        return (self.path, self.line, self.original, self.replacement)

    @property
    def equivalence_reason(self) -> str | None:
        """Return why this mutant is equivalent, or None if it is a genuine gap."""
        return EQUIVALENT_MUTANTS.get(self.key)


@dataclass
class ModuleResult:
    """The outcome of a full sweep over one module."""

    path: str
    killed: int = 0
    survived: int = 0
    survivors: list[Mutant] = field(default_factory=list)

    @property
    def total(self) -> int:
        """Total number of mutants actually run for this module."""
        return self.killed + self.survived

    @property
    def score(self) -> float:
        """Mutation score as a percentage, or 0.0 when no mutants ran."""
        return 100.0 * self.killed / self.total if self.total else 0.0


def generate_mutants(path: Path, *, repo_root: Path, source: str | None = None) -> list[Mutant]:
    """
    Generate every compilable single-token mutant of a source file.

    Tokenizing rather than parsing is deliberate: ``ast`` does not attach source
    positions to operator nodes (``ast.Lt`` and friends carry no ``col_offset``), so
    an AST walk cannot say *where* to rewrite a ``<``. The token stream gives an exact
    span for every token, and comments and strings arrive as their own token types so
    they are skipped for free - a mutation inside a docstring would be a guaranteed
    survivor and pure noise.

    Args:
        path (Path): The module to mutate.
        repo_root (Path): Repository root used to make mutant paths relative.
        source (str | None): Already anchored source text. Direct callers may omit it
            to read ``path`` normally.

    Returns:
        list[Mutant]: One mutant per mutable token whose mutated form still compiles.
            Mutants that do not compile (``*args`` -> ``/args``) are dropped rather
            than counted, since a syntax error is not a defect a test could miss.
    """
    if source is None:
        source = path.read_text(encoding="utf-8")
    lines = source.splitlines(keepends=True)
    line_starts: list[int] = []
    offset = 0
    for line in lines:
        line_starts.append(offset)
        offset += len(line)

    def absolute(row: int, col: int) -> int:
        """Convert a 1-based (row, col) token position to a string offset."""
        return line_starts[row - 1] + col

    mutants: list[Mutant] = []
    readline = io.StringIO(source).readline
    for token in tokenize.generate_tokens(readline):
        if token.type == tokenize.OP:
            replacement = OPERATOR_MUTATIONS.get(token.string)
        elif token.type == tokenize.NAME:
            replacement = KEYWORD_MUTATIONS.get(token.string)
        elif token.type == tokenize.NUMBER:
            replacement = mutate_number(token.string)
        else:
            continue
        if replacement is None:
            continue

        start = absolute(token.start[0], token.start[1])
        end = absolute(token.end[0], token.end[1])
        mutated = source[:start] + replacement + source[end:]

        # A mutant that will not compile is not a defect any test could have caught.
        try:
            compile(mutated, str(path), "exec")
        except SyntaxError:
            continue

        mutants.append(
            Mutant(
                path=str(path.relative_to(repo_root)),
                line=token.start[0],
                original=token.string,
                replacement=replacement,
                source=mutated,
            )
        )
    return mutants


def clear_bytecode_caches(root: Path) -> int:
    """
    Delete every ``__pycache__`` directory under ``root``.

    Half of the defence described in the module-level ``.. warning::``. Setting
    ``PYTHONDONTWRITEBYTECODE`` stops *new* caches being written; it does nothing about
    a ``.pyc`` left behind by an earlier run, which CPython will happily load if its
    recorded (mtime, size) still matches the source it is shadowing. Deleting them is
    the only way to be sure the interpreter compiles the mutant in front of it.

    Args:
        root (Path): Directory to clean recursively.

    Returns:
        int: Number of cache directories removed.
    """
    removed = 0
    for cache in root.rglob("__pycache__"):
        shutil.rmtree(cache, ignore_errors=True)
        removed += 1
    return removed


@dataclass(frozen=True)
class TestRun:
    """
    The outcome of one pytest invocation, including what it printed.

    The sweep only needs :attr:`passed`, but the preflight needs :attr:`output`: a
    refusal that says "the baseline is red" without saying *how* is a refusal nobody can
    act on, and discarding the output is what let a broken environment look like a
    perfect score for as long as it did.
    """

    passed: bool
    output: str
    returncode: int | None
    timed_out: bool


def run_pytest(
    test_paths: tuple[str, ...],
    timeout: int,
    parallel: bool = False,
    *,
    repo_root: WorkspaceRoot,
) -> TestRun:
    """
    Run pytest over the given paths and return the outcome and the captured output.

    Args:
        test_paths (tuple[str, ...]): Test files or directories to run.
        timeout (int): Seconds before the run is abandoned.
        parallel (bool): Distribute across cores with xdist. Used for the full-tier
            escalation, which is the expensive path; the scoped runs are short enough
            that xdist's worker startup would cost more than it saves.
        repo_root (WorkspaceRoot): Disposable repository root capability used as the
            child CWD and first import path.

    Returns:
        TestRun: ``passed`` is True only if pytest exited 0. A timeout counts as a
            failure - a mutant that sends the suite into an infinite loop has been
            detected, not missed - and reports that fact in ``output``.

    Note:
        ``-B`` and ``PYTHONDONTWRITEBYTECODE`` are the other half of the ``.pyc``
        defence, and either alone would disable bytecode writing; both are set so that
        the invariant does not depend on how the harness was invoked.
        ``-p no:cacheprovider`` keeps pytest from writing ``.pytest_cache``
        entries that would make ``--lf`` behave differently between mutants, and
        ``-o log_cli=false`` suppresses the live log pytest.ini turns on, which is pure
        overhead during the sweep itself.
    """
    for node in test_paths:
        path_part = node.split("::", maxsplit=1)[0]
        candidate = Path(path_part)
        if not path_part or candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError(f"pytest node must be workspace-relative: {node}")

    with using_root_capability(repo_root) as capability:
        child_root = git_capability_path(capability)
        env = dict(os.environ)
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        env.pop("PYTEST_ADDOPTS", None)
        env.pop("PYTEST_PLUGINS", None)
        resolved_real_root = REAL_REPO_ROOT.resolve(strict=True)
        retained_python_path: list[str] = []
        for entry in env.get("PYTHONPATH", "").split(os.pathsep):
            if not entry:
                continue
            resolved_entry = Path(entry).resolve(strict=False)
            if resolved_entry == resolved_real_root or resolved_entry.is_relative_to(resolved_real_root):
                continue
            retained_python_path.append(entry)
        env["PYTHONPATH"] = os.pathsep.join((str(child_root), *retained_python_path))
        command = [
            sys.executable,
            "-B",
            "-m",
            "pytest",
            "-m",
            "unit",
            "-x",
            "-q",
            "--no-header",
            "-p",
            "no:cacheprovider",
            "-p",
            "pytest_timeout",
            "-o",
            "log_cli=false",
            "--timeout=120",
        ]
        if parallel:
            # Only whether ANY test failed matters here, so xdist's non-deterministic
            # interaction with -x is irrelevant - it still reports a non-zero exit.
            command += ["-p", "xdist.plugin", "-n", "4", "--dist", "loadfile"]
        command += [*test_paths]
        try:
            completed = subprocess.run(
                command,
                cwd=child_root,
                env=env,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
                pass_fds=(capability.descriptor,),
            )
        except subprocess.TimeoutExpired:
            return TestRun(
                passed=False,
                output=f"pytest exceeded the {timeout}s timeout and was abandoned.",
                returncode=None,
                timed_out=True,
            )
    # `or ""` because a subprocess double may leave these unset; a missing stream is not
    # a reason for the preflight to crash instead of reporting.
    output = (completed.stdout or "") + (completed.stderr or "")
    return TestRun(
        passed=completed.returncode == 0,
        output=output,
        returncode=completed.returncode,
        timed_out=False,
    )


def run_tests(
    test_paths: tuple[str, ...],
    timeout: int,
    parallel: bool = False,
    *,
    repo_root: WorkspaceRoot,
) -> bool:
    """
    Run pytest and report only whether it passed.

    The sweep's kill check is exactly this bit, so it gets its own name; everything that
    needs the output calls :func:`run_pytest` directly.

    Args:
        test_paths (tuple[str, ...]): Test files or directories to run.
        timeout (int): Seconds before the run is abandoned.
        parallel (bool): Distribute across cores with xdist.
        repo_root (WorkspaceRoot): Disposable repository root used by the pytest child.

    Returns:
        bool: True if every test passed.
    """
    return run_pytest(test_paths, timeout, parallel, repo_root=repo_root).passed


def format_baseline_refusal(label: str, test_paths: tuple[str, ...], output: str) -> str:
    """
    Render the refusal printed when the unmutated tree does not pass its own tests.

    Args:
        label (str): What was being checked - a target path, or the escalation tier.
        test_paths (tuple[str, ...]): The pytest arguments that failed.
        output (str): Everything pytest printed, verbatim.

    Returns:
        str: The message to print on stderr before exiting non-zero.
    """
    return "\n".join(
        [
            "REFUSING TO SWEEP: the tests failed on the UNMUTATED tree.",
            "",
            f"  checking: {label}",
            f"  pytest:   {' '.join(test_paths)}",
            "",
            "A mutant is recorded as KILLED whenever pytest exits non-zero, and pytest",
            "exits non-zero for an unrelated failure, a collection error, a missing",
            "dependency or a timeout just as readily as for a detected mutant. Sweeping",
            "now would score those as kills and publish the result - a broken checkout",
            "scores 100%.",
            "",
            "Fix the failure below, then re-run. Nothing has been mutated or written.",
            "",
            "-" * 70,
            output.rstrip() or "(pytest produced no output)",
            "-" * 70,
        ]
    )


def baseline_refusal(
    targets: dict[str, tuple[str, ...]],
    timeout: int,
    *,
    repo_root: WorkspaceRoot,
) -> str | None:
    """
    Run the sweep's own test invocations against the unmutated tree.

    Both invocations the sweep uses are checked, because both feed the score and they
    fail independently:

    * the scoped run per target decides ``killed``;
    * the parallel full-tier run decides whether a scoped survivor is *recorded* as a
      survivor, so if it is broken every genuine gap silently becomes a kill.

    Args:
        targets (dict[str, tuple[str, ...]]): The targets this run will sweep, mapped to
            their scoped test files. Narrowed by ``--module`` exactly as the sweep is.
        timeout (int): Per-pytest-run timeout in seconds.
        repo_root (WorkspaceRoot): Disposable repository root used for every baseline
            child.

    Returns:
        str | None: The refusal to print, or ``None`` when every check passed. Returned
            rather than raised so ``main()`` keeps its single exit point.
    """
    checks: list[tuple[str, tuple[str, ...], bool]] = [
        (f"{path_str} (scoped kill check)", scoped, False) for path_str, scoped in targets.items()
    ]
    # Last, because it is the slow one and a scoped failure is the cheaper diagnosis.
    checks.append(("full unit tier (survivor escalation)", ("tests/unit",), True))

    for label, test_paths, parallel in checks:
        print(f"  baseline: {label} ...", end="", flush=True)
        run = run_pytest(test_paths, timeout, parallel=parallel, repo_root=repo_root)
        print(" ok" if run.passed else " FAILED", flush=True)
        if not run.passed:
            return format_baseline_refusal(label, test_paths, run.output)
    return None


def sweep_module(
    path_str: str,
    scoped_tests: tuple[str, ...],
    timeout: int,
    verbose: bool,
    *,
    repo_root: WorkspaceRoot,
) -> ModuleResult:
    """
    Run every mutant of one module and record which survived.

    The original file is restored in a ``finally`` block, so an exception or a
    ``KeyboardInterrupt`` mid-sweep cannot leave a mutated module on disk.

    "Original" here means *the text on disk when this call started*, not the committed
    text - the harness has no way to tell the two apart. ``main()`` therefore refuses to
    begin unless every target is clean (see ``mutation_guard``), and that precondition is
    what makes the restore, the score and the end-of-run check mean what they say.

    Args:
        path_str (str): Repo-relative path of the module to mutate.
        scoped_tests (tuple[str, ...]): Test files that exercise this module.
        timeout (int): Per-pytest-run timeout in seconds.
        verbose (bool): Print each mutant's outcome as it is decided.
        repo_root (WorkspaceRoot): Disposable repository root containing the mutation
            target.

    Returns:
        ModuleResult: Kill/survive counts and the surviving mutants.
    """
    with using_root_capability(repo_root) as capability:
        child_root = git_capability_path(capability)
        path = child_root / path_str
        original_bytes = read_source_bytes(capability, path_str)
        original = original_bytes.decode("utf-8")
        mutants = generate_mutants(path, repo_root=child_root, source=original)
        result = ModuleResult(path=path_str)

        # flush= throughout: a sweep is long enough that someone will watch the log, and an
        # interrupted run should still show how far it got.
        print(f"\n{path_str}: {len(mutants)} mutants", flush=True)
        try:
            for index, mutant in enumerate(mutants, start=1):
                atomic_replace_source(capability, path_str, mutant.source.encode("utf-8"))
                # Between writing the mutant and running its tests, every time - not once at
                # the start. This is half the `.pyc` invariant in the module docstring: the
                # child must not be able to load bytecode built from the previous mutant, or
                # from the original, both of which can share this file's (mtime, size).
                # tests/unit/test_mutation_test.py::test_no_bytecode_cache_survives_into_any_mutants_test_run
                # is the only thing that notices if this line goes.
                clear_bytecode_caches(child_root / "vntyper")

                if not run_tests(scoped_tests, timeout, repo_root=capability):
                    mutant.killed, mutant.killed_by = True, "scoped"
                # Survived the scoped tests - confirm against the whole tier before
                # believing it, so the score is not biased down by the scoping.
                elif not run_tests(("tests/unit",), timeout, parallel=True, repo_root=capability):
                    mutant.killed, mutant.killed_by = True, "full tier"
                else:
                    mutant.killed = False

                if mutant.killed:
                    result.killed += 1
                else:
                    result.survived += 1
                    result.survivors.append(mutant)

                if verbose:
                    status = f"killed ({mutant.killed_by})" if mutant.killed else "SURVIVED"
                    print(f"  [{index}/{len(mutants)}] {mutant.describe()} -> {status}", flush=True)
                else:
                    print(f"  [{index}/{len(mutants)}] {'.' if mutant.killed else 'S'}", end="", flush=True)
        finally:
            atomic_replace_source(capability, path_str, original_bytes)
            clear_bytecode_caches(child_root / "vntyper")

    if not verbose:
        print()
    print(f"  {result.killed}/{result.total} killed ({result.score:.1f}%)", flush=True)
    return result


def canary_refusal(timeout: int, *, repo_root: WorkspaceRoot) -> str | None:
    """Prove one exact scoring mutant is killed for the expected reason.

    Args:
        timeout (int): Maximum seconds for the canary's scoped pytest run.
        repo_root (WorkspaceRoot): Disposable repository root containing the canary
            target.

    Returns:
        str | None: A refusal diagnostic, or ``None`` only when pytest exits exactly 1
            and prints the expected node and assertion evidence.
    """
    path_str, line, original_token, replacement = CANARY_KEY
    with using_root_capability(repo_root) as capability:
        child_root = git_capability_path(capability)
        path = child_root / path_str
        original_bytes = read_source_bytes(capability, path_str)
        original = original_bytes.decode("utf-8")
        mutants = generate_mutants(path, repo_root=child_root, source=original)
        canary = next((mutant for mutant in mutants if mutant.key == CANARY_KEY), None)
        if canary is None:
            return (
                "REFUSING TO SWEEP: canary mutant is missing: "
                f"path={path_str!r}, line={line}, original={original_token!r}, replacement={replacement!r}"
            )

        try:
            atomic_replace_source(capability, path_str, canary.source.encode("utf-8"))
            clear_bytecode_caches(child_root / "vntyper")
            run = run_pytest(TARGETS[path_str], timeout, repo_root=capability)
        finally:
            atomic_replace_source(capability, path_str, original_bytes)
            clear_bytecode_caches(child_root / "vntyper")

    if run.timed_out:
        return f"REFUSING TO SWEEP: canary timed out.\n{run.output}"
    output_lines = run.output.splitlines()
    failure_summaries = [line for line in output_lines if line.startswith("FAILED ")]
    expected_failure_summary = f"FAILED {CANARY_EXPECTED_NODE}"
    expected_evidence = (
        CANARY_EXPECTED_NODE in run.output
        and CANARY_EXPECTED_ASSERTION in run.output
        and len(failure_summaries) == 1
        and (
            failure_summaries[0] == expected_failure_summary
            or failure_summaries[0].startswith(f"{expected_failure_summary} ")
        )
        and any(line.strip("= ").startswith("1 failed") for line in output_lines)
        and not any(marker in run.output for marker in CANARY_INFRASTRUCTURE_MARKERS)
    )
    if run.returncode == 1 and expected_evidence:
        return None
    if run.returncode == 0:
        return f"REFUSING TO SWEEP: canary survived.\n{run.output}"
    return (
        f"REFUSING TO SWEEP: canary infrastructure failure (pytest exit {run.returncode}); "
        f"expected {CANARY_EXPECTED_NODE!r} and {CANARY_EXPECTED_ASSERTION!r}.\n{run.output}"
    )


def format_report(results: list[ModuleResult], elapsed: float) -> str:
    """
    Render the raw console report for a completed sweep.

    This exact text is embedded verbatim in the committed Markdown page, so it is the
    evidence artefact rather than a summary of one.

    Args:
        results (list[ModuleResult]): One entry per module swept.
        elapsed (float): Wall-clock duration of the sweep in seconds.

    Returns:
        str: A plain-text report naming every surviving mutant.
    """
    killed = sum(r.killed for r in results)
    total = sum(r.total for r in results)
    score = 100.0 * killed / total if total else 0.0

    out: list[str] = []
    out.append("VNtyper mutation testing - advisory score")
    out.append("=" * 60)
    out.append("")
    out.append("Command:  make mutation")
    out.append(f"Total:    {total} mutants, {killed} killed, {total - killed} survived")
    out.append(f"Score:    {score:.1f}%")
    out.append(f"Duration: {elapsed / 60:.1f} min")
    out.append("")
    out.append("Per module")
    out.append("-" * 60)
    out.extend(f"{r.score:6.1f}%  {r.killed:3d}/{r.total:3d}  {r.path}" for r in sorted(results, key=lambda r: r.score))
    out.append("")
    out.append("Surviving mutants  [E] = hand-classified equivalent, [ ] = genuine gap")
    out.append("-" * 60)
    any_survivors = False
    for r in sorted(results, key=lambda r: r.path):
        if not r.survivors:
            continue
        any_survivors = True
        out.append(f"{r.path}")
        for mutant in r.survivors:
            reason = mutant.equivalence_reason
            marker = "E" if reason else " "
            out.append(f"  [{marker}] line {mutant.line:4d}  {mutant.original!r} -> {mutant.replacement!r}")
            if reason:
                out.append(f"          equivalent: {reason}")
        out.append("")
    if not any_survivors:
        out.append("None. Every mutant was killed.")
        out.append("")

    stale = stale_equivalence_keys(results)
    if stale:
        out.append("STALE EQUIVALENCE ENTRIES (no longer match a surviving mutant)")
        out.append("-" * 60)
        out.extend(f"    {key}" for key in stale)
        out.append("")
    return "\n".join(out)


def all_survivors(results: list[ModuleResult]) -> list[Mutant]:
    """Return every surviving mutant across all modules, ordered by path and line."""
    survivors = [m for r in results for m in r.survivors]
    return sorted(survivors, key=lambda m: (m.path, m.line))


def stale_equivalence_keys(results: list[ModuleResult]) -> list[tuple[str, int, str, str]]:
    """
    Return equivalence entries that no longer correspond to a surviving mutant.

    An entry goes stale when the module is edited and its line numbers shift, or when a
    new test finally kills the mutant. Either way the entry is now excusing something
    that is not there, so it is surfaced rather than silently ignored - a silent stale
    entry is how a real gap ends up permanently classified away.

    Args:
        results (list[ModuleResult]): One entry per module swept.

    Returns:
        list[tuple[str, int, str, str]]: The keys with no matching survivor. Only
            modules actually swept are considered, so a scoped ``--module`` run does
            not report every other module's entries as stale.
    """
    swept_paths = {r.path for r in results}
    live = {m.key for m in all_survivors(results)}
    return sorted(key for key in EQUIVALENT_MUTANTS if key[0] in swept_paths and key not in live)


def format_markdown(results: list[ModuleResult], elapsed: float) -> str:
    """
    Render the committed docs page, with the raw console output embedded verbatim.

    Args:
        results (list[ModuleResult]): One entry per module swept.
        elapsed (float): Wall-clock duration of the sweep in seconds.

    Returns:
        str: A Markdown page suitable for ``docs/development/``.
    """
    killed = sum(r.killed for r in results)
    total = sum(r.total for r in results)
    score = 100.0 * killed / total if total else 0.0

    survivors = all_survivors(results)
    equivalent = [m for m in survivors if m.equivalence_reason]
    genuine = [m for m in survivors if not m.equivalence_reason]
    adjusted_total = total - len(equivalent)
    adjusted = 100.0 * killed / adjusted_total if adjusted_total else 0.0

    out: list[str] = []
    out.append("# Mutation Testing")
    out.append("")
    out.append('!!! warning "Advisory only - nothing gates on this number"')
    out.append("")
    out.append("    This score is **not** a pass/fail threshold and CI does not bind against")
    out.append("    it. Read it as a map of which decisions are untested, not as a grade.")
    out.append("")
    out.append("Line coverage answers *did a test execute this line?*. It does not answer")
    out.append("*would a test have noticed if this line were wrong?* - and for VNtyper the")
    out.append("second question is the one that matters, because the characteristic failure")
    out.append("is a silently wrong genotype call rather than a crash.")
    out.append("")
    out.append("Mutation testing answers the second question directly: it introduces a")
    out.append("deliberate defect, runs the tests, and records whether anything failed. A")
    out.append("**surviving** mutant is a defect the suite cannot see.")
    out.append("")
    out.append("## Result")
    out.append("")
    out.append(f"**{killed} of {total} mutants killed - a raw mutation score of {score:.1f}%.**")
    out.append("")
    out.append(f"Of the {len(survivors)} survivors, {len(equivalent)} are hand-classified as")
    out.append("*equivalent* (the mutation cannot change observable behaviour, so no test could")
    out.append(f"ever kill it) and {len(genuine)} are genuine gaps. Excluding the equivalent")
    out.append(f"mutants the score is **{adjusted:.1f}%** ({killed}/{adjusted_total}).")
    out.append("")
    out.append("Both numbers are given because neither alone is honest: the raw score")
    out.append("understates the suite by counting unkillable mutants against it, and the")
    out.append("adjusted score depends on classifications that are a human judgement call.")
    out.append("Every classification is listed below with its reason so it can be checked.")
    out.append("")
    out.append("| Module | Killed | Total | Raw score |")
    out.append("| --- | ---: | ---: | ---: |")
    out.extend(
        f"| `{r.path}` | {r.killed} | {r.total} | {r.score:.1f}% |" for r in sorted(results, key=lambda r: r.score)
    )
    out.append("")
    out.append("## Surviving mutants")
    out.append("")
    out.append("### Genuine gaps")
    out.append("")
    if genuine:
        out.append("Each of these is a change to the source that the **entire** unit tier does")
        out.append("not notice. They are the actionable output of this exercise: a test that")
        out.append("kills one is a test that would have caught a real defect of that shape.")
        out.append("")
        out.append("| Module | Line | Mutation |")
        out.append("| --- | ---: | --- |")
        out.extend(
            f"| `{m.path}` | {m.line} | `{m.original}` &rarr; `{m.replacement or '(deleted)'}` |" for m in genuine
        )
    else:
        out.append("None - every survivor is classified equivalent below.")
    out.append("")
    out.append("### Classified equivalent")
    out.append("")
    if equivalent:
        out.append("These mutations cannot change behaviour that any test could legitimately")
        out.append("observe, so they are **not** gaps in the suite. Each reason below is a claim")
        out.append("you can check against the source; if one turns out to be wrong the entry")
        out.append("should be deleted, not the score explained away.")
        out.append("")
        out.append("Most of them are `.get()` defaults on `kestrel_config.json` keys that the")
        out.append("shipped config always supplies, which makes the default value dead code.")
        out.append("Being precise about the scope of that claim: a `--config-path` omitting the")
        out.append("key *would* reach the default, so these are unreachable **with the shipped")
        out.append("configuration** rather than unreachable in principle. That is the right")
        out.append("standard here - `AGENTS.md` trap 2 records that `--config-path` replaces the")
        out.append("whole config rather than merging it, and that a partial config already fails")
        out.append("with `KeyError` elsewhere in the pipeline, so a config missing these keys is")
        out.append("not a supported input.")
        out.append("")
        out.append("| Module | Line | Mutation | Why it cannot be killed |")
        out.append("| --- | ---: | --- | --- |")
        out.extend(
            f"| `{m.path}` | {m.line} | `{m.original}` &rarr; `{m.replacement or '(deleted)'}` "
            f"| {m.equivalence_reason} |"
            for m in equivalent
        )
    else:
        out.append("None classified yet.")
    out.append("")
    out.append("## How this compares to the 43.5% baseline")
    out.append("")
    out.append("The experiment that motivated this work scored **43.5%** (27 of 62 mutants")
    out.append("killed) across the eight highest-coverage modules. That harness was never")
    out.append("committed, which is why this one exists.")
    out.append("")
    out.append('!!! note "The two totals are not directly comparable"')
    out.append("")
    out.append("    Different mutant population, different modules: 62 mutants over eight")
    out.append(f"    modules then, {total} over five now, generated by a different operator set.")
    out.append("    A higher or lower headline number would not by itself mean the suite has")
    out.append("    improved or regressed. Only per-module figures on the same module carry")
    out.append("    across, and even those only loosely.")
    out.append("")
    out.append("The one comparison that is meaningful is `confidence_assignment.py`, the")
    out.append("module that motivated the whole effort: it had **100% line coverage and a 21%")
    out.append("mutation score**, i.e. four of five deliberate defects in it went undetected by")
    out.append("a fully green build.")
    out.append("")
    for r in results:
        if "confidence_assignment" not in r.path:
            continue
        module_equivalent = [m for m in r.survivors if m.equivalence_reason]
        module_adjusted_total = r.total - len(module_equivalent)
        module_adjusted = 100.0 * r.killed / module_adjusted_total if module_adjusted_total else 0.0
        out.append("| `confidence_assignment.py` | Then | Now |")
        out.append("| --- | ---: | ---: |")
        out.append("| Line coverage | 100% | 100% |")
        out.append(f"| Mutation score | 21% | {r.score:.1f}% raw, {module_adjusted:.1f}% adjusted |")
        out.append("")
    out.append("## Reproducing this")
    out.append("")
    out.append("```bash")
    out.append("make mutation")
    out.append("```")
    out.append("")
    out.append("The harness is `scripts/mutation_test.py`. It mutates one token at a time,")
    out.append("runs the module's own tests first and escalates anything that survives to")
    out.append("the full unit tier before recording it as a survivor, so the score is not")
    out.append("biased by the scoping.")
    out.append("")
    out.append("Before it mutates anything it runs those same pytest invocations against the")
    out.append("**unmutated** tree and aborts unless they pass, printing the failure. A mutant")
    out.append("is counted as killed whenever pytest exits non-zero, and pytest exits non-zero")
    out.append("for an unrelated failure, a collection error or a missing dependency just as")
    out.append("readily - so without that preflight a broken checkout scores 100% and")
    out.append("overwrites this page with the result.")
    out.append("")
    out.append('!!! danger "No child may load bytecode built from a different revision"')
    out.append("")
    out.append("    That is the invariant. CPython validates a cached `.pyc` against the")
    out.append("    source's `(mtime, size)` pair with one-second mtime granularity, so a")
    out.append("    mutant written in the same second as the file it replaces **and of the")
    out.append("    same byte length** (`==` to `!=`, `1` to `2`) is indistinguishable from")
    out.append("    it to the cache validator: the interpreter loads the stale `.pyc`, runs")
    out.append("    the **unmutated** code, every such mutant *survives* and the score is")
    out.append("    fiction. Two sweeps produced exactly that before it was found.")
    out.append("")
    out.append("    Two defences hold it, and both are required. `run_pytest()` passes")
    out.append("    `python -B` and sets `PYTHONDONTWRITEBYTECODE=1` in the child, so no")
    out.append("    `.pyc` is written during the run; and every `__pycache__` under")
    out.append("    `vntyper/` is deleted before the sweep starts and again after each")
    out.append("    mutant is written, so none left by an earlier run or an earlier mutant")
    out.append("    can be loaded. The flags stop new caches, the deletion stops old ones.")
    out.append("")
    out.append("    The `PYTHONDONTWRITEBYTECODE=1` on the `make mutation` recipe is defence")
    out.append("    in depth for the parent process, which never imports a target module -")
    out.append("    it is not what holds the invariant, and the harness is safe without it.")
    out.append("")
    out.append('!!! note "Each measurement runs in an isolated workspace"')
    out.append("")
    out.append("    The harness captures HEAD in a disposable detached worktree and overlays")
    out.append("    the current non-ignored working state, except selected mutation targets")
    out.append("    and requested output paths. Selected targets therefore come from the")
    out.append("    captured commit, while ordinary edits and new tests participate in the")
    out.append("    measurement without being written back.")
    out.append("")
    out.append("    Import provenance is proved against the pinned worktree before testing.")
    out.append("    A green baseline and a known-killed canary must then pass before ordinary")
    out.append("    mutants are measured, and the post-overlay baseline is verified after")
    out.append("    the canary and after every target.")
    out.append("")
    out.append("    Every mutant and bytecode-cache write is confined to that workspace;")
    out.append("    real production source is never mutated. Requested report artifacts are")
    out.append("    built completely and installed atomically in the real checkout.")
    out.append("    The cleanup is best effort: SIGINT, SIGTERM, SIGHUP and SIGQUIT attempt")
    out.append("    the common unwind path, while SIGKILL or a host crash can leave only an")
    out.append("    orphan disposable worktree for later inspection and removal.")
    out.append("")
    out.append("## Related: branch coverage, now enabled")
    out.append("")
    out.append("Mutation testing and branch coverage were investigated together, because both")
    out.append("ask a sharper question than line coverage and they agreed on which modules are")
    out.append("weakest. The branch-coverage half of that work is recorded here so it is not")
    out.append("re-derived from scratch.")
    out.append("")
    out.append("`[tool.coverage.run]` now sets `branch = true`, so an `if` that is entered but")
    out.append("never taken no longer counts as fully covered. It was enabled in **#196**,")
    out.append("measured on `fix/issue-181-197-followups` at `5bb2463`:")
    out.append("")
    out.append("| Measure | Value |")
    out.append("| --- | ---: |")
    out.append("| Line (statement) coverage | 76.60% |")
    out.append("| **Branch-inclusive total** | **74.22%** |")
    out.append("| Branch-only coverage | 66.00% |")
    out.append("| Branch exits never taken | 512 of 1506 |")
    out.append("")
    out.append("`fail_under` was raised **70 &rarr; 74** in the same commit, to the figure")
    out.append("`scripts/coverage_gate.py` printed for that run. **The floor was raised to meet")
    out.append("the measurement; the measurement was not weakened to fit the floor.** That")
    out.append("distinction is the whole point of the ratchet.")
    out.append("")
    out.append('!!! warning "74 is a branch-inclusive floor, and nothing else notices if that changes"')
    out.append("")
    out.append("    Deleting `branch = true` does not fail any gate on its own - it *raises*")
    out.append("    the reported total, because statement-only coverage of the same suite is")
    out.append("    76.60% against the branch-inclusive 74.22%. The build would go green while")
    out.append("    measuring strictly less. `tests/unit/test_coverage_gate.py::test_branch_coverage_is_enabled`")
    out.append("    exists solely to make that edit fail, and it is the only thing that does.")
    out.append("")
    out.append("### Correction: the previously recorded prerequisite was wrong")
    out.append("")
    out.append("This section formerly recorded the opposite conclusion, and it is kept here")
    out.append("rather than deleted, because a document that quietly rewrites its own history")
    out.append("stops being worth trusting.")
    out.append("")
    out.append("The earlier measurement was **63.80%** branch-inclusive (76.60% is the current")
    out.append("line figure; it was 66.82% then) against a `fail_under` of **66**. Enabling")
    out.append("branch coverage at that point really would have failed CI on the enabling")
    out.append("commit, and the decision not to enable it - and specifically not to lower the")
    out.append("floor to admit it - was correct.")
    out.append("")
    out.append("What was wrong was the stated route out. The old text identified")
    out.append("`cohort_summary.py` and `install_references.py` as holding 275 of the 685")
    out.append("untaken exits, both on the oversized-file list in `AGENTS.md`, and concluded:")
    out.append('*"Splitting them is the prerequisite, not writing more tests against them."*')
    out.append("")
    out.append("**It was not a prerequisite.** Branch coverage cleared the floor with both")
    out.append("files still unsplit and still untested. The gap was closed instead by testing")
    out.append("five small, already-testable modules to 100% and a sixth to 98%:")
    out.append("`cross_match.py`, `utils.py`, `file_processing.py`,")
    out.append("`extract_unmapped_from_offset.py`, `variant_parsing.py`, and `docker/app/tasks.py`.")
    out.append("")
    out.append('The generalisable mistake was reading "these two files hold the most untaken')
    out.append('exits" as "these two files are the ones that must be fixed". Concentration of')
    out.append("missing coverage is not the same as cheapness of covering it: the two oversized")
    out.append("modules are expensive precisely because they fuse I/O with logic, while the")
    out.append("same number of units was available across several modules that could be called")
    out.append("directly. Splitting `cohort_summary.py` and `install_references.py` remains")
    out.append("worth doing for the reasons `AGENTS.md` gives - it was simply never a blocker")
    out.append("for this.")
    out.append("")
    out.append("## Raw output")
    out.append("")
    out.append("```text")
    out.append(format_report(results, elapsed))
    out.append("```")
    out.append("")
    return "\n".join(out)


def results_to_dict(results: list[ModuleResult], elapsed: float) -> dict:
    """
    Serialise a completed sweep so the page can be re-rendered without re-measuring.

    A sweep costs 15-30 minutes, but classifying a survivor as equivalent only changes
    how the results are *presented*. Persisting them means adding a classification is a
    second-long re-render rather than a reason to re-run the whole thing - and, more
    importantly, the published page stays a rendering of a real measurement instead of
    something hand-edited afterwards.

    ``Mutant.source`` (the whole mutated file) is deliberately dropped: it is megabytes
    and nothing in the report uses it.

    Args:
        results (list[ModuleResult]): One entry per module swept.
        elapsed (float): Wall-clock duration of the sweep in seconds.

    Returns:
        dict: A JSON-serialisable representation of the sweep.
    """
    return {
        "elapsed": elapsed,
        "modules": [
            {
                "path": r.path,
                "killed": r.killed,
                "survived": r.survived,
                "survivors": [
                    {"line": m.line, "original": m.original, "replacement": m.replacement} for m in r.survivors
                ],
            }
            for r in results
        ],
    }


def results_from_dict(payload: dict) -> tuple[list[ModuleResult], float]:
    """
    Rebuild sweep results from :func:`results_to_dict` output.

    Args:
        payload (dict): A previously serialised sweep.

    Returns:
        tuple[list[ModuleResult], float]: The results and the original elapsed time.
    """
    results = [
        ModuleResult(
            path=module["path"],
            killed=module["killed"],
            survived=module["survived"],
            survivors=[
                Mutant(
                    path=module["path"],
                    line=s["line"],
                    original=s["original"],
                    replacement=s["replacement"],
                    source="",
                    killed=False,
                )
                for s in module["survivors"]
            ],
        )
        for module in payload["modules"]
    ]
    return results, float(payload["elapsed"])


def write_outputs(results: list[ModuleResult], elapsed: float, output: Path | None, results_json: Path | None) -> None:
    """
    Write the report and, optionally, the machine-readable results.

    Args:
        results (list[ModuleResult]): One entry per module swept.
        elapsed (float): Wall-clock duration of the sweep in seconds.
        output (Path | None): Report destination. ``.md`` renders the docs page, any
            other suffix renders the plain-text report.
        results_json (Path | None): Where to persist the raw results, if anywhere.
    """
    validate_disjoint_paths((("report output", output), ("results JSON", results_json)))

    rendered_output = None
    if output:
        rendered_output = (
            format_markdown(results, elapsed) if output.suffix == ".md" else format_report(results, elapsed)
        )
    rendered_json = None
    if results_json:
        rendered_json = json.dumps(results_to_dict(results, elapsed), indent=2) + "\n"

    if results_json and rendered_json is not None:
        results_json.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_text(results_json, rendered_json)
        print(f"Raw results written to {results_json}")
    if output and rendered_output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_text(output, rendered_output)
        print(f"Report written to {output}")


def _refuse_if_dirty(targets: Iterable[str], outputs: Sequence[Path | None]) -> str | None:
    """
    Decide whether this run may write anything at all.

    Args:
        targets (Iterable[str]): Repo-relative sources the run will mutate. Empty on the
            ``--render-only`` path, which mutates nothing.
        outputs (Sequence[Path | None]): Files the run will overwrite, ``None`` for the
            ones it was not asked to write.

    Returns:
        str | None: The refusal to print, or ``None`` when the run may proceed. Both
            refusals are returned rather than raised so the caller keeps its single exit
            point, and both fail closed: a dirty file refuses, and so does an
            indeterminate answer from git.
    """
    guarded = writable_paths(REAL_REPO_ROOT, targets, outputs)
    try:
        dirty = dirty_paths(REAL_REPO_ROOT, guarded)
    except RuntimeError as exc:
        return str(exc)
    return format_dirty_tree_refusal(dirty) if dirty else None


def _real_target_digests(real_root: Path, targets: Mapping[str, object]) -> dict[str, str]:
    """Capture SHA-256 digests for exactly the selected real-checkout targets.

    Args:
        real_root: Immutable real repository root.
        targets: Selected repository-relative mutation targets.

    Returns:
        A target-to-digest mapping captured before workspace creation.

    Raises:
        OSError: If a selected target cannot be read.
        RuntimeError: If a selected target is not a regular file.
    """
    return {relative: hashlib.sha256(read_source_bytes(real_root, relative)).hexdigest() for relative in targets}


def _verify_real_target_digests(real_root: Path, expected: Mapping[str, str]) -> None:
    """Require selected real-checkout targets to retain their startup bytes.

    Args:
        real_root: Immutable real repository root.
        expected: Target digests captured before workspace creation.

    Raises:
        OSError: If a selected target cannot be read.
        RuntimeError: If a target is not regular or its bytes changed.
    """
    for relative, expected_digest in expected.items():
        current_digest = hashlib.sha256(read_source_bytes(real_root, relative)).hexdigest()
        if current_digest != expected_digest:
            raise RuntimeError(f"real mutation target changed during sweep: {relative}")


def _terminate(signum: int, _frame: FrameType | None) -> NoReturn:
    """Convert a catchable termination signal into the common unwind path."""
    raise KeyboardInterrupt(f"terminated by signal {signum}")


def _install_signal_handlers() -> tuple[int, ...]:
    """Install cleanup-aware handlers without replacing Python's SIGINT behavior."""
    registered: list[int] = []
    for name in ("SIGTERM", "SIGHUP", "SIGQUIT"):
        signum = getattr(signal, name, None)
        if signum is not None:
            signal.signal(signum, _terminate)
            registered.append(signum)
    return tuple(registered)


def main() -> int:
    """
    Run the mutation sweep and print the report.

    Returns:
        int: 0 after measurement, output installation, workspace cleanup, and real-source
            digest verification all complete.
            **Non-zero means no usable measurement was produced**, never "the score is
            too low": the refusals are a dirty working tree, a red baseline, a git that
            cannot answer, a disposable baseline that cannot be restored, failed cleanup,
            or a changed real target. The score itself is advisory and is not gated - see
            the module docstring for why gating on it would be wrong.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--module", help="Only mutate targets whose path contains this substring.")
    parser.add_argument("--output", type=Path, help="Also write the report to this file.")
    parser.add_argument("--timeout", type=int, default=600, help="Per-pytest-run timeout in seconds.")
    parser.add_argument("--verbose", action="store_true", help="Print every mutant's outcome.")
    parser.add_argument("--results-json", type=Path, help="Persist the raw results for later re-rendering.")
    parser.add_argument(
        "--render-only",
        type=Path,
        help="Re-render the report from a previous --results-json instead of sweeping. "
        "Use after adding an EQUIVALENT_MUTANTS entry: the measurement is unchanged, "
        "only its presentation.",
    )
    args = parser.parse_args()
    output = None if args.output is None else (REAL_REPO_ROOT / args.output).resolve(strict=False)
    results_json = None if args.results_json is None else (REAL_REPO_ROOT / args.results_json).resolve(strict=False)
    render_input = None if args.render_only is None else (REAL_REPO_ROOT / args.render_only).resolve(strict=False)

    # Re-render path: no mutation, no tests, no source rewritten - but `write_outputs()`
    # still overwrites the docs page wholesale, so the preflight applies here too. It is
    # narrower: only the outputs are at risk, because no target is touched.
    if render_input is not None:
        try:
            validate_disjoint_paths((("render input", render_input), ("report output", output)))
        except ValueError as exc:
            print(f"REFUSING TO WRITE: {exc}", file=sys.stderr)
            return 1
        refusal = _refuse_if_dirty([], [output])
        if refusal is not None:
            print(refusal, file=sys.stderr)
            return 1
        render_results, elapsed = results_from_dict(json.loads(render_input.read_text(encoding="utf-8")))
        print(format_report(render_results, elapsed))
        write_outputs(render_results, elapsed, output, None)
        return 0

    targets = {p: t for p, t in TARGETS.items() if not args.module or args.module in p}
    if not targets:
        print(f"No target matches {args.module!r}. Known targets:", file=sys.stderr)
        for path_str in TARGETS:
            print(f"  {path_str}", file=sys.stderr)
        return 0

    try:
        validate_disjoint_paths(
            (("report output", output), ("results JSON", results_json)),
            tuple((path_str, REAL_REPO_ROOT / path_str) for path_str in targets),
        )
    except ValueError as exc:
        print(f"REFUSING TO WRITE: {exc}", file=sys.stderr)
        return 1

    # Selected targets must identify committed baseline bytes; requested outputs must be
    # safe to replace. `--module` narrows the target side exactly as far as the sweep.
    refusal = _refuse_if_dirty(targets, [output, results_json])
    if refusal is not None:
        print(refusal, file=sys.stderr)
        return 1

    status = 0
    real_digests: dict[str, str] | None = None
    try:
        real_digests = _real_target_digests(REAL_REPO_ROOT, targets)
        _install_signal_handlers()
        excluded_outputs = tuple(path for path in (output, results_json) if path is not None)
        with detached_head_workspace(REAL_REPO_ROOT, tuple(targets), excluded_outputs) as workspace:
            print(f"Captured HEAD: {workspace.head}")
            print(f"Disposable worktree: {workspace.sweep_root}")
            verify_import_provenance(workspace, tuple(targets))
            execution_root = workspace.root_capability or workspace.sweep_root

            with using_root_capability(execution_root) as capability:
                removed = clear_bytecode_caches(git_capability_path(capability) / "vntyper")
            print(f"Cleared {removed} __pycache__ directories before starting.")

            print("Checking the baseline: the tests must pass on the UNMUTATED tree.")
            refusal = baseline_refusal(targets, args.timeout, repo_root=execution_root)
            if refusal is not None:
                raise RuntimeError(refusal)

            refusal = canary_refusal(args.timeout, repo_root=execution_root)
            if refusal is not None:
                raise RuntimeError(refusal)
            print("Known-killed canary: killed (scoped)")
            workspace.verify_baseline()

            start = time.monotonic()
            results: list[ModuleResult] = []
            for path_str, scoped_tests in targets.items():
                results.append(
                    sweep_module(
                        path_str,
                        scoped_tests,
                        args.timeout,
                        args.verbose,
                        repo_root=execution_root,
                    )
                )
                workspace.verify_baseline()
            elapsed = time.monotonic() - start

            print()
            print(format_report(results, elapsed))
            write_outputs(results, elapsed, output, results_json)
    except (RuntimeError, OSError, KeyboardInterrupt) as exc:
        print(f"\n{exc}", file=sys.stderr)
        status = 1
    finally:
        if real_digests is not None:
            try:
                _verify_real_target_digests(REAL_REPO_ROOT, real_digests)
            except (RuntimeError, OSError, KeyboardInterrupt) as exc:
                print(f"\n{exc}", file=sys.stderr)
                status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
