#!/usr/bin/env python3
"""Reproduce the golden-cohort gate (#179): one entry point, one instrument.

``docs/development/golden-cohort-gate.md`` records a series of before-versus-after runs
over the whole local test cohort. The scripts that produced the first three were never
committed, so the instrument backing every genotype claim on this project had to be
reconstructed from prose each time. This is that instrument. (The page grows a section per
run, so nothing here counts them - a docstring that names a total is wrong on the day the
next run lands, and two in this package were.)

What it does that the first three recorded runs did:

* derives the 50 base cases from what is in ``tests/data`` rather than from a list, and
  self-checks the whole matrix against the counts the page records. The five non-fast ids,
  the three adVNTR ids and the three probes are **not** derived: they are declared policy,
  resolved against the derived set. See :mod:`golden_cohort.matrix`;
* launches **every** run through a wrapper that prints its resolved ``vntyper.__file__``
  and marker-module state as its first line and refuses to start unless both agree with
  its side, because the editable finder otherwise makes a baseline checkout run candidate
  code;
* compares the Kestrel row set, the pre-filter frame, the adVNTR frame, the coverage
  summary, the screening sentence and its computed emphasis, the recorded pipeline steps,
  the executed shell commands and the exit code.

What it adds: **cohort mode**, uncovered by runs 1-3 and refactored across five new modules
by ``1c6c9d6``.

What it refuses, all of which it used to permit (see :mod:`golden_cohort.admissibility`):

* a case that did not exit as its matrix entry declared, or that exited zero and wrote none
  of the artefacts it must produce - without this, two sides that both die before writing
  a genotype artefact compare *equal* and the gate exits 0;
* a comparison of a run root against itself, of two identically-labelled sides, of two
  identical trees or commits, or of two sides that expected the same marker state;
* an unfiltered attestation run whose derived matrix deviates from the documented contract,
  and any matrix with no cases at all.

What it records that it used to assert: the ``HEAD`` each side ran and whether its
genotype-bearing paths were dirty, into ``side.json``. ``--expect-before-sha`` /
``--expect-after-sha`` turn that record into a check.

Typical use (the baseline worktree is created by whoever runs the gate, not by this
script - it takes the path):

    python scripts/golden_cohort_gate.py probe --tree /path/to/before \\
        --marker vntyper.scripts.pipeline_guards
    python scripts/golden_cohort_gate.py run --side before --tree /path/to/before \\
        --run-root /scratch/gate/before --marker vntyper.scripts.pipeline_guards \\
        --expect-marker absent
    python scripts/golden_cohort_gate.py run --side after --tree "$PWD" \\
        --run-root /scratch/gate/after --marker vntyper.scripts.pipeline_guards \\
        --expect-marker present
    python scripts/golden_cohort_gate.py compare --before-root /scratch/gate/before \\
        --after-root /scratch/gate/after --json /scratch/gate/result.json \\
        --text /scratch/gate/result.md --expect-after-sha ec67fff

Exit codes:
    0: the command succeeded. For ``compare`` that means the verdict is ``IDENTICAL`` and
       nothing weaker - a clean run over a reduced matrix is ``REDUCED`` and exits 1, on
       purpose, so that `gate compare && echo PASS` cannot print PASS for a run that
       measured less than the attestation matrix.
    1: a delta, an unmet case expectation, an unverified run, a blocked cohort case, a
       reduced matrix, a revision mismatch, two sides that are not opposed, or a usage
       error.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from golden_cohort import (  # noqa: E402
    HARNESS_VERSION,
    admissibility,
    compare,
    launcher,
    matrix,
    normalise,
    runner,
)

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _split_launch_argv(argv: list[str]) -> tuple[list[str], list[str]]:
    """Split the harness's own arguments from the ``vntyper`` ones at the first ``--``.

    ``argparse.REMAINDER`` mis-handles a trailing argument list that itself starts with a
    flag, which every ``vntyper`` invocation does.

    Args:
        argv: The process arguments, without the program name.

    Returns:
        tuple[list[str], list[str]]: The harness's arguments and the ``vntyper`` ones.
    """
    if "--" not in argv:
        return argv, []
    index = argv.index("--")
    return argv[:index], argv[index + 1 :]


def build_parser() -> argparse.ArgumentParser:
    """Build the harness's argument parser.

    Returns:
        argparse.ArgumentParser: The parser, with ``matrix``, ``probe``, ``launch``,
        ``run`` and ``compare`` registered under ``dest="command"``.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"golden-cohort-gate {HARNESS_VERSION}")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_matrix_options(target: argparse.ArgumentParser) -> None:
        target.add_argument("--data-dir", type=Path, default=REPO_ROOT / "tests" / "data")
        target.add_argument("--case", action="append", default=[], help="Keep only cases whose id contains this.")
        target.add_argument("--no-probes", action="store_true")
        target.add_argument("--no-cohort", action="store_true")
        target.add_argument("--advntr-max-coverage", type=int, default=300)
        target.add_argument(
            "--non-fast-cases",
            default=None,
            help="Comma-separated base case ids that repeat without --fast-mode. See matrix.py.",
        )
        target.add_argument("--advntr-cases", default=None, help="Comma-separated base case ids that run adVNTR.")
        target.add_argument(
            "--allow-matrix-drift",
            action="store_true",
            help=(
                "Build the matrix even though it deviates from the documented contract. The run is then "
                "recorded as not attestation-grade and a clean comparison of it reads REDUCED, not IDENTICAL."
            ),
        )

    parser_matrix = subparsers.add_parser("matrix", help="Derive and print the case matrix.")
    add_matrix_options(parser_matrix)
    parser_matrix.add_argument("--out", type=Path, default=None)

    parser_probe = subparsers.add_parser(
        "probe", help="Demonstrate that a bare import resolves to the editable install, not the CWD."
    )
    parser_probe.add_argument("--tree", type=Path, required=True)
    parser_probe.add_argument("--marker", required=True)

    parser_launch = subparsers.add_parser("launch", help="Internal: run one vntyper invocation under the wrapper.")
    parser_launch.add_argument("--tree", type=Path, required=True)
    parser_launch.add_argument("--side", required=True)
    parser_launch.add_argument("--marker", required=True)
    parser_launch.add_argument("--expect-marker", required=True, choices=["present", "absent"])
    parser_launch.add_argument("--commands-log", type=Path, default=None)

    parser_run = subparsers.add_parser("run", help="Run one side of the comparison.")
    add_matrix_options(parser_run)
    parser_run.add_argument("--side", required=True, choices=["before", "after"])
    parser_run.add_argument("--tree", type=Path, required=True)
    parser_run.add_argument("--run-root", type=Path, required=True)
    parser_run.add_argument("--marker", required=True)
    parser_run.add_argument("--expect-marker", required=True, choices=["present", "absent"])
    parser_run.add_argument("--threads", type=int, default=4)
    parser_run.add_argument("--advntr-threads", type=int, default=8)
    parser_run.add_argument("--jobs", type=int, default=1)
    parser_run.add_argument("--timeout", type=int, default=5400)

    parser_compare = subparsers.add_parser("compare", help="Compare two completed sides.")
    parser_compare.add_argument("--before-root", type=Path, required=True)
    parser_compare.add_argument("--after-root", type=Path, required=True)
    parser_compare.add_argument("--json", dest="json_out", type=Path, default=None)
    parser_compare.add_argument("--text", dest="text_out", type=Path, default=None)
    parser_compare.add_argument(
        "--expect-before-sha", default=None, help="Refuse unless the baseline side recorded this commit."
    )
    parser_compare.add_argument(
        "--expect-after-sha", default=None, help="Refuse unless the candidate side recorded this commit."
    )
    parser_compare.add_argument(
        "--require-clean",
        action="store_true",
        help="Refuse if either side's vntyper/, docker/ or scripts/ had uncommitted changes when it ran.",
    )

    return parser


def _matrix_from_args(args: argparse.Namespace) -> dict[str, Any]:
    """Build the matrix from the shared matrix options.

    Args:
        args: The parsed arguments.

    Returns:
        dict[str, Any]: The matrix.
    """
    non_fast = (
        tuple(part.strip() for part in args.non_fast_cases.split(",") if part.strip())
        if args.non_fast_cases
        else matrix.NON_FAST_CASE_IDS
    )
    advntr = (
        tuple(part.strip() for part in args.advntr_cases.split(",") if part.strip())
        if args.advntr_cases
        else matrix.ADVNTR_CASE_IDS
    )
    return matrix.build_matrix(
        args.data_dir.resolve(),
        non_fast_ids=non_fast,
        advntr_ids=advntr,
        advntr_max_coverage=args.advntr_max_coverage,
        case_filter=args.case or None,
        include_probes=not args.no_probes,
        include_cohort=not args.no_cohort,
        strict=not args.allow_matrix_drift,
    )


def cmd_matrix(args: argparse.Namespace) -> int:
    """Derive the matrix and write or print it.

    Args:
        args: The parsed arguments.

    Returns:
        int: 0 if the matrix could be built, 1 if it deviates from the documented contract
        or has no cases. This used to return 0 unconditionally, which meant a silently
        reduced matrix passed its own self-check.
    """
    try:
        built = _matrix_from_args(args)
    except ValueError as exc:
        logger.error(str(exc))
        return 1
    if not built.get("cases") and not built.get("probes"):
        logger.error("The matrix contains no pipeline cases or probes.")
        return 1
    text = json.dumps(built, indent=2)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
        logger.info(f"matrix written to {args.out}")
    else:
        print(text)
    return 0


def cmd_probe(args: argparse.Namespace) -> int:
    """Show what ``import vntyper`` resolves to, three ways, from inside one tree.

    This is the check the gate page describes as "the part that could have made the whole
    exercise worthless". The three modes are not decoration; they are three different
    ``sys.path[0]`` values and only one of them is the one the ``vntyper`` console script
    actually has:

    * ``dash_c`` - ``python -c`` puts the **current directory** on ``sys.path``, so a bare
      import from inside the tree finds the tree. This is the reassuring measurement, and
      it is not the one that matters.
    * ``script_outside_tree`` - a script file living somewhere else puts **its own
      directory** on ``sys.path``, and the current directory nowhere. That is exactly the
      console script's situation, and it is where the editable finder - appended to
      ``sys.meta_path`` and pointing at whichever worktree the editable install was made
      from - answers instead. A leak here means a baseline checkout would have executed
      candidate code.
    * ``pythonpath_pinned`` - what :mod:`golden_cohort.launcher` does, and what every gate
      run therefore does.

    Args:
        args: The parsed arguments.

    Returns:
        int: 0 if the pinned probe resolved into the tree, 1 otherwise.
    """
    tree = args.tree.resolve()
    snippet = (
        "import importlib.util, json, sys, vntyper\n"
        "print(json.dumps({'file': vntyper.__file__, 'sys_path_0': sys.path[0],"
        f" 'marker': importlib.util.find_spec({args.marker!r}) is not None}}))\n"
    )

    results: dict[str, Any] = {}
    with tempfile.TemporaryDirectory(prefix="gate_probe_") as tmp:
        script = Path(tmp) / "resolve_probe.py"
        script.write_text(snippet, encoding="utf-8")
        modes: tuple[tuple[str, list[str], bool], ...] = (
            ("dash_c", [sys.executable, "-c", snippet], False),
            ("script_outside_tree", [sys.executable, str(script)], False),
            ("pythonpath_pinned", [sys.executable, str(script)], True),
        )
        for label, command, pinned in modes:
            env = dict(os.environ)
            env.pop("PYTHONPATH", None)
            if pinned:
                env["PYTHONPATH"] = str(tree)
            completed = subprocess.run(command, cwd=str(tree), env=env, capture_output=True, text=True, check=False)
            try:
                results[label] = json.loads(completed.stdout.strip().splitlines()[-1])
            except (ValueError, IndexError):
                results[label] = {"error": completed.stderr.strip()[-500:]}

    for label, info in results.items():
        logger.info(f"probe {label}: {info}")

    def resolved_in_tree(label: str) -> bool:
        found = results[label].get("file")
        return isinstance(found, str) and str(tree) in found

    ok = resolved_in_tree("pythonpath_pinned")
    leaked = not resolved_in_tree("script_outside_tree")
    print(
        json.dumps(
            {"tree": str(tree), "marker": args.marker, "results": results, "unpinned_script_leaked": leaked}, indent=2
        )
    )
    if leaked:
        logger.warning(
            "An unpinned script run from inside this tree resolved `vntyper` somewhere else. That is the "
            "failure the wrapper exists to prevent, and it is why every gate run goes through `launch`."
        )
    else:
        logger.info(
            "No leak observed here, which does not make the wrapper optional: whether a bare import leaks "
            "depends on where the editable install points, which is a property of the machine, not of the run."
        )
    return 0 if ok else 1


def cmd_launch(args: argparse.Namespace, vntyper_argv: list[str]) -> int:
    """Run one ``vntyper`` invocation under the resolution wrapper.

    Args:
        args: The parsed arguments.
        vntyper_argv: Everything after ``--``.

    Returns:
        int: The wrapper's exit code.
    """
    return launcher.launch(
        tree=args.tree,
        side=args.side,
        marker=args.marker,
        expect_marker=args.expect_marker == "present",
        commands_log=args.commands_log,
        argv=vntyper_argv,
    )


def cmd_run(args: argparse.Namespace) -> int:
    """Run one side of the comparison.

    Args:
        args: The parsed arguments.

    Returns:
        int: 0 if every run verified its package resolution, every case did what the matrix
        declared it would, and no cohort case was refused. The matrix is built **before**
        anything launches, so a drifted or empty matrix fails here rather than after an
        hour of pipeline runs.
    """
    try:
        built = _matrix_from_args(args)
    except ValueError as exc:
        logger.error(str(exc))
        return 1
    record = runner.run_side(
        matrix=built,
        tree=args.tree.resolve(),
        run_root=args.run_root.resolve(),
        side=args.side,
        marker=args.marker,
        expect_marker=args.expect_marker == "present",
        threads=args.threads,
        advntr_threads=args.advntr_threads,
        jobs=args.jobs,
        timeout=args.timeout,
        skip_cohort=args.no_cohort,
    )
    blocked = [case_id for case_id, item in record["cohort_results"].items() if item.get("blocked")]
    if blocked:
        logger.error(f"cohort cases refused for missing inputs: {blocked}")
    return 0 if record["launch_verified"] and record["expectations_met"] and not blocked else 1


def cmd_compare(args: argparse.Namespace) -> int:
    """Compare two completed sides and write both result documents.

    Args:
        args: The parsed arguments.

    Returns:
        int: 0 when the verdict is ``IDENTICAL``, 1 otherwise. Two sides that are not
        genuinely opposed, or whose recorded revisions disagree with what the caller
        expected, are refused before any artefact is read - a comparison of a tree against
        itself is not a passing gate, it is the absence of one.
    """
    before_root = args.before_root.resolve()
    after_root = args.after_root.resolve()
    before_side = runner.load_side(before_root)
    after_side = runner.load_side(after_root)

    problems = admissibility.check_sides_opposed(before_root, after_root, before_side, after_side)
    problems.extend(admissibility.verify_revision(before_side.get("revision"), args.expect_before_sha, side="before"))
    problems.extend(admissibility.verify_revision(after_side.get("revision"), args.expect_after_sha, side="after"))
    if args.require_clean:
        for label, side in (("before", before_side), ("after", after_side)):
            revision = side.get("revision") or {}
            if revision.get("dirty_relevant"):
                problems.append(
                    f"{label}: {len(revision.get('dirty_paths', []))} uncommitted change(s) under "
                    f"{revision.get('revision_paths')} when it ran, so it attests no commit"
                )
    for label, side in (("before", before_side), ("after", after_side)):
        revision = side.get("revision") or {}
        if not revision.get("head"):
            logger.warning(
                f"{label}: no revision was recorded for this side, so this comparison names a path and not a "
                "commit. Re-run the side with a harness that records one before writing an attestation."
            )
        elif revision.get("dirty_relevant"):
            logger.warning(
                f"{label}: ran {revision['head']} with {len(revision.get('dirty_paths', []))} uncommitted "
                f"change(s) under {revision.get('revision_paths')}, so it attests that commit plus edits."
            )
    if problems:
        for problem in problems:
            logger.error(f"refusing to compare: {problem}")
        return 1

    # Both sides read the same `tests/data` by absolute path, exactly as run 2 and run 3
    # did, so the data root is shared and must normalise to one placeholder on both sides.
    data_roots: set[str] = set()
    for root in (before_root, after_root):
        declared = (compare.artifacts.read_json(root / "matrix.json") or {}).get("data_dir")
        if isinstance(declared, str):
            data_roots.add(declared)
    if len(data_roots) > 1:
        logger.error(f"The two sides read different test-data directories: {sorted(data_roots)}")
        return 1
    data_root = Path(next(iter(data_roots))) if data_roots else None

    before_rules = normalise.build_rules(
        source_root=Path(before_side["tree"]), run_root=before_root, data_root=data_root
    )
    after_rules = normalise.build_rules(source_root=Path(after_side["tree"]), run_root=after_root, data_root=data_root)

    result = compare.compare_sides(
        before_root,
        after_root,
        before_side,
        after_side,
        normalise.manifest(after_rules),
        before_rules,
        after_rules,
    )
    text = compare.render_text(result)

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        logger.info(f"machine-readable result written to {args.json_out}")
    if args.text_out:
        args.text_out.parent.mkdir(parents=True, exist_ok=True)
        args.text_out.write_text(text, encoding="utf-8")
        logger.info(f"human summary written to {args.text_out}")
    print(text)
    return 0 if result["verdict"] == "IDENTICAL" else 1


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and dispatch.

    Args:
        argv: The arguments, without the program name. Defaults to ``sys.argv[1:]``.

    Returns:
        int: The process exit code.
    """
    raw = list(sys.argv[1:] if argv is None else argv)
    own, vntyper_argv = _split_launch_argv(raw)
    args = build_parser().parse_args(own)

    if args.command == "launch":
        # Deliberately no logging configuration here: `vntyper.cli.main` is the sole place
        # the pipeline's logging is set up, and configuring the root logger first would put
        # the harness between a run and its own log file.
        return cmd_launch(args, vntyper_argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    if args.command == "matrix":
        return cmd_matrix(args)
    if args.command == "probe":
        return cmd_probe(args)
    if args.command == "run":
        return cmd_run(args)
    return cmd_compare(args)


if __name__ == "__main__":
    sys.exit(main())
