"""Execute one side of the comparison.

One side is one source tree plus one expectation about the marker module. Every case is a
separate process, launched through :mod:`golden_cohort.launcher` so that which package it
resolved is recorded before it does any work, and every process writes into its own
directory under the side's run root.

Ordering is explicit and is not an implementation detail: ``vntyper cohort`` reads
``pipeline_summary.json`` out of directories the per-sample runs create, so the per-sample
cases run first and a cohort case whose inputs are missing is **refused** rather than run.
A cohort over nothing produces a report over nothing and compares clean, which is the worst
possible failure for a gate.

Every case is also judged against what its matrix entry declared - its ``expect_exit`` and,
for a case expected to exit zero, the artefacts it must have written. See
:mod:`golden_cohort.admissibility` for why that is not optional: without it, two sides that
both die before writing anything compare *equal*.

The side record names the revision each side ran, not merely the path. ``side.json``
carries ``revision`` with the tree's ``HEAD``, its branch and whether the genotype-bearing
paths were dirty, so "this run attests commit X" is something the instrument recorded
rather than something its operator asserted afterwards.
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from golden_cohort import HARNESS_VERSION, admissibility, launcher
from golden_cohort.artifacts import read_json

logger = logging.getLogger(__name__)

#: Where the harness's own entry point lives, so a launched process can re-enter it.
GATE_SCRIPT = Path(__file__).resolve().parents[1] / "golden_cohort_gate.py"


def pipeline_argv(case: dict[str, Any], output_dir: Path, *, threads: int, advntr_threads: int) -> list[str]:
    """Build the ``vntyper pipeline`` argument list for one case.

    Args:
        case: The case from the matrix.
        output_dir: Where the run writes.
        threads: ``--threads`` for an ordinary case.
        advntr_threads: ``--threads`` for an adVNTR case, which the page ran wider.

    Returns:
        list[str]: The argument list, without the program name.
    """
    argv = [
        "pipeline",
        "--bam",
        case["bam"],
        "--reference-assembly",
        case["assembly"],
        "--output-dir",
        str(output_dir),
        "--threads",
        str(advntr_threads if case.get("advntr") else threads),
    ]
    if case.get("fast_mode", True):
        argv.append("--fast-mode")
    if case.get("advntr"):
        argv.extend(["--extra-modules", "advntr"])
        if case.get("advntr_max_coverage"):
            argv.extend(["--advntr-max-coverage", str(case["advntr_max_coverage"])])
    return argv


def cohort_argv(case: dict[str, Any], output_dir: Path, input_dirs: list[Path]) -> list[str]:
    """Build the ``vntyper cohort`` argument list for one cohort case.

    The flags are those ``vntyper/scripts/cli_parser.py`` actually registers: one of
    ``-i/--input-dirs`` or ``--input-file`` is required, ``-o/--output-dir`` is required,
    and ``--summary-file``, ``--summary-formats`` and ``--pseudonymize-samples`` are
    optional. ``--pseudonymize-samples`` takes an optional basename and defaults to
    ``sample_`` when given bare.

    Args:
        case: The cohort case from the matrix.
        output_dir: Where the cohort report is written.
        input_dirs: The per-sample directories to aggregate.

    Returns:
        list[str]: The argument list, without the program name.
    """
    argv = ["cohort", "--output-dir", str(output_dir)]
    argv.append("--input-dirs")
    argv.extend(str(path) for path in input_dirs)
    if case.get("summary_formats"):
        argv.extend(["--summary-formats", case["summary_formats"]])
    if case.get("pseudonymize"):
        argv.extend(["--pseudonymize-samples", case["pseudonymize"]])
    return argv


def _run_one(
    *,
    case: dict[str, Any],
    argv: list[str],
    tree: Path,
    side: str,
    marker: str,
    expect_marker: bool,
    output_dir: Path,
    log_dir: Path,
    timeout: int,
) -> dict[str, Any]:
    """Launch one case through the wrapper and record everything about the launch.

    Args:
        case: The matrix entry, which supplies the case id and the expectations.
        argv: The ``vntyper`` argument list.
        tree: The side's source tree.
        side: ``before`` or ``after``.
        marker: The marker module name.
        expect_marker: Whether the marker must be present on this side.
        output_dir: Where the case writes, so its required artefacts can be checked.
        log_dir: Where to write the logs, the command record and ``result.json``.
        timeout: Seconds before the case is killed.

    Returns:
        dict[str, Any]: The record also written to ``<log_dir>/result.json``, including
        the ``expectation_met`` / ``expectation_problems`` judgement.
    """
    case_id = case["case_id"]
    log_dir.mkdir(parents=True, exist_ok=True)
    commands_log = log_dir / "commands.jsonl"
    command = [
        sys.executable,
        str(GATE_SCRIPT),
        "launch",
        "--tree",
        str(tree),
        "--side",
        side,
        "--marker",
        marker,
        "--expect-marker",
        "present" if expect_marker else "absent",
        "--commands-log",
        str(commands_log),
        "--",
        *argv,
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(tree), *([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])])

    started = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            command,
            cwd=str(tree),
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        stdout, stderr, exit_code = completed.stdout, completed.stderr, completed.returncode
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        exit_code = 99

    (log_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (log_dir / "stderr.log").write_text(stderr, encoding="utf-8")

    launch_line = ""
    for line in stdout.splitlines():
        if line.startswith(launcher.LAUNCH_PREFIX):
            launch_line = line
            break
    aborted = any(line.startswith(launcher.ABORT_PREFIX) for line in stdout.splitlines())

    record = {
        "case_id": case_id,
        "argv": argv,
        "exit_code": exit_code,
        "launch_line": launch_line,
        "launch_verified": bool(launch_line) and not aborted and exit_code != launcher.EXIT_ABORT,
        "aborted": aborted,
        "timed_out": timed_out,
        "seconds": round(time.monotonic() - started, 2),
    }
    record.update(admissibility.check_case(case, record, output_dir))
    (log_dir / "result.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    status = "ABORTED" if aborted else f"exit {exit_code}"
    logger.info(f"[{side}] {case_id}: {status} in {record['seconds']}s")
    for problem in record["expectation_problems"]:
        logger.error(f"[{side}] {case_id} did not do what the matrix declared: {problem}")
    return record


def run_side(
    *,
    matrix: dict[str, Any],
    tree: Path,
    run_root: Path,
    side: str,
    marker: str,
    expect_marker: bool,
    threads: int,
    advntr_threads: int,
    jobs: int,
    timeout: int,
    skip_cohort: bool = False,
) -> dict[str, Any]:
    """Run every case of the matrix on one side.

    Args:
        matrix: The matrix from :func:`golden_cohort.matrix.build_matrix`.
        tree: The side's source tree.
        run_root: Where this side's outputs and logs go.
        side: ``before`` or ``after``.
        marker: The marker module name.
        expect_marker: Whether the marker must be present on this side.
        threads: ``--threads`` for ordinary cases.
        advntr_threads: ``--threads`` for adVNTR cases.
        jobs: How many per-sample cases to run concurrently.
        timeout: Per-case timeout in seconds.
        skip_cohort: Run only the per-sample cases.

    Returns:
        dict[str, Any]: The side record, also written to ``<run_root>/side.json``.
    """
    run_root.mkdir(parents=True, exist_ok=True)
    cases_root = run_root / "cases"
    logs_root = run_root / "logs"
    (run_root / "matrix.json").write_text(json.dumps(matrix, indent=2), encoding="utf-8")

    pipeline_cases = [*matrix["cases"], *matrix["probes"]]
    results: dict[str, Any] = {}

    def launch_case(case: dict[str, Any]) -> tuple[str, dict[str, Any]]:
        output_dir = cases_root / case["case_id"]
        argv = pipeline_argv(case, output_dir, threads=threads, advntr_threads=advntr_threads)
        record = _run_one(
            case=case,
            argv=argv,
            tree=tree,
            side=side,
            marker=marker,
            expect_marker=expect_marker,
            output_dir=output_dir,
            log_dir=logs_root / case["case_id"],
            timeout=timeout,
        )
        return case["case_id"], record

    logger.info(f"[{side}] running {len(pipeline_cases)} per-sample cases with {jobs} concurrent job(s)")
    if jobs <= 1:
        for case in pipeline_cases:
            case_id, record = launch_case(case)
            results[case_id] = record
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            results = dict(pool.map(launch_case, pipeline_cases))

    cohort_results: dict[str, Any] = {}
    if not skip_cohort:
        cohort_results = _run_cohorts(
            matrix=matrix,
            cases_root=cases_root,
            run_root=run_root,
            logs_root=logs_root,
            tree=tree,
            side=side,
            marker=marker,
            expect_marker=expect_marker,
            timeout=timeout,
        )

    # A refused cohort case never started, so it has no resolution to verify. Counting it as
    # unverified would report the wrong failure - the failure is the missing inputs, which
    # is reported on its own.
    launched = {
        case_id: item for case_id, item in {**results, **cohort_results}.items() if not item.get("blocked", False)
    }
    # `all()` over an empty mapping is True, so a side that launched nothing used to report
    # itself verified. It is not verified; it is unmeasured. `build_matrix` refuses a
    # zero-case matrix, and this is the second lock on the same door.
    launch_verified = bool(launched) and all(item["launch_verified"] for item in launched.values())
    unmet = sorted(case_id for case_id, item in launched.items() if not item.get("expectation_met", True))
    revision = admissibility.describe_tree(tree)

    record = {
        "harness_version": HARNESS_VERSION,
        "side": side,
        "tree": str(tree.resolve()),
        "revision": revision,
        "run_root": str(run_root.resolve()),
        "marker": marker,
        "expect_marker": "present" if expect_marker else "absent",
        "threads": threads,
        "advntr_threads": advntr_threads,
        "jobs": jobs,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "pipeline_results": results,
        "cohort_results": cohort_results,
        "launch_verified": launch_verified,
        "cases_launched": len(launched),
        "expectations_unmet": unmet,
        "expectations_met": not unmet,
    }
    (run_root / "side.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    unverified = [case_id for case_id, item in launched.items() if not item["launch_verified"]]
    if not launched:
        logger.error(f"[{side}] launched no cases at all, so nothing about this side is verified")
    elif unverified:
        logger.error(f"[{side}] {len(unverified)} case(s) did not verify their package resolution: {unverified}")
    else:
        logger.info(f"[{side}] all {len(launched)} runs verified their package resolution")

    if unmet:
        logger.error(f"[{side}] {len(unmet)} case(s) did not do what the matrix declared they would: {unmet}")
    if revision.get("head"):
        dirty = " (DIRTY)" if revision.get("dirty_relevant") else ""
        logger.info(f"[{side}] ran {revision['head']} on {revision.get('branch')}{dirty}")
    else:
        logger.warning(f"[{side}] could not record a revision for {tree}: {revision.get('error')}")
    return record


def _run_cohorts(
    *,
    matrix: dict[str, Any],
    cases_root: Path,
    run_root: Path,
    logs_root: Path,
    tree: Path,
    side: str,
    marker: str,
    expect_marker: bool,
    timeout: int,
) -> dict[str, Any]:
    """Run the cohort cases, refusing any whose per-sample inputs are not on disk.

    Args:
        matrix: The matrix.
        cases_root: Where the per-sample runs wrote.
        run_root: The side's run root.
        logs_root: Where per-case logs go.
        tree: The side's source tree.
        side: ``before`` or ``after``.
        marker: The marker module name.
        expect_marker: Whether the marker must be present on this side.
        timeout: Per-case timeout in seconds.

    Returns:
        dict[str, Any]: One record per cohort case. A refused case carries
        ``"blocked": True`` and the inputs that were missing.
    """
    cohorts_root = run_root / "cohorts"
    results: dict[str, Any] = {}

    for case in matrix.get("cohort_cases", []):
        output_dir = cohorts_root / case["case_id"]
        output_dir.mkdir(parents=True, exist_ok=True)

        reason = ""
        if case.get("empty_input_dir"):
            empty_dir = run_root / "empty_cohort_input"
            empty_dir.mkdir(parents=True, exist_ok=True)
            input_dirs = [empty_dir]
            missing: list[str] = []
        else:
            input_dirs = [cases_root / case_id for case_id in case["inputs"]]
            missing = [str(path) for path in input_dirs if not (path / "pipeline_summary.json").is_file()]
            if not input_dirs:
                # `--input-dirs` takes nargs="+", so an empty list is not a small cohort -
                # it is an argparse usage error dressed up as one. Refuse it here instead,
                # where the reason ("the matrix produced no per-sample cases") survives.
                missing = ["<none>"]
                reason = "the matrix produced no per-sample cases for this cohort to aggregate"
            elif missing:
                reason = (
                    f"{len(missing)} of {len(input_dirs)} input directories have no pipeline_summary.json "
                    f"({', '.join(missing[:3])})"
                )

        if missing and not case.get("allow_missing_inputs"):
            logger.error(
                f"[{side}] cohort case {case['case_id']} refused: {reason}. "
                "Running it anyway would aggregate fewer samples than intended and compare clean."
            )
            results[case["case_id"]] = {
                "case_id": case["case_id"],
                "blocked": True,
                "blocked_reason": reason,
                "missing_inputs": missing,
                "exit_code": None,
                "launch_line": "",
                "launch_verified": False,
                "aborted": False,
                "timed_out": False,
                "seconds": 0.0,
                # A refused case has no exit code and wrote nothing, so it has no
                # expectation to meet or miss. It is reported as blocked, which is a
                # stronger finding than an unmet expectation and is never silent.
                "expectation_met": True,
                "expectation_problems": [],
                "missing_artifacts": [],
            }
            continue

        argv = cohort_argv(case, output_dir, input_dirs)
        record = _run_one(
            case=case,
            argv=argv,
            tree=tree,
            side=side,
            marker=marker,
            expect_marker=expect_marker,
            output_dir=output_dir,
            log_dir=logs_root / case["case_id"],
            timeout=timeout,
        )
        record["blocked"] = False
        record["missing_inputs"] = missing
        record["input_count"] = len(input_dirs)
        results[case["case_id"]] = record

    return results


def load_side(run_root: Path) -> dict[str, Any]:
    """Read a side's ``side.json`` back.

    Args:
        run_root: The side's run root.

    Returns:
        dict[str, Any]: The side record.

    Raises:
        ValueError: If the side record is not there, which means that side never ran.
    """
    record = read_json(run_root / "side.json")
    if record is None:
        msg = f"No side.json under {run_root}. That side has not been run, so there is nothing to compare."
        logger.error(msg)
        raise ValueError(msg)
    return record
