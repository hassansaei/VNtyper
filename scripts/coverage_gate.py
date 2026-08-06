#!/usr/bin/env python3
"""
Report unit-test coverage against a hard floor and a soft target.

Two thresholds, deliberately different:

* **Hard floor** (`fail_under` in ``pyproject.toml``) - CI fails below it. It exists to
  stop regressions, and it is a ratchet: it only ever moves up.
* **Soft target** (``--target``, 80% by default) - what the project is aiming for.
  Falling short is a warning, never a failure.

The floor is what makes the target reachable. Every time coverage climbs meaningfully
above the floor, this script prints the exact edit to raise it, so the floor follows
coverage upward instead of sitting at whatever it was set to years ago.

Usage:
    python scripts/coverage_gate.py [--target 80] [--ratchet-margin 0.25]

Run it after a coverage run, i.e. once a .coverage data file exists.

Exit codes:
    0: at or above the hard floor (a warning is printed if below target)
    1: below the hard floor, or coverage data could not be read
"""

from __future__ import annotations

import argparse
import io
import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = REPO_ROOT / "pyproject.toml"


def read_floor() -> float:
    """Read the hard floor from ``[tool.coverage.report] fail_under``.

    Returns:
        float: The configured floor.

    Raises:
        SystemExit: If the floor cannot be read. A gate that cannot find its
            threshold must fail rather than silently pass everything.
    """
    try:
        text = PYPROJECT.read_text(encoding="utf-8")
    except OSError as exc:
        raise SystemExit(f"Could not read {PYPROJECT}: {exc}. Coverage gate cannot run.") from exc
    match = re.search(r"^fail_under\s*=\s*([0-9.]+)", text, re.MULTILINE)
    if match is None:
        raise SystemExit(f"No `fail_under` found in {PYPROJECT} [tool.coverage.report]. Coverage gate cannot run.")
    return float(match.group(1))


def read_total() -> float | None:
    """Read total coverage percentage from the existing coverage data file.

    Uses coverage.py's Python API rather than shelling out to `coverage report`. That
    is faster, does not depend on the CLI being on PATH, and avoids spawning a
    subprocess altogether.

    Returns:
        Optional[float]: The total percentage, or None if coverage could not be read.
    """
    try:
        import coverage
    except ImportError:
        return None

    try:
        cov = coverage.Coverage(config_file=str(PYPROJECT))
        cov.load()
        # report() returns the total percentage; the text report itself is discarded.
        # fail_under is enforced by this script, so ignore coverage's own exception.
        return float(cov.report(file=io.StringIO()))
    except Exception:
        # A missing/corrupt .coverage file, or no measured data. Callers treat None as
        # "could not read" and exit non-zero with an actionable message.
        return None


def emit(line: str, *, gh_command: str = "") -> None:
    """Print a line, additionally annotating it for GitHub Actions when running in CI.

    Args:
        line: The human-readable message.
        gh_command: One of "warning", "notice", "error", or "" for a plain line.
    """
    if gh_command and os.environ.get("GITHUB_ACTIONS") == "true":
        print(f"::{gh_command}::{line}")
    else:
        print(line)


def write_summary(total: float, floor: float, target: float) -> None:
    """Append a progress table to the GitHub Actions step summary, if available.

    Args:
        total: Measured coverage percentage.
        floor: Hard floor percentage.
        target: Soft target percentage.
    """
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary:
        return
    filled = int(total / 100 * 40)
    bar = "█" * filled + "░" * (40 - filled)
    status = "✅ at target" if total >= target else "🎯 below target"
    with open(summary, "a", encoding="utf-8") as handle:
        handle.write(
            f"### Unit test coverage\n\n"
            f"`{bar}` **{total:.2f}%**\n\n"
            f"| | % |\n| --- | --- |\n"
            f"| Measured | **{total:.2f}** |\n"
            f"| Hard floor (CI fails below) | {floor:.0f} |\n"
            f"| Target | {target:.0f} — {status} |\n\n"
        )


def main() -> int:
    """Compare measured coverage against the floor and the target.

    Returns:
        int: Process exit code.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        type=float,
        default=80.0,
        help="Soft target percentage to strive for (default: 80).",
    )
    parser.add_argument(
        "--ratchet-margin",
        type=float,
        default=0.25,
        help="Suggest raising the floor once coverage exceeds it by this many points.",
    )
    args = parser.parse_args()

    floor = read_floor()
    total = read_total()

    if total is None:
        emit(
            "Could not read coverage data. Run a coverage command first, e.g. "
            "`make test-unit-cov`.",
            gh_command="error",
        )
        return 1

    print(f"coverage {total:.2f}%  (floor {floor:.0f}%, target {args.target:.0f}%)")
    write_summary(total, floor, args.target)

    if total + 1e-9 < floor:
        emit(
            f"Coverage {total:.2f}% is below the hard floor of {floor:.0f}%. "
            "Add tests for the code you changed; do not lower the floor.",
            gh_command="error",
        )
        return 1

    # The ratchet: lock in gains so the floor tracks reality and the target stays live.
    # `suggested > floor` matters as much as the margin: the floor is an integer, so at
    # 36.81% against a floor of 36 there is nothing to lock in, and printing
    # `set fail_under = 36` when it is already 36 trains readers to ignore the notice.
    suggested = int(total)
    if total >= floor + args.ratchet_margin and suggested > floor:
        emit(
            f"Coverage {total:.2f}% now exceeds the floor of {floor:.0f}%. "
            f"Raise it: set `fail_under = {suggested}` in pyproject.toml "
            "[tool.coverage.report].",
            gh_command="notice",
        )

    if total < args.target:
        emit(
            f"Coverage {total:.2f}% is below the {args.target:.0f}% target "
            f"({args.target - total:.2f} points to go). Not a failure - but any file you "
            "touch should leave this higher than you found it.",
            gh_command="warning",
        )
    else:
        emit(f"Coverage {total:.2f}% meets the {args.target:.0f}% target.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
