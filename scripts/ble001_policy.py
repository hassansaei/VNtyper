"""Measure and validate the repository's reviewed BLE001 exception policy."""

from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, order=True)
class Diagnostic:
    """One normalized Ruff BLE001 diagnostic."""

    path: str
    row: int
    column: int
    code: str
    message: str


@dataclass(frozen=True)
class Measurement:
    """A Ruff version and the diagnostics it emitted."""

    ruff_version: str
    diagnostics: tuple[Diagnostic, ...]


@dataclass(frozen=True)
class HandlerPolicy:
    """Reviewed classification for the broad handlers in one symbol."""

    path: str
    symbol: str
    normal_count: int
    all_count: int
    category: str
    rationale: str


@dataclass(frozen=True)
class FailOpenPolicy:
    """Reviewed behavior evidence for a category-C handler."""

    path: str
    symbol: str
    disposition: str
    outcome: str
    rationale: str
    behavior_test: str
    observable_via: str


@dataclass(frozen=True)
class Policy:
    """Complete reviewed BLE001 policy data."""

    reviewed_ruff_version: str
    expected_normal: int
    expected_all: int
    handlers: tuple[HandlerPolicy, ...]
    fail_open: tuple[FailOpenPolicy, ...]


def read_ruff_paths(makefile: Path) -> tuple[str, ...]:
    """Read the canonical Ruff path scope from a Makefile.

    Args:
        makefile: Makefile containing the single authoritative ``RUFF_PATHS`` assignment.

    Returns:
        The ordered path tokens from the assignment.

    Raises:
        ValueError: If the Makefile is unreadable or does not contain exactly one non-empty assignment.
    """
    try:
        source = makefile.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Could not read Makefile {makefile}: {exc}") from exc
    assignments = re.findall(r"^RUFF_PATHS\s*:?=\s*(.*?)\s*$", source, flags=re.MULTILINE)
    if len(assignments) != 1:
        raise ValueError(f"Expected exactly one RUFF_PATHS assignment in {makefile}; found {len(assignments)}")
    scope = tuple(assignments[0].split())
    if not scope:
        raise ValueError(f"Expected exactly one RUFF_PATHS assignment with at least one path in {makefile}")
    return scope


def normalize_diagnostics(payload: str, repo_root: Path) -> tuple[Diagnostic, ...]:
    """Normalize Ruff JSON diagnostics to deterministic repository-relative records.

    Args:
        payload: Ruff JSON output.
        repo_root: Repository root against which diagnostic paths are resolved.

    Returns:
        Sorted immutable diagnostic records.

    Raises:
        ValueError: If JSON, record fields, locations, or paths are malformed.
    """
    try:
        decoded: Any = json.loads(payload)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError(f"Malformed Ruff JSON: {exc}") from exc
    if not isinstance(decoded, list):
        raise ValueError("Malformed Ruff JSON: top-level value must be an array")

    root = repo_root.resolve()
    diagnostics: list[Diagnostic] = []
    for index, record in enumerate(decoded):
        if not isinstance(record, dict):
            raise ValueError(f"Malformed Ruff diagnostic at index {index}: expected an object")
        filename = record.get("filename")
        code = record.get("code")
        message = record.get("message")
        if not all(isinstance(value, str) and value for value in (filename, code, message)):
            raise ValueError(f"Malformed Ruff diagnostic at index {index}: invalid filename, code, or message")
        assert isinstance(filename, str)
        assert isinstance(code, str)
        assert isinstance(message, str)
        diagnostic_path = Path(filename)
        resolved = diagnostic_path.resolve() if diagnostic_path.is_absolute() else (root / diagnostic_path).resolve()
        try:
            relative_path = resolved.relative_to(root).as_posix()
        except ValueError as exc:
            raise ValueError(f"Ruff diagnostic path escapes repository root: {filename}") from exc

        location = record.get("location")
        if not isinstance(location, dict):
            raise ValueError(f"Malformed Ruff diagnostic at index {index}: missing location")
        position: object = location
        cell = location.get("cell")
        if isinstance(cell, dict) and "start" in cell:
            position = cell["start"]
        if not isinstance(position, dict):
            raise ValueError(f"Malformed Ruff diagnostic at index {index}: invalid location")
        row = position.get("row")
        column = position.get("column")
        if type(row) is not int or type(column) is not int or row <= 0 or column <= 0:
            raise ValueError(f"Malformed Ruff diagnostic at index {index}: invalid row or column")
        diagnostics.append(Diagnostic(relative_path, row, column, code, message))
    return tuple(sorted(diagnostics))


def measure_ble001(
    repo_root: Path,
    scope: tuple[str, ...],
    *,
    ignore_noqa: bool,
    ruff_executable: str = "ruff",
) -> Measurement:
    """Run a read-only no-cache BLE001 measurement.

    Args:
        repo_root: Repository root used as Ruff's working directory.
        scope: Ordered repository path scope to scan.
        ignore_noqa: Whether Ruff should include explicitly suppressed diagnostics.
        ruff_executable: Ruff executable name or path.

    Returns:
        The actual Ruff version and normalized diagnostics.

    Raises:
        RuntimeError: If Ruff is unavailable, its version is unusable, or its check fails operationally.
        ValueError: If Ruff emits malformed diagnostic JSON.
    """
    try:
        version_result = subprocess.run(
            [ruff_executable, "--version"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise RuntimeError(f"Could not execute Ruff {ruff_executable}: {exc}") from exc
    if version_result.returncode != 0:
        detail = version_result.stderr.strip() or version_result.stdout.strip() or "no diagnostic output"
        raise RuntimeError(f"Could not determine Ruff version from {ruff_executable}: {detail}")
    ruff_version = version_result.stdout.strip()
    if re.fullmatch(r"ruff [0-9]+\.[0-9]+\.[0-9]+(?:[^\s]*)?", ruff_version) is None:
        raise RuntimeError(f"Ruff returned malformed or empty version output: {ruff_version!r}")

    command = [ruff_executable, "check", "--no-cache", "--select", "BLE001", "--output-format", "json"]
    if ignore_noqa:
        command.append("--ignore-noqa")
    command.extend(scope)
    try:
        check_result = subprocess.run(
            command,
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise RuntimeError(f"Could not execute Ruff {ruff_executable} ({ruff_version}): {exc}") from exc
    if check_result.returncode not in (0, 1):
        detail = check_result.stderr.strip() or check_result.stdout.strip() or "no diagnostic output"
        raise RuntimeError(f"Ruff check failed under {ruff_version} with status {check_result.returncode}: {detail}")
    return Measurement(ruff_version, normalize_diagnostics(check_result.stdout, repo_root))


def enclosing_symbol(source: str, row: int) -> str:
    """Return the deepest qualified AST symbol containing a source row.

    Args:
        source: Python source to parse without importing or executing it.
        row: One-based source row for a Ruff diagnostic.

    Returns:
        A qualified class/function path, or ``<module>`` for module-level rows.

    Raises:
        ValueError: If the row is invalid or the source cannot be parsed.
    """
    if row <= 0:
        raise ValueError("Diagnostic row must be positive")
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        raise ValueError(f"Could not parse Python source: {exc}") from exc

    candidates: list[tuple[int, str]] = []

    def visit(node: ast.AST, parents: tuple[str, ...]) -> None:
        next_parents = parents
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            end_lineno = node.end_lineno
            if end_lineno is not None and node.lineno <= row <= end_lineno:
                next_parents = (*parents, node.name)
                candidates.append((len(next_parents), ".".join(next_parents)))
        for child in ast.iter_child_nodes(node):
            visit(child, next_parents)

    visit(tree, ())
    if not candidates:
        return "<module>"
    return max(candidates)[1]


def load_policy(path: Path) -> Policy:
    """Load and validate a frozen BLE001 policy artifact."""
    raise NotImplementedError("not implemented")


def validate_policy(
    repo_root: Path,
    normal: Measurement,
    all_handlers: Measurement,
    policy: Policy,
) -> list[str]:
    """Return all mismatches between live measurements and reviewed policy."""
    raise NotImplementedError("not implemented")


def main(argv: list[str] | None = None) -> int:
    """Measure and validate the repository's BLE001 policy.

    Args:
        argv: Optional command arguments excluding the program name.

    Returns:
        Zero when policy validation succeeds, otherwise one.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--policy", type=Path, default=Path("scripts/ble001_policy.json"))
    parser.add_argument("--ruff", default="ruff")
    args = parser.parse_args(argv)

    repo_root = args.repo_root.resolve()
    policy_path = args.policy if args.policy.is_absolute() else repo_root / args.policy
    try:
        scope = read_ruff_paths(repo_root / "Makefile")
        policy = load_policy(policy_path)
        normal = measure_ble001(repo_root, scope, ignore_noqa=False, ruff_executable=args.ruff)
        all_handlers = measure_ble001(repo_root, scope, ignore_noqa=True, ruff_executable=args.ruff)
        errors = validate_policy(repo_root, normal, all_handlers, policy)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"BLE001 policy validation failed: {exc}")
        return 1

    if errors:
        for error in errors:
            print(error)
        return 1

    category_counts = {
        category: sum(handler.all_count for handler in policy.handlers if handler.category == category)
        for category in ("A", "B", "C")
    }
    print(f"reviewed Ruff: {policy.reviewed_ruff_version}")
    print(f"actual Ruff: {normal.ruff_version}")
    print(f"normal/all: {policy.expected_normal}/{policy.expected_all}")
    print(f"categories A/B/C: {category_counts['A']}/{category_counts['B']}/{category_counts['C']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
