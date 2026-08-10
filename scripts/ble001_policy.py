"""Measure and validate the repository's reviewed BLE001 exception policy."""

from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from scripts.ble001_policy_validation import (
        behavior_node_error,
        require_count,
        require_exact_keys,
        require_non_empty_string,
        require_relative_path,
        validate_scope_paths,
    )
elif __package__:
    from .ble001_policy_validation import (
        behavior_node_error,
        require_count,
        require_exact_keys,
        require_non_empty_string,
        require_relative_path,
        validate_scope_paths,
    )
else:
    from ble001_policy_validation import (
        behavior_node_error,
        require_count,
        require_exact_keys,
        require_non_empty_string,
        require_relative_path,
        validate_scope_paths,
    )

TEST_NODE_ID = re.compile(
    r"^tests/unit/(?:[A-Za-z0-9_.-]+/)*test_[A-Za-z0-9_]+\.py"
    r"(?:::[A-Za-z_][A-Za-z0-9_]*)*::test_[A-Za-z_][A-Za-z0-9_]*"
    r"(?:\[[^\]\r\n]+\])?$"
)
RUFF_VERSION = re.compile(r"^ruff [0-9]+\.[0-9]+\.[0-9]+(?:[^\s]*)?$")
DISPOSITIONS = {
    "preserved-contract",
    "conformed-to-existing-contract",
    "preserved-no-authorized-alternative",
}


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
    expected_identities: int
    expected_categories: tuple[int, int, int]
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
    validate_scope_paths(makefile.parent, scope)
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
    validate_scope_paths(repo_root, scope)
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
    if RUFF_VERSION.fullmatch(ruff_version) is None:
        raise RuntimeError(f"Ruff returned malformed or empty version output: {ruff_version!r}")

    command = [ruff_executable, "check", "--no-cache", "--select", "BLE001", "--output-format", "json"]
    if ignore_noqa:
        command.append("--ignore-noqa")
    command.append("--")
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
    """Load and validate a frozen BLE001 policy artifact.

    Args:
        path: JSON policy artifact path.

    Returns:
        Immutable, fully validated policy data.

    Raises:
        ValueError: If the file is unreadable or any schema invariant is violated.
    """
    try:
        source = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Could not read policy {path}: {exc}") from exc
    try:
        payload: Any = json.loads(source)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Could not parse policy JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Policy must be a JSON object")
    require_exact_keys(
        payload,
        {"reviewed_ruff_version", "expected_counts", "handlers", "fail_open"},
        "policy",
    )

    reviewed_version = require_non_empty_string(payload, "reviewed_ruff_version", "policy")
    if RUFF_VERSION.fullmatch(reviewed_version) is None:
        raise ValueError(f"policy.reviewed_ruff_version is malformed: {reviewed_version!r}")

    expected_counts = payload["expected_counts"]
    if not isinstance(expected_counts, dict):
        raise ValueError("policy.expected_counts must be an object")
    require_exact_keys(
        expected_counts,
        {"normal", "ignore_noqa", "identities", "categories"},
        "policy.expected_counts",
    )
    expected_normal = require_count(expected_counts, "normal", "policy.expected_counts")
    expected_all = require_count(expected_counts, "ignore_noqa", "policy.expected_counts")
    expected_identities = require_count(expected_counts, "identities", "policy.expected_counts")
    raw_categories = expected_counts["categories"]
    if not isinstance(raw_categories, dict):
        raise ValueError("policy.expected_counts.categories must be an object")
    require_exact_keys(raw_categories, {"A", "B", "C"}, "policy.expected_counts.categories")
    expected_categories = (
        require_count(raw_categories, "A", "policy.expected_counts.categories"),
        require_count(raw_categories, "B", "policy.expected_counts.categories"),
        require_count(raw_categories, "C", "policy.expected_counts.categories"),
    )
    if expected_normal > expected_all:
        raise ValueError("policy expected normal count cannot exceed ignore_noqa count")

    raw_handlers = payload["handlers"]
    if not isinstance(raw_handlers, list):
        raise ValueError("policy.handlers must be an array")
    if not raw_handlers:
        raise ValueError("policy.handlers must contain at least one handler")
    if not expected_normal or not expected_all or not expected_identities:
        raise ValueError("policy expected normal, ignore_noqa, and identity counts must be strictly positive")
    handlers: list[HandlerPolicy] = []
    handler_keys: set[tuple[str, str]] = set()
    for index, raw_handler in enumerate(raw_handlers):
        context = f"policy.handlers[{index}]"
        if not isinstance(raw_handler, dict):
            raise ValueError(f"{context} must be an object")
        require_exact_keys(
            raw_handler,
            {"path", "symbol", "normal_count", "all_count", "category", "rationale"},
            context,
        )
        handler_path = require_relative_path(raw_handler, "path", context)
        symbol = require_non_empty_string(raw_handler, "symbol", context)
        normal_count = require_count(raw_handler, "normal_count", context)
        all_count = require_count(raw_handler, "all_count", context)
        if not all_count:
            raise ValueError(f"{context}.all_count must be strictly positive")
        if normal_count > all_count:
            raise ValueError(f"{context}.normal_count cannot exceed all_count")
        category = require_non_empty_string(raw_handler, "category", context)
        if category not in {"A", "B", "C"}:
            raise ValueError(f"{context}.category must be A, B, or C")
        rationale = require_non_empty_string(raw_handler, "rationale", context)
        key = (handler_path, symbol)
        if key in handler_keys:
            raise ValueError(f"duplicate handler policy for {handler_path}::{symbol}")
        handler_keys.add(key)
        handlers.append(HandlerPolicy(handler_path, symbol, normal_count, all_count, category, rationale))

    normal_total = sum(handler.normal_count for handler in handlers)
    all_total = sum(handler.all_count for handler in handlers)
    if (normal_total, all_total) != (expected_normal, expected_all):
        raise ValueError(
            "policy handler totals do not match expected_counts: "
            f"handlers {normal_total}/{all_total}, expected {expected_normal}/{expected_all}"
        )
    if len(handlers) != expected_identities:
        raise ValueError(
            f"policy handler identity count does not match expected_counts: handlers {len(handlers)}, "
            f"expected {expected_identities}"
        )
    actual_categories = tuple(sum(handler.category == category for handler in handlers) for category in ("A", "B", "C"))
    if actual_categories != expected_categories:
        raise ValueError(
            "policy category identity counts do not match expected_counts: "
            f"handlers {actual_categories}, expected {expected_categories}"
        )

    raw_fail_open = payload["fail_open"]
    if not isinstance(raw_fail_open, list):
        raise ValueError("policy.fail_open must be an array")
    fail_open: list[FailOpenPolicy] = []
    fail_open_keys: set[tuple[str, str]] = set()
    for index, raw_record in enumerate(raw_fail_open):
        context = f"policy.fail_open[{index}]"
        if not isinstance(raw_record, dict):
            raise ValueError(f"{context} must be an object")
        require_exact_keys(
            raw_record,
            {"path", "symbol", "disposition", "outcome", "rationale", "behavior_test", "observable_via"},
            context,
        )
        record_path = require_relative_path(raw_record, "path", context)
        symbol = require_non_empty_string(raw_record, "symbol", context)
        disposition = require_non_empty_string(raw_record, "disposition", context)
        if disposition not in DISPOSITIONS:
            raise ValueError(f"{context}.disposition is not an allowed policy disposition")
        outcome = require_non_empty_string(raw_record, "outcome", context)
        rationale = require_non_empty_string(raw_record, "rationale", context)
        behavior_test = require_non_empty_string(raw_record, "behavior_test", context)
        if TEST_NODE_ID.fullmatch(behavior_test) is None or any(
            part in {".", ".."} for part in behavior_test.split("::", maxsplit=1)[0].split("/")
        ):
            raise ValueError(f"{context}.behavior_test is not a valid unit-test node ID: {behavior_test!r}")
        observable_via = require_non_empty_string(raw_record, "observable_via", context)
        key = (record_path, symbol)
        if key in fail_open_keys:
            raise ValueError(f"duplicate fail_open policy for {record_path}::{symbol}")
        fail_open_keys.add(key)
        fail_open.append(
            FailOpenPolicy(
                record_path,
                symbol,
                disposition,
                outcome,
                rationale,
                behavior_test,
                observable_via,
            )
        )

    category_c_keys = {(handler.path, handler.symbol) for handler in handlers if handler.category == "C"}
    missing_fail_open = sorted(category_c_keys - fail_open_keys)
    if missing_fail_open:
        rendered = ", ".join(f"{item_path}::{symbol}" for item_path, symbol in missing_fail_open)
        raise ValueError(f"category C handlers require exactly one fail_open record; missing: {rendered}")
    orphan_fail_open = sorted(fail_open_keys - category_c_keys)
    if orphan_fail_open:
        rendered = ", ".join(f"{item_path}::{symbol}" for item_path, symbol in orphan_fail_open)
        raise ValueError(f"fail_open record exists for non-category-C handler: {rendered}")

    return Policy(
        reviewed_version,
        expected_normal,
        expected_all,
        expected_identities,
        expected_categories,
        tuple(handlers),
        tuple(fail_open),
    )


def validate_policy(
    repo_root: Path,
    normal: Measurement,
    all_handlers: Measurement,
    policy: Policy,
) -> list[str]:
    """Return all mismatches between live measurements and reviewed policy.

    Args:
        repo_root: Repository root containing every diagnostic source path.
        normal: Ruff measurement using normal suppression semantics.
        all_handlers: Ruff measurement including explicit suppressions.
        policy: Reviewed policy to compare with live discovery.

    Returns:
        Actionable mismatch descriptions. An empty list means exact agreement.
    """
    errors: list[str] = []
    if normal.ruff_version != all_handlers.ruff_version:
        errors.append(
            f"Measurement version mismatch: normal Ruff {normal.ruff_version}; "
            f"ignore-noqa Ruff {all_handlers.ruff_version}. Rerun both modes with one executable."
        )
        return errors
    inventory_errors: list[str] = []
    if len(normal.diagnostics) != policy.expected_normal:
        inventory_errors.append(
            f"normal count mismatch: expected {policy.expected_normal}, actual {len(normal.diagnostics)} "
            f"under reviewed {policy.reviewed_ruff_version} and actual {normal.ruff_version}"
        )
    if len(all_handlers.diagnostics) != policy.expected_all:
        inventory_errors.append(
            f"ignore-noqa count mismatch: expected {policy.expected_all}, actual {len(all_handlers.diagnostics)} "
            f"under reviewed {policy.reviewed_ruff_version} and actual {all_handlers.ruff_version}"
        )

    normal_records = Counter(normal.diagnostics)
    all_records = Counter(all_handlers.diagnostics)
    if normal_records - all_records:
        inventory_errors.append("normal diagnostics are not a subset of ignore-noqa diagnostics")

    normal_counts, normal_rows, normal_mapping_errors = _symbol_counts(repo_root, normal.diagnostics)
    all_counts, all_rows, all_mapping_errors = _symbol_counts(repo_root, all_handlers.diagnostics)
    inventory_errors.extend(normal_mapping_errors)
    inventory_errors.extend(error for error in all_mapping_errors if error not in normal_mapping_errors)

    expected_normal = {
        (handler.path, handler.symbol): handler.normal_count for handler in policy.handlers if handler.normal_count
    }
    expected_all = {
        (handler.path, handler.symbol): handler.all_count for handler in policy.handlers if handler.all_count
    }
    inventory_errors.extend(_count_mismatches("normal", expected_normal, normal_counts, normal_rows))
    inventory_errors.extend(_count_mismatches("ignore-noqa", expected_all, all_counts, all_rows))

    actual_categories = tuple(
        sum(handler.category == category for handler in policy.handlers) for category in ("A", "B", "C")
    )
    if len(policy.handlers) != policy.expected_identities:
        inventory_errors.append(
            f"handler identity count mismatch: expected {policy.expected_identities}, actual {len(policy.handlers)}"
        )
    if actual_categories != policy.expected_categories:
        inventory_errors.append(
            f"category identity counts mismatch: expected {policy.expected_categories}, actual {actual_categories}"
        )
    if inventory_errors and normal.ruff_version != policy.reviewed_ruff_version:
        errors.append(
            f"Inventory drift measured under reviewed Ruff {policy.reviewed_ruff_version} and actual Ruff "
            f"{normal.ruff_version}; review and reclassify the changed identities."
        )
    errors.extend(inventory_errors)

    category_c_keys = {(handler.path, handler.symbol) for handler in policy.handlers if handler.category == "C"}
    fail_open_keys = {(record.path, record.symbol) for record in policy.fail_open}
    for item_path, symbol in sorted(category_c_keys - fail_open_keys):
        errors.append(f"category C handler lacks fail_open policy: {item_path}::{symbol}")
    for item_path, symbol in sorted(fail_open_keys - category_c_keys):
        errors.append(f"stale fail_open policy for non-category-C handler: {item_path}::{symbol}")
    for record in policy.fail_open:
        node_error = behavior_node_error(repo_root, record.behavior_test)
        if node_error is not None:
            errors.append(f"unresolved behavior-test node {record.behavior_test}: {node_error}")
    return errors


def _symbol_counts(
    repo_root: Path, diagnostics: tuple[Diagnostic, ...]
) -> tuple[dict[tuple[str, str], int], dict[tuple[str, str], list[int]], list[str]]:
    """Map normalized diagnostics to stable symbol counts and observed rows."""
    counts: Counter[tuple[str, str]] = Counter()
    rows: dict[tuple[str, str], list[int]] = {}
    errors: list[str] = []
    source_cache: dict[str, str | None] = {}
    for diagnostic in diagnostics:
        if diagnostic.path not in source_cache:
            try:
                source_cache[diagnostic.path] = (repo_root / diagnostic.path).read_text(encoding="utf-8")
            except OSError as exc:
                source_cache[diagnostic.path] = None
                errors.append(f"Could not read diagnostic source {diagnostic.path}: {exc}")
        source = source_cache[diagnostic.path]
        if source is None:
            continue
        try:
            symbol = enclosing_symbol(source, diagnostic.row)
        except ValueError as exc:
            errors.append(f"Could not map diagnostic source {diagnostic.path} row {diagnostic.row}: {exc}")
            continue
        identity = (diagnostic.path, symbol)
        counts[identity] += 1
        rows.setdefault(identity, []).append(diagnostic.row)
    return dict(counts), rows, errors


def _count_mismatches(
    mode: str,
    expected: dict[tuple[str, str], int],
    actual: dict[tuple[str, str], int],
    rows: dict[tuple[str, str], list[int]],
) -> list[str]:
    """Describe exact added, removed, and changed stable identity counts."""
    errors: list[str] = []
    for identity in sorted(set(expected) | set(actual)):
        expected_count = expected.get(identity, 0)
        actual_count = actual.get(identity, 0)
        if expected_count == actual_count:
            continue
        item_path, symbol = identity
        if actual_count == 0:
            errors.append(
                f"{mode} removed reviewed identity {item_path}::{symbol}: expected {expected_count}, actual 0"
            )
        else:
            errors.append(
                f"{mode} identity {item_path}::{symbol}: expected {expected_count}, actual {actual_count}; "
                f"observed rows {rows.get(identity, [])}"
            )
    return errors


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
    print(f"suppression delta: {policy.expected_all - policy.expected_normal}")
    print(f"categories A/B/C: {category_counts['A']}/{category_counts['B']}/{category_counts['C']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
