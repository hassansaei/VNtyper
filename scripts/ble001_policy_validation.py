"""Pure path and test-node validation for the BLE001 policy adapter."""

from __future__ import annotations

import ast
from pathlib import Path
from subprocess import run as run_process
from typing import Any


def validate_scope_paths(repo_root: Path, scope: tuple[str, ...]) -> None:
    """Require every Ruff scope token to be ordinary, contained, and Git-tracked.

    Args:
        repo_root: Repository root against which path operands are interpreted.
        scope: Ruff path operands to validate without invoking Ruff.

    Raises:
        ValueError: If a token is malformed, missing, a symlink, outside the repository, or untracked.
    """
    if not scope:
        raise ValueError("RUFF_PATHS scope token list must not be empty")
    repository_root = repo_root.resolve()
    for token in scope:
        candidate = Path(token)
        if token.startswith("-") or candidate.is_absolute() or "\\" in token:
            raise ValueError(f"RUFF_PATHS scope token must be a repository-relative ordinary path: {token!r}")
        if any(part in {"", ".", ".."} for part in token.rstrip("/").split("/")):
            raise ValueError(f"RUFF_PATHS scope token must be a normalized repository-relative path: {token!r}")
        _reject_symlink_components(repository_root, candidate, f"RUFF_PATHS scope token {token!r}")
        resolved = (repository_root / candidate).resolve()
        try:
            resolved.relative_to(repository_root)
        except ValueError as exc:
            raise ValueError(f"RUFF_PATHS scope token escapes repository root: {token!r}") from exc
        if not resolved.is_file() and not resolved.is_dir():
            raise ValueError(f"RUFF_PATHS scope token is not an existing ordinary path: {token!r}")
        tracked_error = _tracked_path_error(repository_root, candidate)
        if tracked_error is not None:
            raise ValueError(f"RUFF_PATHS scope token is not tracked evidence: {token!r}: {tracked_error}")


def behavior_node_error(repo_root: Path, node_id: str) -> str | None:
    """Return why a unit-test node cannot be resolved from tracked ordinary source.

    Args:
        repo_root: Repository root containing the reviewed unit-test source.
        node_id: Pytest node ID with a repository-relative source path.

    Returns:
        ``None`` when the node resolves uniquely, otherwise an actionable error string.
    """
    path_text, *qualifiers = node_id.split("::")
    if not qualifiers:
        return f"behavior-test node must name a test function after '::': {path_text}"
    source_path = Path(path_text)
    if source_path.is_absolute() or any(part in {"", ".", ".."} for part in source_path.parts):
        return f"source path is not normalized and repository-relative: {path_text}"
    try:
        _reject_symlink_components(repo_root.resolve(), source_path, f"behavior-test source {path_text}")
    except ValueError as exc:
        return str(exc)
    try:
        source = (repo_root / path_text).read_text(encoding="utf-8")
    except OSError as exc:
        return f"could not read {path_text}: {exc}"
    tracked_error = _tracked_path_error(repo_root.resolve(), source_path)
    if tracked_error is not None:
        return f"behavior-test source is not tracked evidence: {path_text}: {tracked_error}"
    try:
        body = ast.parse(source).body
    except SyntaxError as exc:
        return f"could not parse {path_text}: {exc}"
    qualifiers[-1] = qualifiers[-1].split("[", maxsplit=1)[0]
    for index, name in enumerate(qualifiers):
        allowed = (ast.FunctionDef, ast.AsyncFunctionDef) if index == len(qualifiers) - 1 else (ast.ClassDef,)
        matches = [node for node in body if isinstance(node, allowed) and node.name == name]
        if len(matches) != 1:
            return f"could not resolve {'::'.join(qualifiers[: index + 1])} in {path_text}"
        body = matches[0].body
    return None


def require_exact_keys(record: dict[str, Any], expected: set[str], context: str) -> None:
    """Reject missing or unknown JSON object keys.

    Args:
        record: JSON object to validate.
        expected: Exact accepted key set.
        context: Human-readable object location for errors.

    Raises:
        ValueError: If any expected key is missing or any unknown key is present.
    """
    missing = sorted(expected - set(record))
    unknown = sorted(set(record) - expected)
    if missing or unknown:
        raise ValueError(f"{context} has missing keys {missing} and unknown keys {unknown}")


def require_non_empty_string(record: dict[str, Any], key: str, context: str) -> str:
    """Return one required non-empty string field.

    Args:
        record: JSON object containing the field.
        key: Required field name.
        context: Human-readable object location for errors.

    Returns:
        The validated single-line string.

    Raises:
        ValueError: If the field is not a non-empty single-line string.
    """
    value = record[key]
    if not isinstance(value, str) or not value.strip() or "\n" in value or "\r" in value:
        raise ValueError(f"{context}.{key} must be a non-empty single-line string")
    return value


def require_count(record: dict[str, Any], key: str, context: str) -> int:
    """Return one required non-negative integer count.

    Args:
        record: JSON object containing the count.
        key: Required count field name.
        context: Human-readable object location for errors.

    Returns:
        The validated count.

    Raises:
        ValueError: If the field is not a non-negative integer.
    """
    value = record[key]
    if type(value) is not int or value < 0:
        raise ValueError(f"{context}.{key} must be a non-negative integer")
    return value


def require_relative_path(record: dict[str, Any], key: str, context: str) -> str:
    """Return one normalized repository-relative POSIX path.

    Args:
        record: JSON object containing the path.
        key: Required path field name.
        context: Human-readable object location for errors.

    Returns:
        The validated repository-relative path.

    Raises:
        ValueError: If the field is absolute, non-normal, or not a non-empty string.
    """
    value = require_non_empty_string(record, key, context)
    if Path(value).is_absolute():
        raise ValueError(f"{context}.{key} must be relative")
    if "\\" in value or any(part in {".", "..", ""} for part in value.split("/")):
        raise ValueError(f"{context}.{key} contains parent traversal or a non-normal path component")
    return value


def _reject_symlink_components(repo_root: Path, candidate: Path, context: str) -> None:
    """Reject a symlink at any component of a repository-relative path."""
    current = repo_root
    for part in candidate.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"{context} must not contain a symlink component: {current}")


def _tracked_path_error(repo_root: Path, candidate: Path) -> str | None:
    """Return why Git cannot verify one literal path as tracked evidence."""
    try:
        result = run_process(
            [
                "git",
                "-C",
                str(repo_root),
                "--literal-pathspecs",
                "ls-files",
                "--error-unmatch",
                "--",
                candidate.as_posix(),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        return f"could not run git ls-files: {exc}"
    if result.returncode != 0 or not result.stdout.strip():
        return result.stderr.strip() or "git ls-files returned no tracked path"
    return None
