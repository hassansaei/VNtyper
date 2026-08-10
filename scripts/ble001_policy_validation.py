"""Pure path and test-node validation for the BLE001 policy adapter."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any


def validate_scope_paths(makefile: Path, scope: tuple[str, ...]) -> None:
    """Require every Ruff scope token to name an ordinary in-repository path."""
    repository_root = makefile.parent.resolve()
    for token in scope:
        candidate = Path(token)
        if token.startswith("-") or candidate.is_absolute() or "\\" in token:
            raise ValueError(f"RUFF_PATHS token must be a repository-relative ordinary path: {token!r}")
        if any(part in {"", ".", ".."} for part in token.rstrip("/").split("/")):
            raise ValueError(f"RUFF_PATHS token must be a normalized repository-relative path: {token!r}")
        resolved = (repository_root / candidate).resolve()
        try:
            resolved.relative_to(repository_root)
        except ValueError as exc:
            raise ValueError(f"RUFF_PATHS token escapes repository root: {token!r}") from exc
        if not resolved.is_file() and not resolved.is_dir():
            raise ValueError(f"RUFF_PATHS token is not an existing ordinary path: {token!r}")


def behavior_node_error(repo_root: Path, node_id: str) -> str | None:
    """Return why a unit-test node cannot be found through static source inspection."""
    path_text, *qualifiers = node_id.split("::")
    try:
        source = (repo_root / path_text).read_text(encoding="utf-8")
    except OSError as exc:
        return f"could not read {path_text}: {exc}"
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
    """Reject missing or unknown JSON object keys."""
    missing = sorted(expected - set(record))
    unknown = sorted(set(record) - expected)
    if missing or unknown:
        raise ValueError(f"{context} has missing keys {missing} and unknown keys {unknown}")


def require_non_empty_string(record: dict[str, Any], key: str, context: str) -> str:
    """Return one required non-empty string field."""
    value = record[key]
    if not isinstance(value, str) or not value.strip() or "\n" in value or "\r" in value:
        raise ValueError(f"{context}.{key} must be a non-empty single-line string")
    return value


def require_count(record: dict[str, Any], key: str, context: str) -> int:
    """Return one required non-negative integer count."""
    value = record[key]
    if type(value) is not int or value < 0:
        raise ValueError(f"{context}.{key} must be a non-negative integer")
    return value


def require_relative_path(record: dict[str, Any], key: str, context: str) -> str:
    """Return one normalized repository-relative POSIX path."""
    value = require_non_empty_string(record, key, context)
    if Path(value).is_absolute():
        raise ValueError(f"{context}.{key} must be relative")
    if "\\" in value or any(part in {".", "..", ""} for part in value.split("/")):
        raise ValueError(f"{context}.{key} contains parent traversal or a non-normal path component")
    return value
