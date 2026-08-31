"""Architectural invariants for the non-executable comparator rule boundary."""

from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import pytest

import vntyper.scripts.comparator_rules as comparator_rules_module
import vntyper.scripts.cross_match as cross_match_module
import vntyper.scripts.flagging as flagging_module

pytestmark = pytest.mark.unit

_BOUNDARY_MODULES = (comparator_rules_module, flagging_module, cross_match_module)
_FORBIDDEN_CALL_TERMINALS = frozenset({"eval", "exec", "compile", "literal_eval", "safe_eval", "simple_eval"})
_FORBIDDEN_CALL_PATHS = frozenset(
    {
        "ast.parse",
        "asteval.Interpreter",
        "jinja2.Environment",
        "jinja2.Template",
        "numexpr.evaluate",
        "simpleeval.EvalWithCompoundTypes",
        "simpleeval.SimpleEval",
    }
)
_FORBIDDEN_IMPORT_ROOTS = frozenset({"ast", "asteval", "jinja2", "numexpr", "simpleeval"})


def _source_tree(module: ModuleType) -> ast.Module:
    source_path = Path(module.__file__ or "")
    assert source_path.is_file(), f"source file missing for {module.__name__}"
    return ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))


def _dotted_name(node: ast.expr) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _dotted_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return None


def test_rule_boundary_has_no_executable_expression_or_template_entry_point() -> None:
    violations: list[str] = []

    for module in _BOUNDARY_MODULES:
        for node in ast.walk(_source_tree(module)):
            if isinstance(node, ast.Call):
                call_name = _dotted_name(node.func)
                if call_name is not None and (
                    call_name.rsplit(".", 1)[-1] in _FORBIDDEN_CALL_TERMINALS or call_name in _FORBIDDEN_CALL_PATHS
                ):
                    violations.append(f"{module.__name__}:{node.lineno}: call {call_name}")
            elif isinstance(node, ast.Import):
                violations.extend(
                    f"{module.__name__}:{node.lineno}: import {alias.name}"
                    for alias in node.names
                    if alias.name.split(".", 1)[0] in _FORBIDDEN_IMPORT_ROOTS
                )
            elif isinstance(node, ast.ImportFrom):
                imported_root = (node.module or "").split(".", 1)[0]
                if imported_root in _FORBIDDEN_IMPORT_ROOTS:
                    violations.append(f"{module.__name__}:{node.lineno}: import from {node.module}")

    assert violations == []


def test_flagging_does_not_export_retired_executable_string_helpers() -> None:
    assert not hasattr(flagging_module, "evaluate_condition")
    assert not hasattr(flagging_module, "regex_match")
