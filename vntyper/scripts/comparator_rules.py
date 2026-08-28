"""Validated comparator rules for configuration-driven row decisions."""

from __future__ import annotations

import copy
import logging
import math
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from numbers import Real
from typing import Literal, NoReturn, TypeAlias, cast

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

JsonScalar: TypeAlias = str | int | float | bool | None
Operator: TypeAlias = Literal["eq", "lt", "in", "casefold_eq"]
_OPERATORS = frozenset({"eq", "lt", "in", "casefold_eq"})
_NULL = object()


@dataclass(frozen=True)
class ColumnOperand:
    """An immutable reference to one allowlisted row column."""

    name: str


@dataclass(frozen=True)
class LiteralOperand:
    """An immutable JSON scalar or membership collection."""

    value: JsonScalar | tuple[JsonScalar, ...]


Operand: TypeAlias = ColumnOperand | LiteralOperand


@dataclass(frozen=True)
class Predicate:
    """One validated binary comparison."""

    left: Operand
    operator: Operator
    right: Operand


@dataclass(frozen=True)
class CompiledRule:
    """A non-empty immutable conjunction of validated predicates."""

    predicates: tuple[Predicate, ...]


def _invalid(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _scalar_family(value: JsonScalar) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, str):
        return "string"
    return "number"


def _validate_scalar(value: object, *, path: str) -> JsonScalar:
    if isinstance(value, float) and not math.isfinite(value):
        _invalid(f"{path} must be a finite JSON number")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    _invalid(f"{path} must be a JSON scalar")


def _validate_literal(value: object, *, path: str) -> JsonScalar | tuple[JsonScalar, ...]:
    if not isinstance(value, list):
        return _validate_scalar(value, path=path)
    if not value:
        _invalid(f"{path} must be a non-empty homogeneous JSON-scalar list")
    validated = tuple(_validate_scalar(item, path=f"{path}[{index}]") for index, item in enumerate(value))
    first_family = _scalar_family(validated[0])
    if any(_scalar_family(item) != first_family for item in validated[1:]):
        _invalid(f"{path} must be a homogeneous JSON-scalar list")
    return validated


def _validate_operand(value: object, *, allowed_columns: frozenset[str], path: str) -> Operand:
    if not isinstance(value, dict) or set(value) not in ({"column"}, {"literal"}):
        _invalid(f"{path} must contain exactly one of 'column' or 'literal'")
    if "column" in value:
        name = value["column"]
        if not isinstance(name, str) or not name:
            _invalid(f"{path}.column must be a non-empty string")
        if name not in allowed_columns:
            _invalid(f"{path}.column names unknown column {name!r}; allowed columns are {sorted(allowed_columns)!r}")
        return ColumnOperand(name=name)
    return LiteralOperand(value=_validate_literal(value["literal"], path=f"{path}.literal"))


def _reject_collection(operand: Operand, *, path: str) -> None:
    if isinstance(operand, LiteralOperand) and isinstance(operand.value, tuple):
        _invalid(f"{path} cannot be a literal list for this operator")


def _validate_literal_number(operand: Operand, *, path: str) -> None:
    if not isinstance(operand, LiteralOperand):
        return
    value = operand.value
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _invalid(f"{path}.literal must be a real number other than a boolean for operator 'lt'")


def _validate_literal_string(operand: Operand, *, path: str) -> None:
    if isinstance(operand, LiteralOperand) and not isinstance(operand.value, str):
        _invalid(f"{path}.literal must be a string for operator 'casefold_eq'")


def _validate_predicate_compatibility(predicate: Predicate, *, path: str) -> None:
    if predicate.operator == "in":
        if not isinstance(predicate.right, LiteralOperand) or not isinstance(predicate.right.value, tuple):
            _invalid(f"{path}.right must be a non-empty literal list for operator 'in'")
        if isinstance(predicate.right.value[0], bool):
            _invalid(f"{path}.right.literal does not support boolean values for operator 'in'")
        _reject_collection(predicate.left, path=f"{path}.left")
        if isinstance(predicate.left, LiteralOperand) and predicate.left.value is not None:
            if isinstance(predicate.left.value, bool):
                _invalid(f"{path}.left.literal does not support boolean values for operator 'in'")
            right_family = _scalar_family(predicate.right.value[0])
            if _scalar_family(cast(JsonScalar, predicate.left.value)) != right_family:
                _invalid(f"{path} has incompatible literal operand families")
        return

    _reject_collection(predicate.left, path=f"{path}.left")
    _reject_collection(predicate.right, path=f"{path}.right")
    if predicate.operator == "lt":
        _validate_literal_number(predicate.left, path=f"{path}.left")
        _validate_literal_number(predicate.right, path=f"{path}.right")
    elif predicate.operator == "casefold_eq":
        _validate_literal_string(predicate.left, path=f"{path}.left")
        _validate_literal_string(predicate.right, path=f"{path}.right")
    elif (
        isinstance(predicate.left, LiteralOperand)
        and isinstance(predicate.right, LiteralOperand)
        and predicate.left.value is not None
        and predicate.right.value is not None
        and _scalar_family(cast(JsonScalar, predicate.left.value))
        != _scalar_family(cast(JsonScalar, predicate.right.value))
    ):
        _invalid(f"{path} has incompatible literal operand families")


def validate_rule(rule: object, *, allowed_columns: Collection[str], context: str) -> CompiledRule:
    """Validate untrusted JSON data and compile it into immutable comparator values.

    Args:
        rule: Untrusted value loaded from configuration JSON.
        allowed_columns: Column names the consumer permits this rule to read.
        context: Configuration path included in every diagnostic.

    Returns:
        The immutable compiled rule.

    Raises:
        ValueError: If any part of the rule is malformed or unsupported.
    """
    if not isinstance(rule, dict) or set(rule) != {"all"}:
        _invalid(f"{context} must contain exactly the key 'all'")
    configured_predicates = rule["all"]
    if not isinstance(configured_predicates, list):
        _invalid(f"{context}.all must be a non-empty list")
    if not configured_predicates:
        _invalid(f"{context}.all must be a non-empty list")

    allowed = frozenset(allowed_columns)
    predicates: list[Predicate] = []
    for index, configured in enumerate(configured_predicates):
        path = f"{context}.all[{index}]"
        if not isinstance(configured, dict) or set(configured) != {"left", "operator", "right"}:
            _invalid(f"{path} must contain exactly 'left', 'operator', and 'right'")
        configured_operator = configured["operator"]
        if not isinstance(configured_operator, str) or configured_operator not in _OPERATORS:
            _invalid(f"{path}.operator must be one of {sorted(_OPERATORS)!r}")
        predicate = Predicate(
            left=_validate_operand(configured["left"], allowed_columns=allowed, path=f"{path}.left"),
            operator=cast(Operator, configured_operator),
            right=_validate_operand(configured["right"], allowed_columns=allowed, path=f"{path}.right"),
        )
        _validate_predicate_compatibility(predicate, path=path)
        predicates.append(predicate)
    return CompiledRule(predicates=tuple(predicates))


def _operand_value(operand: Operand, row: Mapping[str, object], *, context: str) -> object:
    if isinstance(operand, LiteralOperand):
        return operand.value
    if operand.name not in row:
        _invalid(f"{context} cannot evaluate missing row column {operand.name!r}")
    return row[operand.name]


def _normalize_null(value: object) -> object:
    if value is None or value is pd.NA:
        return _NULL
    if isinstance(value, (float, np.floating)) and math.isnan(value):
        return _NULL
    return value


def _runtime_family(value: object, *, context: str) -> str:
    if isinstance(value, (bool, np.bool_)):
        return "boolean"
    if isinstance(value, str):
        return "string"
    if isinstance(value, Real):
        if isinstance(value, (float, np.floating)) and not math.isfinite(value):
            _invalid(f"{context} requires finite real-number row values")
        return "number"
    _invalid(f"{context} requires JSON-scalar row values")


def _require_same_family(left: object, right: object, *, context: str) -> str:
    left_family = _runtime_family(left, context=context)
    right_family = _runtime_family(right, context=context)
    if left_family != right_family:
        _invalid(f"{context} requires compatible families")
    return left_family


def _evaluate_predicate(predicate: Predicate, row: Mapping[str, object], *, context: str) -> bool:
    left = _normalize_null(_operand_value(predicate.left, row, context=context))
    right = _operand_value(predicate.right, row, context=context)
    if not isinstance(right, tuple):
        right = _normalize_null(right)
    if predicate.operator == "eq":
        if left is not _NULL:
            _runtime_family(left, context=context)
        if right is not _NULL:
            _runtime_family(right, context=context)
        if left is _NULL or right is _NULL:
            return False
        _require_same_family(left, right, context=context)
        return bool(left == right)
    if predicate.operator == "lt":
        if (left is not _NULL and _runtime_family(left, context=context) != "number") or (
            right is not _NULL and _runtime_family(right, context=context) != "number"
        ):
            _invalid(f"{context} requires real numbers other than booleans for operator 'lt'")
        if left is _NULL or right is _NULL:
            return False
        return bool(left < right)  # type: ignore[operator]
    if predicate.operator == "in":
        if not isinstance(right, tuple):  # pragma: no cover - prevented by validation and immutable types
            _invalid(f"{context} requires a literal list on the right for operator 'in'")
        if left is not _NULL:
            left_family = _runtime_family(left, context=context)
            if left_family == "boolean":
                _invalid(f"{context} does not support boolean operands for operator 'in'")
        normalized_members = tuple(_normalize_null(member) for member in right)
        if left is _NULL or any(member is _NULL for member in normalized_members):
            return False
        for member in normalized_members:
            _require_same_family(left, member, context=context)
        return bool(left in right)  # type: ignore[operator]
    if (left is not _NULL and _runtime_family(left, context=context) != "string") or (
        right is not _NULL and _runtime_family(right, context=context) != "string"
    ):
        _invalid(f"{context} requires strings for operator 'casefold_eq'")
    if left is _NULL or right is _NULL:
        return False
    return cast(str, left).casefold() == cast(str, right).casefold()


def evaluate_rule(rule: CompiledRule, row: Mapping[str, object], *, context: str) -> bool:
    """Evaluate every predicate in a compiled rule against one row.

    Args:
        rule: A rule returned by :func:`validate_rule`.
        row: Column values for one row or comparison record.
        context: Configuration path included in every diagnostic.

    Returns:
        Whether every predicate evaluates true.

    Raises:
        ValueError: If a required row value is missing or incompatible.
    """
    results: list[bool] = []
    for index, predicate in enumerate(rule.predicates):
        results.append(_evaluate_predicate(predicate, row, context=f"{context}.all[{index}]"))
    return all(results)


def adapt_legacy_rule(rule: object, *, exact_rules: Mapping[str, object], context: str) -> object:
    """Migrate only an exact expression shipped immediately before issue #286.

    Args:
        rule: Configured structured rule or historical expression string.
        exact_rules: Exact historical strings mapped to fixed structured rules.
        context: Configuration path included in rejection diagnostics.

    Returns:
        Non-string input unchanged, or a deep copy of an exact migration target.

    Raises:
        ValueError: If a string is not an exact historical expression.
    """
    if not isinstance(rule, str):
        return rule
    if rule in exact_rules:
        return copy.deepcopy(exact_rules[rule])
    _invalid(
        f"{context} uses an unsupported legacy expression. Only the exact expression shipped by VNtyper "
        "immediately before issue #286 can be migrated; use {'all': [{'left': {'column': '...'}, "
        "'operator': 'eq', 'right': {'literal': '...'}}]}."
    )
