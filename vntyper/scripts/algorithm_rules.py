"""
vntyper/scripts/algorithm_rules.py

Module Purpose:
---------------
The one interpreter for the ``algorithm_logic`` section of ``report_config.json``.

The per-sample screening state and cohort categorisation previously used separate
interpreters. Both treated an absent condition column as a failed rule and a ``NaN``
cell as the string ``"nan"``. That made the same result shape negative per sample but
a flagged positive after cohort concatenation supplied a ``Flag = NaN`` cell.

The tri-state verdict is opt-in via configuration, mirroring
``screening_summary.supports_subthreshold`` (#266). It exists only when the block
declares :data:`UNESTABLISHED_RESULT` under ``non_finding_results``. Configurations
written before this change keep the legacy absent-column, NaN, malformed-condition,
and unsupported-operator behaviour.

Functions:
    supports_unestablished: Whether this block can describe an unestablished result.
    evaluate_rule: Evaluate one row against one rule's conditions.
    compute_algorithm_result: Reduce one result frame and rule block to a state value.
"""

from __future__ import annotations

import logging
from typing import Any, Final

import pandas as pd

logger = logging.getLogger(__name__)

#: Used when a block declares no ``default``.
FALLBACK_ALGORITHM_RESULT: Final[str] = "none"

#: Used when one or more opted-in rules could not be evaluated and none matched.
UNESTABLISHED_RESULT: Final[str] = "unestablished"

#: The three verdicts one rule can reach.
RULE_MATCHED: Final[str] = "matched"
RULE_FAILED: Final[str] = "failed"
RULE_UNEVALUABLE: Final[str] = "unevaluable"

#: The three states one condition can reach.
CONDITION_TRUE: Final[str] = "true"
CONDITION_FALSE: Final[str] = "false"
CONDITION_UNEVALUABLE: Final[str] = "unevaluable"

_SUPPORTED_OPERATORS: Final[frozenset[str]] = frozenset({"==", "!=", "in", "not in"})


def supports_unestablished(logic_config: dict[str, Any]) -> bool:
    """Return whether an algorithm block opts in to tri-state evaluation.

    Args:
        logic_config: One block of ``algorithm_logic``.

    Returns:
        bool: Whether ``unestablished`` is a declared non-finding result.
    """
    return UNESTABLISHED_RESULT in tuple(logic_config.get("non_finding_results", ()))


def _condition_holds(actual: str, expected: Any) -> bool:
    """Compare one stringified cell value using the legacy condition language."""
    if isinstance(expected, dict):
        operator = expected.get("operator")
        expected_value = expected.get("value")
        if operator == "==":
            return actual == str(expected_value).strip()
        if operator == "!=":
            return actual != str(expected_value).strip()
        if operator in ("in", "not in"):
            options = expected_value if isinstance(expected_value, list) else [expected_value]
            return (actual in options) if operator == "in" else (actual not in options)
        logger.debug("Unsupported or malformed operator %r; condition fails under legacy semantics.", operator)
        return False
    if isinstance(expected, list):
        return actual in expected
    return actual == str(expected)


def _condition_is_evaluable(expected: Any) -> bool:
    """Return whether an opted-in condition has a complete, non-missing expectation."""
    if not isinstance(expected, dict):
        return isinstance(expected, list) or not _value_is_missing(expected)
    operator = expected.get("operator")
    if not isinstance(operator, str) or operator not in _SUPPORTED_OPERATORS or "value" not in expected:
        return False
    expected_value = expected["value"]
    return isinstance(expected_value, list) or not _value_is_missing(expected_value)


def _value_is_missing(value: Any) -> bool:
    """Return whether a cell has no scalar value that can be compared."""
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _condition_state(row: pd.Series, column: str, expected: Any) -> str:
    """Evaluate one condition against one row using three-valued semantics."""
    if column not in row or not _condition_is_evaluable(expected):
        return CONDITION_UNEVALUABLE
    value = row[column]
    if _value_is_missing(value):
        return CONDITION_UNEVALUABLE
    actual = str(value).strip()
    return CONDITION_TRUE if _condition_holds(actual, expected) else CONDITION_FALSE


def evaluate_rule(row: pd.Series, conditions: dict[str, Any]) -> str:
    """Evaluate one rule over one row without condition-order dependence.

    Every condition is evaluated. A definitive false condition outranks an
    unevaluable one, so a negative placeholder stays negative even if another column
    is absent. If no condition is false, any absent value, malformed condition, or
    unsupported operator makes the rule unevaluable.

    Args:
        row: One result row.
        conditions: The rule's condition mapping.

    Returns:
        str: ``matched``, ``failed``, or ``unevaluable``.
    """
    if not isinstance(conditions, dict):
        return RULE_UNEVALUABLE

    states = [_condition_state(row, column, expected) for column, expected in conditions.items()]
    if CONDITION_FALSE in states:
        return RULE_FAILED
    if CONDITION_UNEVALUABLE in states:
        return RULE_UNEVALUABLE
    return RULE_MATCHED


def _legacy_rule_verdict(row: pd.Series, conditions: Any) -> str:
    """Evaluate a rule with the exact pre-change fail-on-invalid semantics."""
    for column, expected in conditions.items():
        if column not in row:
            logger.debug("Column %r not found; rule fails under legacy semantics.", column)
            return RULE_FAILED
        actual = str(row.get(column, "")).strip()
        if not _condition_holds(actual, expected):
            logger.debug("Condition on %r not met under legacy semantics (actual=%r).", column, actual)
            return RULE_FAILED
    return RULE_MATCHED


def compute_algorithm_result(df: pd.DataFrame, logic_config: dict[str, Any]) -> Any:
    """Reduce one algorithm's results frame to a single state value.

    The first matching ordered rule wins. In opted-in configurations, an unevaluable
    rule produces ``unestablished`` only after all later rules have had a chance to
    match. An empty frame and configurations without an unevaluable rule return the
    configured default, preserving legacy behaviour.

    Args:
        df: The results frame; only its first row is evaluated.
        logic_config: One block of ``report_config.json``'s ``algorithm_logic``.

    Returns:
        Any: A matching rule result, ``unestablished``, or the configured default.
    """
    default = logic_config.get("default", FALLBACK_ALGORITHM_RESULT)
    if df.empty:
        logger.debug("DataFrame is empty; returning default result %r.", default)
        return default

    row = df.iloc[0]
    logger.debug("Data row for evaluation: %s", row.to_dict())
    tri_state = supports_unestablished(logic_config)
    any_unevaluable = False

    for index, rule in enumerate(logic_config.get("rules", [])):
        conditions = rule.get("conditions", {})
        if tri_state and "conditions" not in rule:
            verdict = RULE_UNEVALUABLE
        else:
            verdict = evaluate_rule(row, conditions) if tri_state else _legacy_rule_verdict(row, conditions)
        if verdict == RULE_MATCHED:
            result = rule.get("result")
            logger.debug("Rule %s matched; returning result %r.", index, result)
            return result
        if verdict == RULE_UNEVALUABLE:
            logger.debug("Rule %s was unevaluable.", index)
            any_unevaluable = True

    if tri_state and any_unevaluable:
        logger.debug("No rule matched and at least one was unevaluable; returning %r.", UNESTABLISHED_RESULT)
        return UNESTABLISHED_RESULT
    logger.debug("No rule matched; returning default result %r.", default)
    return default
