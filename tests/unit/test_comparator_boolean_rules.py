"""Tests for recursively nested boolean comparator-rule nodes."""

from __future__ import annotations

import logging
from collections.abc import Mapping

import pytest

from vntyper.scripts.comparator_rules import evaluate_rule, validate_rule

pytestmark = pytest.mark.unit


def predicate(column: str, operator: str, literal: object) -> dict[str, object]:
    """Return one configured predicate leaf."""
    return {
        "left": {"column": column},
        "operator": operator,
        "right": {"literal": literal},
    }


@pytest.mark.parametrize(
    ("configured", "row", "expected"),
    [
        (
            {
                "all": [
                    predicate("REF", "eq", "C"),
                    {
                        "any": [
                            predicate("Depth", "lt", 10),
                            {"not": predicate("Kind", "casefold_eq", "deletion")},
                        ]
                    },
                ]
            },
            {"REF": "C", "Depth": 12, "Kind": "Insertion"},
            True,
        ),
        (
            {"any": [predicate("REF", "eq", "A"), {"all": [predicate("Depth", "lt", 10)]}]},
            {"REF": "C", "Depth": 9},
            True,
        ),
        (
            {"not": {"any": [predicate("REF", "eq", "A"), predicate("REF", "eq", "C")]}},
            {"REF": "C"},
            False,
        ),
    ],
)
def test_nested_all_any_not_semantics(configured: object, row: Mapping[str, object], expected: bool) -> None:
    compiled = validate_rule(configured, allowed_columns=row, context="test.rule")

    assert evaluate_rule(compiled, row, context="test.rule") is expected


@pytest.mark.parametrize(
    ("configured", "expected_path"),
    [
        ({"any": []}, "test.rule.any"),
        ({"all": []}, "test.rule.all"),
        ({"not": []}, "test.rule.not"),
        ({"not": {}}, "test.rule.not"),
        ({"not": predicate("REF", "eq", "C"), "extra": True}, "test.rule"),
        ({"all": [predicate("REF", "eq", "C")], "any": []}, "test.rule"),
        ({"xor": [predicate("REF", "eq", "C")]}, "test.rule"),
        (predicate("REF", "eq", "C"), "test.rule"),
        ({"all": [{"left": {"column": "REF"}, "operator": "eq"}]}, "test.rule.all[0]"),
    ],
)
def test_malformed_boolean_nodes_fail_closed_with_paths(
    configured: object,
    expected_path: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.comparator_rules")

    with pytest.raises(ValueError) as raised:
        validate_rule(configured, allowed_columns=["REF"], context="test.rule")

    assert expected_path in str(raised.value)
    assert caplog.messages[-1] == str(raised.value)


def test_complete_nested_tree_is_validated_before_evaluation() -> None:
    configured = {
        "any": [
            predicate("REF", "eq", "C"),
            {"not": predicate("UNDECLARED", "eq", "unsafe")},
        ]
    }

    with pytest.raises(ValueError) as raised:
        validate_rule(configured, allowed_columns=["REF"], context="test.rule")

    assert str(raised.value).startswith("test.rule.any[1].not.left.column names unknown column 'UNDECLARED'")


@pytest.mark.parametrize(
    "configured",
    [
        {"all": [predicate("A", "eq", "no"), predicate("B", "eq", "yes")]},
        {"any": [predicate("A", "eq", "yes"), predicate("B", "eq", "yes")]},
    ],
)
def test_boolean_nodes_evaluate_every_child_without_short_circuiting(configured: object) -> None:
    compiled = validate_rule(configured, allowed_columns=["A", "B"], context="test.rule")

    with pytest.raises(ValueError) as raised:
        evaluate_rule(compiled, {"A": "yes"}, context="test.rule")

    assert "cannot evaluate missing row column 'B'" in str(raised.value)


def test_nested_configuration_is_compiled_to_immutable_values() -> None:
    configured = {"not": {"any": [predicate("REF", "eq", "C")]}}
    compiled = validate_rule(configured, allowed_columns=["REF"], context="test.rule")

    configured["not"]["any"].clear()  # type: ignore[index,union-attr]

    assert evaluate_rule(compiled, {"REF": "A"}, context="test.rule") is True


def test_excessive_boolean_nesting_is_rejected_deterministically() -> None:
    configured: object = predicate("REF", "eq", "C")
    for _index in range(34):
        configured = {"not": configured}

    with pytest.raises(ValueError) as raised:
        validate_rule(configured, allowed_columns=["REF"], context="test.rule")

    assert str(raised.value).endswith("exceeds maximum boolean nesting of 32")


def test_maximum_boolean_nesting_is_accepted() -> None:
    configured: object = predicate("REF", "eq", "C")
    for _index in range(32):
        configured = {"not": configured}

    compiled = validate_rule(configured, allowed_columns=["REF"], context="test.rule")

    assert evaluate_rule(compiled, {"REF": "C"}, context="test.rule") is True
