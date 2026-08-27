"""Tests for the validated, non-executable comparator rule language."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import TypeAlias, cast

import numpy as np
import pandas as pd
import pytest

from vntyper.scripts.comparator_rules import adapt_legacy_rule, evaluate_rule, validate_rule

pytestmark = pytest.mark.unit

RuleConfig: TypeAlias = dict[str, list[dict[str, object]]]


def column(name: str) -> dict[str, str]:
    return {"column": name}


def literal(value: object) -> dict[str, object]:
    return {"literal": value}


def rule(left: Mapping[str, object], operator: str, right: Mapping[str, object]) -> RuleConfig:
    return {"all": [{"left": dict(left), "operator": operator, "right": dict(right)}]}


@pytest.mark.parametrize(
    ("configured", "row", "expected"),
    [
        (rule(column("REF"), "eq", literal("C")), {"REF": "C"}, True),
        (rule(literal("C"), "eq", column("REF")), {"REF": "A"}, False),
        (rule(column("DP"), "lt", literal(10)), {"DP": 9}, True),
        (rule(literal(10), "lt", column("DP")), {"DP": 11}, True),
        (rule(column("Motif"), "in", literal(["1", "2"])), {"Motif": "2"}, True),
        (rule(column("Kind"), "casefold_eq", literal("insertion")), {"Kind": "Insertion"}, True),
        (rule(literal("INSERTION"), "casefold_eq", column("Kind")), {"Kind": "deletion"}, False),
        (rule(literal("2"), "in", literal(["1", "2"])), {}, True),
        (rule(column("Boolean"), "in", literal([True, False])), {"Boolean": True}, True),
        (rule(literal(True), "eq", literal(True)), {}, True),
        (rule(literal(1), "eq", literal(1.0)), {}, True),
        (rule(literal(None), "eq", literal("x")), {}, False),
        (rule(literal("x"), "eq", literal(None)), {}, False),
    ],
)
def test_supported_operators(configured: object, row: dict[str, object], expected: bool) -> None:
    compiled = validate_rule(configured, allowed_columns=row, context="test.rule")
    assert evaluate_rule(compiled, row, context="test.rule") is expected


@pytest.mark.parametrize(
    ("configured", "allowed_columns", "expected_path"),
    [
        (None, ["REF"], "flagging_rules.Low_Coverage"),
        ({"all": [], "extra": True}, ["REF"], "flagging_rules.Low_Coverage"),
        ({"all": ()}, ["REF"], "flagging_rules.Low_Coverage.all"),
        ({"all": []}, ["REF"], "flagging_rules.Low_Coverage.all"),
        ({"all": ["REF == 'C'"]}, ["REF"], "flagging_rules.Low_Coverage.all[0]"),
        (
            {"all": [{"left": column("REF"), "right": literal("C")}]},
            ["REF"],
            "flagging_rules.Low_Coverage.all[0]",
        ),
        (
            {"all": [{"left": column("REF"), "operator": "eq", "right": literal("C"), "extra": 1}]},
            ["REF"],
            "flagging_rules.Low_Coverage.all[0]",
        ),
        ({"all": [{"all": []}]}, ["REF"], "flagging_rules.Low_Coverage.all[0]"),
        ({"all": [{"any": []}]}, ["REF"], "flagging_rules.Low_Coverage.all[0]"),
        ({"all": [{"or": []}]}, ["REF"], "flagging_rules.Low_Coverage.all[0]"),
        ({"all": [{"not": {}}]}, ["REF"], "flagging_rules.Low_Coverage.all[0]"),
        (rule({}, "eq", literal("C")), ["REF"], "flagging_rules.Low_Coverage.all[0].left"),
        (
            rule({"column": "REF", "literal": "C"}, "eq", literal("C")),
            ["REF"],
            "flagging_rules.Low_Coverage.all[0].left",
        ),
        (
            rule({"column": "REF", "extra": True}, "eq", literal("C")),
            ["REF"],
            "flagging_rules.Low_Coverage.all[0].left",
        ),
        (rule(column(""), "eq", literal("C")), ["REF"], "flagging_rules.Low_Coverage.all[0].left.column"),
        (
            rule({"column": 7}, "eq", literal("C")),
            ["REF"],
            "flagging_rules.Low_Coverage.all[0].left.column",
        ),
        (
            rule(column("MISSING"), "eq", literal("C")),
            ["REF", "ALT"],
            "flagging_rules.Low_Coverage.all[0].left.column",
        ),
        (rule(column("REF"), "==", literal("C")), ["REF"], "flagging_rules.Low_Coverage.all[0].operator"),
        (
            {"all": [{"left": column("REF"), "operator": 7, "right": literal("C")}]},
            ["REF"],
            "flagging_rules.Low_Coverage.all[0].operator",
        ),
        (rule(column("Motif"), "in", literal("1")), ["Motif"], "flagging_rules.Low_Coverage.all[0].right"),
        (
            rule(column("Motif"), "in", column("Allowed")),
            ["Allowed", "Motif"],
            "flagging_rules.Low_Coverage.all[0].right",
        ),
        (rule(column("Motif"), "in", literal([])), ["Motif"], "flagging_rules.Low_Coverage.all[0].right.literal"),
        (
            rule(column("Motif"), "in", literal(["1", 2])),
            ["Motif"],
            "flagging_rules.Low_Coverage.all[0].right.literal",
        ),
        (
            rule(column("Motif"), "in", literal([["1"]])),
            ["Motif"],
            "flagging_rules.Low_Coverage.all[0].right.literal[0]",
        ),
        (rule(column("REF"), "eq", literal(["C"])), ["REF"], "flagging_rules.Low_Coverage.all[0].right"),
        (rule(literal(True), "lt", column("DP")), ["DP"], "flagging_rules.Low_Coverage.all[0].left.literal"),
        (rule(literal("C"), "eq", literal(3)), [], "flagging_rules.Low_Coverage.all[0]"),
        (rule(literal("1"), "in", literal([1, 2])), [], "flagging_rules.Low_Coverage.all[0]"),
        (rule(literal("9"), "lt", literal(10)), [], "flagging_rules.Low_Coverage.all[0].left.literal"),
        (rule(literal(9), "casefold_eq", literal("9")), [], "flagging_rules.Low_Coverage.all[0].left.literal"),
    ],
)
def test_invalid_schema_is_rejected_with_logged_context_and_path(
    configured: object,
    allowed_columns: list[str],
    expected_path: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.comparator_rules")

    with pytest.raises(ValueError) as raised:
        validate_rule(configured, allowed_columns=allowed_columns, context="flagging_rules.Low_Coverage")

    message = str(raised.value)
    assert expected_path in message
    assert caplog.messages[-1] == message


def test_unknown_column_message_sorts_the_allowlist() -> None:
    with pytest.raises(ValueError) as raised:
        validate_rule(
            rule(column("MISSING"), "eq", literal("C")),
            allowed_columns=["REF", "ALT"],
            context="test.rule",
        )

    assert "['ALT', 'REF']" in str(raised.value)


@pytest.mark.parametrize(
    "payload",
    [
        "__import__('os').system('touch {marker}')",
        "obj.attr",
        "obj[0]",
        "[x for x in y]",
        "(lambda: 1)()",
        "open('/etc/passwd')",
    ],
)
def test_code_shaped_literal_is_only_data(payload: str, tmp_path: Path) -> None:
    marker = tmp_path / "executed"
    value = payload.format(marker=marker)
    compiled = validate_rule(
        rule(column("Value"), "eq", literal(value)), allowed_columns=["Value"], context="test.rule"
    )

    assert evaluate_rule(compiled, {"Value": value}, context="test.rule") is True
    assert not marker.exists()


@pytest.mark.parametrize(
    "payload",
    [
        "__import__('os').system('touch {marker}')",
        "obj.attr",
        "obj[0]",
        "[x for x in y]",
        "(lambda: 1)()",
        "open('/etc/passwd')",
    ],
)
def test_top_level_code_shaped_string_is_rejected_without_execution(payload: str, tmp_path: Path) -> None:
    marker = tmp_path / "executed"

    with pytest.raises(ValueError, match="test.rule"):
        validate_rule(payload.format(marker=marker), allowed_columns=["Value"], context="test.rule")

    assert not marker.exists()


def test_legacy_adapter_returns_non_strings_unchanged() -> None:
    configured = rule(column("REF"), "eq", literal("C"))

    assert adapt_legacy_rule(configured, exact_rules={}, context="test.rule") is configured


def test_legacy_adapter_deep_copies_an_exact_shipped_expression() -> None:
    migrated = rule(column("REF"), "eq", literal("C"))
    exact_rules = {"REF == 'C'": migrated}

    adapted = adapt_legacy_rule("REF == 'C'", exact_rules=exact_rules, context="test.rule")

    assert adapted == migrated
    assert adapted is not migrated
    assert isinstance(adapted, dict)
    adapted_predicates = cast(list[dict[str, object]], adapted["all"])
    adapted_right = cast(dict[str, object], adapted_predicates[0]["right"])
    adapted_right["literal"] = "A"
    migrated_right = cast(dict[str, object], migrated["all"][0]["right"])
    assert migrated_right["literal"] == "C"


@pytest.mark.parametrize("configured", [" REF == 'C'", "REF == 'C' ", "REF  == 'C'", "REF == 'A'"])
def test_legacy_adapter_rejects_every_non_exact_string(configured: str, caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.comparator_rules")

    with pytest.raises(ValueError) as raised:
        adapt_legacy_rule(configured, exact_rules={"REF == 'C'": {}}, context="flagging_rules.Reference")

    message = str(raised.value)
    assert message == (
        "flagging_rules.Reference uses an unsupported legacy expression. Only the exact expression shipped by "
        "VNtyper immediately before issue #286 can be migrated; use {'all': [{'left': {'column': '...'}, "
        "'operator': 'eq', 'right': {'literal': '...'}}]}."
    )
    assert caplog.messages[-1] == message


@pytest.mark.parametrize("null_value", [None, pd.NA, float("nan"), np.float64("nan")])
@pytest.mark.parametrize(
    ("operator", "right"),
    [
        ("eq", literal("x")),
        ("lt", literal(10)),
        ("in", literal(["x", "y"])),
        ("casefold_eq", literal("x")),
    ],
)
def test_null_row_values_make_every_operator_false(null_value: object, operator: str, right: dict[str, object]) -> None:
    compiled = validate_rule(rule(column("Value"), operator, right), allowed_columns=["Value"], context="test.rule")

    assert evaluate_rule(compiled, {"Value": null_value}, context="test.rule") is False


@pytest.mark.parametrize("null_value", [None, pd.NA, float("nan"), np.float64("nan")])
def test_null_row_values_are_false_on_the_right_side(null_value: object) -> None:
    compiled = validate_rule(rule(literal("x"), "eq", column("Value")), allowed_columns=["Value"], context="test.rule")

    assert evaluate_rule(compiled, {"Value": null_value}, context="test.rule") is False


@pytest.mark.parametrize(
    ("configured", "row", "message_fragment"),
    [
        (rule(column("DP"), "lt", literal(10)), {"DP": "9"}, "real numbers"),
        (rule(column("DP"), "lt", literal(10)), {"DP": True}, "real numbers"),
        (rule(column("Kind"), "casefold_eq", literal("insertion")), {"Kind": 1}, "strings"),
        (rule(column("Value"), "eq", literal("1")), {"Value": 1}, "compatible families"),
        (rule(column("Value"), "in", literal(["1", "2"])), {"Value": 1}, "compatible families"),
        (rule(column("Value"), "eq", literal("1")), {"Value": ["1"]}, "JSON-scalar row values"),
    ],
)
def test_incompatible_non_null_row_values_raise_logged_value_error(
    configured: object,
    row: dict[str, object],
    message_fragment: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.comparator_rules")
    compiled = validate_rule(configured, allowed_columns=row, context="test.rule")

    with pytest.raises(ValueError) as raised:
        evaluate_rule(compiled, row, context="test.rule")

    message = str(raised.value)
    assert "test.rule.all[0]" in message
    assert message_fragment in message
    assert caplog.messages[-1] == message


def test_missing_row_column_raises_even_when_it_was_compile_time_allowlisted(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.comparator_rules")
    compiled = validate_rule(rule(column("REF"), "eq", literal("C")), allowed_columns=["REF"], context="test.rule")

    with pytest.raises(ValueError) as raised:
        evaluate_rule(compiled, {}, context="test.rule")

    message = str(raised.value)
    assert message == "test.rule.all[0] cannot evaluate missing row column 'REF'"
    assert caplog.messages[-1] == message


def test_a_false_predicate_does_not_hide_a_later_incompatible_value() -> None:
    configured = {
        "all": [
            {"left": column("REF"), "operator": "eq", "right": literal("C")},
            {"left": column("DP"), "operator": "lt", "right": literal(10)},
        ]
    }
    compiled = validate_rule(configured, allowed_columns=["REF", "DP"], context="test.rule")

    with pytest.raises(ValueError, match=r"test\.rule\.all\[1\].*real numbers"):
        evaluate_rule(compiled, {"REF": "A", "DP": "9"}, context="test.rule")


def test_all_requires_every_predicate_to_be_true() -> None:
    configured = {
        "all": [
            {"left": column("REF"), "operator": "eq", "right": literal("A")},
            {"left": column("ALT"), "operator": "eq", "right": literal("T")},
        ]
    }
    compiled = validate_rule(configured, allowed_columns=["REF", "ALT"], context="test.rule")

    assert evaluate_rule(compiled, {"REF": "C", "ALT": "T"}, context="test.rule") is False


@pytest.mark.parametrize(
    ("configured", "row"),
    [
        (rule(column("Kind"), "casefold_eq", literal("INSERTION")), {"Kind": "Insertion"}),
        (rule(literal("Insertion"), "casefold_eq", column("Kind")), {"Kind": "INSERTION"}),
    ],
)
def test_casefold_is_applied_to_both_operands(configured: object, row: dict[str, object]) -> None:
    compiled = validate_rule(configured, allowed_columns=row, context="test.rule")

    assert evaluate_rule(compiled, row, context="test.rule") is True


def test_null_membership_collection_makes_the_predicate_false() -> None:
    compiled = validate_rule(
        rule(column("Value"), "in", literal([None])), allowed_columns=["Value"], context="test.rule"
    )

    assert evaluate_rule(compiled, {"Value": "anything"}, context="test.rule") is False


def test_compilation_copies_mutable_configuration_into_immutable_values() -> None:
    configured_values = ["1", "2"]
    configured = rule(column("Motif"), "in", literal(configured_values))
    compiled = validate_rule(configured, allowed_columns=["Motif"], context="test.rule")

    configured_values.append("3")
    configured_left = cast(dict[str, object], configured["all"][0]["left"])
    configured_left["column"] = "Changed"
    configured["all"].clear()

    assert evaluate_rule(compiled, {"Motif": "2"}, context="test.rule") is True
    assert evaluate_rule(compiled, {"Motif": "3"}, context="test.rule") is False


def test_validation_reaches_an_unknown_column_in_a_later_predicate() -> None:
    configured = {
        "all": [
            {"left": column("REF"), "operator": "eq", "right": literal("A")},
            {"left": column("MISSING"), "operator": "eq", "right": literal("T")},
        ]
    }

    with pytest.raises(ValueError, match=r"test\.rule\.all\[1\]\.left\.column"):
        validate_rule(configured, allowed_columns=["REF"], context="test.rule")
