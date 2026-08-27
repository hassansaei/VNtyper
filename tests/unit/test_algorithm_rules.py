"""The one interpreter for ``algorithm_logic``, and its tri-state verdict (spec §2.1).

Legacy behaviour (a config that does NOT declare ``unestablished`` under
``non_finding_results``) is pinned exactly: an absent column fails its rule, a NaN cell
is compared as the string ``"nan"``, an unsupported operator fails its rule. The
tri-state semantics exist only when the configuration opts in, mirroring
``screening_summary.supports_subthreshold`` (#266).
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from vntyper.scripts.algorithm_rules import (
    CONDITION_FALSE,
    CONDITION_TRUE,
    CONDITION_UNEVALUABLE,
    FALLBACK_ALGORITHM_RESULT,
    RULE_FAILED,
    RULE_MATCHED,
    RULE_UNEVALUABLE,
    UNESTABLISHED_RESULT,
    compute_algorithm_result,
    evaluate_rule,
    supports_unestablished,
)

pytestmark = pytest.mark.unit

#: A shape-faithful miniature of the shipped Kestrel block: Confidence listed before
#: Flag, exactly the order that made the legacy interpreter's verdict depend on JSON
#: key order.
LEGACY_LOGIC = {
    "rules": [
        {
            "conditions": {
                "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                "Flag": {"operator": "==", "value": "Not flagged"},
            },
            "result": "High_Precision",
        },
        {
            "conditions": {
                "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                "Flag": {"operator": "!=", "value": "Not flagged"},
            },
            "result": "High_Precision_flagged",
        },
    ],
    "default": "negative",
}

OPTED_IN = {**LEGACY_LOGIC, "non_finding_results": [UNESTABLISHED_RESULT]}


def _row(**cells: Any) -> pd.DataFrame:
    return pd.DataFrame([cells])


# ---- the public constants and opt-in predicate ---------------------------------


def test_the_public_tokens_have_the_declared_values() -> None:
    assert FALLBACK_ALGORITHM_RESULT == "none"
    assert UNESTABLISHED_RESULT == "unestablished"
    assert (RULE_MATCHED, RULE_FAILED, RULE_UNEVALUABLE) == ("matched", "failed", "unevaluable")
    assert (CONDITION_TRUE, CONDITION_FALSE, CONDITION_UNEVALUABLE) == ("true", "false", "unevaluable")


def test_a_config_written_before_this_change_does_not_support_unestablished() -> None:
    assert supports_unestablished(LEGACY_LOGIC) is False
    assert supports_unestablished({}) is False


def test_the_declaration_is_the_opt_in() -> None:
    assert supports_unestablished(OPTED_IN) is True


# ---- legacy semantics, pinned exactly -------------------------------------------


def test_legacy_an_absent_column_fails_the_rule_and_falls_to_the_default() -> None:
    assert compute_algorithm_result(_row(Confidence="High_Precision"), LEGACY_LOGIC) == "negative"


def test_legacy_a_nan_cell_is_compared_as_the_string_nan() -> None:
    frame = _row(Confidence="High_Precision", Flag=float("nan"))
    assert compute_algorithm_result(frame, LEGACY_LOGIC) == "High_Precision_flagged"


def test_legacy_an_empty_frame_returns_the_default_or_fallback() -> None:
    assert compute_algorithm_result(pd.DataFrame(), LEGACY_LOGIC) == "negative"
    assert compute_algorithm_result(pd.DataFrame(), {}) == FALLBACK_ALGORITHM_RESULT


@pytest.mark.parametrize("condition", [{"operator": "~=", "value": "x"}, {"value": "x"}])
def test_legacy_an_unsupported_or_malformed_operator_fails_the_rule(condition: dict[str, Any]) -> None:
    logic = {"rules": [{"conditions": {"C": condition}, "result": "hit"}], "default": "miss"}
    assert compute_algorithm_result(_row(C="x"), logic) == "miss"


def test_legacy_a_malformed_conditions_container_keeps_raising() -> None:
    """The opted-in guard must not change the old interpreter's error surface."""
    logic = {"rules": [{"conditions": ["not", "a", "mapping"], "result": "hit"}], "default": "miss"}
    with pytest.raises(AttributeError, match="items"):
        compute_algorithm_result(_row(C="x"), logic)


@pytest.mark.parametrize(
    "condition,actual,expected",
    [
        ({"operator": "==", "value": "x"}, " x ", "hit"),
        ({"operator": "==", "value": "x"}, "y", "miss"),
        ({"operator": "!=", "value": "x"}, "y", "hit"),
        ({"operator": "!=", "value": "x"}, "x", "miss"),
        ({"operator": "in", "value": ["x", "y"]}, "y", "hit"),
        ({"operator": "in", "value": "x"}, "y", "miss"),
        ({"operator": "not in", "value": ["x", "y"]}, "z", "hit"),
        ({"operator": "not in", "value": "x"}, "x", "miss"),
        (["x", "y"], "y", "hit"),
        (5, 5, "hit"),
    ],
)
def test_legacy_supported_condition_spellings_keep_their_results(
    condition: Any,
    actual: Any,
    expected: str,
) -> None:
    logic = {"rules": [{"conditions": {"C": condition}, "result": "hit"}], "default": "miss"}
    assert compute_algorithm_result(_row(C=actual), logic) == expected


def test_legacy_the_first_matching_rule_and_first_frame_row_win() -> None:
    logic = {
        "rules": [
            {"conditions": {"C": "x"}, "result": "first"},
            {"conditions": {"C": "x"}, "result": "second"},
        ],
        "default": "miss",
    }
    frame = pd.DataFrame([{"C": "x"}, {"C": "not-x"}])
    assert compute_algorithm_result(frame, logic) == "first"


# ---- tri-state semantics, under the opted-in config -----------------------------


def test_an_absent_condition_column_is_unestablished_not_negative() -> None:
    """§2.1 half (a): a High_Precision row with no Flag column."""
    assert compute_algorithm_result(_row(Confidence="High_Precision"), OPTED_IN) == UNESTABLISHED_RESULT


@pytest.mark.parametrize("missing", [float("nan"), None, pd.NA, "   "])
def test_a_missing_cell_is_unestablished_not_a_flagged_positive(missing: Any) -> None:
    """§2.1 half (b): the mixed-cohort union frame's missing Flag false positive."""
    frame = _row(Confidence="High_Precision", Flag=missing)
    assert compute_algorithm_result(frame, OPTED_IN) == UNESTABLISHED_RESULT


def test_the_negative_placeholder_stays_negative_because_confidence_is_definitively_false() -> None:
    """The placeholder's Confidence is definitively false for every rule."""
    assert compute_algorithm_result(_row(Confidence="Negative"), OPTED_IN) == "negative"


def test_the_verdict_does_not_depend_on_condition_order() -> None:
    """No short-circuit: Flag before Confidence must classify identically."""
    reversed_logic = {
        "rules": [
            {
                "conditions": {
                    "Flag": {"operator": "==", "value": "Not flagged"},
                    "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                },
                "result": "High_Precision",
            }
        ],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert compute_algorithm_result(_row(Confidence="Negative"), reversed_logic) == "negative"
    assert compute_algorithm_result(_row(Confidence="Negative"), OPTED_IN) == "negative"


def test_failed_outranks_unevaluable_within_one_rule() -> None:
    logic = {
        "rules": [{"conditions": {"C": "x", "Missing": "y"}, "result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert compute_algorithm_result(_row(C="not-x"), logic) == "negative"


def test_a_later_matched_rule_wins_over_an_earlier_unevaluable_one() -> None:
    logic = {
        "rules": [
            {"conditions": {"Missing": "x"}, "result": "first"},
            {"conditions": {"C": "x"}, "result": "second"},
        ],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert compute_algorithm_result(_row(C="x"), logic) == "second"


@pytest.mark.parametrize("condition", [None, "", "   "])
def test_opted_in_missing_shorthand_expected_value_is_unevaluable(condition: Any) -> None:
    """A missing expectation cannot establish either equality result."""
    logic = {
        "rules": [{"conditions": {"C": condition}, "result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    actual = "None" if condition is None else "x"

    assert evaluate_rule(pd.Series({"C": actual}), {"C": condition}) == RULE_UNEVALUABLE
    assert compute_algorithm_result(_row(C=actual), logic) == UNESTABLISHED_RESULT


@pytest.mark.parametrize(
    "operator,expected,actual",
    [
        ("==", None, "None"),
        ("==", "", "x"),
        ("!=", None, "x"),
        ("!=", "   ", "x"),
        ("in", None, "x"),
        ("in", "", "x"),
        ("not in", None, "x"),
        ("not in", "   ", "x"),
    ],
)
def test_opted_in_explicit_operator_with_missing_expected_value_is_unevaluable(
    operator: str,
    expected: Any,
    actual: str,
) -> None:
    condition = {"operator": operator, "value": expected}
    logic = {
        "rules": [{"conditions": {"C": condition}, "result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }

    assert evaluate_rule(pd.Series({"C": actual}), {"C": condition}) == RULE_UNEVALUABLE
    assert compute_algorithm_result(_row(C=actual), logic) == UNESTABLISHED_RESULT


def test_opted_in_rule_without_conditions_is_unestablished() -> None:
    logic = {
        "rules": [{"result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert compute_algorithm_result(_row(C="x"), logic) == UNESTABLISHED_RESULT


def test_explicit_empty_conditions_mapping_remains_an_unconditional_match() -> None:
    logic = {
        "rules": [{"conditions": {}, "result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert compute_algorithm_result(_row(C="x"), logic) == "hit"


def test_legacy_rule_without_conditions_remains_an_unconditional_match() -> None:
    logic = {"rules": [{"result": "hit"}], "default": "negative"}
    assert compute_algorithm_result(_row(C="x"), logic) == "hit"


def test_an_empty_frame_is_the_default_even_under_the_opted_in_config() -> None:
    """A stage that produced no row established nothing new, preserving today's path."""
    assert compute_algorithm_result(pd.DataFrame(), OPTED_IN) == "negative"


@pytest.mark.parametrize(
    "condition",
    [
        {"operator": "~=", "value": "x"},
        {"value": "x"},
        {"operator": "=="},
        {"operator": None, "value": "x"},
        {"operator": ["=="], "value": "x"},
    ],
)
def test_opted_in_unsupported_or_malformed_conditions_are_unestablished(condition: dict[str, Any]) -> None:
    """A configuration defect must not become a definite negative in opted-in mode."""
    logic = {
        "rules": [{"conditions": {"C": condition}, "result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert evaluate_rule(pd.Series({"C": "x"}), {"C": condition}) == RULE_UNEVALUABLE
    assert compute_algorithm_result(_row(C="x"), logic) == UNESTABLISHED_RESULT


def test_a_malformed_conditions_container_is_unevaluable() -> None:
    row = pd.Series({"C": "x"})
    assert evaluate_rule(row, ["not", "a", "mapping"]) == RULE_UNEVALUABLE  # type: ignore[arg-type]

    logic = {
        "rules": [{"conditions": ["not", "a", "mapping"], "result": "hit"}],
        "default": "negative",
        "non_finding_results": [UNESTABLISHED_RESULT],
    }
    assert compute_algorithm_result(_row(C="x"), logic) == UNESTABLISHED_RESULT


def test_evaluate_rule_returns_the_three_verdicts() -> None:
    row = pd.Series({"C": "x", "E": ""})
    assert evaluate_rule(row, {"C": "x"}) == RULE_MATCHED
    assert evaluate_rule(row, {"C": "y"}) == RULE_FAILED
    assert evaluate_rule(row, {"Missing": "x"}) == RULE_UNEVALUABLE
    assert evaluate_rule(row, {"E": "x"}) == RULE_UNEVALUABLE
    assert evaluate_rule(row, {"C": "y", "Missing": "x"}) == RULE_FAILED


def test_an_empty_conditions_mapping_matches_and_a_non_scalar_cell_can_be_compared() -> None:
    assert evaluate_rule(pd.Series({"C": ["x", "y"]}), {}) == RULE_MATCHED
    assert evaluate_rule(pd.Series({"C": ["x", "y"]}), {"C": "['x', 'y']"}) == RULE_MATCHED


# ---- the deduplication itself is pinned -----------------------------------------


def test_the_two_call_sites_share_one_interpreter() -> None:
    """AGENTS.md extraction rule + spec §2.1: BOTH call sites use the one module."""
    from vntyper.scripts import algorithm_rules, cohort_categories, screening_summary

    assert screening_summary.compute_algorithm_result is algorithm_rules.compute_algorithm_result
    assert cohort_categories.compute_algorithm_result is algorithm_rules.compute_algorithm_result


def test_the_duplicate_cohort_rules_module_is_gone() -> None:
    import importlib.util

    assert importlib.util.find_spec("vntyper.scripts.cohort_rules") is None
