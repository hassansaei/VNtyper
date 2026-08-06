"""The config-driven rule engine that turns one result row into an algorithm verdict.

`report_config.json` states, as data, what makes a Kestrel or adVNTR row positive:
a list of rules, each a set of column conditions and the result to return when they all
hold. `compute_algorithm_result` is the interpreter for that little language, and it was
the single largest uncovered region of `cohort_summary.py` - 130 lines carrying four
operators, of which only the `in`/`==` pair the shipped config happens to use was ever
executed by a test.

Everything here is **characterisation**: it records what the interpreter does today and
claims none of it is intended. The rules it interprets decide whether a sample is
reported positive, so a silent change in any operator is a changed diagnosis; that is
what these pin.
"""

from __future__ import annotations

import pandas as pd
import pytest

from vntyper.scripts.cohort_rules import compute_algorithm_result

pytestmark = pytest.mark.unit


def _row(**values: object) -> pd.DataFrame:
    """Build the single-row frame the interpreter is always called with.

    Args:
        **values: The row's columns.

    Returns:
        pd.DataFrame: A one-row frame.
    """
    return pd.DataFrame([values])


# ---------------------------------------------------------------------------
# No rule to apply
# ---------------------------------------------------------------------------


def test_an_empty_frame_returns_the_configured_default() -> None:
    assert compute_algorithm_result(pd.DataFrame(), {"default": "negative"}) == "negative"


def test_an_empty_frame_with_no_configured_default_returns_none_as_a_string() -> None:
    """The fallback default is the string `"none"`, not `None`."""
    assert compute_algorithm_result(pd.DataFrame(), {}) == "none"


def test_a_row_matching_no_rule_returns_the_configured_default() -> None:
    logic = {"rules": [{"conditions": {"Flag": "Not flagged"}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(Flag="flagged"), logic) == "negative"


def test_a_config_with_no_rules_at_all_returns_the_default() -> None:
    assert compute_algorithm_result(_row(Flag="x"), {"default": "negative"}) == "negative"


# ---------------------------------------------------------------------------
# Which row, and which rule
# ---------------------------------------------------------------------------


def test_only_the_first_row_is_evaluated() -> None:
    """The interpreter is documented as taking "generally a single row" and takes
    `iloc[0]` unconditionally; rows after the first are invisible to it."""
    frame = pd.DataFrame([{"Flag": "no"}, {"Flag": "yes"}])
    logic = {"rules": [{"conditions": {"Flag": "yes"}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(frame, logic) == "negative"


def test_the_first_matching_rule_wins() -> None:
    logic = {
        "rules": [
            {"conditions": {"Flag": "x"}, "result": "first"},
            {"conditions": {"Flag": "x"}, "result": "second"},
        ],
        "default": "negative",
    }

    assert compute_algorithm_result(_row(Flag="x"), logic) == "first"


def test_a_matching_rule_with_no_result_key_returns_none() -> None:
    """A rule that matches but declares no `result` yields Python `None`, which then
    flows on into the report as a category of its own."""
    logic = {"rules": [{"conditions": {"Flag": "x"}}], "default": "negative"}

    assert compute_algorithm_result(_row(Flag="x"), logic) is None


def test_a_rule_with_no_conditions_matches_everything() -> None:
    """An empty condition set is vacuously true, so such a rule is an unconditional
    result and shadows every rule after it."""
    logic = {"rules": [{"conditions": {}, "result": "always"}], "default": "negative"}

    assert compute_algorithm_result(_row(Flag="anything"), logic) == "always"


def test_a_condition_on_an_absent_column_fails_the_rule() -> None:
    """A renamed upstream column silently stops matching rather than raising."""
    logic = {"rules": [{"conditions": {"Missing": "x"}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(Flag="x"), logic) == "negative"


def test_every_condition_of_a_rule_must_hold() -> None:
    logic = {
        "rules": [{"conditions": {"Flag": "x", "Confidence": "High"}, "result": "positive"}],
        "default": "negative",
    }

    assert compute_algorithm_result(_row(Flag="x", Confidence="Low"), logic) == "negative"
    assert compute_algorithm_result(_row(Flag="x", Confidence="High"), logic) == "positive"


# ---------------------------------------------------------------------------
# The four operators
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "operator,expected,actual,matches",
    [
        ("==", "Not flagged", "Not flagged", True),
        ("==", "Not flagged", "flagged", False),
        ("!=", "Not flagged", "flagged", True),
        ("!=", "Not flagged", "Not flagged", False),
        ("in", ["High_Precision", "High_Precision*"], "High_Precision*", True),
        ("in", ["High_Precision"], "Low_Precision", False),
        ("not in", ["Not flagged", "None"], "Depth_flagged", True),
        ("not in", ["Not flagged", "None"], "Not flagged", False),
    ],
)
def test_each_operator_decides_a_rule(operator: str, expected: object, actual: str, matches: bool) -> None:
    logic = {
        "rules": [{"conditions": {"C": {"operator": operator, "value": expected}}, "result": "positive"}],
        "default": "negative",
    }

    assert compute_algorithm_result(_row(C=actual), logic) == ("positive" if matches else "negative")


@pytest.mark.parametrize("operator", ["in", "not in"])
def test_a_scalar_value_is_promoted_to_a_one_element_list(operator: str) -> None:
    """`in` against a bare string would otherwise be a substring test."""
    logic = {
        "rules": [{"conditions": {"C": {"operator": operator, "value": "High"}}, "result": "positive"}],
        "default": "negative",
    }

    assert compute_algorithm_result(_row(C="Hi"), logic) == ("negative" if operator == "in" else "positive")


def test_an_unsupported_operator_fails_the_rule_rather_than_raising() -> None:
    """A typo such as `"eq"` disables the rule silently - the failure mode AGENTS.md
    trap 3 describes for the `eval`-based flagging expressions, in a second place."""
    logic = {
        "rules": [{"conditions": {"C": {"operator": "eq", "value": "x"}}, "result": "positive"}],
        "default": "negative",
    }

    assert compute_algorithm_result(_row(C="x"), logic) == "negative"


def test_a_condition_dict_with_no_operator_fails_the_rule() -> None:
    logic = {"rules": [{"conditions": {"C": {"value": "x"}}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(C="x"), logic) == "negative"


# ---------------------------------------------------------------------------
# Bare (non-dict) condition values
# ---------------------------------------------------------------------------


def test_a_bare_string_condition_is_an_equality_test() -> None:
    logic = {"rules": [{"conditions": {"C": "x"}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(C="x"), logic) == "positive"
    assert compute_algorithm_result(_row(C="y"), logic) == "negative"


def test_a_bare_list_condition_is_a_membership_test() -> None:
    logic = {"rules": [{"conditions": {"C": ["x", "y"]}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(C="y"), logic) == "positive"
    assert compute_algorithm_result(_row(C="z"), logic) == "negative"


# ---------------------------------------------------------------------------
# How cell values are coerced before comparison
# ---------------------------------------------------------------------------


def test_the_cell_value_is_stringified_and_stripped() -> None:
    """Every `parsed_result` value is a string by the time the cohort reads it
    (AGENTS.md trap 8), but the interpreter coerces anyway, and it strips."""
    logic = {"rules": [{"conditions": {"C": "5"}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(C="  5  "), logic) == "positive"
    assert compute_algorithm_result(_row(C=5), logic) == "positive"


def test_a_missing_value_compares_as_the_string_nan() -> None:
    """`str(NaN)` is `"nan"`, so a null cell is not equal to `""` and not equal to
    anything a config is likely to declare."""
    logic = {"rules": [{"conditions": {"C": "nan"}, "result": "positive"}], "default": "negative"}

    assert compute_algorithm_result(_row(C=float("nan")), logic) == "positive"


def test_the_equality_operators_stringify_the_expected_value_but_the_membership_ones_do_not() -> None:
    """Characterisation of an asymmetry, not of an intention.

    `==` and `!=` compare against `str(expected).strip()`, so a numeric literal in the
    config works. `in` and `not in` compare the stringified cell against the raw list,
    so the same numeric literal never matches. The shipped `report_config.json` uses
    only strings, so nothing hits this today; it is pinned because a future config that
    writes `"value": [5]` would be disabled with no error.
    """
    equality = {"rules": [{"conditions": {"C": {"operator": "==", "value": 5}}, "result": "hit"}], "default": "miss"}
    membership = {
        "rules": [{"conditions": {"C": {"operator": "in", "value": [5]}}, "result": "hit"}],
        "default": "miss",
    }

    assert compute_algorithm_result(_row(C="5"), equality) == "hit"
    assert compute_algorithm_result(_row(C="5"), membership) == "miss"


# ---------------------------------------------------------------------------
# The shipped configuration, end to end
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "confidence,flag,expected",
    [
        ("High_Precision", "Not flagged", "High_Precision"),
        ("High_Precision*", "Not flagged", "High_Precision"),
        ("Low_Precision", "Not flagged", "Low_Precision"),
        ("High_Precision", "Depth_Score_flagged", "High_Precision_flagged"),
        ("Low_Precision", "Depth_Score_flagged", "Low_Precision_flagged"),
        ("Negative", "Not flagged", "negative"),
    ],
)
def test_the_shipped_kestrel_rules_produce_the_recorded_verdicts(confidence: str, flag: str, expected: str) -> None:
    """Read out of `report_config.json` itself so a config edit fails here."""
    from vntyper.scripts.cohort_summary import load_report_config

    logic = load_report_config()["algorithm_logic"]["kestrel"]

    assert compute_algorithm_result(_row(Confidence=confidence, Flag=flag), logic) == expected


@pytest.mark.parametrize(
    "vid,flag,expected",
    [
        ("25561", "Not flagged", "positive"),
        ("Negative", "Not flagged", "negative"),
        ("25561", "Low_Coverage", "positive flagged"),
        ("Negative", "Low_Coverage", "negative"),
        ("Negative", "Not applicable", "negative"),
        ("Negative", "None", "negative"),
    ],
)
def test_the_shipped_advntr_rules_produce_the_recorded_verdicts(vid: str, flag: str, expected: str) -> None:
    """The verdicts the shipped rules produce, read out of `report_config.json` itself.

    The fourth case is the one to watch. Rule 2 (`"positive flagged"`) once guarded on
    the flag alone while rule 1 (`"positive"`) guarded on both `VID != "Negative"` and
    the flag, so on paper a `VID == "Negative"` row carrying any flag outside the three
    excused values came out `"positive flagged"`. Rule 2 now carries the same `VID`
    guard, and that case is `"negative"`.

    The asymmetry was never reachable: both sentinel-row constructions set
    `Flag = "Not applicable"`, which rule 2 excludes; neither is routed through
    `add_flags`; and the shipped flagging rules cannot fire on sentinel values anyway.
    So closing it moved no verdict the pipeline can produce - `tests/unit/test_advntr_rule_symmetry.py`
    owns that property and enumerates all five reachable `(VID, Flag)` pairs. What this
    parametrisation owns is the interpreter's output over the shipped table, case by
    case. See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-advntr-vid-column-overloaded.md`.

    The flag strings here are real ones: `Low_Coverage` is a rule name from
    `advntr_config.json`'s `flagging_rules`, which is what `add_flags` writes into the
    column.
    """
    from vntyper.scripts.cohort_summary import load_report_config

    logic = load_report_config()["algorithm_logic"]["advntr"]

    assert compute_algorithm_result(_row(VID=vid, Flag=flag), logic) == expected
