"""
vntyper/scripts/cohort_rules.py

Module Purpose:
---------------
The interpreter for the ``algorithm_logic`` section of ``report_config.json``.

That section states, as data rather than as code, what makes a Kestrel or an adVNTR
result row positive: a list of rules, each pairing a set of per-column conditions with
the result to return when all of them hold. This module evaluates one result row against
one such rule list and returns the verdict.

Extracted verbatim from ``cohort_summary.py`` (Task 22 of the #181-#197 follow-ups),
which had grown to 911 lines with this - its single largest region - reachable only
through a filesystem walk and a Jinja2 render. Nothing about the evaluation changed in
the move; ``tests/unit/test_cohort_rules.py`` characterises it and
``tests/unit/test_cohort_summary_oracle.py`` pins the report it feeds.

A condition that cannot be evaluated - an absent column, an operator this module does
not implement - fails its rule rather than raising. That is the same silent-disable
failure mode AGENTS.md trap 3 records for the ``eval``-based flagging expressions, and
it is characterised rather than changed.
"""

import logging
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


def compute_algorithm_result(df: pd.DataFrame, logic_config: dict[str, Any]) -> Any:
    """
    Computes the algorithm result (for Kestrel or adVNTR) based on the provided logic configuration.
    Iterates over each rule in logic_config["rules"]. For each condition, compares the plain text value
    from the DataFrame (using the column name) with the expected value.

    Supported operators are: "==", "!=", "in", and "not in". If expected is a list, membership is checked.
    Returns the rule's "result" if matched; otherwise returns logic_config["default"].

    Args:
        df (pandas.DataFrame): DataFrame containing the results (generally a single row).
        logic_config (dict): Configuration dictionary with rules.

    Returns:
        str: The computed algorithm result (e.g., 'High_Precision', 'positive', etc.).
    """
    if df.empty:
        logger.debug("DataFrame is empty; returning default result.")
        return logic_config.get("default", "none")
    row = df.iloc[0]
    logger.debug("Data row for evaluation: %s", row.to_dict())
    logger.debug("Logic configuration: %s", logic_config)
    for idx, rule in enumerate(logic_config.get("rules", [])):
        logger.debug("Evaluating rule %s: %s", idx, rule)
        conditions = rule.get("conditions", {})
        rule_matches = True
        for col, expected in conditions.items():
            if col not in row:
                logger.debug("Rule %s: Column '%s' not found; rule fails.", idx, col)
                rule_matches = False
                break
            actual = str(row.get(col, "")).strip()
            logger.debug(
                "Rule %s, column '%s': actual='%s', expected='%s'",
                idx,
                col,
                actual,
                expected,
            )
            if isinstance(expected, dict):
                op = expected.get("operator")
                exp_val = expected.get("value")
                if op == "==":
                    if actual != str(exp_val).strip():
                        logger.debug(
                            "Rule %s: Condition '%s == %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                elif op == "!=":
                    if actual == str(exp_val).strip():
                        logger.debug(
                            "Rule %s: Condition '%s != %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                elif op == "in":
                    if not isinstance(exp_val, list):
                        exp_val = [exp_val]
                    if actual not in exp_val:
                        logger.debug(
                            "Rule %s: Condition '%s in %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                elif op == "not in":
                    if not isinstance(exp_val, list):
                        exp_val = [exp_val]
                    if actual in exp_val:
                        logger.debug(
                            "Rule %s: Condition '%s not in %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                else:
                    logger.debug(
                        "Rule %s: Unsupported operator '%s' for column '%s'.",
                        idx,
                        op,
                        col,
                    )
                    rule_matches = False
                    break
            else:
                if isinstance(expected, list):
                    if actual not in expected:
                        logger.debug(
                            "Rule %s: Condition '%s in %s' not met (actual='%s').",
                            idx,
                            col,
                            expected,
                            actual,
                        )
                        rule_matches = False
                        break
                else:
                    if actual != str(expected):
                        logger.debug(
                            "Rule %s: Condition '%s == %s' not met (actual='%s').",
                            idx,
                            col,
                            expected,
                            actual,
                        )
                        rule_matches = False
                        break
        if rule_matches:
            result = rule.get("result")
            logger.debug("Rule %s PASSED; returning result: %s", idx, result)
            return result
        else:
            logger.debug("Rule %s did not pass.", idx)
    logger.debug("No rule matched; returning default result.")
    return logic_config.get("default", "none")
