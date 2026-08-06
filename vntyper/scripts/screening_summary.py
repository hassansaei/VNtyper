"""
screening_summary.py

Module Purpose:
---------------
Turns the Kestrel and adVNTR result frames plus one quality gate into the
screening state, and looks that state up in ``report_config.json`` to get the
sentence the report prints.

Extracted from ``generate_report.py`` (861 LOC, 4% covered) under AGENTS.md rule
3. It is pure: two DataFrames, two numbers and a config mapping in, a
:class:`ScreeningSummary` out. No filesystem beyond reading the packaged config,
no rendering.

**The clinical text is config-driven and must stay that way** (AGENTS.md). Nothing
here composes, rewords or conditions a sentence; it computes a state and reads the
message for that state out of the configuration. A state with no configured
message falls back to ``screening_summary_default`` -- which is why
``ScreeningSummary.matched_rule`` exists: a report can otherwise announce a
negative screening for a sample two algorithms called positive, with nothing in
the output to show a rule was missing.

The state has three axes:

* ``kestrel_result`` -- one of the ``result`` values under
  ``algorithm_logic.kestrel``, or that block's ``default``.
* ``advntr_result`` -- likewise for ``algorithm_logic.advntr``, plus
  :data:`NOT_PERFORMED` when the adVNTR stage did not run.
* ``quality_metrics_pass`` -- mean VNTR coverage against its threshold.

``is_positive`` is derived from the first two by comparing against each block's
declared ``default``, not by looking for a word in the rendered sentence. The
template used to test ``'negative' not in summary_text`` to decide whether to
style the box as a positive finding, and **not one** of the configured messages
contains that word -- only the fallback default does -- so every configured
message, including "No variant detected by either genotyping method", rendered as
a positive finding.

Functions:
    load_report_config: Read the packaged ``report_config.json``
    compute_algorithm_result: One results frame + one rule block to a state value
    build_screening_summary: The three state axes to the configured message
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

#: The ``advntr_result`` value recorded when the adVNTR stage did not run at all.
#: Distinct from that block's ``default`` ("negative"), which means adVNTR ran and
#: found nothing.
NOT_PERFORMED = "none"

#: Used when ``report_config.json`` declares no default for an algorithm block.
FALLBACK_ALGORITHM_RESULT = "none"

#: Used when ``report_config.json`` declares no ``screening_summary_default``.
FALLBACK_SUMMARY_MESSAGE = "The screening was negative (no valid Kestrel or adVNTR data)."

#: What the report says when the summary could not be computed at all.
UNAVAILABLE_SUMMARY_MESSAGE = "No summary available."


@dataclass(frozen=True)
class ScreeningSummary:
    """The screening state and the sentence configured for it.

    Attributes:
        text: The message to render.
        is_positive: Whether either algorithm reported a finding. Drives the
            template's emphasis; computed from the state, never from ``text``.
        kestrel_result: The computed Kestrel state value.
        advntr_result: The computed adVNTR state value.
        quality_metrics_pass: Whether mean VNTR coverage met its threshold.
        matched_rule: Whether a configured rule matched. False means ``text`` is
            the fallback default and the state has no message of its own.
    """

    text: str
    is_positive: bool
    kestrel_result: str
    advntr_result: str
    quality_metrics_pass: bool
    matched_rule: bool


def load_report_config() -> dict[str, Any]:
    """Load ``report_config.json`` from beside this module.

    Returns:
        dict[str, Any]: The report configuration, or ``{}`` if it cannot be read.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, "report_config.json")
    try:
        with open(config_path) as f:
            report_config = json.load(f)
        logger.info("Loaded report config from %s", config_path)
        return report_config
    except Exception as e:
        logger.error("Failed to load report config: %s", e)
        return {}


def _condition_holds(actual: str, expected: Any) -> bool:
    """Evaluate one configured condition against one cell value.

    Two spellings are supported. A mapping is an explicit
    ``{"operator": ..., "value": ...}`` with operators ``==``, ``!=``, ``in`` and
    ``not in``; anything else is an implicit equality, or membership when it is a
    list.

    Args:
        actual: The cell value, already stringified and stripped.
        expected: The configured expectation.

    Returns:
        bool: Whether the condition holds. An unsupported operator is False.
    """
    if isinstance(expected, dict):
        op = expected.get("operator")
        exp_val = expected.get("value")
        if op == "==":
            return actual == str(exp_val).strip()
        if op == "!=":
            return actual != str(exp_val).strip()
        if op in ("in", "not in"):
            options = exp_val if isinstance(exp_val, list) else [exp_val]
            return (actual in options) if op == "in" else (actual not in options)
        logger.debug("Unsupported operator %r; condition fails.", op)
        return False
    if isinstance(expected, list):
        return actual in expected
    return actual == str(expected)


def compute_algorithm_result(df: pd.DataFrame, logic_config: dict[str, Any]) -> str:
    """Reduce one algorithm's results frame to a single state value.

    Only the first row is examined: the frame is sorted so the most confident
    call leads, and the screening summary describes the sample, not each variant.

    Args:
        df: The results frame. An empty frame yields the configured default.
        logic_config: One block of ``report_config.json``'s ``algorithm_logic``.

    Returns:
        str: The ``result`` of the first matching rule, else the block's
        ``default``.
    """
    default = logic_config.get("default", FALLBACK_ALGORITHM_RESULT)
    if df.empty:
        logger.debug("DataFrame is empty; returning default result %r.", default)
        return default

    row = df.iloc[0]
    logger.debug("Data row for evaluation: %s", row.to_dict())

    for idx, rule in enumerate(logic_config.get("rules", [])):
        conditions = rule.get("conditions", {})
        for col, expected in conditions.items():
            if col not in row:
                logger.debug("Rule %s: column %r not found; rule fails.", idx, col)
                break
            actual = str(row.get(col, "")).strip()
            if not _condition_holds(actual, expected):
                logger.debug("Rule %s: condition on %r not met (actual=%r).", idx, col, actual)
                break
        else:
            result = rule.get("result")
            logger.debug("Rule %s PASSED; returning result: %s", idx, result)
            return result

    logger.debug("No rule matched; returning default result %r.", default)
    return default


def _is_finding(result: str, default: str) -> bool:
    """Whether one algorithm's state value represents a finding.

    Derived from the configured ``default`` rather than a hardcoded list, so a
    renamed ``result`` in ``report_config.json`` cannot silently invert this.
    ``tests/unit/test_screening_summary.py`` asserts every shipped ``result``
    classifies as a finding and every ``default`` does not.

    Args:
        result: The computed state value.
        default: The block's configured ``default``.

    Returns:
        bool: True when the algorithm reported something.
    """
    return result not in (default, NOT_PERFORMED, "")


def _rule_matches(current: dict[str, Any], conditions: dict[str, Any]) -> bool:
    """Whether every condition of one screening rule holds.

    Args:
        current: The three state axes.
        conditions: The rule's ``conditions`` mapping.

    Returns:
        bool: True when every named axis matches.
    """
    for key, rule_value in conditions.items():
        actual = current.get(key)
        if isinstance(rule_value, list):
            if actual not in rule_value:
                return False
        elif actual != rule_value:
            return False
    return True


def build_screening_summary(
    kestrel_df: pd.DataFrame,
    advntr_df: pd.DataFrame,
    advntr_available: bool,
    mean_vntr_coverage: float | None,
    mean_vntr_cov_threshold: float,
    report_config: dict[str, Any],
) -> ScreeningSummary:
    """Compute the screening state and look up its configured message.

    Args:
        kestrel_df: Kestrel results, unformatted (no HTML in the cells).
        advntr_df: adVNTR results.
        advntr_available: Whether the adVNTR stage ran.
        mean_vntr_coverage: Mean coverage over the VNTR region, or None.
        mean_vntr_cov_threshold: The coverage threshold.
        report_config: The parsed ``report_config.json``.

    Returns:
        ScreeningSummary: The state and its message.
    """
    default_message = report_config.get("screening_summary_default", FALLBACK_SUMMARY_MESSAGE)
    try:
        algorithm_logic = report_config.get("algorithm_logic", {})
        kestrel_logic = algorithm_logic.get("kestrel", {})
        advntr_logic = algorithm_logic.get("advntr", {})

        kestrel_result = compute_algorithm_result(kestrel_df, kestrel_logic)
        advntr_result = compute_algorithm_result(advntr_df, advntr_logic) if advntr_available else NOT_PERFORMED
        logger.debug("Computed Kestrel result: %s; adVNTR result: %s", kestrel_result, advntr_result)

        quality_metrics_pass = not (mean_vntr_coverage is not None and mean_vntr_coverage < mean_vntr_cov_threshold)

        current = {
            "kestrel_result": kestrel_result,
            "advntr_result": advntr_result,
            "quality_metrics_pass": quality_metrics_pass,
        }
        logger.debug("Unified screening conditions: %s", current)

        text = ""
        for rule in report_config.get("screening_summary_rules", []):
            if _rule_matches(current, rule.get("conditions", {})):
                text = rule.get("message", "")
                logger.debug("Unified rule matched: %s", rule.get("conditions"))
                break

        matched_rule = bool(text)
        if not matched_rule:
            text = default_message
            logger.warning(
                "No screening rule covers kestrel_result=%r, advntr_result=%r, quality_metrics_pass=%r; "
                "falling back to the default message.",
                kestrel_result,
                advntr_result,
                quality_metrics_pass,
            )

        is_positive = _is_finding(
            kestrel_result, kestrel_logic.get("default", FALLBACK_ALGORITHM_RESULT)
        ) or _is_finding(advntr_result, advntr_logic.get("default", FALLBACK_ALGORITHM_RESULT))
    except Exception as ex:
        logger.error("Exception in build_screening_summary: %s", ex)
        return ScreeningSummary(
            text=UNAVAILABLE_SUMMARY_MESSAGE,
            is_positive=False,
            kestrel_result="",
            advntr_result="",
            quality_metrics_pass=False,
            matched_rule=False,
        )

    logger.debug("Final screening summary: %s", text)
    return ScreeningSummary(
        text=text,
        is_positive=is_positive,
        kestrel_result=kestrel_result,
        advntr_result=advntr_result,
        quality_metrics_pass=quality_metrics_pass,
        matched_rule=matched_rule,
    )
