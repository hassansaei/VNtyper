"""Pure side-specific outcome declarations for the golden-cohort matrix."""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)


def without_side_expectations(case: dict[str, Any]) -> dict[str, Any]:
    """Return a repeat whose processing mode has its own ordinary success contract."""
    independent = dict(case)
    independent.pop("side_expectations", None)
    return independent


def materialize_side_expectation(case: dict[str, Any], side: str) -> dict[str, Any]:
    """Resolve a differential declaration into the legacy admissibility fields.

    Args:
        case: Matrix case, optionally carrying ``side_expectations``.
        side: ``before`` or ``after``.

    Returns:
        A shallow copy with the selected ``expect_exit`` and ``required_artifacts``.

    Raises:
        ValueError: If a differential case does not declare the requested side.
    """
    runtime = dict(case)
    side_expectations = case.get("side_expectations")
    if side_expectations is None:
        return runtime
    if not isinstance(side_expectations, dict) or side not in side_expectations:
        msg = f"Case {case['case_id']} has no {side!r} expectation in side_expectations"
        logger.error(msg)
        raise ValueError(msg)
    selected = side_expectations[side]
    if not isinstance(selected, dict):
        msg = f"Case {case['case_id']} has a malformed {side!r} expectation"
        logger.error(msg)
        raise ValueError(msg)
    if selected.get("expect_exit") not in {"zero", "nonzero"}:
        msg = f"Case {case['case_id']} has an invalid {side!r} expect_exit"
        logger.error(msg)
        raise ValueError(msg)
    required_artifacts = selected.get("required_artifacts")
    if not isinstance(required_artifacts, list) or not all(isinstance(item, str) for item in required_artifacts):
        msg = f"Case {case['case_id']} has invalid {side!r} required_artifacts"
        logger.error(msg)
        raise ValueError(msg)
    if "expected_stderr_contains" in selected:
        expected_stderr = selected["expected_stderr_contains"]
        if not isinstance(expected_stderr, str) or not expected_stderr.strip():
            msg = f"Case {case['case_id']} has invalid {side!r} expected_stderr_contains"
            logger.error(msg)
            raise ValueError(msg)
    runtime.update(selected)
    return runtime
