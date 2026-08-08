"""Pure side-specific outcome declarations for the golden-cohort matrix."""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)

# Measured by the milestone-four full 50-BAM candidate run. The baseline exits zero for
# these inputs only because samtools's stranded FASTQ was ignored; the candidate names and
# rejects every produced file instead. Keeping the ids explicit makes this a reviewable
# gate policy rather than an inference from sample names.
MIXED_LAYOUT_BASE_CASE_IDS: frozenset[str] = frozenset(
    {
        "40cf_hg38_subset",
        "6449_hg19_subset",
        "66bf_GRCh37_bwa",
        "66bf_GRCh38_bwa",
        "66bf_hg19_bwa",
        "66bf_hg19_ensembl_bwa",
        "66bf_hg19_subset",
        "66bf_hg38_bwa",
        "66bf_hg38_ensembl_bwa",
        "6c28_hg19_subset",
        "7a61_GRCh37_bwa",
        "7a61_GRCh38_bwa",
        "7a61_hg19_bwa",
        "7a61_hg19_ensembl_bwa",
        "7a61_hg19_subset",
        "7a61_hg38_bwa",
        "7a61_hg38_ensembl_bwa",
        "a5c1_GRCh37_bwa",
        "a5c1_GRCh38_bwa",
        "a5c1_hg19_bwa",
        "a5c1_hg19_ensembl_bwa",
        "a5c1_hg19_subset",
        "a5c1_hg38_bwa",
        "a5c1_hg38_ensembl_bwa",
        "b178_hg19_subset",
        "dfc3_GRCh37_bwa",
        "dfc3_GRCh38_bwa",
        "dfc3_hg19_bwa",
        "dfc3_hg19_ensembl_bwa",
        "dfc3_hg19_subset",
        "dfc3_hg38_bwa",
        "dfc3_hg38_ensembl_bwa",
    }
)


def declare_mixed_layout_outcome(case: dict[str, Any]) -> dict[str, Any]:
    """Attach differential baseline/candidate expectations to a measured mixed case.

    Args:
        case: A derived base case.

    Returns:
        A shallow copy carrying side-specific expectations when its id is measured mixed.
    """
    declared = dict(case)
    if declared["case_id"] in MIXED_LAYOUT_BASE_CASE_IDS:
        required = list(declared["required_artifacts"])
        declared["side_expectations"] = {
            "before": {"expect_exit": "zero", "required_artifacts": required},
            "after": {"expect_exit": "nonzero", "required_artifacts": []},
        }
    return declared


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
    runtime.update(selected)
    return runtime
