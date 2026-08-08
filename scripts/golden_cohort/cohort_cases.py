"""Pure cohort-case declarations for the golden-cohort matrix."""

from __future__ import annotations

from typing import Any

from golden_cohort.admissibility import COHORT_REQUIRED_ARTIFACTS


def build_cohort_cases(pipeline_case_ids: list[str]) -> list[dict[str, Any]]:
    """Build cohort-mode cases over candidate-success per-sample outputs.

    Args:
        pipeline_case_ids: Per-sample case ids expected to write candidate summaries.

    Returns:
        The four cohort cases, in run order.
    """
    single = pipeline_case_ids[:1]
    return [
        {
            "case_id": "cohort_multi",
            "kind": "cohort",
            "group": "cohort",
            "inputs": list(pipeline_case_ids),
            "summary_formats": "csv,tsv,json",
            "pseudonymize": None,
            "expect_exit": "zero",
            "allow_missing_inputs": False,
            "required_artifacts": list(COHORT_REQUIRED_ARTIFACTS),
        },
        {
            "case_id": "cohort_multi_pseudonymized",
            "kind": "cohort",
            "group": "cohort",
            "inputs": list(pipeline_case_ids),
            "summary_formats": "csv,tsv,json",
            "pseudonymize": "sample_",
            "expect_exit": "zero",
            "allow_missing_inputs": False,
            "required_artifacts": [*COHORT_REQUIRED_ARTIFACTS, "pseudonymization_table.tsv"],
        },
        {
            "case_id": "cohort_single",
            "kind": "cohort",
            "group": "cohort",
            "inputs": single,
            "summary_formats": "csv,tsv,json",
            "pseudonymize": None,
            "expect_exit": "zero",
            "allow_missing_inputs": False,
            "required_artifacts": list(COHORT_REQUIRED_ARTIFACTS),
        },
        {
            "case_id": "cohort_empty",
            "kind": "cohort",
            "group": "cohort",
            "inputs": [],
            "empty_input_dir": True,
            "summary_formats": "csv,tsv,json",
            "pseudonymize": None,
            "expect_exit": "zero",
            "allow_missing_inputs": True,
            "required_artifacts": [],
        },
    ]
