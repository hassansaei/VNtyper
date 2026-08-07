"""
vntyper/scripts/cohort_categories.py

Module Purpose:
---------------
Reduce row-level algorithm verdicts to the per-sample categories the cohort report
counts: ``Positive``, ``Positive_Flagged`` and ``Negative``.

``cohort_rules.compute_algorithm_result`` answers a question about one result *row*. A
cohort report answers a question about one *sample*, and a sample can contribute several
rows. This module is the reduction in between - map each verdict onto a category, then
take the highest category any of the sample's rows reached - together with the count of
each category that feeds the donut charts.

Kestrel and adVNTR name their verdicts differently (``High_Precision_flagged`` against
``positive flagged``), so the two mappings are separate functions rather than one shared
table. Extracted verbatim from ``cohort_summary.py`` in Task 22 of the #181-#197
follow-ups; ``tests/unit/test_cohort_categories.py`` characterises it.
"""

import logging
from collections.abc import Callable
from typing import Any

import pandas as pd

from vntyper.scripts.cohort_rules import compute_algorithm_result

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# New helper functions to correctly aggregate row-level results into
# sample-level categories (Positive, Positive_Flagged, Negative).
# --------------------------------------------------------------------------
def unify_kestrel_result(row_result: Any) -> str:
    """
    Convert a row-level Kestrel result (e.g. 'High_Precision', 'Low_Precision_flagged')
    into a broader category: 'Positive', 'Positive_Flagged', or 'Negative'.
    """
    if row_result in ["High_Precision", "Low_Precision"]:
        return "Positive"
    elif row_result in ["High_Precision_flagged", "Low_Precision_flagged"]:
        return "Positive_Flagged"
    else:
        return "Negative"


def unify_advntr_result(row_result: Any) -> str:
    """
    Convert a row-level adVNTR result (e.g. 'positive', 'positive flagged')
    into a broader category: 'Positive', 'Positive_Flagged', or 'Negative'.
    """
    if row_result == "positive":
        return "Positive"
    elif row_result == "positive flagged":
        return "Positive_Flagged"
    else:
        return "Negative"


def aggregate_sample_category(results: list[str]) -> str:
    """
    Given a list of final row-level categories for a sample (each in
    {'Positive', 'Positive_Flagged', 'Negative'}), pick the highest category
    following the rule:
      - If there's at least one 'Positive' => 'Positive'
      - Else if there's at least one 'Positive_Flagged' => 'Positive_Flagged'
      - Else => 'Negative'
    """
    if any(r == "Positive" for r in results):
        return "Positive"
    elif any(r == "Positive_Flagged" for r in results):
        return "Positive_Flagged"
    else:
        return "Negative"


def sample_categories(
    frame: pd.DataFrame,
    logic_config: dict[str, Any],
    unify: Callable[[Any], str],
) -> pd.Series:
    """Reduce a results frame to one category per sample.

    Each row is evaluated against ``logic_config``, the verdict is mapped through
    ``unify``, and the rows of each sample are aggregated with
    :func:`aggregate_sample_category`.

    ``frame`` is **not modified**. The two working columns ``__row_result`` and
    ``__unified`` are written onto a copy, because the caller goes on to publish the
    frame: ``aggregate_cohort`` renders the report and then exports the very same frame
    as CSV/TSV/JSON, so annotating it in place - which this did until the #181-#197
    follow-ups - put both internal columns into every machine-readable output. Pinned by
    ``test_cohort_categories.py::test_the_reduction_leaves_the_caller_s_frame_untouched``.

    Args:
        frame: One row per result. Not modified. A frame that is empty, or that carries
            no ``Sample`` column, yields an empty series.
        logic_config: The ``algorithm_logic`` entry for this algorithm.
        unify: :func:`unify_kestrel_result` or :func:`unify_advntr_result`.

    Returns:
        pd.Series: One category per sample, indexed by sample name.
    """
    if not frame.empty and "Sample" in frame.columns:
        annotated = frame.copy()
        annotated["__row_result"] = annotated.apply(
            lambda row: compute_algorithm_result(pd.DataFrame([row]), logic_config),
            axis=1,
        )
        annotated["__unified"] = annotated["__row_result"].apply(unify)
        return annotated.groupby("Sample")["__unified"].apply(list).apply(aggregate_sample_category)
    return pd.Series(dtype=str)


def category_counts(categories: pd.Series) -> tuple[int, int, int, int]:
    """Count each category, in the order the donut chart draws them.

    The total is the sum of the three known categories rather than the length of
    ``categories``, so a label the reduction never produces is left out of the chart
    instead of being attributed to one of the three segments.

    Args:
        categories: One category per sample, as returned by :func:`sample_categories`.

    Returns:
        tuple[int, int, int, int]: ``(positive, positive_flagged, negative, total)``.
    """
    counts = categories.value_counts()
    positive = counts.get("Positive", 0)
    flagged = counts.get("Positive_Flagged", 0)
    negative = counts.get("Negative", 0)
    return positive, flagged, negative, positive + flagged + negative
