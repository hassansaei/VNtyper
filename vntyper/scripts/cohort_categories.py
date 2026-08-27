"""
vntyper/scripts/cohort_categories.py

Module Purpose:
---------------
Reduce row-level algorithm verdicts to the per-sample categories the cohort report
counts: ``Positive``, ``Positive_Flagged``, ``Negative`` and ``Unestablished``.

``algorithm_rules.compute_algorithm_result`` answers a question about one result *row*. A
cohort report answers a question about one *sample*, and a sample can contribute several
rows. This module is the reduction in between - map each verdict onto a category, then
take the highest category any of the sample's rows reached - together with the count of
each category that feeds the donut charts.

Kestrel and adVNTR name their verdicts differently (``High_Precision_flagged`` against
``positive flagged``), so the two mappings are separate functions rather than one shared
table. The shared interpreter lives in ``algorithm_rules.py``; this reduction was
extracted from ``cohort_summary.py`` in Task 22 of the #181-#197 follow-ups and is
characterised by ``tests/unit/test_cohort_categories.py``.
"""

import logging
from collections.abc import Callable, Sequence
from typing import Any

import pandas as pd

from vntyper.scripts.algorithm_rules import UNESTABLISHED_RESULT, compute_algorithm_result

logger = logging.getLogger(__name__)

#: The sample-level category for a result the run did not establish. It covers
#: unevaluable rows and named samples that contributed no rows; neither is Negative.
UNESTABLISHED_CATEGORY = "Unestablished"


# --------------------------------------------------------------------------
# Aggregate row-level results into sample-level categories.
# --------------------------------------------------------------------------
def unify_kestrel_result(row_result: Any) -> str:
    """Map one row-level Kestrel result to its cohort category.

    Args:
        row_result: The result emitted by the configured Kestrel rule interpreter.

    Returns:
        str: Positive, Positive_Flagged, Unestablished, or Negative.
    """
    if row_result in ["High_Precision", "Low_Precision"]:
        return "Positive"
    elif row_result in ["High_Precision_flagged", "Low_Precision_flagged"]:
        return "Positive_Flagged"
    elif row_result == UNESTABLISHED_RESULT:
        return UNESTABLISHED_CATEGORY
    else:
        return "Negative"


def unify_advntr_result(row_result: Any) -> str:
    """Map one row-level adVNTR result to its cohort category.

    Args:
        row_result: The result emitted by the configured adVNTR rule interpreter.

    Returns:
        str: Positive, Positive_Flagged, Unestablished, or Negative.
    """
    if row_result == "positive":
        return "Positive"
    elif row_result == "positive flagged":
        return "Positive_Flagged"
    elif row_result == UNESTABLISHED_RESULT:
        return UNESTABLISHED_CATEGORY
    else:
        return "Negative"


def aggregate_sample_category(results: list[str]) -> str:
    """Choose one sample category from its row categories.

    Precedence is Positive, Positive_Flagged, Unestablished, then Negative. An
    unevaluable row prevents a negative sample assertion unless another row establishes
    a positive result.

    Args:
        results: The mapped categories of every row for one sample.

    Returns:
        str: The highest-precedence category present, or Negative for an empty list.
    """
    if any(r == "Positive" for r in results):
        return "Positive"
    elif any(r == "Positive_Flagged" for r in results):
        return "Positive_Flagged"
    elif any(r == UNESTABLISHED_CATEGORY for r in results):
        return UNESTABLISHED_CATEGORY
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

    ``frame`` is **not modified**. The ``Sample`` values are normalised to strings and
    the two working columns ``__row_result`` and ``__unified`` are written onto a copy,
    because the caller goes on to publish the frame: ``aggregate_cohort`` renders the
    report and then exports the very same frame as CSV/TSV/JSON. Normalising the source
    frame would change an export; leaving the copied sample values numeric would make a
    direct caller's string roster count one sample twice. Pinned by
    ``test_cohort_categories.py::test_the_reduction_leaves_the_caller_s_frame_untouched``
    and the numeric direct-render regression in ``test_cohort_summary_oracle.py``.

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
        annotated["Sample"] = annotated["Sample"].astype(str)
        annotated["__row_result"] = annotated.apply(
            lambda row: compute_algorithm_result(pd.DataFrame([row]), logic_config),
            axis=1,
        )
        annotated["__unified"] = annotated["__row_result"].apply(unify)
        return annotated.groupby("Sample")["__unified"].apply(list).apply(aggregate_sample_category)
    return pd.Series(dtype=str)


def complete_sample_categories(categories: pd.Series, sample_names: Sequence[str]) -> pd.Series:
    """Give every named sample a category; one with no rows is Unestablished.

    Args:
        categories: One category per sample that contributed rows.
        sample_names: Every sample the cohort knows about.

    Returns:
        pd.Series: One category per named sample. The input is not modified.
    """
    completed = categories.copy()
    for sample in sample_names:
        if sample not in completed.index:
            completed.loc[sample] = UNESTABLISHED_CATEGORY
    return completed


def samples_without_rows(frame: pd.DataFrame, sample_names: Sequence[str]) -> list[str]:
    """Return named samples absent from a result frame, in roster order.

    Args:
        frame: Aggregated result rows; it may lack a ``Sample`` column.
        sample_names: Every sample the cohort knows about.

    Returns:
        list[str]: Samples counted as Unestablished because they contributed no rows.
    """
    present = set(frame["Sample"].astype(str)) if "Sample" in frame.columns else set()
    return [sample for sample in sample_names if sample not in present]


def category_counts(categories: pd.Series) -> tuple[int, int, int, int, int]:
    """Count each category, in the order the donut chart draws them.

    The total is the sum of the four known categories rather than the length of
    ``categories``, so a label the reduction never produces is left out of the chart
    instead of being attributed to one of the four segments.

    Args:
        categories: One category per sample, as returned by :func:`sample_categories`.

    Returns:
        tuple[int, int, int, int, int]: ``(positive, positive_flagged, negative,
        unestablished, total)``.
    """
    counts = categories.value_counts()
    positive = counts.get("Positive", 0)
    flagged = counts.get("Positive_Flagged", 0)
    negative = counts.get("Negative", 0)
    unestablished = counts.get(UNESTABLISHED_CATEGORY, 0)
    return positive, flagged, negative, unestablished, positive + flagged + negative + unestablished
