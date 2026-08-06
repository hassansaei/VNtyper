"""
vntyper/scripts/cohort_exports.py

Module Purpose:
---------------
Write the machine-readable cohort outputs that accompany the HTML report.

``vntyper cohort --output-format csv,tsv,json`` writes the aggregated Kestrel and adVNTR
frames beside ``cohort_summary.html``. Those are the files a downstream analysis reads,
so the file names, the separators and the log lines announcing them are a contract.
:func:`write_pseudonymization_table` writes the mapping that makes a pseudonymised
cohort resolvable again.

Note that nothing here strips columns: whatever the frame carries is exported. The
report is rendered before these are written and rendering annotates the frames with two
working columns, so the exports currently carry ``__row_result`` and ``__unified`` as
well. That is a defect, characterised rather than fixed - see
``tests/unit/test_cohort_exports.py::test_a_frame_carrying_the_reductions_working_columns_exports_them_today``.

Extracted from ``cohort_summary.py`` in Task 22 of the #181-#197 follow-ups.
"""

import logging
import os
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def parse_output_formats(additional_formats: str) -> list[str]:
    """Split the ``--output-format`` argument into the formats it names.

    Nothing validates the result: a name this module does not implement is carried
    through and simply never matches, so an unsupported format produces the HTML report
    and no error.

    Args:
        additional_formats: A comma-separated list, e.g. ``"csv, TSV"``.

    Returns:
        list[str]: The named formats, lower-cased and stripped, empties dropped.
    """
    return [fmt.strip().lower() for fmt in additional_formats.split(",") if fmt.strip()]


def write_cohort_frame(
    frame: pd.DataFrame,
    output_dir: str | os.PathLike[str],
    stem: str,
    label: str,
    formats: list[str],
) -> None:
    """Write one aggregated frame in each of the requested formats.

    An empty frame writes nothing at all, so a cohort in which no sample produced
    adVNTR results gets no ``cohort_advntr.csv`` rather than a header-only one.

    Args:
        frame: The aggregated results. Not modified.
        output_dir: The cohort's output directory.
        stem: The file name without its extension, e.g. ``"cohort_kestrel"``.
        label: How the algorithm is named in the log line, e.g. ``"Kestrel"``.
        formats: The formats to write, as :func:`parse_output_formats` returns them.
    """
    if frame.empty:
        return
    if "csv" in formats:
        csv_path = Path(output_dir) / f"{stem}.csv"
        frame.to_csv(csv_path, index=False)
        logger.info(f"Cohort {label} CSV written to: {csv_path}")
    if "tsv" in formats:
        tsv_path = Path(output_dir) / f"{stem}.tsv"
        frame.to_csv(tsv_path, sep="\t", index=False)
        logger.info(f"Cohort {label} TSV written to: {tsv_path}")
    if "json" in formats:
        json_path = Path(output_dir) / f"{stem}.json"
        frame.to_json(json_path, orient="records", indent=4)
        logger.info(f"Cohort {label} JSON written to: {json_path}")


def write_pseudonymization_table(output_dir: str | os.PathLike[str], sample_mapping: dict[str, str]) -> None:
    """Write the pseudonym-to-original mapping as a two-column TSV.

    This runs after the report is already on disk, so a failure is logged rather than
    raised: a completed cohort must not become a traceback.

    The rows go out through ``writelines`` rather than the ``write``-per-row loop this
    was extracted from. That is the one non-verbatim change in the extraction and it is
    ruff's own fix for FURB122, which the loop violates and which ruff did not report
    while the block was buried in ``aggregate_cohort``. The bytes are unchanged -
    ``writelines`` inserts no separator of its own - and
    ``test_cohort_exports.py::test_the_pseudonymization_table_is_a_two_column_tsv``
    pins them.

    Args:
        output_dir: The cohort's output directory.
        sample_mapping: Pseudonym -> original sample name.
    """
    pseudonym_table_path = Path(output_dir) / "pseudonymization_table.tsv"
    try:
        with open(pseudonym_table_path, "w") as pt:
            pt.write("Pseudonym\tOriginal\n")
            pt.writelines(f"{pseudonym}\t{original}\n" for pseudonym, original in sample_mapping.items())
        logger.info(f"Pseudonymization table written to: {pseudonym_table_path}")
    except Exception as e:
        logger.error(f"Failed to write pseudonymization table: {e}")
