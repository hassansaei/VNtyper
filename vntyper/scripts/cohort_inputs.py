"""
vntyper/scripts/cohort_inputs.py

Module Purpose:
---------------
Locate the samples of a cohort and read what the pipeline recorded for each of them.

``vntyper cohort`` is pointed at whatever the caller has: a directory that is itself one
sample's output, a directory holding many, a zip file of either, or something that turns
out to be none of those. :func:`discover_sample_directories` reduces all of that to the
set of directories holding a ``pipeline_summary.json``, extracting any zip files on the
way; :func:`load_pipeline_summary_for_sample` then reads one such directory into the
three pieces the report needs.

Two failure modes are deliberate and both are characterised rather than changed. A
sample that cannot be read - a missing summary, malformed JSON, an unparseable timestamp
- is logged and dropped so one bad sample cannot abort a 40-sample cohort. And the
discovery returns a **set**, so the order samples appear in the report is not
reproducible between processes; that one is a defect, recorded in
``tests/unit/test_cohort_inputs.py::test_discovery_returns_an_unordered_set_today``.

The four step names this module matches are compared by exact string against what
``pipeline.py`` writes. A typo does not fail - it silently drops a section of the report
(AGENTS.md trap 5) - so they are imported from ``summary_steps``, never spelled out.

Extracted from ``cohort_summary.py`` in Task 22 of the #181-#197 follow-ups.
"""

import hashlib
import json
import logging
import shutil
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any

from vntyper.scripts.summary_steps import (
    STEP_ADVNTR,
    STEP_BAM_HEADER,
    STEP_COVERAGE,
    STEP_KESTREL,
)

logger = logging.getLogger(__name__)

#: The file every sample directory is identified by.
PIPELINE_SUMMARY_FILENAME = "pipeline_summary.json"


def discover_sample_directories(input_paths: list[str]) -> tuple[set[Path], list[str]]:
    """Resolve the cohort's input paths to the sample directories they contain.

    A directory holding a ``pipeline_summary.json`` is taken as one sample and is not
    searched further; any other directory is searched recursively. A zip file is
    extracted to a temporary directory and then treated the same way. Anything else -
    a path that does not exist, a file that is not a zip - is logged and skipped, so a
    cohort of forty runs is not lost to one bad argument.

    Args:
        input_paths: Directories or zip files, as ``vntyper cohort`` received them.

    Returns:
        tuple[set[Path], list[str]]: The sample directories, and the temporary
        directories the caller must pass to :func:`cleanup_temp_dirs` once the report
        has been written. The directories come back as a **set**, so their iteration
        order is not stable between processes.
    """
    temp_dirs: list[str] = []
    processed_dirs: set[Path] = set()  # use a set to avoid duplicate directories

    # Identify valid directories/zip files (no changes here)
    for path_str in input_paths:
        path = Path(path_str)
        if not path.exists():
            logger.warning(f"Input path does not exist and will be skipped: {path}")
            continue
        if path.is_dir():
            if (path / PIPELINE_SUMMARY_FILENAME).exists():
                logger.info(f"Found {PIPELINE_SUMMARY_FILENAME} in {path}")
                processed_dirs.add(path)
            else:
                found = False
                for summary_file_path in path.rglob(PIPELINE_SUMMARY_FILENAME):
                    sample_dir = summary_file_path.parent
                    logger.info(f"Found {PIPELINE_SUMMARY_FILENAME} in {sample_dir}")
                    processed_dirs.add(sample_dir)
                    found = True
                if not found:
                    logger.warning(f"No {PIPELINE_SUMMARY_FILENAME} found in directory {path}")
        elif zipfile.is_zipfile(path):
            logger.info(f"Extracting zip file: {path}")
            temp_dir = tempfile.mkdtemp(prefix="cohort_zip_")
            try:
                with zipfile.ZipFile(path, "r") as zip_ref:
                    zip_ref.extractall(temp_dir)
                temp_path = Path(temp_dir)
                if (temp_path / PIPELINE_SUMMARY_FILENAME).exists():
                    logger.info(f"Found {PIPELINE_SUMMARY_FILENAME} in {temp_path}")
                    processed_dirs.add(temp_path)
                else:
                    found = False
                    for summary_file_path in temp_path.rglob(PIPELINE_SUMMARY_FILENAME):
                        sample_dir = summary_file_path.parent
                        logger.info(f"Found {PIPELINE_SUMMARY_FILENAME} in {sample_dir}")
                        processed_dirs.add(sample_dir)
                        found = True
                    if not found:
                        logger.warning(f"No {PIPELINE_SUMMARY_FILENAME} found in extracted zip file: {path}")
                temp_dirs.append(temp_dir)
            except zipfile.BadZipFile as e:
                logger.error(f"Bad zip file {path}: {e}")
                shutil.rmtree(temp_dir)
            except Exception as e:
                logger.error(f"Error extracting zip file {path}: {e}")
                shutil.rmtree(temp_dir)
        else:
            logger.warning(f"Unsupported file type (not a directory or zip): {path}")

    return processed_dirs, temp_dirs


def cleanup_temp_dirs(temp_dirs: list[str]) -> None:
    """Remove the directories the zip inputs were extracted into.

    Called after the report has been written, so a directory that will not delete is
    logged rather than raised - the cohort's output is already on disk by then.

    Args:
        temp_dirs: The temporary directories from :func:`discover_sample_directories`.
    """
    for temp_dir in temp_dirs:
        try:
            shutil.rmtree(temp_dir)
            logger.debug(f"Cleaned up temporary directory: {temp_dir}")
        except Exception as e:
            logger.error(f"Failed to remove temporary directory {temp_dir}: {e}")


def pseudonymized_sample_name(prefix: Any, original_sample: str) -> str:
    """Build the pseudonym a sample is reported under.

    The pseudonym is the caller's prefix followed by the first five hex digits of the
    MD5 of the original name, so it is stable across runs and the pseudonymization
    table written beside the report stays meaningful.

    Args:
        prefix: The value ``--pseudonymize`` supplied. Interpolated rather than
            concatenated, so a non-string does not raise.
        original_sample: The sample directory's name.

    Returns:
        str: The pseudonym.
    """
    # Compute MD5 hash of the original sample name and take first 5 characters.
    hash_suffix = hashlib.md5(original_sample.encode()).hexdigest()[:5]
    return f"{prefix}{hash_suffix}"


def parse_pipeline_summary(summary: dict[str, Any]) -> tuple[list[dict], list[dict], dict[str, Any]]:
    """Extract the cohort's three inputs from one parsed ``pipeline_summary.json``.

    Steps this cohort does not consume are ignored, and a step recorded more than once
    leaves the last occurrence in place.

    Args:
        summary: A parsed ``pipeline_summary.json`` mapping.

    Returns:
        tuple[list[dict], list[dict], dict[str, Any]]: The Kestrel rows, the adVNTR
        rows, and the per-sample statistics (``runtime``, ``version``, ``assembly``,
        ``pipeline``, ``coverage``).

    Raises:
        ValueError: If ``pipeline_start`` or ``pipeline_end`` is present but is not an
            ISO-8601 timestamp. :func:`load_pipeline_summary_for_sample` catches this
            and drops the sample.
    """
    kestrel_data: list[dict] = []
    advntr_data: list[dict] = []
    additional_stats: dict[str, Any] = {}

    # Compute runtime from top-level timestamps if available
    pipeline_start = summary.get("pipeline_start")
    pipeline_end = summary.get("pipeline_end")
    if pipeline_start and pipeline_end:
        start_dt = datetime.fromisoformat(pipeline_start)
        end_dt = datetime.fromisoformat(pipeline_end)
        runtime_sec = (end_dt - start_dt).total_seconds()
        additional_stats["runtime"] = f"{runtime_sec:.2f} seconds"
    else:
        additional_stats["runtime"] = "N/A"

    # Pipeline version from top-level field
    additional_stats["version"] = summary.get("version", "N/A")

    # Initialize defaults for assembly, pipeline and coverage
    additional_stats["assembly"] = "N/A"
    additional_stats["pipeline"] = "N/A"
    additional_stats["coverage"] = {}

    for step in summary.get("steps", []):
        if step.get("step") == STEP_KESTREL:
            kestrel_data = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == STEP_ADVNTR:
            advntr_data = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == STEP_BAM_HEADER:
            parsed = step.get("parsed_result", {})
            additional_stats["assembly"] = parsed.get("assembly_text", "N/A")
            additional_stats["pipeline"] = parsed.get("alignment_pipeline", "N/A")
        elif step.get("step") == STEP_COVERAGE:
            parsed = step.get("parsed_result", {})
            data_list = parsed.get("data", [])
            if data_list:
                additional_stats["coverage"] = data_list[0]
    return kestrel_data, advntr_data, additional_stats


def load_pipeline_summary_for_sample(sample_dir) -> tuple[list[dict], list[dict], dict[str, Any]]:
    """
    Load the pipeline_summary.json from a sample directory and extract Kestrel,
    adVNTR data and additional statistics (runtime, coverage, version, assembly, pipeline).

    For the adVNTR step, the algorithm result will later be computed based on the logic
    defined in report_config.json.

    Anything that goes wrong - the file is absent, the JSON does not parse, a timestamp
    does not - is logged and yields three empties, so one unreadable sample does not
    abort the cohort.

    Parameters
    ----------
    sample_dir : str or Path
        Directory containing the pipeline_summary.json file.

    Returns
    -------
    tuple
        Three elements: (kestrel_data, advntr_data, additional_stats).
        additional_stats is a dict with keys:
          - runtime: pipeline run duration (in seconds)
          - version: pipeline version
          - assembly: assembly text from BAM Header Parsing
          - pipeline: alignment pipeline from BAM Header Parsing
          - coverage: coverage metrics dict (mean, median, stdev, min, max)
    """
    sample_dir = Path(sample_dir)
    summary_path = sample_dir / PIPELINE_SUMMARY_FILENAME
    if not summary_path.exists():
        logger.warning(f"Pipeline summary file not found in {sample_dir}")
        return [], [], {}
    try:
        with open(summary_path) as f:
            summary = json.load(f)
        return parse_pipeline_summary(summary)
    except Exception as e:
        logger.error(f"Error loading pipeline summary from {sample_dir}: {e}")
        return [], [], {}
