"""
vntyper/scripts/resume.py

Module Purpose:
---------------
Pure logic and checkpoint inspection for pipeline resumption (--resume, #20).

Provides:
- load_prior_summary: Load and validate an existing pipeline_summary.json.
- resume_refusals: Enforce fatal run-identity invariant checks.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Final

from vntyper.scripts import summary_steps

logger = logging.getLogger(__name__)

#: Table of required extra sibling files per reusable stage
STEP_OUTPUT_SIBLINGS: Final[dict[str, tuple[str, ...]]] = {
    summary_steps.STEP_KESTREL: (
        "output.vcf",
        "output.bam",
        "kestrel_pre_result.tsv",
    ),
    summary_steps.STEP_ADVNTR: (),
    summary_steps.STEP_BAM_TO_FASTQ: ("output_sliced.bam",),
    summary_steps.STEP_CRAM_TO_FASTQ: ("output_sliced.bam",),
    summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT: ("output_sliced.bam",),
    summary_steps.STEP_FASTQ_ALIGNMENT: (),
}

_REUSABLE_STEPS: Final[frozenset[str]] = frozenset(STEP_OUTPUT_SIBLINGS.keys())


def _compute_md5(path: Path) -> str | None:
    """Compute the MD5 checksum of a file in 64 KiB chunks."""
    if not path.is_file():
        return None
    md5 = hashlib.md5()
    try:
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                md5.update(chunk)
        return md5.hexdigest()
    except OSError as exc:
        logger.warning("Error calculating MD5 for %s: %s", path, exc)
        return None


def load_prior_summary(path: str | Path) -> dict[str, Any] | None:
    """Load a prior pipeline summary from a file or directory path.

    Args:
        path: Path to pipeline_summary.json or the run directory containing it.

    Returns:
        dict[str, Any] | None: The parsed JSON dictionary, or None if missing/invalid.
    """
    p = Path(path)
    if p.is_dir():
        p = p / "pipeline_summary.json"

    if not p.is_file():
        logger.debug("Prior summary not found at %s", p)
        return None

    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            logger.warning("Prior summary at %s is not a dictionary", p)
            return None
        return data
    except (json.JSONDecodeError, OSError) as exc:
        logger.warning("Failed to parse prior summary at %s: %s", p, exc)
        return None


def resume_refusals(
    prior: dict[str, Any],
    *,
    version: str,
    input_files: dict[str, Any],
    sample_name: str,
    reference_key_used: str | None,
    decision_profile_sha256: str,
) -> list[str]:
    """Validate that run identity invariants match the prior summary.

    Refusals concern run identity only: version, input files, sample name,
    reference key, or decision-profile digest.

    Args:
        prior: Prior summary dictionary.
        version: Pipeline version for the current run.
        input_files: Input file paths mapping for the current run.
        sample_name: Resolved sample name for the current run.
        reference_key_used: Reference key used for the current run.
        decision_profile_sha256: SHA256 digest of the resolved decision profile.

    Returns:
        list[str]: Descriptions of all detected mismatches; empty if allowed.
    """
    refusals: list[str] = []

    prior_version = prior.get("version")
    if prior_version != version:
        refusals.append(f"version differs (prior: {prior_version!r}, current: {version!r})")

    prior_inputs = prior.get("input_files")
    if prior_inputs != input_files:
        refusals.append(f"input files differ (prior: {prior_inputs!r}, current: {input_files!r})")

    prior_sample = prior.get("sample_name")
    if prior_sample != sample_name:
        refusals.append(f"sample name differs (prior: {prior_sample!r}, current: {sample_name!r})")

    prior_ref_key = prior.get("reference_key_used")
    if prior_ref_key != reference_key_used:
        refusals.append(f"reference key differs (prior: {prior_ref_key!r}, current: {reference_key_used!r})")

    prior_profile_sha = prior.get("decision_profile_sha256")
    if prior_profile_sha != decision_profile_sha256:
        refusals.append(
            f"decision profile digest differs (prior: {prior_profile_sha!r}, current: {decision_profile_sha256!r})"
        )

    return refusals


def step_is_reusable(
    prior: dict[str, Any],
    step_name: str,
    output_root: str | Path,
) -> bool:
    """Determine whether a step from a prior run is valid for reuse.

    Checks:
    1. The stage is in the reusable set (conversion/alignment, Kestrel, adVNTR).
    2. The step was recorded in the prior summary without result_file_missing.
    3. The recorded result_file exists on disk.
    4. The file's current MD5 matches the recorded checksum.
    5. Every required sibling output file exists on disk.

    Args:
        prior: Prior summary dictionary.
        step_name: Step name to evaluate.
        output_root: Root directory of the output.

    Returns:
        bool: True if all artifacts exist and are uncorrupted; False otherwise.
    """
    if step_name not in _REUSABLE_STEPS:
        return False

    prior_steps = prior.get("steps", [])
    step_record = next((s for s in prior_steps if isinstance(s, dict) and s.get("step") == step_name), None)
    if step_record is None:
        logger.debug("Step %r was not recorded in prior summary", step_name)
        return False

    if step_record.get("result_file_missing"):
        logger.debug("Step %r recorded a missing result file in prior summary", step_name)
        return False

    raw_result_file = step_record.get("result_file")
    if not raw_result_file:
        logger.debug("Step %r has no result_file recorded", step_name)
        return False

    recorded_md5 = step_record.get("md5sum")
    if not recorded_md5:
        logger.debug("Step %r has no md5sum recorded", step_name)
        return False

    result_path = Path(raw_result_file)
    if not result_path.is_file():
        # Fallback relative to output_root if directory was relocated
        candidate = Path(output_root) / result_path.name
        if candidate.is_file():
            result_path = candidate
        else:
            logger.debug("Result file %s for step %r does not exist", result_path, step_name)
            return False

    actual_md5 = _compute_md5(result_path)
    if actual_md5 != recorded_md5:
        logger.warning(
            "MD5 mismatch for step %r (%s): recorded %s, current %s",
            step_name,
            result_path,
            recorded_md5,
            actual_md5,
        )
        return False

    # Check required sibling files
    stage_dir = result_path.parent
    for sibling in STEP_OUTPUT_SIBLINGS.get(step_name, ()):
        sibling_path = stage_dir / sibling
        if not sibling_path.is_file():
            logger.debug("Required sibling %s for step %r does not exist", sibling_path, step_name)
            return False

    # Check mate FASTQ for BAM/CRAM conversion if R1 is present
    if "_R1" in result_path.name:
        mate_name = result_path.name.replace("_R1", "_R2")
        mate_path = stage_dir / mate_name
        # If mate was created in the stage, it must still exist
        if not mate_path.is_file():
            logger.debug("Expected mate FASTQ %s for step %r does not exist", mate_path, step_name)
            return False

    return True


def make_reused_step_record(
    prior_step: dict[str, Any],
    prior_pipeline_start: str | None,
) -> dict[str, Any]:
    """Construct a carry-forward step record with reused_from provenance.

    Args:
        prior_step: The exact step record from the prior summary.
        prior_pipeline_start: The pipeline_start timestamp of the donor run.

    Returns:
        dict[str, Any]: Verbatim copy with reused_from set.
    """
    record = dict(prior_step)
    record["reused_from"] = prior_pipeline_start
    return record
