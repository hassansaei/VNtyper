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

import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


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
    except Exception as exc:
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
