"""Run-local admission of approved caller profiles during report verification."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import NoReturn

from vntyper.scripts.artifact_names import DECISION_PROFILE_SNAPSHOT_RELATIVE
from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)

logger = logging.getLogger(__name__)


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def resolve_recorded_explicit_profile(
    raw: bytes, snapshot_path: Path, summary: Mapping[str, object] | None
) -> ResolvedDecisionProfile:
    """Resolve an explicit snapshot, requiring full approval for caller schema v2.

    Args:
        raw: Already hash-checked canonical snapshot bytes.
        snapshot_path: The run's fixed provenance/decision_profile.json location.
        summary: Recorded analysis settings, required for caller calibration.

    Returns:
        A strictly resolved profile whose caller approval matches its actual bytes.

    Raises:
        ValueError: If a caller snapshot lacks its recorded, approved portable bundle.
    """
    document = load_strict_json_object(raw)
    metadata = document.get("generated_metadata")
    if not (
        document.get("schema_version") == 2
        and isinstance(metadata, Mapping)
        and metadata.get("generation_target") == "callers"
    ):
        return parse_decision_profile(raw, packaged_document=load_packaged_decision_profile().document)
    from vntyper.scripts.calibration_caller_bundle import load_caller_model_bundle

    if tuple(snapshot_path.parts[-2:]) != DECISION_PROFILE_SNAPSHOT_RELATIVE.parts:
        _fail("caller profile snapshot must use its fixed run-local provenance path")
    settings = None if summary is None else summary.get("analysis_settings")
    if not isinstance(settings, Mapping):
        _fail("caller profile snapshot requires recorded caller calibration analysis settings")
    bundle = load_caller_model_bundle(snapshot_path.parent.parent / "caller_calibration")
    if settings.get("caller_calibration_bundle_sha256") != bundle.sha256:
        _fail("recorded caller calibration bundle identity differs from its verified snapshot")
    if raw != bundle.profile.canonical_bytes:
        _fail("recorded decision profile differs from the approved caller bundle")
    return bundle.profile
