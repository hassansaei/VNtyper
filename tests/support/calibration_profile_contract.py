"""Release-independent fitted content with separately verified release metadata."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from vntyper.scripts.calibration_profiles import build_generated_profile
from vntyper.scripts.canonical_json import canonical_sha256, load_strict_json_object
from vntyper.scripts.decision_profile import ResolvedDecisionProfile

# Only generator_version and its derived profile_id vary with a package release.
FITTED_CONTENT_SHA256 = "0c241ebb2d1176b4a3c76d38c42676ea3aa7a80a4751dba70c3eaad7e2dec672"


def assert_fitted_profile(profile: ResolvedDecisionProfile, *, generator_version: str) -> None:
    """Pin the fitted decisions and provenance, then verify the release metadata.

    Args:
        profile: The resolved generated profile from the golden fit.
        generator_version: The installed package version expected for this fit.

    Raises:
        AssertionError: If fitted content, generator version or derived identity differs.
    """
    document = load_strict_json_object(profile.canonical_bytes)
    metadata = document["generated_metadata"]
    assert metadata["generator_version"] == generator_version, "unexpected generator version"
    derived = build_generated_profile(
        cast(Mapping[str, object], profile.components["dominance"]),
        dataset_manifest_hash=metadata["dataset_manifest_hash"],
        partition_manifest_hash=metadata["partition_manifest_hash"],
        seed=metadata["seed"],
        objective=metadata["objective"],
        generator_version=metadata["generator_version"],
    )
    assert profile.profile_id == derived.profile_id, "incorrect derived profile identity"
    del document["profile_id"]
    del metadata["generator_version"]
    assert canonical_sha256(document) == FITTED_CONTENT_SHA256, "fitted content changed"
