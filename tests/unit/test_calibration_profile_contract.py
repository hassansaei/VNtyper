"""Unit proof that the fitted golden contract ignores only release metadata."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from tests.support.calibration_profile_contract import FITTED_CONTENT_SHA256, assert_fitted_profile
from vntyper.scripts.calibration_profiles import build_generated_profile
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)

pytestmark = pytest.mark.unit


def _profile(version: str = "2.0.42", minimum_share: float = 0.5) -> ResolvedDecisionProfile:
    return build_generated_profile(
        {
            "abstain_on_inadmissible_advntr": False,
            "enabled": True,
            "minimum_record_count_margin": 1,
            "minimum_record_share": minimum_share,
            "minimum_record_share_margin": 0.0,
            "xd_veto": "disabled",
        },
        dataset_manifest_hash="fcbcf50db55b2b915296ca44d41628fd1599ba81cdf555cf55da918fd7562e9f",
        partition_manifest_hash="2a520f9c8cd989def64bca26d2886e2beddc277e2dd780d522554349d210a758",
        seed=295,
        objective="lexicographic-safety-v1",
        generator_version=version,
    )


def test_release_bump_changes_real_digest_but_preserves_fitted_contract() -> None:
    before = _profile()
    after = _profile("2.0.43")

    assert before.digest != after.digest
    assert before.profile_id != after.profile_id
    assert_fitted_profile(before, generator_version="2.0.42")
    assert_fitted_profile(after, generator_version="2.0.43")


def test_fitted_contract_does_not_mutate_profile() -> None:
    profile = _profile()
    before = deepcopy(load_strict_json_object(profile.canonical_bytes))

    assert_fitted_profile(profile, generator_version="2.0.42")

    assert load_strict_json_object(profile.canonical_bytes) == before


def test_fitted_value_change_fails_even_with_correct_derived_identity() -> None:
    changed = _profile(minimum_share=0.75)

    with pytest.raises(AssertionError, match="fitted content changed"):
        assert_fitted_profile(changed, generator_version="2.0.42")


def test_stale_generator_version_fails() -> None:
    with pytest.raises(AssertionError, match="generator version"):
        assert_fitted_profile(_profile(), generator_version="2.0.43")


def test_wrong_derived_profile_id_fails() -> None:
    document = load_strict_json_object(_profile().canonical_bytes)
    document["profile_id"] = "vntyper-generated-wrong"
    profile = parse_decision_profile(
        canonical_json_bytes(document), packaged_document=load_packaged_decision_profile().document
    )

    with pytest.raises(AssertionError, match="derived profile identity"):
        assert_fitted_profile(profile, generator_version="2.0.42")


def test_golden_documentation_records_content_pin_and_conditional_gate() -> None:
    page = (Path(__file__).resolve().parents[2] / "docs/development/calibration-validation.md").read_text(
        encoding="utf-8"
    )

    assert FITTED_CONTENT_SHA256 in page
    assert "excluding only `profile_id` and `generated_metadata.generator_version`" in page
    assert "make test-golden" in page
    assert "`make check-all` also runs this tier when both roots are set" in page
