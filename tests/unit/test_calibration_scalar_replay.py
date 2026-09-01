"""Scalar dominance replay preserves closed abstention reasons."""

import pytest

from vntyper.scripts.calibration_scalar_replay import replay_scalar_dominance

pytestmark = pytest.mark.unit


def _component(**updates):
    value = {
        "enabled": True,
        "minimum_record_count_margin": 1,
        "minimum_record_share": 0.5,
        "minimum_record_share_margin": 0.0,
        "xd_veto": "disabled",
        "abstain_on_inadmissible_advntr": False,
    }
    value.update(updates)
    return value


def _features(**updates):
    value = {
        "canonical_identity": "identity-a",
        "haplotype_record_count_margin": 2,
        "haplotype_record_share": 0.75,
        "haplotype_record_share_margin": 0.5,
        "haplotype_record_tie": False,
        "xd_availability_fraction": 1.0,
        "advntr_evidence_disposition": "admissible",
    }
    value.update(updates)
    return value


def test_scalar_replay_selects_or_reports_each_replayable_abstention() -> None:
    assert replay_scalar_dominance(_features(), _component()).selected_identity == "identity-a"
    assert replay_scalar_dominance(_features(haplotype_record_tie=True), _component()).abstention_reason == "record-tie"
    assert (
        replay_scalar_dominance(_features(haplotype_record_count_margin=0), _component()).abstention_reason
        == "insufficient-dominance"
    )
    assert (
        replay_scalar_dominance(
            _features(xd_availability_fraction=0.5), _component(xd_veto="missingness")
        ).abstention_reason
        == "xd-missingness"
    )
    assert (
        replay_scalar_dominance(
            _features(advntr_evidence_disposition="identity-insufficient"),
            _component(abstain_on_inadmissible_advntr=True),
        ).abstention_reason
        == "inadmissible-advntr"
    )


def test_scalar_replay_refuses_vetoes_that_require_per_record_xd() -> None:
    with pytest.raises(ValueError, match="concentration"):
        replay_scalar_dominance(_features(), _component(xd_veto="concentration"))


def test_missing_identity_or_dominance_fields_are_not_applicable() -> None:
    assert not replay_scalar_dominance(_features(canonical_identity=None), _component()).applicable
    incomplete = _features()
    incomplete.pop("haplotype_record_share")
    assert not replay_scalar_dominance(incomplete, _component()).applicable
