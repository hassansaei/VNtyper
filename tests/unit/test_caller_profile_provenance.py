"""Reports require the same approved caller bundle as the recorded run profile."""

import pytest

from tests.unit.test_calibration_caller_bundle import bundle
from tests.unit.test_profile_provenance import _schema_three_summary
from vntyper.scripts.calibration_caller_bundle import load_caller_model_bundle
from vntyper.scripts.profile_provenance import resolve_summary_profile, snapshot_decision_profile

pytestmark = pytest.mark.unit


def prepared(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    portable = bundle(source)
    run = tmp_path / "run"
    run.mkdir()
    portable.rename(run / "caller_calibration")
    loaded = load_caller_model_bundle(run / "caller_calibration")
    snapshot_decision_profile(loaded.profile, run / "provenance" / "decision_profile.json")
    summary = _schema_three_summary(loaded.profile)
    summary["analysis_settings"] = {"caller_calibration_bundle_sha256": loaded.sha256}
    return run, loaded, summary


def test_report_admits_only_verified_approved_caller_snapshot(tmp_path):
    run, loaded, summary = prepared(tmp_path)
    verified = resolve_summary_profile(summary, run)
    assert verified.sha256 == loaded.profile.digest
    assert verified.profile_kind == "generated"
    assert verified.source == "explicit-cli"


@pytest.mark.parametrize(
    "change", ["bundle-missing", "approval-changed", "unrecorded", "wrong-digest", "bundle-symlink"]
)
def test_report_cannot_admit_standalone_or_changed_caller_profile(tmp_path, change):
    run, _, summary = prepared(tmp_path)
    root = run / "caller_calibration"
    if change == "bundle-missing":
        root.rename(run / "removed")
    elif change == "approval-changed":
        (root / "portable-approval.json").write_text("{}")
    elif change == "unrecorded":
        summary.pop("analysis_settings")
    elif change == "wrong-digest":
        summary["analysis_settings"]["caller_calibration_bundle_sha256"] = "0" * 64
    else:
        real = run / "outside"
        root.rename(real)
        root.symlink_to(real, target_is_directory=True)
    with pytest.raises(ValueError):
        resolve_summary_profile(summary, run)


def test_profile_snapshot_cannot_use_a_different_relative_path(tmp_path):
    from vntyper.scripts.profile_provenance import DecisionProfileProvenance, verify_profile_snapshot

    run, loaded, summary = prepared(tmp_path)
    # Even an approved bundle cannot authorize arbitrary profile bytes at another snapshot path.
    alternate = run / "elsewhere" / "decision_profile.json"
    snapshot_decision_profile(loaded.profile, alternate)
    provenance = DecisionProfileProvenance(
        loaded.profile.profile_id,
        loaded.profile.profile_revision,
        "generated",
        "explicit-cli",
        loaded.profile.digest,
        "provenance/decision_profile.json",
    )
    with pytest.raises(ValueError):
        verify_profile_snapshot(provenance, alternate, schema_three_summary=summary)


def test_approved_bundle_cannot_authorize_different_profile_bytes(tmp_path):
    from vntyper.scripts.profile_provenance import profile_summary_fields

    run, _, summary = prepared(tmp_path)
    other = tmp_path / "other"
    other.mkdir()
    different = load_caller_model_bundle(bundle(other, mode="legacy"))
    snapshot_decision_profile(different.profile, run / "provenance" / "decision_profile.json")
    summary.update(profile_summary_fields(different.profile))
    with pytest.raises(ValueError, match="differs from the approved caller bundle"):
        resolve_summary_profile(summary, run)


def test_direct_caller_snapshot_verification_requires_recorded_settings(tmp_path):
    from vntyper.scripts.caller_profile_provenance import resolve_recorded_explicit_profile

    run, loaded, _ = prepared(tmp_path)
    with pytest.raises(ValueError, match="recorded"):
        resolve_recorded_explicit_profile(
            loaded.profile.canonical_bytes, run / "provenance" / "decision_profile.json", None
        )
