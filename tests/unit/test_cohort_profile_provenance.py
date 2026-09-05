"""Decision-profile provenance across cohort loading, grouping, and exports."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.cli import load_config
from vntyper.scripts import cohort_summary, summary
from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.cohort_inputs import load_pipeline_summary_for_sample
from vntyper.scripts.cohort_profiles import annotate_profile_rows, group_decision_profiles, profile_export_fields
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)
from vntyper.scripts.profile_provenance import snapshot_decision_profile
from vntyper.scripts.summary_steps import STEP_KESTREL

pytestmark = pytest.mark.unit

IDENTITY = {
    "Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
    "Molecular_Identity_Status": "unique",
    "Equivalent_Representation_Count": 1,
    "Identity_Hypothesis_Count": 1,
}


def _explicit_profile() -> ResolvedDecisionProfile:
    packaged = load_packaged_decision_profile()
    document = json.loads(packaged.canonical_bytes)
    document["profile_id"] = "cohort-explicit"
    document["profile_revision"] = "study-7"
    document["profile_kind"] = "explicit-custom"
    inventory = document["inventory"]
    assert isinstance(inventory, dict)
    inventory["/components/kestrel/duplicate_flagging/flag_name"]["value"] = "Study_Duplicate"
    return parse_decision_profile(canonical_json_bytes(document), packaged_document=packaged.document)


def _write_profile_sample(
    directory: Path,
    profile: ResolvedDecisionProfile,
    *,
    confidence: str,
    flag: str,
) -> Path:
    directory.mkdir(parents=True)
    payload = summary.start_summary(version="9.9.9", decision_profile=profile)
    payload["steps"] = [
        {
            "step": STEP_KESTREL,
            "parsed_result": {
                "data": [
                    {
                        "Motif": "5",
                        "Confidence": confidence,
                        "Flag": flag,
                        **IDENTITY,
                    }
                ]
            },
        }
    ]
    summary.write_summary(payload, directory / "pipeline_summary.json")
    snapshot_decision_profile(profile, directory / "provenance" / "decision_profile.json")
    return directory


def test_cohort_loader_verifies_and_returns_the_samples_recorded_profile(tmp_path: Path) -> None:
    profile = _explicit_profile()
    sample = _write_profile_sample(
        tmp_path / "sample",
        profile,
        confidence="High_Precision",
        flag="Not flagged",
    )

    rows, _, stats = load_pipeline_summary_for_sample(sample)

    expected = {
        "Decision_Profile_ID": profile.profile_id,
        "Decision_Profile_Revision": profile.profile_revision,
        "Decision_Profile_SHA256": profile.digest,
    }
    assert {key: rows[0][key] for key in expected} == expected
    assert {key: stats[key] for key in expected} == expected


def test_cohort_loader_rejects_a_tampered_schema_three_snapshot(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    profile = load_packaged_decision_profile()
    sample = _write_profile_sample(
        tmp_path / "sample",
        profile,
        confidence="High_Precision",
        flag="Not flagged",
    )
    snapshot = sample / "provenance" / "decision_profile.json"
    snapshot.write_bytes(profile.canonical_bytes.replace(b'"profile_revision":"2"', b'"profile_revision":"9"'))

    with pytest.raises(ValueError, match="decision profile"):
        load_pipeline_summary_for_sample(sample)
    assert "decision profile" in caplog.text.lower()


def test_mixed_profile_cohort_exports_hashes_and_separates_decision_metrics(tmp_path: Path) -> None:
    root = tmp_path / "cohort"
    packaged = load_packaged_decision_profile()
    explicit = _explicit_profile()
    _write_profile_sample(
        root / "packaged_sample",
        packaged,
        confidence="High_Precision",
        flag="Not flagged",
    )
    _write_profile_sample(
        root / "explicit_sample",
        explicit,
        confidence="Low_Precision",
        flag="flagged",
    )
    output = tmp_path / "out"

    cohort_summary.aggregate_cohort(
        input_paths=[str(root)],
        output_dir=str(output),
        summary_file="cohort_summary.html",
        config=load_config(None),
        additional_formats="csv",
    )

    results = pd.read_csv(output / "cohort_kestrel.csv", keep_default_na=False)
    stats = pd.read_csv(output / "cohort_stats.csv", keep_default_na=False)
    profile_columns = ["Decision_Profile_ID", "Decision_Profile_Revision", "Decision_Profile_SHA256"]
    assert all(column in results.columns for column in profile_columns)
    assert all(column in stats.columns for column in profile_columns)
    assert set(results["Decision_Profile_SHA256"]) == {packaged.digest, explicit.digest}
    assert set(stats["Decision_Profile_SHA256"]) == {packaged.digest, explicit.digest}

    rendered = (output / "cohort_summary.html").read_text(encoding="utf-8")
    assert "Decision profile groups" in rendered
    assert packaged.digest in rendered
    assert explicit.digest in rendered
    assert "packaged_sample" in rendered
    assert "explicit_sample" in rendered
    assert "Pooled decision-performance plots are suppressed across decision profile hashes." in rendered
    assert rendered.count("Profile-specific Kestrel summary") == 2
    assert rendered.count("Profile-specific adVNTR summary") == 2


def test_legacy_cohort_profile_fields_never_infer_the_current_package(tmp_path: Path) -> None:
    sample = tmp_path / "legacy"
    sample.mkdir()
    (sample / "pipeline_summary.json").write_text(json.dumps({"schema_version": 2, "steps": []}), encoding="utf-8")

    _, _, stats = load_pipeline_summary_for_sample(sample)

    legacy = "decision profile not recorded by legacy run"
    assert stats["Decision_Profile_ID"] == legacy
    assert stats["Decision_Profile_Revision"] == legacy
    assert stats["Decision_Profile_SHA256"] == legacy
    assert load_packaged_decision_profile().digest not in stats.values()


def test_profile_groups_preserve_sha_order_and_sample_membership() -> None:
    entries = [
        {
            "Sample": "one",
            "Decision_Profile_ID": "profile-a",
            "Decision_Profile_Revision": "1",
            "Decision_Profile_SHA256": "a" * 64,
        },
        {
            "Sample": "two",
            "Decision_Profile_ID": "profile-b",
            "Decision_Profile_Revision": "2",
            "Decision_Profile_SHA256": "b" * 64,
        },
        {
            "Sample": "three",
            "Decision_Profile_ID": "profile-a",
            "Decision_Profile_Revision": "1",
            "Decision_Profile_SHA256": "a" * 64,
        },
    ]

    groups = group_decision_profiles(entries)

    assert [group.sha256 for group in groups] == ["a" * 64, "b" * 64]
    assert groups[0].samples == ("one", "three")
    assert groups[1].samples == ("two",)


def test_profile_group_rejects_one_hash_with_conflicting_metadata() -> None:
    entries = [
        {
            "Sample": "one",
            "Decision_Profile_ID": "profile-a",
            "Decision_Profile_Revision": "1",
            "Decision_Profile_SHA256": "a" * 64,
        },
        {
            "Sample": "two",
            "Decision_Profile_ID": "different-id",
            "Decision_Profile_Revision": "1",
            "Decision_Profile_SHA256": "a" * 64,
        },
    ]

    with pytest.raises(ValueError, match="conflicting metadata"):
        group_decision_profiles(entries)


def test_profile_export_fields_rejects_an_unverified_value() -> None:
    with pytest.raises(ValueError, match="verified"):
        profile_export_fields(object())  # type: ignore[arg-type]


def test_profile_row_projection_copies_rows_and_rejects_incomplete_fields() -> None:
    row = {"Motif": "5"}
    fields = {
        "Decision_Profile_ID": "profile-a",
        "Decision_Profile_Revision": "1",
        "Decision_Profile_SHA256": "a" * 64,
    }

    projected = annotate_profile_rows([row], fields)

    assert projected == [{"Motif": "5", **fields}]
    assert row == {"Motif": "5"}
    with pytest.raises(ValueError, match="incomplete"):
        annotate_profile_rows([row], {"Decision_Profile_ID": "profile-a"})


def test_profile_group_requires_complete_nonempty_sample_fields() -> None:
    with pytest.raises(ValueError, match="Sample"):
        group_decision_profiles(
            [
                {
                    "Decision_Profile_ID": "profile-a",
                    "Decision_Profile_Revision": "1",
                    "Decision_Profile_SHA256": "a" * 64,
                }
            ]
        )

    with pytest.raises(ValueError, match="lowercase hexadecimal"):
        group_decision_profiles(
            [
                {
                    "Sample": "one",
                    "Decision_Profile_ID": "profile-a",
                    "Decision_Profile_Revision": "1",
                    "Decision_Profile_SHA256": "not-a-digest",
                }
            ]
        )
