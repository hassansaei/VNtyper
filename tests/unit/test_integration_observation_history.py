"""Immutable released observations and the 2.0.43 to 2.0.46 successors."""

import copy
import json
import subprocess
from pathlib import Path

import pytest

from scripts.integration_compatibility import validate_manifest
from scripts.integration_compatibility_observations import effective_contracts

pytestmark = pytest.mark.unit

#: The #266 message for a subthreshold Kestrel negative without adVNTR, as shipped in 2.0.24.
ISSUE_266_SUBTHRESHOLD_MESSAGE = (
    "No variant called by Kestrel. A Kestrel candidate below the reporting floor was identified and filtered out."
    "<br>The subthreshold candidate is not a call; it is reported so that this result can be distinguished from a "
    "sample in which nothing was found at all.<br>Note: adVNTR genotyping was not performed."
)


def test_issue_293_preserves_history_and_pins_the_shipped_report_observation() -> None:
    """Bind the exact two 40cf observations to unchanged history and report configuration."""
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())
    base = json.loads(
        subprocess.run(
            ["git", "show", "f9e57f73e10d88d0c27cc4c4e8501c892594f0db:tests/compatibility/real_success_baseline.json"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )
    report_config = json.loads(Path("vntyper/scripts/report_config.json").read_text())
    shipped = [
        rule
        for rule in report_config["screening_summary_rules"]
        if rule["conditions"] == {"kestrel_result": "negative_subthreshold", "advntr_result": "none"}
    ]
    observations = manifest["observation_sets"]
    overrides = observations[0]["report_overrides"]
    live = json.loads(Path("tests/test_data_config.json").read_text())
    live_by_identity = {
        (suite, case["test_name"]): case for suite, cases in live["integration_tests"].items() for case in cases
    }
    expected_identities = {
        ("bam_tests", "example_40cf_hg38_subset_fast_gdp_guard"),
        ("bam_tests", "example_40cf_hg38_subset_default"),
    }

    assert json.dumps(manifest["contracts"], sort_keys=True, separators=(",", ":")) == json.dumps(
        base["contracts"], sort_keys=True, separators=(",", ":")
    )
    assert observations[:-4] == [
        {
            "version": "2.0.24",
            "provenance_commit": "f9e57f73e10d88d0c27cc4c4e8501c892594f0db",
            "extends": None,
            "report_overrides": overrides,
        },
        {
            "version": "2.0.25",
            "provenance_commit": "ec62d8f4e02212634b63d399275ba40d9e24fc1d",
            "extends": "2.0.24",
            "report_overrides": [],
        },
        {
            "version": "2.0.26",
            "provenance_commit": "c503c186e55edfc6b1bd140ca4ae9101551254e3",
            "extends": "2.0.25",
            "report_overrides": [],
        },
        {
            "version": "2.0.27",
            "provenance_commit": "e71e1ef4b98dfb6250ade623914cfc9eeccd8843",
            "extends": "2.0.26",
            "report_overrides": [],
        },
        {
            "version": "2.0.28",
            "provenance_commit": "65707249532efa47e09b6f91e2ce31b162f02627",
            "extends": "2.0.27",
            "report_overrides": [],
        },
        {
            "version": "2.0.29",
            "provenance_commit": "39c2afd7c9b03203dfd833773811d5e9cdb31f00",
            "extends": "2.0.28",
            "report_overrides": [],
        },
        {
            "version": "2.0.30",
            "provenance_commit": "48e876cdd0bcb5299339fe06db5faea9d7214d61",
            "extends": "2.0.29",
            "report_overrides": [],
        },
        {
            "version": "2.0.31",
            "provenance_commit": "2bb526a2c058f0ee7cd20a6f82f8b7459393be5b",
            "extends": "2.0.30",
            "report_overrides": [],
        },
        {
            "version": "2.0.32",
            "provenance_commit": "44bb06954a4738b7739c81ab3a0bdb8e8fa01396",
            "extends": "2.0.31",
            "report_overrides": [],
        },
        {
            "version": "2.0.33",
            "provenance_commit": "07d66ada02dcd80a702a18577bab609e6b8eaaae",
            "extends": "2.0.32",
            "report_overrides": [],
            "kestrel_overrides": [
                {
                    "suite": "fastq_tests",
                    "test_name": "example_6449_hg19_subset_fastq_shark",
                    "kestrel": {
                        "Estimated_Depth_AlternateVariant": {
                            "value": 268,
                            "tolerance": {
                                "kind": "percentage",
                                "value": 0,
                            },
                        },
                        "Estimated_Depth_Variant_ActiveRegion": {
                            "value": 16583,
                            "tolerance": {
                                "kind": "percentage",
                                "value": 0,
                            },
                        },
                        "Depth_Score": {
                            "value": 0.016161128866911897,
                            "tolerance": {
                                "kind": "percentage",
                                "value": 0,
                            },
                        },
                        "Confidence": {
                            "value": "High_Precision*",
                            "tolerance": None,
                        },
                    },
                },
                {
                    "suite": "fastq_tests",
                    "test_name": "example_6449_hg19_subset_single_fastq",
                    "kestrel": {
                        "Estimated_Depth_AlternateVariant": {
                            "value": 262,
                            "tolerance": {
                                "kind": "percentage",
                                "value": 0,
                            },
                        },
                        "Estimated_Depth_Variant_ActiveRegion": {
                            "value": 11484,
                            "tolerance": {
                                "kind": "percentage",
                                "value": 0,
                            },
                        },
                        "Depth_Score": {
                            "value": 0.022814350400557296,
                            "tolerance": {
                                "kind": "percentage",
                                "value": 0,
                            },
                        },
                        "Confidence": {
                            "value": "High_Precision*",
                            "tolerance": None,
                        },
                    },
                },
            ],
        },
        {
            "version": "2.0.34",
            "provenance_commit": "72176945d341b78d4361ddf0fe54eb076de2ca01",
            "extends": "2.0.33",
            "report_overrides": [],
        },
        {
            "version": "2.0.35",
            "provenance_commit": "7c584c3cc41bd6d2f66d1c5a5cd332896e5a54dd",
            "extends": "2.0.34",
            "report_overrides": [],
        },
        {
            "version": "2.0.36",
            "provenance_commit": "681f1b77209693ce58ab0851c55847a6daaac935",
            "extends": "2.0.35",
            "report_overrides": [],
        },
        {
            "version": "2.0.37",
            "provenance_commit": "5844282379551b9f3666c420843d6bbb4a7d33fa",
            "extends": "2.0.36",
            "report_overrides": [],
        },
        {
            "version": "2.0.38",
            "provenance_commit": "0ce7e755884fc8853da60220074deafd3a2aa0e7",
            "extends": "2.0.37",
            "report_overrides": [],
        },
        {
            "version": "2.0.39",
            "provenance_commit": "6bbd3e01e926ec63638f10abbbef32db3631a8b3",
            "extends": "2.0.38",
            "report_overrides": [],
        },
        {
            "version": "2.0.40",
            "provenance_commit": "bad094027f4fc11d3cef5440413e88f1f223c23b",
            "extends": "2.0.39",
            "report_overrides": [],
        },
        {
            "version": "2.0.41",
            "provenance_commit": "b69203eb51bdda46b9f67407b0c2fc7696d6cdf4",
            "extends": "2.0.40",
            "report_overrides": [],
        },
        {
            "version": "2.0.42",
            "provenance_commit": "22ef08be4b675835fdda26d49bdc4537d94069f8",
            "extends": "2.0.41",
            "report_overrides": observations[-5]["report_overrides"],
            "coverage_qc_overrides": [
                {"suite": "bam_tests", "test_name": "example_7a61_hg19_subset_fast", "coverage_qc": "REDUCED"}
            ],
        },
    ]
    assert len(shipped) == 1
    assert len(overrides) == 2
    assert {(row["suite"], row["test_name"]) for row in overrides} == expected_identities
    # History is immutable: 2.0.24 keeps the message as it was shipped then.
    assert all(row["report"] == [ISSUE_266_SUBTHRESHOLD_MESSAGE] for row in overrides)
    # 2.0.42 reworded every message for a call, the two #266 subthreshold messages and the
    # plain Kestrel-only negative, so each real success that asserts one carries an
    # override, and each override is what the live declaration now asserts.
    reworded = observations[-5]["report_overrides"]
    reworded_identities = {(row["suite"], row["test_name"]) for row in reworded}
    assert len(reworded) == 24
    assert expected_identities <= reworded_identities
    assert all(
        row["report"] == [shipped[0]["message"]]
        for row in reworded
        if (row["suite"], row["test_name"]) in expected_identities
    )
    assert (
        sum(
            row["report"][0].startswith("Kestrel called a frameshift variant with high precision.<br>")
            for row in reworded
        )
        == 21
    )
    assert all(
        live_by_identity[(row["suite"], row["test_name"])]["report_assertions"] == row["report"] for row in reworded
    )


def test_release_2043_inherits_all_2042_outcomes_without_overrides() -> None:
    """A tooling-only release preserves every effective real-data expectation."""
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())
    resources = json.loads(Path("tests/test_data_config.json").read_text())
    assert manifest["observation_sets"][-4] == {
        "version": "2.0.43",
        "provenance_commit": "19f3bc6edd447683f5d0f8eea6a76aeda4ec67dc",
        "extends": "2.0.42",
        "report_overrides": [],
    }
    current = copy.deepcopy(manifest)
    del current["observation_sets"][-3:]
    previous = copy.deepcopy(current)
    previous["observation_sets"].pop()
    assert effective_contracts(current, validate_manifest(current, resources), "2.0.43") == effective_contracts(
        previous, validate_manifest(previous, resources), "2.0.42"
    )


def test_release_2044_inherits_all_2043_outcomes_without_overrides() -> None:
    """A report-wording and logging release preserves every effective real-data expectation."""
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())
    resources = json.loads(Path("tests/test_data_config.json").read_text())
    assert manifest["observation_sets"][-3] == {
        "version": "2.0.44",
        "provenance_commit": "7f035392f3efbd6055b39b80dddcd0cafde87ead",
        "extends": "2.0.43",
        "report_overrides": [],
    }
    current = copy.deepcopy(manifest)
    del current["observation_sets"][-2:]
    previous = copy.deepcopy(current)
    previous["observation_sets"].pop()
    assert effective_contracts(current, validate_manifest(current, resources), "2.0.44") == effective_contracts(
        previous, validate_manifest(previous, resources), "2.0.43"
    )


def test_release_2045_inherits_all_2044_outcomes_without_overrides() -> None:
    """A nomenclature-reconciliation release preserves every effective real-data expectation."""
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())
    resources = json.loads(Path("tests/test_data_config.json").read_text())
    assert manifest["observation_sets"][-2] == {
        "version": "2.0.45",
        "provenance_commit": "d89b021e9e1b764bf6d99089e88ba3895c7c619c",
        "extends": "2.0.44",
        "report_overrides": [],
    }
    current = copy.deepcopy(manifest)
    current["observation_sets"].pop()
    previous = copy.deepcopy(current)
    previous["observation_sets"].pop()
    assert effective_contracts(current, validate_manifest(current, resources), "2.0.45") == effective_contracts(
        previous, validate_manifest(previous, resources), "2.0.44"
    )


def test_release_2046_changes_only_the_cram_step_lists() -> None:
    """#342 adds the header step to CRAM summaries; every other expectation is inherited."""
    manifest = json.loads(Path("tests/compatibility/real_success_baseline.json").read_text())
    resources = json.loads(Path("tests/test_data_config.json").read_text())
    steps = ["BAM Header Parsing", "CRAM to FASTQ Conversion", "Coverage Calculation", "Kestrel Genotyping"]
    latest = manifest["observation_sets"][-1]
    assert {key: latest[key] for key in ("version", "provenance_commit", "extends", "report_overrides")} == {
        "version": "2.0.46",
        "provenance_commit": "b160e0833d87655f9b7f18654719a1ddc7ec264c",
        "extends": "2.0.45",
        "report_overrides": [],
    }
    assert {row["suite"] for row in latest["summary_steps_overrides"]} == {"cram_tests"}
    assert all(row["steps"] == steps for row in latest["summary_steps_overrides"])
    previous = copy.deepcopy(manifest)
    previous["observation_sets"].pop()
    after = effective_contracts(manifest, validate_manifest(manifest, resources), "2.0.46")
    before = effective_contracts(previous, validate_manifest(previous, resources), "2.0.45")
    changed = {identity for identity in after if after[identity] != before[identity]}
    assert changed == {("cram_tests", row["test_name"]) for row in latest["summary_steps_overrides"]}
    for identity in changed:
        restored = copy.deepcopy(after[identity])
        restored["outcomes"]["summary"]["steps"] = steps[1:]
        assert restored == before[identity]
