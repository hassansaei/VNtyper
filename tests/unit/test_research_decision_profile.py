"""Runtime admission of a derived research caller profile.

A calibration run that derives cutoffs is useless unless the pipeline can consume
the result. These tests pin the supported research path and, just as importantly,
prove that admitting it does not weaken the approved-bundle or explicit-profile
gates that already exist.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

import pytest

from vntyper.scripts.calibration_caller_policy import (
    KESTREL_CALLER_POLICY_POINTERS,
    decode_caller_policy_values,
)
from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile
from vntyper.scripts.decision_profile import (
    load_packaged_decision_profile,
    resolve_decision_profile,
    resolve_research_decision_profile,
)
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit

LINKED_POINTERS = (
    "/components/kestrel/alt_filtering/gg_depth_score_threshold",
    "/components/kestrel/confidence_assignment/depth_score_thresholds/low",
    "/components/kestrel/confidence_assignment/reporting_floor",
)


def _research_profile_bytes(value: float = 0.002) -> bytes:
    """Build the exact bytes a calibration export would emit for a lowered floor."""
    packaged = load_packaged_decision_profile()
    inventory = json.loads(packaged.canonical_bytes)["inventory"]
    values = {pointer: inventory[pointer]["value"] for pointer in KESTREL_CALLER_POLICY_POINTERS}
    for pointer in LINKED_POINTERS:
        values[pointer] = value
    policy = decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": ["kestrel"],
            "values": values,
        }
    )
    profile = build_caller_generated_profile(
        policy,
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=20260915,
        generator_version="unit-test",
    )
    return profile.canonical_bytes


def test_research_profile_threads_every_linked_cutoff_into_the_kestrel_component(tmp_path: Path) -> None:
    path = tmp_path / "research-profile.json"
    path.write_bytes(_research_profile_bytes(0.002))

    resolved = resolve_research_decision_profile(path)
    kestrel = cast(Mapping[str, Any], resolved.components["kestrel"])
    confidence = cast(Mapping[str, Any], kestrel["confidence_assignment"])
    depth_scores = cast(Mapping[str, Any], confidence["depth_score_thresholds"])

    assert resolved.profile_kind == "generated"
    assert confidence["reporting_floor"] == 0.002
    assert depth_scores["low"] == 0.002
    assert cast(Mapping[str, Any], kestrel["alt_filtering"])["gg_depth_score_threshold"] == 0.002
    # The band edge is a separate axis and must not move as a side effect.
    assert depth_scores["high"] == 0.00515


def test_the_ordinary_decision_profile_flag_still_refuses_a_caller_generated_profile(tmp_path: Path) -> None:
    path = tmp_path / "research-profile.json"
    path.write_bytes(_research_profile_bytes())

    with pytest.raises(ValueError, match="calibration bundle"):
        resolve_decision_profile(path)


def test_research_admission_refuses_a_profile_that_is_not_caller_generated(tmp_path: Path) -> None:
    packaged = load_packaged_decision_profile()
    document = json.loads(packaged.canonical_bytes)
    document["profile_id"] = "unit-test-explicit"
    document["profile_revision"] = "test-1"
    document["profile_kind"] = "explicit-custom"
    path = tmp_path / "explicit.json"
    path.write_text(json.dumps(document))

    with pytest.raises(ValueError, match="research"):
        resolve_research_decision_profile(path)


def test_research_admission_reports_an_unreadable_path(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="cannot read"):
        resolve_research_decision_profile(tmp_path / "absent.json")


def test_run_configuration_accepts_the_research_profile_and_records_its_identity(tmp_path: Path) -> None:
    path = tmp_path / "research-profile.json"
    path.write_bytes(_research_profile_bytes(0.0031))

    configuration = resolve_run_configuration(research_profile=path)

    confidence = cast(Mapping[str, Any], configuration.kestrel["confidence_assignment"])

    assert configuration.decision_profile.profile_kind == "generated"
    assert confidence["reporting_floor"] == 0.0031
    assert configuration.caller_calibration is None


def test_research_profile_is_exclusive_with_the_explicit_profile(tmp_path: Path) -> None:
    path = tmp_path / "research-profile.json"
    path.write_bytes(_research_profile_bytes())

    with pytest.raises(ValueError, match="exclusive"):
        resolve_run_configuration(path, research_profile=path)


def test_research_profile_is_exclusive_with_an_approved_bundle(tmp_path: Path) -> None:
    path = tmp_path / "research-profile.json"
    path.write_bytes(_research_profile_bytes())

    with pytest.raises(ValueError, match="exclusive"):
        resolve_run_configuration(
            research_profile=path,
            calibration_bundle=tmp_path / "bundle",
            calibration_context=tmp_path / "context.json",
        )


def test_pipeline_parser_exposes_the_research_profile_flag() -> None:
    from vntyper.scripts.cli_parser import build_parser

    parser = build_parser()
    arguments = parser.parse_args(
        ["pipeline", "-o", "out", "--bam", "in.bam", "--research-decision-profile", "research.json"]
    )

    assert arguments.research_decision_profile == Path("research.json")
