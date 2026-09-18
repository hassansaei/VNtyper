"""Tests for the strict adVNTR 2.4 calibration ABI contracts."""

from __future__ import annotations

import json
import subprocess
from collections.abc import Mapping
from dataclasses import replace

import pytest

from vntyper.modules.advntr import advntr_calibration_policy as policy
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values

pytestmark = pytest.mark.unit


def capabilities() -> dict[str, object]:
    return {
        "schema_version": "advntr-capabilities-v1",
        "package_version": "2.4.0",
        "build_id": "a" * 64,
        "source_revision": "b" * 40,
        "capabilities": [
            "fit-background-installed-v1",
            "frameshift-calibration-capture-v1",
            "frameshift-run-local-policy-v1",
            "frameshift-calibration-capture-v2",
            "frameshift-replay-v1",
        ],
        "capture_schema_versions": [1, 2],
        "policy_schema_versions": [
            "advntr-frameshift-policy-v1",
            "advntr-frameshift-replay-policy-v1",
        ],
        "background_recipe_ids": ["recipe-v1"],
    }


def capture_policy() -> dict[str, object]:
    return {
        "schema_version": "advntr-runtime-capture-policy-v1",
        "parameters": {
            "platform": "illumina",
            "frameshift_mode": True,
            "is_haploid": False,
            "caller_mode": "legacy",
            "threads": 2,
            "minimum_read_length": None,
            "prune_reverse": False,
            "filter_adapter_readthrough": False,
            "minimum_read_match_ratio": 0.6,
            "minimum_relative_ru_coverage": None,
            "use_reference_alignment": True,
            "fully_covered_ru_only": False,
            "maximum_error_rate": 0.05,
            "legacy_error_rate": 0.01,
            "mapq_cutoff": 0,
            "base_quality_cutoff": 20,
            "maximum_low_quality_fraction": 0.1,
            "enhanced_hmm": True,
            "trained_hmms": False,
        },
    }


def caller_values(*, mode: str = "legacy") -> dict[str, object]:
    values: dict[str, object] = {
        "/components/kestrel/alt_filtering/gg_depth_score_threshold": 0.5,
        "/components/kestrel/confidence_assignment/reporting_floor": 0.5,
        "/components/kestrel/confidence_assignment/var_active_region_threshold": 1,
        "/components/kestrel/confidence_assignment/depth_score_thresholds/low": 0.2,
        "/components/kestrel/confidence_assignment/depth_score_thresholds/high": 0.8,
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/low": 1,
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_low": 2,
        "/components/kestrel/confidence_assignment/alt_depth_thresholds/mid_high": 3,
        "/components/advntr/calibrated_calling/mode": mode,
        "/components/advntr/calibrated_calling/cutoff": 0.002,
        "/components/advntr/calibrated_calling/minimum_read_support": 4,
        "/components/advntr/calibrated_calling/rare_unit_fraction": 0.2,
        "/components/advntr/calibrated_calling/adapter_filter": True,
        "/components/advntr/calibrated_calling/minimum_read_match_ratio": 0.7,
        "/components/advntr/calibrated_calling/prune_reverse": True,
    }
    return {
        "schema_version": "calibration-caller-policy-values-v1",
        "required_callers": ["advntr", "kestrel"],
        "values": values,
    }


def test_decodes_and_revalidates_exact_published_capabilities() -> None:
    decoded = policy.decode_advntr_capabilities(capabilities())
    pin = policy.AdvntrToolPin("2.4.0", "a" * 64, "b" * 40)

    assert policy.advntr_capabilities_document(decoded) == capabilities()
    assert policy.require_advntr_capabilities(decoded, pin) == decoded

    with pytest.raises(ValueError, match="canonical content"):
        policy.advntr_capabilities_document(replace(decoded, package_version="2.4.1"))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("extra", 1),
        ("capture_schema_versions", [2, 1]),
        ("capture_schema_versions", [1, True]),
        ("capabilities", ["frameshift-replay-v1"]),
        ("source_revision", "not-a-revision"),
    ],
)
def test_capability_contract_rejects_unknown_or_nonpublished_inventory(field: str, value: object) -> None:
    document = capabilities()
    document[field] = value
    with pytest.raises(ValueError):
        policy.decode_advntr_capabilities(document)


def test_probe_uses_shell_free_argv_and_refuses_failed_or_ambiguous_output() -> None:
    calls: list[tuple[tuple[str, ...], dict[str, object]]] = []

    def runner(argv: tuple[str, ...], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0, json.dumps(capabilities()) + "\n", "")

    pin = policy.AdvntrToolPin("2.4.0", "a" * 64, "b" * 40)
    result = policy.probe_advntr_capabilities(("/tool dir/python", "-m", "advntr"), pin, runner=runner)

    assert result.package_version == "2.4.0"
    assert calls == [
        (
            ("/tool dir/python", "-m", "advntr", "capabilities", "--json"),
            {"capture_output": True, "text": True, "check": False},
        )
    ]

    def failed(argv: tuple[str, ...], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(argv, 2, "", "private stderr")

    with pytest.raises(RuntimeError, match="capability probe failed"):
        policy.probe_advntr_capabilities(("advntr",), pin, runner=failed)


def test_capture_policy_projection_is_complete_immutable_and_strictly_typed() -> None:
    decoded = policy.decode_capture_policy(capture_policy())
    assert policy.capture_policy_document(decoded) == capture_policy()
    assert decoded.minimum_read_match_ratio == 0.6

    for field, value in (
        ("threads", True),
        ("minimum_read_match_ratio", 1),
        ("maximum_error_rate", float("inf")),
        ("mapq_cutoff", False),
    ):
        document = capture_policy()
        document["parameters"][field] = value  # type: ignore[index]
        with pytest.raises(ValueError):
            policy.decode_capture_policy(document)

    raw = capture_policy()
    raw["parameters"]["unknown"] = 1  # type: ignore[index]
    with pytest.raises(ValueError, match="every field"):
        policy.decode_capture_policy(raw)


def test_projects_calibrated_values_into_capture_and_replay_policy() -> None:
    baseline = policy.decode_capture_policy(capture_policy())
    caller = decode_caller_policy_values(caller_values(mode="exact"))

    projected = policy.capture_policy_for_caller(baseline, caller)
    replay = policy.replay_policy_document(projected, caller)

    assert projected.caller_mode == "exact"
    assert projected.minimum_relative_ru_coverage == 0.2
    assert projected.minimum_read_match_ratio == 0.7
    assert replay["capture_policy"] == policy.capture_policy_document(projected)
    assert replay["caller_policy"] == {
        "schema_version": "advntr-frameshift-policy-v1",
        "mode": "exact",
        "cutoff": 0.002,
        "minimum_read_support": 4,
    }


def test_calibrated_projection_requires_advntr_and_explicit_finite_ratio() -> None:
    baseline_document = capture_policy()
    baseline_document["parameters"]["minimum_read_match_ratio"] = None  # type: ignore[index]
    baseline = policy.decode_capture_policy(baseline_document)
    kestrel = caller_values()
    kestrel["required_callers"] = ["kestrel"]
    original_values = kestrel["values"]
    assert isinstance(original_values, Mapping)
    kestrel["values"] = {key: value for key, value in original_values.items() if "/advntr/" not in key}

    with pytest.raises(ValueError, match="requires adVNTR"):
        policy.capture_policy_for_caller(baseline, decode_caller_policy_values(kestrel))

    assert policy.capture_policy_document(baseline)["parameters"]["minimum_read_match_ratio"] is None  # type: ignore[index]
