from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import MappingProxyType
from unittest.mock import patch

import pytest

from tests.builders import kestrel_config
from tests.unit.test_calibration_caller_policy import policy_document, policy_values
from tests.unit.test_calibration_kestrel_replay import _GG, _candidate, _capture, _raw
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, decode_caller_policy_values
from vntyper.scripts.calibration_cutoff_kestrel import (
    kestrel_grid_replay_document,
    replay_kestrel_grid,
)
from vntyper.scripts.calibration_kestrel_capture import kestrel_capture_document
from vntyper.scripts.calibration_kestrel_replay import kestrel_replay_selected_frame, replay_kestrel_capture
from vntyper.scripts.canonical_json import canonical_json_bytes

pytestmark = pytest.mark.unit


def _write_capture(
    path: Path,
    *,
    empty: bool = False,
    config: dict[str, object] | None = None,
    baseline: CallerPolicyValues | None = None,
) -> CallerPolicyValues:
    actual = config or kestrel_config()
    frame = _raw().iloc[0:0] if empty else _raw()
    capture = _capture(frame, actual, policy=baseline)
    path.write_bytes(canonical_json_bytes(kestrel_capture_document(capture)))
    return capture.baseline_policy


def _native_negative(path: Path) -> None:
    path.write_text(
        "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
        "Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
        "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n",
        encoding="utf-8",
    )


def _native_positive(path: Path, capture_path: Path) -> None:
    from vntyper.scripts.calibration_kestrel_capture import decode_kestrel_capture
    from vntyper.scripts.canonical_json import load_strict_json_object

    capture = decode_kestrel_capture(load_strict_json_object(capture_path.read_bytes()))
    replay = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    kestrel_replay_selected_frame(replay).to_csv(path, sep="\t", index=False)


def test_grid_replays_baseline_first_and_deduplicates_identical_kestrel_parameters(tmp_path: Path) -> None:
    capture_path = tmp_path / "capture.json"
    base_kestrel = _capture(_raw(), kestrel_config()).baseline_policy
    baseline_raw = policy_document()
    policy_values(baseline_raw).update(base_kestrel.values)
    baseline = decode_caller_policy_values(baseline_raw)
    _write_capture(capture_path, baseline=baseline)
    changed = _candidate(baseline, **{_GG: 0.02})
    alias_raw = baseline_raw
    alias_raw = {**alias_raw, "values": dict(policy_values(alias_raw))}
    policy_values(alias_raw)["/components/advntr/calibrated_calling/cutoff"] = 0.002
    alias = decode_caller_policy_values(alias_raw)

    with patch(
        "vntyper.scripts.calibration_cutoff_kestrel.replay_kestrel_capture",
        wraps=replay_kestrel_capture,
    ) as replay:
        result = replay_kestrel_grid(
            {"case-a": capture_path},
            {"ad-only": alias, "changed": changed},
        )

    assert result.policy_ids == ("baseline", "ad-only", "changed")
    assert result.unique_kestrel_parameter_count == 2
    assert replay.call_count == 2
    assert (
        result.observations["baseline"]["case-a"].replay_sha256
        == result.observations["ad-only"]["case-a"].replay_sha256
    )
    assert result.observations["baseline"]["case-a"].called_positive is True
    assert result.observations["baseline"]["case-a"].canonical_identity is not None
    assert type(result.observations).__name__ == type(MappingProxyType({})).__name__
    assert type(result.observations["baseline"]).__name__ == type(MappingProxyType({})).__name__


def test_complete_filtered_no_call_is_negative_but_empty_capture_without_native_is_unavailable(tmp_path: Path) -> None:
    filtered_path = tmp_path / "filtered.json"
    strict = kestrel_config(**{"confidence_assignment.reporting_floor": 1.0})
    _write_capture(filtered_path, config=strict)
    empty_path = tmp_path / "empty.json"
    _write_capture(empty_path, empty=True)

    filtered = replay_kestrel_grid({"case": filtered_path}, {})
    empty = replay_kestrel_grid({"case": empty_path}, {})

    assert filtered.observations["baseline"]["case"].disposition == "no-call"
    assert filtered.observations["baseline"]["case"].called_positive is False
    assert empty.observations["baseline"]["case"].disposition == "unassessable-no-candidates"
    assert empty.observations["baseline"]["case"].called_positive is None


def test_exact_native_negative_proves_empty_capture_negative_for_every_policy(tmp_path: Path) -> None:
    capture_path = tmp_path / "capture.json"
    baseline = _write_capture(capture_path, empty=True)
    native = tmp_path / "native.tsv"
    _native_negative(native)
    changed = _candidate(baseline, **{_GG: 0.02})

    result = replay_kestrel_grid(
        {"case": capture_path},
        {"changed": changed},
        native_paths={"case": native},
    )

    assert result.baseline_parity == {"case": "native-exact"}
    assert all(result.observations[policy]["case"].called_positive is False for policy in result.policy_ids)


def test_native_positive_requires_exact_baseline_selected_fields(tmp_path: Path) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    native = tmp_path / "native.tsv"
    _native_positive(native, capture_path)
    result = replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native})
    assert result.baseline_parity["case"] == "native-exact"

    raw = native.read_text(encoding="utf-8").replace("Low_Precision", "High_Precision")
    native.write_text(raw, encoding="utf-8")
    with pytest.raises(ValueError, match="native Kestrel baseline differs"):
        replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native})


def test_grid_binds_exact_file_inputs_and_revalidates_immutable_projection(tmp_path: Path) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    result = replay_kestrel_grid({"case": capture_path}, {})
    document = kestrel_grid_replay_document(result)

    assert document["schema_version"] == "calibration-kestrel-grid-replay-v1"
    assert document["input_sha256"] == result.input_sha256
    assert len(result.capture_file_sha256["case"]) == 64
    with pytest.raises(ValueError, match="canonical content or digest"):
        kestrel_grid_replay_document(replace(result, sha256="0" * 64))


def test_process_path_replays_each_case_and_returns_deterministic_result(tmp_path: Path) -> None:
    paths = {}
    for key in ("a", "b"):
        path = tmp_path / f"{key}.json"
        _write_capture(path)
        paths[key] = path

    sequential = replay_kestrel_grid(paths, {}, workers=1)
    parallel = replay_kestrel_grid(paths, {}, workers=2)
    assert kestrel_grid_replay_document(parallel) == kestrel_grid_replay_document(sequential)


@pytest.mark.parametrize("workers", [True, 0, -1, 1.5])
def test_grid_rejects_invalid_worker_counts(tmp_path: Path, workers: object) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    with pytest.raises(ValueError, match="workers"):
        replay_kestrel_grid({"case": capture_path}, {}, workers=workers)  # type: ignore[arg-type]


def test_grid_rejects_missing_malformed_or_unbound_inputs(tmp_path: Path) -> None:
    missing = tmp_path / "missing.json"
    with pytest.raises(ValueError, match="capture"):
        replay_kestrel_grid({"case": missing}, {})

    malformed = tmp_path / "malformed.json"
    malformed.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="capture"):
        replay_kestrel_grid({"case": malformed}, {})

    valid = tmp_path / "valid.json"
    baseline = _write_capture(valid)
    with pytest.raises(ValueError, match="exact capture roster"):
        replay_kestrel_grid({"case": valid}, {}, native_paths={"other": None})
    with pytest.raises(ValueError, match="reserved"):
        replay_kestrel_grid({"case": valid}, {"baseline": baseline})
    with pytest.raises(ValueError, match="CallerPolicyValues"):
        replay_kestrel_grid({"case": valid}, {"bad": object()})  # type: ignore[dict-item]


def test_grid_rejects_different_capture_baselines(tmp_path: Path) -> None:
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    _write_capture(first)
    _write_capture(second, config=kestrel_config(**{"confidence_assignment.reporting_floor": 0.01}))

    with pytest.raises(ValueError, match="baseline policies differ"):
        replay_kestrel_grid({"a": first, "b": second}, {})
