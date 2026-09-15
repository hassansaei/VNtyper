from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType
from unittest.mock import patch

import pandas as pd
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


# The published production column set, taken from the header a real 2.0.35 run wrote to
# kestrel/kestrel_result.tsv. It carries every __Identity_* column but no __Calibration_*
# column, and it adds nomenclature and resolved-identity columns the replay never emits.
PRODUCTION_KESTREL_COLUMNS: tuple[str, ...] = (
    "Motifs",
    "POS",
    "REF",
    "ALT",
    "Sample",
    "Motif_sequence",
    "Variant",
    "Del",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "ref_len",
    "alt_len",
    "Frame_Score",
    "is_frameshift",
    "direction",
    "frameshift_amount",
    "is_valid_frameshift",
    "Depth_Score",
    "Confidence",
    "depth_confidence_pass",
    "__Identity_Raw_Representation_Key",
    "__Identity_Molecular_Identity",
    "__Identity_Translation_Status",
    "__Identity_Translation_Failure",
    "__Identity_Context_Diverges",
    "__Identity_Observation_Ordinal",
    "haplo_count",
    "alt_filter_pass",
    "motif_filter_pass",
    "Motif_fasta",
    "POS_fasta",
    "Motif",
    "Flag",
    "flag_filter_pass",
    "__Identity_Selected_Raw_Representation_Key",
    "__Identity_Equivalent_Representation_Count",
    "__Identity_Hypothesis_Count",
    "__Identity_Group_Blocking_Gates",
    "__Identity_Group_Flags",
    "__Identity_Selected_Observation_Ordinal",
    "__Identity_Group_Context_Diverges",
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
    "Molecular_Identity",
    "Molecular_Identity_Status",
    "Equivalent_Representation_Count",
    "Identity_Hypothesis_Count",
)


def _baseline_frame(capture_path: Path) -> pd.DataFrame:
    from vntyper.scripts.calibration_kestrel_capture import decode_kestrel_capture
    from vntyper.scripts.canonical_json import load_strict_json_object

    capture = decode_kestrel_capture(load_strict_json_object(capture_path.read_bytes()))
    replay = replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )
    return kestrel_replay_selected_frame(replay)


def _native_production(
    path: Path,
    capture_path: Path,
    *,
    changes: Mapping[str, str] | None = None,
    drop: Sequence[str] = (),
) -> None:
    """Write a native TSV that has the production shape, not the replay frame shape."""
    columns = [name for name in PRODUCTION_KESTREL_COLUMNS if name not in set(drop)]
    frame = _baseline_frame(capture_path)
    native = frame.drop(columns=[name for name in frame.columns if name not in set(columns)])
    for name in columns:
        if name not in native.columns:
            # Production-only annotation columns carry synthetic, never private, values.
            native[name] = f"synthetic-{name}"
    for name, value in (changes or {}).items():
        native[name] = value
    body = native[columns].to_csv(sep="\t", index=False)
    path.write_text(f"## VNtyper Kestrel result\n## VNtyper Version: 0.0.0-test\n{body}", encoding="utf-8")


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


def test_native_production_column_set_is_accepted_and_private_columns_are_excluded(tmp_path: Path) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    native = tmp_path / "native.tsv"
    _native_production(native, capture_path)

    header = native.read_text(encoding="utf-8").splitlines()[2].split("\t")
    assert header == list(PRODUCTION_KESTREL_COLUMNS)
    # The replay frame carries a calibration-private column production never publishes.
    assert "__Calibration_Source_Row_Ordinal" in _baseline_frame(capture_path).columns
    assert "__Calibration_Source_Row_Ordinal" not in header

    result = replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native})
    assert result.baseline_parity["case"] == "native-exact"
    assert result.observations["baseline"]["case"].called_positive is True


@pytest.mark.parametrize(
    "column,value", [("Depth_Score", "0.9999"), ("Confidence", "High_Precision"), ("POS", "424242")]
)
def test_native_row_disagreeing_on_a_decision_column_is_rejected(tmp_path: Path, column: str, value: str) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    native = tmp_path / "native.tsv"
    _native_production(native, capture_path)
    assert replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native}).baseline_parity == {
        "case": "native-exact"
    }

    _native_production(native, capture_path, changes={column: value})
    with pytest.raises(ValueError, match="native Kestrel baseline differs"):
        replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native})


@pytest.mark.parametrize("column", ["POS", "REF", "ALT", "Depth_Score", "Confidence"])
def test_native_without_the_required_decision_columns_is_rejected(tmp_path: Path, column: str) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    native = tmp_path / "native.tsv"
    _native_production(native, capture_path, drop=(column,))

    with pytest.raises(ValueError, match="missing required comparable columns"):
        replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native})


@pytest.mark.parametrize("column", ["Motifs", "Estimated_Depth_AlternateVariant"])
def test_native_lacking_an_optional_shared_column_is_still_compared(tmp_path: Path, column: str) -> None:
    capture_path = tmp_path / "capture.json"
    _write_capture(capture_path)
    native = tmp_path / "native.tsv"
    _native_production(native, capture_path, drop=(column,))
    assert replay_kestrel_grid({"case": capture_path}, {}, native_paths={"case": native}).baseline_parity == {
        "case": "native-exact"
    }

    _native_production(native, capture_path, drop=(column,), changes={"Depth_Score": "0.9999"})
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
