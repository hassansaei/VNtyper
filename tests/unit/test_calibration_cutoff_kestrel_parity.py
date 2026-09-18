"""Unit tests for native Kestrel baseline result parity validation."""

from __future__ import annotations

import pytest

from tests.builders import kestrel_config
from tests.unit.test_calibration_kestrel_replay import _capture, _raw
from vntyper.scripts.calibration_cutoff_kestrel_parity import (
    NATIVE_PARITY_REQUIRED,
    _cell_text,
    validate_native_kestrel_baseline_parity,
)
from vntyper.scripts.calibration_kestrel_replay import replay_kestrel_capture

pytestmark = pytest.mark.unit


def _called_replay_result():
    config = kestrel_config()
    capture = _capture(_raw(), config)
    return replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )


def _negative_replay_result():
    config = kestrel_config()
    capture = _capture(_raw().iloc[0:0], config)
    return replay_kestrel_capture(
        capture,
        capture.baseline_policy,
        capture_policy_sha256=capture.provenance.capture_policy_sha256,
    )


def test_cell_text_handling() -> None:
    assert _cell_text(None) == ""
    assert _cell_text(float("nan")) == ""
    assert _cell_text(42) == "42"
    assert _cell_text("val") == "val"
    assert _cell_text(0.5) == "0.5"


def test_validate_native_kestrel_baseline_parity_requires_bytes() -> None:
    replay = _negative_replay_result()
    with pytest.raises(ValueError, match="must be exact bytes"):
        validate_native_kestrel_baseline_parity("not bytes", replay)  # type: ignore[arg-type]


def test_validate_native_kestrel_baseline_parity_negative_match() -> None:
    replay = _negative_replay_result()
    # Valid TSV with single negative placeholder row
    negative_tsv = (
        b"Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
        b"Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
        b"None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
    )
    verdict = validate_native_kestrel_baseline_parity(negative_tsv, replay)
    assert verdict == "negative"


def test_validate_native_kestrel_baseline_parity_negative_mismatch_called_in_native() -> None:
    replay = _negative_replay_result()
    called_tsv = b"POS\tREF\tALT\tDepth_Score\tConfidence\n67\tG\tGG\t0.5\tHigh\n"
    with pytest.raises(ValueError, match="differs from capture replay"):
        validate_native_kestrel_baseline_parity(called_tsv, replay)


def test_validate_native_kestrel_baseline_parity_called_mismatch_negative_in_native() -> None:
    replay = _called_replay_result()
    negative_tsv = (
        b"Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
        b"Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
        b"None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n"
    )
    with pytest.raises(ValueError, match="differs from capture replay"):
        validate_native_kestrel_baseline_parity(negative_tsv, replay)


def test_validate_native_kestrel_baseline_parity_missing_required_columns() -> None:
    replay = _called_replay_result()
    # Missing Depth_Score and Confidence in native TSV
    tsv = b"POS\tREF\tALT\n67\tG\tGG\n"
    with pytest.raises(ValueError, match="missing required comparable columns"):
        validate_native_kestrel_baseline_parity(tsv, replay)


def test_validate_native_kestrel_baseline_parity_field_mismatch() -> None:
    replay = _called_replay_result()
    # Different Depth_Score
    tsv = b"POS\tREF\tALT\tDepth_Score\tConfidence\n67\tG\tGG\t0.9999\tHigh\n"
    with pytest.raises(ValueError, match="differs from capture replay selected fields"):
        validate_native_kestrel_baseline_parity(tsv, replay)


def test_validate_native_kestrel_baseline_parity_success_called() -> None:
    from vntyper.scripts.calibration_kestrel_replay import kestrel_replay_selected_frame

    replay = _called_replay_result()
    selected_row = kestrel_replay_selected_frame(replay).iloc[0]

    cols = list(NATIVE_PARITY_REQUIRED)
    header = "\t".join(cols)
    values = "\t".join(str(selected_row[c]) for c in cols)
    tsv = f"{header}\n{values}\n".encode()

    verdict = validate_native_kestrel_baseline_parity(tsv, replay)
    assert verdict == "called"
