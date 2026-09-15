"""Cohort orchestration accepts independent truth targets without promotion claims."""

import json
from argparse import Namespace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

pytestmark = pytest.mark.unit


def inputs(tmp_path, *, target="auto", labels=True):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample_id\tbam\tassembly\tgenotype\tallele_1\tallele_2\n"
        + "".join(
            f"s{i}\ts{i}.bam\tGRCh38\t{'positive' if i % 2 else 'negative'}\t{str(i + 10) if labels else ''}\t{str(i + 20) if labels else ''}\n"
            for i in range(8)
        )
    )
    output = tmp_path / "out"
    output.mkdir()
    return Namespace(
        manifest=manifest,
        reference=tmp_path / "ref.fa",
        target=target,
        folds=4,
        seed=7,
        count_convention="source-reported",
        caller_runs=None,
        caller_policies=None,
    ), output


def test_both_targets_produce_held_out_length_and_explicit_unavailable_callers(tmp_path):
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path)
    training_sets = []

    def read(path, _reference, *, assembly):
        return SimpleNamespace(values=(float(path.stem[1:]),), reasons=(), sha256=path.stem)

    def fit(measurements, targets, *, count_convention):
        training_sets.append(tuple(row.values[0] for row in measurements))
        x = np.array([[1, row.values[0]] for row in measurements])
        return np.linalg.lstsq(x, targets, rcond=None)[0]

    def predict(measurement, model):
        return SimpleNamespace(estimated_repeat_count=float(model[0] + model[1] * measurement.values[0]), reasons=())

    with patch.object(
        module,
        "_length_services",
        return_value=(
            read,
            fit,
            predict,
            lambda model: {"coefficients": model.tolist()},
            lambda measurement: measurement.reasons,
        ),
    ):
        assert module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())
    assert result["length"]["metrics"]["mae"] == pytest.approx(0, abs=1e-10)
    assert result["length"]["metrics"]["baseline_mae_paired"] > 0
    assert len(training_sets) == 5  # Four training-only fits, then final all-data research model.
    assert all(len(rows) == 6 for rows in training_sets[:-1])
    assert result["callers"]["baseline"]["no_calls"] == 8
    assert result["callers"]["status"] == "unavailable"
    assert not (output / "portable-approval.json").exists()
    assert (output / "length-model.json").is_file()
    assert "exploratory" in (output / "report.html").read_text().lower()


def test_mutation_only_needs_no_reference_and_never_opens_length_reads(tmp_path):
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path, labels=False)
    args.reference = None
    with patch.object(module, "_length_services", side_effect=AssertionError("length opened")):
        assert not module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())
    assert result["targets"] == ["callers"]
    assert not (output / "length-model.json").exists()


def test_held_out_prediction_warnings_are_reported_without_excluding_predictions(tmp_path):
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path, target="length")
    warnings = ("feature_A_outside_training_range", "synthetic_<range>")
    with patch.object(
        module,
        "_length_services",
        return_value=(
            lambda *args, **kwargs: SimpleNamespace(reasons=()),
            lambda *args, **kwargs: {},
            lambda *args, **kwargs: SimpleNamespace(estimated_repeat_count=40.0, reasons=(), warnings=warnings),
            lambda model: model,
            lambda measurement: (),
        ),
    ):
        assert module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())["length"]
    assert result["metrics"]["predicted"] == 8
    assert result["measurement_reasons"] == {}
    assert result["prediction_warning_counts"] == dict.fromkeys(warnings, 8)
    assert all(row["prediction_warnings"] == list(warnings) and row["prediction"] == 40 for row in result["rows"])
    report = (output / "report.html").read_text()
    assert "feature_A_outside_training_range" in report
    assert "synthetic_&lt;range&gt;" in report
    assert "synthetic_<range>" not in report
    assert "Predictions with warnings: 8" in report


def test_length_target_requires_reference_before_feature_io(tmp_path):
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path, target="length")
    args.reference = None
    with (
        patch.object(module, "_length_services", side_effect=AssertionError("opened")),
        pytest.raises(ValueError, match="reference"),
    ):
        module.run_cohort_calibration(args, output)
    assert not tuple(output.iterdir())


def test_explicit_length_without_truth_writes_unavailable_result(tmp_path):
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path, target="length", labels=False)
    args.reference = None
    assert not module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())
    assert result["length"]["status"] == "unavailable"


def test_real_shared_length_fit_and_predictor_on_invented_measurements(tmp_path):
    from tests.unit.test_calibration_standard_length import _measurement
    from vntyper.scripts import calibration_cohort as module
    from vntyper.scripts.length_standard_model import decode_standard_length_model

    args, output = inputs(tmp_path, target="length")
    measurements = {f"s{i}": _measurement(10 + i * 5) for i in range(8)}
    with patch(
        "vntyper.scripts.length_standard_io.read_standard_length_features",
        side_effect=lambda path, ref, assembly: measurements[path.stem],
    ):
        assert module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())
    assert result["length"]["metrics"]["predicted"] == 8
    assert result["length"]["metrics"]["mae"] < 0.01
    model = decode_standard_length_model(json.loads((output / "length-model.json").read_text()))
    assert model.count_convention == "source-reported"
    assert len(model.coefficients) == 13


def test_report_escapes_sample_identifiers_and_retains_failed_measurements(tmp_path):
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path, target="length")
    text = args.manifest.read_text().replace("s0\t", "<script>alert(1)</script>\t")
    args.manifest.write_text(text)

    def read(path, reference, assembly):
        raise ValueError("invented unavailable evidence")

    def fit(*args, **kwargs):
        raise ValueError("no measurements")

    with patch.object(module, "_length_services", return_value=(read, fit, object(), object(), object())):
        assert not module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())
    assert result["length"]["metrics"]["unavailable"] == 8
    assert "<script>alert(1)</script>" not in (output / "report.html").read_text()
    assert "&lt;script&gt;" in (output / "report.html").read_text()
    assert not (output / "length-model.json").exists()


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("manifest", "bad"),
        ("target", "bad"),
        ("count_convention", "bad"),
        ("folds", 1),
        ("seed", True),
        ("caller_runs", "bad"),
    ],
)
def test_invalid_arguments_fail_before_output(tmp_path, field, value):
    from vntyper.scripts.calibration_cohort import run_cohort_calibration

    args, output = inputs(tmp_path)
    setattr(args, field, value)
    with pytest.raises(ValueError):
        run_cohort_calibration(args, output)
    assert not tuple(output.iterdir())


def test_no_clobber_staging_guard(tmp_path):
    from vntyper.scripts.calibration_cohort import run_cohort_calibration

    args, output = inputs(tmp_path)
    (output / "keep").write_text("kept")
    with pytest.raises(ValueError, match="empty staged"):
        run_cohort_calibration(args, output)
    assert (output / "keep").read_text() == "kept"


def test_low_qc_row_does_not_poison_other_samples_training_folds(tmp_path):
    from tests.unit.test_calibration_standard_length import _measurement
    from vntyper.scripts import calibration_cohort as module

    args, output = inputs(tmp_path, target="length")
    from vntyper.scripts.length_standard_features import (
        decode_standard_length_measurement,
        encode_standard_length_measurement,
    )

    measurements = {f"s{i}": _measurement(10 + i * 5) for i in range(8)}
    raw = encode_standard_length_measurement(measurements["s0"])
    raw["qc"]["invariant_supporting_fragments"] = 0
    from vntyper.scripts.canonical_json import canonical_sha256

    raw["sha256"] = canonical_sha256({key: value for key, value in raw.items() if key != "sha256"})
    measurements["s0"] = decode_standard_length_measurement(raw)
    with patch(
        "vntyper.scripts.length_standard_io.read_standard_length_features",
        side_effect=lambda path, ref, assembly: measurements[path.stem],
    ):
        assert module.run_cohort_calibration(args, output)
    result = json.loads((output / "metrics.json").read_text())["length"]
    assert result["metrics"]["eligible"] == 8
    assert result["metrics"]["predicted"] == 7
    assert result["measurement_reasons"]["s0"]
    assert all("s0" not in fold["fitting"] for fold in result["folds"])
    assert all(fold["status"] == "fitted" for fold in result["folds"])
