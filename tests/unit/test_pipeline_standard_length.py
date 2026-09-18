"""Standard research length inference stays independent from mutation outcomes."""

import json
from dataclasses import replace
from pathlib import Path
from unittest.mock import Mock

import pytest

from vntyper.scripts import pipeline_standard_length as subject

pytestmark = pytest.mark.unit


_DEFAULT_CONFIG = {"reference_data": {"standard_length_model_grch38": "reference/grch38-standard-length-model-v1.json"}}


def _model(path=None, *, expected_source=None):
    from tests.unit.test_length_standard_model import _model_document
    from vntyper.scripts.length_standard_model import decode_standard_length_model

    document = _model_document()
    document["model_source"] = (
        "packaged-research"
        if path is None or path == Path("reference/grch38-standard-length-model-v1.json")
        else "local-research"
    )
    model = decode_standard_length_model(document)
    if expected_source is not None and model.model_source != expected_source:
        raise ValueError(f"standard length model source must be {expected_source}")
    return model


@pytest.mark.parametrize("change", ["model", "source", "sha256"])
def test_replaced_configuration_cannot_reuse_resume_identity(monkeypatch: pytest.MonkeyPatch, change: str) -> None:
    from vntyper.scripts.length_standard_model import decode_standard_length_model, encode_standard_length_model

    monkeypatch.setattr(subject, "_load_model", _model)
    config = subject.resolve_standard_length_configuration(
        _DEFAULT_CONFIG, enabled=True, model_path=None, approved_enabled=False
    )
    if change == "model":
        assert config.model is not None
        document = encode_standard_length_model(config.model)
        document["intercept"] = 99.0
        forged = replace(config, model=decode_standard_length_model(document))
    elif change == "source":
        forged = replace(config, source="local-research")
    else:
        forged = replace(config, sha256="f" * 64)
    with pytest.raises(ValueError, match="configuration"):
        subject.standard_length_summary(forged, None, None, ("unsupported-assembly",))


def test_missing_configuration_preserves_legacy_disabled_behavior() -> None:
    result = subject.resolve_standard_length_configuration({}, enabled=None, model_path=None, approved_enabled=False)
    assert result.enabled is False
    assert result.model is None


@pytest.mark.parametrize("override", [None, False])
def test_approved_path_prevents_automatic_duplicate_measurement(override: bool | None) -> None:
    result = subject.resolve_standard_length_configuration(
        {"length_estimation": {"enabled": True}}, enabled=override, model_path=None, approved_enabled=True
    )
    assert result.enabled is False


@pytest.mark.parametrize(
    ("enabled", "model_path", "approved"),
    [(True, None, True), (None, Path("model.json"), True), (False, Path("model.json"), False)],
)
def test_conflicting_explicit_model_requests_fail_before_loading(
    monkeypatch: pytest.MonkeyPatch, enabled: bool | None, model_path: Path | None, approved: bool
) -> None:
    loader = Mock(side_effect=AssertionError("must not read model"))
    monkeypatch.setattr(subject, "_load_model", loader)
    with pytest.raises(ValueError, match="conflict"):
        subject.resolve_standard_length_configuration(
            {}, enabled=enabled, model_path=model_path, approved_enabled=approved
        )
    loader.assert_not_called()


@pytest.mark.parametrize("invalid", ["true", 1, {}, None])
def test_malformed_config_enabled_is_rejected(invalid: object) -> None:
    with pytest.raises(ValueError, match="boolean"):
        subject.resolve_standard_length_configuration(
            {"length_estimation": {"enabled": invalid}}, enabled=None, model_path=None, approved_enabled=False
        )


def test_configuration_identity_binds_model_and_source(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subject, "_load_model", _model)
    packaged = subject.resolve_standard_length_configuration(
        {"length_estimation": {"enabled": True}, **_DEFAULT_CONFIG},
        enabled=None,
        model_path=None,
        approved_enabled=False,
    )
    local = subject.resolve_standard_length_configuration(
        {}, enabled=None, model_path=Path("model.json"), approved_enabled=False
    )
    assert packaged.source == "packaged-research"
    assert local.source == "local-research"
    assert packaged.sha256 != local.sha256


def test_unavailable_summary_keeps_approval_absent(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subject, "_load_model", _model)
    configuration = subject.resolve_standard_length_configuration(
        _DEFAULT_CONFIG, enabled=True, model_path=None, approved_enabled=False
    )
    result = subject.standard_length_summary(configuration, None, None, ("unsupported-assembly",))
    assert result["length_estimation_status"] == "unavailable"
    assert result["estimated_total_repeat_count"] is None
    assert result["length_portable_approval_sha256"] is None
    assert result["length_calibration_id"] is None
    assert result["length_model_source"] == "packaged-research"
    assert result["length_count_convention"] == "source-reported"


def test_unsupported_assembly_skips_reader(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subject, "_load_model", _model)
    configuration = subject.resolve_standard_length_configuration(
        _DEFAULT_CONFIG, enabled=True, model_path=None, approved_enabled=False
    )
    runner = subject.StandardLengthRunner(configuration, "hg19", None, Path("/work"))
    runner(Mock())
    assert runner.summary["length_estimation_reasons"] == ["unsupported-assembly"]
    with pytest.raises(ValueError, match="once"):
        runner(Mock())


@pytest.mark.parametrize(
    ("bwa_reference", "explicit_reference", "expected_reference"),
    [
        ("/refs/hg38.fa", None, "/refs/hg38.fa"),
        (None, "/refs/explicit.fa", "/refs/explicit.fa"),
        ("/refs/hg38.fa", "/refs/explicit.fa", "/refs/explicit.fa"),
    ],
)
def test_actual_pipeline_records_standard_model_inside_retained_lifetime(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    bwa_reference: str | None,
    explicit_reference: str | None,
    expected_reference: str,
) -> None:
    from tests.support.pipeline_harness import run_pipeline_under_harness
    from tests.unit.test_length_standard_model import _measurement, _model_document
    from vntyper.scripts.length_standard_model import decode_standard_length_model

    model = decode_standard_length_model(_model_document())
    measured = _measurement()
    monkeypatch.setattr(subject, "_load_model", lambda _path, **_kw: model)
    configuration = subject.resolve_standard_length_configuration(
        {}, enabled=True, model_path=tmp_path / "model.json", approved_enabled=False
    )
    calls = []

    def read(input_path, reference, *, assembly, index_path):
        calls.append((input_path, reference, assembly, index_path))
        return measured

    monkeypatch.setattr("vntyper.scripts.length_standard_io.read_standard_length_features", read)
    harness = run_pipeline_under_harness(
        tmp_path / "out",
        reference_assembly="hg38",
        standard_length_configuration=configuration,
        header="@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n",
        bwa_reference=bwa_reference,
        reference_fasta=explicit_reference,
    )
    assert len(calls) == 1
    assert calls[0][1] == Path(expected_reference)
    assert calls[0][2] == "hg38"
    assert calls[0][3] is not None
    recorded = json.loads((harness.output_dir / "pipeline_summary.json").read_text())
    assert recorded["estimated_total_repeat_count"] == pytest.approx(110)
    assert recorded["length_model_source"] == "local-research"
    assert recorded["length_model_sha256"] == model.sha256
    assert recorded["length_portable_approval_sha256"] is None
    assert recorded["analysis_settings"]["length_configuration_sha256"] == configuration.sha256


@pytest.mark.parametrize("error", [OSError("missing"), RuntimeError("changed"), ValueError("wrong reference")])
def test_optional_reader_failure_yields_unavailable_length(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error: Exception
) -> None:
    from tests.unit.test_pipeline_length_execution import _plan

    monkeypatch.setattr(subject, "_load_model", _model)
    config = subject.resolve_standard_length_configuration(
        _DEFAULT_CONFIG, enabled=True, model_path=None, approved_enabled=False
    )
    monkeypatch.setattr("vntyper.scripts.length_standard_io.read_standard_length_features", Mock(side_effect=error))
    runner = subject.StandardLengthRunner(config, "hg38", tmp_path / "reference.fa", tmp_path)
    runner(_plan())
    assert runner.summary["estimated_total_repeat_count"] is None
    assert runner.summary["length_estimation_reasons"] == ["measurement-unavailable"]


def test_missing_reference_model_file_fails_with_path(tmp_path: Path) -> None:
    nonexistent = tmp_path / "missing_model.json"
    with pytest.raises(ValueError, match=f"standard length model path must be an existing file: '{nonexistent}'"):
        subject.resolve_standard_length_configuration(
            {"reference_data": {"standard_length_model_grch38": str(nonexistent)}},
            enabled=True,
            model_path=None,
            approved_enabled=False,
        )


def test_missing_reference_configuration_fails_when_enabled() -> None:
    with pytest.raises(ValueError, match="standard length model reference is not configured in reference_data"):
        subject.resolve_standard_length_configuration({}, enabled=True, model_path=None, approved_enabled=False)


def test_explicit_packaged_model_path_resolves_packaged_source(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subject, "_load_model", _model)
    config = subject.resolve_standard_length_configuration(
        {},
        enabled=None,
        model_path=Path("reference/grch38-standard-length-model-v1.json"),
        approved_enabled=False,
    )
    assert config.enabled is True
    assert config.source == "packaged-research"
    assert config.model is not None
    assert config.model.model_source == "packaged-research"
    encoded = subject.encode_standard_length_configuration(config)
    assert encoded["source"] == "packaged-research"


def test_packaged_reference_model_enforces_packaged_source(monkeypatch: pytest.MonkeyPatch) -> None:
    from tests.unit.test_length_standard_model import _model_document
    from vntyper.scripts.length_standard_model import decode_standard_length_model

    doc = _model_document()
    doc["model_source"] = "local-research"
    wrong_model = decode_standard_length_model(doc)
    monkeypatch.setattr(
        subject,
        "_load_model",
        lambda _p, *, expected_source=None: (
            wrong_model
            if expected_source is None
            else (_ for _ in ()).throw(ValueError(f"standard length model source must be {expected_source}"))
        ),
    )
    with pytest.raises(ValueError, match="source must be packaged-research"):
        subject.resolve_standard_length_configuration(
            {"length_estimation": {"enabled": True}, **_DEFAULT_CONFIG},
            enabled=None,
            model_path=None,
            approved_enabled=False,
        )


def test_resolve_configuration_default_and_explicit_warning_threshold(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subject, "_load_model", _model)
    # Default from config.json or absent
    config_default = subject.resolve_standard_length_configuration(
        _DEFAULT_CONFIG, enabled=True, model_path=None, approved_enabled=False
    )
    assert config_default.warning_threshold == 110.0

    # Config-specified threshold
    config_specified = subject.resolve_standard_length_configuration(
        {"length_estimation": {"warning_threshold": 125.5}, **_DEFAULT_CONFIG},
        enabled=True,
        model_path=None,
        approved_enabled=False,
    )
    assert config_specified.warning_threshold == 125.5

    # CLI explicit override takes precedence over config
    config_override = subject.resolve_standard_length_configuration(
        {"length_estimation": {"warning_threshold": 125.5}, **_DEFAULT_CONFIG},
        enabled=True,
        model_path=None,
        approved_enabled=False,
        warning_threshold=95.0,
    )
    assert config_override.warning_threshold == 95.0


@pytest.mark.parametrize("invalid", [0, -10.0, True, False, float("nan"), float("inf"), "110"])
def test_resolve_configuration_rejects_invalid_warning_threshold(invalid: object) -> None:
    with pytest.raises(ValueError, match="standard length warning threshold must be a positive number"):
        subject.resolve_standard_length_configuration(
            _DEFAULT_CONFIG,
            enabled=True,
            model_path=None,
            approved_enabled=False,
            warning_threshold=invalid,  # type: ignore[arg-type]
        )


def test_standard_length_summary_records_warning_threshold_and_warnings(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(subject, "_load_model", _model)
    configuration = subject.resolve_standard_length_configuration(
        _DEFAULT_CONFIG, enabled=True, model_path=None, approved_enabled=False, warning_threshold=105.0
    )
    from tests.unit.test_length_standard_model import _measurement
    from vntyper.scripts.length_standard_model import predict_standard_length

    assert configuration.model is not None
    measurement = _measurement()
    prediction = predict_standard_length(
        measurement,
        configuration.model,
        warning_threshold=configuration.warning_threshold,
    )
    summary = subject.standard_length_summary(configuration, measurement, prediction)
    assert summary["length_warning_threshold"] == 105.0
    assert summary["length_estimation_warnings"] == ["vntr_length_exceeds_sensitivity_cutoff"]


def test_encode_configuration_rejects_invalid_warning_threshold() -> None:
    config = subject.StandardLengthConfiguration(
        enabled=False,
        model=None,
        source=None,
        sha256="0" * 64,
        warning_threshold=-5.0,
    )
    with pytest.raises(ValueError, match="standard length warning threshold must be a positive number"):
        subject.encode_standard_length_configuration(config)
