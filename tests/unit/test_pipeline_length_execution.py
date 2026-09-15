"""Length measurement executes inside the retained alignment-plan lifetime."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.unit.test_length_estimation import features as _features
from tests.unit.test_length_features import _annotation, _context
from tests.unit.test_pipeline_length import _write_applicable_bundle
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.canonical_json import canonical_sha256
from vntyper.scripts.length_feature_provenance import encode_length_feature_context
from vntyper.scripts.pipeline_length import (
    build_length_pipeline_configuration,
    decode_length_pipeline_context,
    encode_length_pipeline_configuration,
)

pytestmark = pytest.mark.unit


def _plan(*, file_format: str = "bam", reference_path: str | None = None) -> AlignmentPlan:
    return AlignmentPlan(
        input_path="/input/original.bam",
        view_path="/run/bound.bam",
        file_format=file_format,
        index_path="/run/bound.bam.bai",
        reference_path=reference_path,
        reference_source="synthetic",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )


def _measurement_only_configuration():
    annotation = _annotation()
    context = _context(annotation)
    return build_length_pipeline_configuration(
        annotation=annotation,
        pipeline_context=decode_length_pipeline_context(
            {
                "schema_version": "length-pipeline-context-v1",
                "evidence_domain": "synthetic",
                "measurement_context": encode_length_feature_context(context),
            }
        ),
        bundle=None,
    )


def test_measurement_uses_bound_view_observed_reference_and_pinned_tool(monkeypatch: pytest.MonkeyPatch) -> None:
    from vntyper.scripts import pipeline_length_execution

    configuration = _measurement_only_configuration()
    features = _features()
    observed: dict[str, object] = {}

    def fake_depth(*args: object) -> tuple[str, ...]:
        observed["depth_args"] = args
        return ("depth",)

    def fake_features(*args: object) -> object:
        observed["feature_args"] = args
        return features

    monkeypatch.setattr(pipeline_length_execution, "read_length_depth", fake_depth)
    monkeypatch.setattr(pipeline_length_execution, "extract_length_features", fake_features)

    result = pipeline_length_execution.measure_pipeline_length(
        plan=_plan(),
        reference_path=Path("/reference/synthetic.fa"),
        samtools_path=Path("/tools/samtools"),
        configuration=configuration,
    )

    assert observed["depth_args"] == (
        Path("/run/bound.bam"),
        Path("/reference/synthetic.fa"),
        configuration.annotation,
        configuration.measurement_context,
        Path("/tools/samtools"),
    )
    assert observed["feature_args"] == (("depth",), configuration.annotation, configuration.measurement_context)
    assert result.features == features
    assert result.estimate.status == "measured-only"
    assert result.configuration_sha256 == configuration.sha256


def test_estimated_result_uses_only_approved_configuration_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from vntyper.scripts import pipeline_length_execution
    from vntyper.scripts.pipeline_length import build_length_pipeline_configuration, decode_length_pipeline_context

    bundle, sidecar = _write_applicable_bundle(tmp_path / "bundle")
    configuration = build_length_pipeline_configuration(
        annotation=bundle.annotation,
        pipeline_context=decode_length_pipeline_context(sidecar),
        bundle=bundle,
    )
    features = _features()
    # The execution seam delegates scientific availability and prediction to the
    # already-tested estimator; this test observes the exact approved model/domain.
    observed: dict[str, object] = {}
    monkeypatch.setattr(pipeline_length_execution, "read_length_depth", lambda *args: ())
    monkeypatch.setattr(pipeline_length_execution, "extract_length_features", lambda *args: features)

    def fake_state(features_value: object, model_value: object, **kwargs: object):
        observed.update(features=features_value, model=model_value, kwargs=kwargs)
        return __import__("vntyper.scripts.length_estimation", fromlist=["LengthEstimate"]).LengthEstimate(
            "estimated", 110.0, (), bundle.model.sha256, features.sha256
        )

    monkeypatch.setattr(pipeline_length_execution, "length_estimation_state", fake_state)
    result = pipeline_length_execution.measure_pipeline_length(
        plan=_plan(),
        reference_path=Path("/reference/synthetic.fa"),
        samtools_path=Path("/tools/samtools"),
        configuration=configuration,
    )

    assert result.estimate.estimated_total_repeat_count == 110.0
    assert observed == {
        "features": features,
        "model": bundle.model,
        "kwargs": {"measurement_enabled": True, "evidence_domain": "synthetic"},
    }


def test_execution_rejects_disabled_config_missing_reference_and_cram_reference_drift() -> None:
    from dataclasses import replace

    from vntyper.scripts.pipeline_length import resolve_length_pipeline_configuration
    from vntyper.scripts.pipeline_length_execution import measure_pipeline_length

    disabled = resolve_length_pipeline_configuration(
        measurement_enabled=False, model_path=None, annotation_path=None, context_path=None
    )
    with pytest.raises(ValueError, match="enabled"):
        measure_pipeline_length(
            plan=_plan(),
            reference_path=Path("/reference.fa"),
            samtools_path=Path("/samtools"),
            configuration=disabled,
        )
    enabled = _measurement_only_configuration()
    with pytest.raises(ValueError, match="reference"):
        measure_pipeline_length(
            plan=_plan(), reference_path=None, samtools_path=Path("/samtools"), configuration=enabled
        )
    with pytest.raises(ValueError, match="CRAM reference"):
        measure_pipeline_length(
            plan=_plan(file_format="cram", reference_path="/proven/reference.fa"),
            reference_path=Path("/different/reference.fa"),
            samtools_path=Path("/samtools"),
            configuration=enabled,
        )
    with pytest.raises(ValueError, match="digest"):
        measure_pipeline_length(
            plan=_plan(),
            reference_path=Path("/reference.fa"),
            samtools_path=Path("/samtools"),
            configuration=replace(enabled, sha256="f" * 64),
        )


def test_summary_fields_have_stable_disabled_measured_and_estimated_shapes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from vntyper.scripts import pipeline_length_execution
    from vntyper.scripts.pipeline_length import resolve_length_pipeline_configuration

    disabled = resolve_length_pipeline_configuration(
        measurement_enabled=False, model_path=None, annotation_path=None, context_path=None
    )
    disabled_fields = pipeline_length_execution.length_summary_fields(disabled, None)
    assert disabled_fields["length_estimation_status"] == "disabled"
    assert disabled_fields["estimated_total_repeat_count"] is None
    assert disabled_fields["length_features"] is None
    assert disabled_fields["length_configuration_sha256"] == disabled.sha256

    configuration = _measurement_only_configuration()
    features = _features()
    monkeypatch.setattr(pipeline_length_execution, "read_length_depth", lambda *args: ())
    monkeypatch.setattr(pipeline_length_execution, "extract_length_features", lambda *args: features)
    measured = pipeline_length_execution.measure_pipeline_length(
        plan=_plan(),
        reference_path=Path("/reference.fa"),
        samtools_path=Path("/samtools"),
        configuration=configuration,
    )
    measured_fields = pipeline_length_execution.length_summary_fields(configuration, measured)
    assert measured_fields["length_estimation_status"] == "measured-only"
    assert measured_fields["length_features_sha256"] == features.sha256
    assert measured_fields["length_features"]["rows"][0]["A"] == features.a  # type: ignore[index]
    assert canonical_sha256(encode_length_pipeline_configuration(configuration)) == configuration.sha256


def test_summary_rejects_missing_or_cross_configuration_measurement() -> None:
    from dataclasses import replace

    from vntyper.scripts.pipeline_length_execution import LengthPipelineMeasurement, length_summary_fields

    configuration = _measurement_only_configuration()
    features = _features()
    estimate = __import__(
        "vntyper.scripts.length_estimation", fromlist=["length_estimation_state"]
    ).length_estimation_state(features, None, measurement_enabled=True, evidence_domain="synthetic")
    result = LengthPipelineMeasurement(features, estimate, configuration.sha256)
    with pytest.raises(ValueError, match="requires"):
        length_summary_fields(configuration, None)
    with pytest.raises(ValueError, match="configuration"):
        length_summary_fields(configuration, replace(result, configuration_sha256="f" * 64))

    other_features = _features(manifest_key="other-member")
    other_estimate = __import__(
        "vntyper.scripts.length_estimation", fromlist=["length_estimation_state"]
    ).length_estimation_state(other_features, None, measurement_enabled=True, evidence_domain="synthetic")
    with pytest.raises(ValueError, match="measurement context"):
        length_summary_fields(
            configuration,
            LengthPipelineMeasurement(other_features, other_estimate, configuration.sha256),
        )


def test_single_use_runner_selects_cram_or_bwa_reference_and_retains_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from vntyper.scripts import pipeline_length_execution

    configuration = _measurement_only_configuration()
    sentinel = object()
    observed: list[tuple[Path | None, Path]] = []

    def fake_measure(**kwargs: object) -> object:
        observed.append((kwargs["reference_path"], kwargs["samtools_path"]))  # type: ignore[arg-type]
        return sentinel

    monkeypatch.setattr(pipeline_length_execution, "measure_pipeline_length", fake_measure)
    runner = pipeline_length_execution.LengthMeasurementRunner(
        configuration=configuration,
        bwa_reference="relative/reference.fa",
        project_root="/project",
        samtools_path="/tools/samtools",
    )
    runner(_plan())
    assert runner.result is sentinel
    assert observed == [(Path("/project/relative/reference.fa"), Path("/tools/samtools"))]
    with pytest.raises(ValueError, match="exactly once"):
        runner(_plan())

    cram_runner = pipeline_length_execution.LengthMeasurementRunner(
        configuration=configuration,
        bwa_reference="ignored.fa",
        project_root="/project",
        samtools_path="/tools/samtools",
    )
    cram_runner(_plan(file_format="cram", reference_path="/bound/reference.fa"))
    assert observed[-1][0] == Path("/bound/reference.fa")


def test_reference_selection_rejects_missing_bwa_and_relative_project_root() -> None:
    from vntyper.scripts.pipeline_length_execution import length_reference_path

    with pytest.raises(ValueError, match="reference FASTA"):
        length_reference_path(_plan(), None, "/project")
    with pytest.raises(ValueError, match="project root"):
        length_reference_path(_plan(), "/reference.fa", "relative")
