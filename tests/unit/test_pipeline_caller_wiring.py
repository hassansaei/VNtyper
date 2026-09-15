"""Actual pipeline ordering and calibrated native callback propagation with invented assets."""

import json
from copy import deepcopy

import pytest

from tests.support.pipeline_harness import (
    MINIMAL_CONFIG,
    _minimal_advntr_model,
    advntr_stub,
    run_pipeline_under_harness,
)
from tests.unit.test_advntr_calibration_policy import capabilities
from tests.unit.test_pipeline_caller_configuration import dual_configuration
from vntyper.modules.advntr.advntr_calibration_policy import decode_advntr_capabilities, require_advntr_capabilities
from vntyper.scripts import pipeline_caller_native
from vntyper.scripts.calibration_caller_bundle import load_caller_model_bundle
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def setup(root):
    model = root / "model.db"
    _minimal_advntr_model(model)
    source = root / "source"
    source.mkdir()
    paths = dual_configuration(source, assembly="hg19", model_bytes=model.read_bytes())
    run = resolve_run_configuration(calibration_bundle=paths[0], calibration_context=paths[1])
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = advntr_stub("2.4.0", root)
    config["reference_data"]["advntr_reference_vntr_hg19"] = str(model)
    return run, config


def test_pipeline_checks_native_identity_then_snapshots_and_executes_frozen_policy(tmp_path, monkeypatch):
    run, config = setup(tmp_path)
    observed = []

    def probe(prefix, pin, **_kwargs):
        observed.append((prefix, pin))
        document = capabilities()
        document["source_revision"] = "d" * 40
        return require_advntr_capabilities(decode_advntr_capabilities(document), pin)

    monkeypatch.setattr(pipeline_caller_native, "probe_advntr_capabilities", probe)
    output = tmp_path / "run"
    harness = run_pipeline_under_harness(
        output, config=config, extra_modules=["advntr"], threads=2, run_configuration=run
    )
    assert len(observed) == 2  # before read preparation, again immediately before native execution
    assert load_caller_model_bundle(output / "caller_calibration").sha256 == run.caller_calibration.bundle.sha256
    argv = harness.kwargs("run_advntr")["calibrated_policy_arguments"]
    assert argv[:2] == ("-t", "2")
    assert "--exact-frameshift-caller" in argv
    assert argv[-1] == str(output / "caller_calibration" / "background.json")
    assert (
        json.loads((output / "pipeline_summary.json").read_text())["analysis_settings"][
            "caller_calibration_bundle_sha256"
        ]
        == run.caller_calibration.bundle.sha256
    )
    assert (
        json.loads((output / "pipeline_summary.json").read_text())["analysis_settings"][
            "caller_calibration_context_sha256"
        ]
        == run.caller_calibration.context.sha256
    )


def test_wrong_native_build_fails_before_alignment_and_kestrel(tmp_path, monkeypatch):
    run, config = setup(tmp_path)

    def refuse(*_args, **_kwargs):
        raise ValueError("wrong installed native build")

    monkeypatch.setattr(pipeline_caller_native, "probe_advntr_capabilities", refuse)
    harness = run_pipeline_under_harness(
        tmp_path / "run", config=config, extra_modules=["advntr"], threads=2, run_configuration=run, expect_failure=True
    )
    assert harness.error is not None
    harness.stages["run_kestrel"].assert_not_called()
    harness.stages["process_bam_to_fastq"].assert_not_called()
    harness.stages["run_advntr"].assert_not_called()


def test_missing_explicit_native_stage_is_refused_before_input_fingerprinting(tmp_path, monkeypatch):
    from vntyper.scripts import pipeline

    run, config = setup(tmp_path)
    fingerprint = pytest.MonkeyPatch()
    calls = []
    fingerprint.setattr(pipeline, "build_canonical_inputs_and_fingerprints", lambda *args: calls.append(args))
    try:
        harness = run_pipeline_under_harness(
            tmp_path / "run", config=config, extra_modules=[], threads=2, run_configuration=run, expect_failure=True
        )
    finally:
        fingerprint.undo()
    assert harness.error is not None
    assert calls == []
