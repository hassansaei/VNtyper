"""Production native stage consumes approved settings and rechecks live assets."""

from importlib import import_module
from unittest.mock import Mock

import pytest

from tests.unit.test_pipeline_caller_native import native_case, runner
from vntyper.scripts.pipeline_caller_native import prepare_caller_native_execution
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def module():
    return import_module("vntyper.scripts.pipeline_advntr_execution")


def test_legacy_stage_preserves_explicit_runtime_and_callback_contract(tmp_path):
    from vntyper.scripts.pipeline_advntr_run_context import AdvntrRunContext

    config = resolve_run_configuration()
    native = AdvntrRunContext("source.db", "snapshot.db", {}, (2, 4, 0), {"advntr": "pinned-tool"})
    invoke = Mock(return_value=0)
    assert (
        module().execute_advntr_genotype(
            configuration=config,
            native_context=native,
            config={"tools": {"advntr": "other"}},
            alignment=tmp_path / "input.bam",
            output=tmp_path,
            cwd="cwd",
            threads=2,
            additional_commands="--prune-reverse",
            background=None,
            invoke=invoke,
        )
        == 0
    )
    assert invoke.call_args.args == ("snapshot.db", tmp_path / "input.bam", tmp_path, "output")
    assert invoke.call_args.kwargs["config"]["tools"] == {"advntr": "pinned-tool"}
    assert invoke.call_args.kwargs["runtime_component"]["settings"]["additional_commands"] == "--prune-reverse"
    assert "calibrated_policy_arguments" not in invoke.call_args.kwargs


@pytest.mark.parametrize("mode", ["legacy", "exact"])
def test_calibrated_stage_rechecks_native_assets_and_passes_all_selected_arguments(tmp_path, monkeypatch, mode):
    from dataclasses import replace

    calibration, native = native_case(tmp_path, mode)
    run = resolve_run_configuration()
    run = replace(
        run,
        decision_profile=calibration.bundle.profile,
        advntr=calibration.bundle.profile.components["advntr"],
        caller_calibration=calibration,
    )
    observed = prepare_caller_native_execution(calibration, native, runner=runner)
    verify = Mock(return_value=observed)
    monkeypatch.setattr(module(), "prepare_caller_native_execution", verify)
    invoke = Mock(return_value=0)
    background = None
    if mode == "exact":
        background = tmp_path / "background.json"
        background.write_bytes(observed.background_bytes)
    assert (
        module().execute_advntr_genotype(
            configuration=run,
            native_context=native,
            config={},
            alignment=tmp_path / "input.bam",
            output=tmp_path,
            cwd="cwd",
            threads=2,
            additional_commands=None,
            background=background,
            invoke=invoke,
        )
        == 0
    )
    verify.assert_called_once_with(calibration, native)
    argv = invoke.call_args.kwargs["calibrated_policy_arguments"]
    assert argv[:2] == ("-t", "2")
    assert ("--exact-frameshift-caller" in argv) is (mode == "exact")
    verify.side_effect = ValueError("changed native build")
    invoke.reset_mock()
    with pytest.raises(ValueError, match="changed native build"):
        module().execute_advntr_genotype(
            configuration=run,
            native_context=native,
            config={},
            alignment=tmp_path / "input.bam",
            output=tmp_path,
            cwd="cwd",
            threads=2,
            additional_commands=None,
            background=background,
            invoke=invoke,
        )
    invoke.assert_not_called()
