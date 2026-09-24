"""Research decision profiles run adVNTR with their calibrated legacy policy (#269)."""

from __future__ import annotations

from dataclasses import replace
from importlib import import_module
from pathlib import Path
from unittest.mock import Mock

import pytest

from tests.unit.test_calibration_caller_policy import policy_document
from vntyper.modules.advntr.advntr_command_builder import build_advntr_command
from vntyper.scripts.calibration_caller_policy import KESTREL_CALLER_POLICY_POINTERS, decode_caller_policy_values
from vntyper.scripts.calibration_caller_profile import build_caller_generated_profile
from vntyper.scripts.pipeline_advntr_run_context import AdvntrRunContext
from vntyper.scripts.run_configuration import RunConfiguration, resolve_run_configuration

pytestmark = pytest.mark.unit

_PREFIX = "/components/advntr/calibrated_calling/"
_OLD_FAILURE = "calibrated native settings require their approved explicit policy arguments"


def module():
    return import_module("vntyper.scripts.pipeline_research_advntr")


def _research_configuration(tmp_path: Path, *, mode: str = "legacy", cutoff: float = 0.004) -> RunConfiguration:
    document = policy_document(include_advntr=True)
    values = document["values"]
    assert isinstance(values, dict)
    values[f"{_PREFIX}mode"] = mode
    values[f"{_PREFIX}cutoff"] = cutoff
    values[f"{_PREFIX}minimum_read_match_ratio"] = 0.6
    profile = build_caller_generated_profile(
        decode_caller_policy_values(document),
        dataset_manifest_hash="a" * 64,
        partition_manifest_hash="b" * 64,
        seed=269,
        generator_version="unit-test",
    )
    path = tmp_path / "research_profile.json"
    path.write_bytes(profile.canonical_bytes)
    return resolve_run_configuration(None, research_profile=path)


def test_legacy_research_profile_renders_the_calibrated_policy_argv(tmp_path: Path) -> None:
    argv = module().research_advntr_policy_argv(_research_configuration(tmp_path), 4)

    assert argv is not None
    assert argv[:6] == ("-t", "4", "--frameshift-pvalue-cutoff", "0.004", "--min-frameshift-read-support", "3")
    assert argv[6:8] == ("--min-read-match-ratio", "0.6")
    assert "--prune-reverse" not in argv
    assert "--filter-adapter-readthrough" not in argv
    assert "--exact-frameshift-caller" not in argv


def test_exact_research_profile_is_refused(tmp_path: Path) -> None:
    configuration = _research_configuration(tmp_path, mode="exact")

    with pytest.raises(ValueError, match="exact mode requires an approved calibration bundle"):
        module().research_advntr_policy_argv(configuration, 4)


def test_packaged_profile_has_no_research_policy() -> None:
    assert module().research_advntr_policy_argv(resolve_run_configuration(), 4) is None


def test_bundle_configuration_is_left_to_its_approved_path(tmp_path: Path) -> None:
    configuration = replace(_research_configuration(tmp_path), caller_calibration=Mock())

    assert module().research_advntr_policy_argv(configuration, 4) is None


def test_research_capture_parameters_are_immutable() -> None:
    parameters = module().RESEARCH_CAPTURE_PARAMETERS

    assert "threads" not in parameters
    with pytest.raises(TypeError):
        parameters["caller_mode"] = "exact"


def _execute(configuration: RunConfiguration, tmp_path: Path, invoke: Mock, threads: int) -> None:
    execution = import_module("vntyper.scripts.pipeline_advntr_execution")
    native = AdvntrRunContext("source.db", "snapshot.db", {}, (2, 4, 0), {"advntr": "pinned-tool"})
    execution.execute_advntr_genotype(
        configuration=configuration,
        native_context=native,
        config={},
        alignment=tmp_path / "input.bam",
        output=tmp_path,
        cwd="cwd",
        threads=threads,
        additional_commands=None,
        background=None,
        invoke=invoke,
    )


def test_research_profile_executes_adVNTR_with_its_derived_cutoff(tmp_path: Path) -> None:
    configuration = _research_configuration(tmp_path, cutoff=0.0025)
    invoke = Mock(return_value=0)

    _execute(configuration, tmp_path, invoke, threads=3)

    arguments = invoke.call_args.kwargs["calibrated_policy_arguments"]
    assert arguments[:4] == ("-t", "3", "--frameshift-pvalue-cutoff", "0.0025")
    # Pins the builder guard: a calibrated component without explicit arguments is still refused.
    with pytest.raises(ValueError, match=_OLD_FAILURE):
        build_advntr_command(
            "adVNTR",
            vid=25561,
            alignment="input.bam",
            result="out.vcf",
            model="snapshot.db",
            working_directory="out",
            threads=3,
            additional_commands="",
            calibrated="calibrated_calling" in configuration.advntr,
            calibrated_policy_arguments=None,
        )
    command = build_advntr_command(
        "adVNTR",
        vid=25561,
        alignment="input.bam",
        result="out.vcf",
        model="snapshot.db",
        working_directory="out",
        threads=3,
        additional_commands="",
        calibrated=True,
        calibrated_policy_arguments=arguments,
    )
    assert "--frameshift-pvalue-cutoff 0.0025" in command
    assert "--min-read-match-ratio 0.6" in command


def test_research_threads_follow_the_runtime_resolution(tmp_path: Path) -> None:
    configuration = _research_configuration(tmp_path)
    runtime = {**configuration.advntr_runtime, "settings": {"additional_commands": "", "threads": 7}}
    invoke = Mock(return_value=0)

    _execute(replace(configuration, advntr_runtime=runtime), tmp_path, invoke, threads=2)

    assert invoke.call_args.kwargs["calibrated_policy_arguments"][:2] == ("-t", "7")


def test_packaged_profile_execution_passes_no_calibrated_arguments(tmp_path: Path) -> None:
    invoke = Mock(return_value=0)

    _execute(resolve_run_configuration(), tmp_path, invoke, threads=2)

    assert "calibrated_policy_arguments" not in invoke.call_args.kwargs


def test_projecting_a_pointer_the_components_lack_names_it() -> None:
    configuration = resolve_run_configuration()
    missing = "/components/advntr/calibrated_calling/cutoff"

    with pytest.raises(ValueError, match=f"research decision profile does not expose {missing} at runtime"):
        module().project_caller_policy(
            {"advntr": configuration.advntr, "kestrel": configuration.kestrel},
            (missing, *KESTREL_CALLER_POLICY_POINTERS),
            ("advntr", "kestrel"),
        )


def test_projecting_the_packaged_kestrel_pointers_rebuilds_a_kestrel_policy() -> None:
    configuration = resolve_run_configuration()

    policy = module().project_caller_policy(
        {"kestrel": configuration.kestrel}, KESTREL_CALLER_POLICY_POINTERS, ("kestrel",)
    )

    assert policy.required_callers == ("kestrel",)
    assert tuple(policy.values) == KESTREL_CALLER_POLICY_POINTERS


def test_an_optimize_exported_advntr_profile_runs_adVNTR_at_its_selected_cutoff(tmp_path: Path) -> None:
    from tests.unit.test_calibration_cutoff_optimize import _run_advntr

    successful, document, output = _run_advntr(tmp_path, caller="advntr", min_specificity=1.0)
    assert successful is True
    selected = document["selection"]["value"]
    configuration = resolve_run_configuration(None, research_profile=output / "research-decision-profile.json")
    invoke = Mock(return_value=0)

    _execute(configuration, tmp_path, invoke, threads=2)

    command = build_advntr_command(
        "adVNTR",
        vid=25561,
        alignment="input.bam",
        result="out.vcf",
        model="snapshot.db",
        working_directory="out",
        threads=2,
        additional_commands="",
        calibrated="calibrated_calling" in configuration.advntr,
        calibrated_policy_arguments=invoke.call_args.kwargs["calibrated_policy_arguments"],
    )
    assert f"--frameshift-pvalue-cutoff {selected}" in command
