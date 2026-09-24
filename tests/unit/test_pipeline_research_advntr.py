"""Research decision profiles run adVNTR with their calibrated legacy policy (#269)."""

from __future__ import annotations

from dataclasses import replace
from importlib import import_module
from pathlib import Path
from unittest.mock import Mock

import pytest

from tests.unit.test_calibration_caller_policy import policy_document
from vntyper.modules.advntr.advntr_command_builder import build_advntr_command
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values
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
    # Without the research argv this exact call raised the pre-fix defect.
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
