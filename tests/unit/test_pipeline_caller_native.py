"""Actual selected native model and producer are checked before calibrated calls."""

import subprocess
from dataclasses import replace
from importlib import import_module
from types import MappingProxyType

import pytest

from tests.unit.test_advntr_calibration_policy import capabilities
from tests.unit.test_pipeline_caller_configuration import dual_configuration
from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.pipeline_advntr_run_context import AdvntrRunContext
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def module():
    return import_module("vntyper.scripts.pipeline_caller_native")


def native_case(root, mode="exact"):
    path, context_path = dual_configuration(root, mode)
    config = resolve_run_configuration(calibration_bundle=path, calibration_context=context_path)
    model = root / "model.db"
    model.write_bytes(b"invented model bytes")
    native = AdvntrRunContext(
        str(model),
        str(model),
        MappingProxyType({"sha256": config.caller_calibration.bundle.advntr_policy.model_sha256}),
        (2, 4, 0),
        MappingProxyType({"advntr": "mamba run -n envadvntr advntr"}),
    )
    return config.caller_calibration, native


def runner(argv, **kwargs):
    assert argv == ("mamba", "run", "-n", "envadvntr", "advntr", "capabilities", "--json")
    assert kwargs == {"capture_output": True, "text": True, "check": False}
    document = capabilities()
    document["source_revision"] = "d" * 40
    return subprocess.CompletedProcess(argv, 0, canonical_json_bytes(document).decode(), "")


@pytest.mark.parametrize("mode", ["legacy", "exact"])
def test_native_preflight_binds_observed_build_model_and_exact_policy(tmp_path, mode):
    configuration, native = native_case(tmp_path, mode)
    execution = module().prepare_caller_native_execution(configuration, native, runner=runner)
    assert execution.capabilities.package_version == "2.4.0"
    assert execution.model_sha256 == configuration.bundle.advntr_policy.model_sha256
    assert execution.argv_prefix == ("mamba", "run", "-n", "envadvntr", "advntr")
    assert execution.background_bytes == configuration.bundle.background_bytes
    background = tmp_path / "runtime-background.json" if mode == "exact" else None
    if background is not None:
        background.write_bytes(execution.background_bytes)
    argv = module().caller_native_policy_argv(configuration, execution, background)
    assert argv[:2] == ("-t", "2")
    assert "--min-frameshift-read-support" in argv
    assert ("--exact-frameshift-caller" in argv) is (mode == "exact")
    assert (str(background) in argv) is (mode == "exact")


@pytest.mark.parametrize(
    "change",
    ["model-bytes", "model-claim", "model-symlink", "model-wal", "build", "version", "revision", "failed-process"],
)
def test_native_preflight_cannot_be_satisfied_by_claimed_assets(tmp_path, change):
    from pathlib import Path

    configuration, native = native_case(tmp_path)
    model = Path(native.model_snapshot)
    if change == "model-bytes":
        model.write_bytes(b"changed model")
    elif change == "model-claim":
        native = replace(native, model={"sha256": "0" * 64})
    elif change == "model-symlink":
        target = tmp_path / "other.db"
        model.rename(target)
        model.symlink_to(target)
    elif change == "model-wal":
        Path(str(model) + "-wal").write_bytes(b"unbound query state")

    def changed_runner(argv, **kwargs):
        result = runner(argv, **kwargs)
        document = capabilities()
        document["source_revision"] = "d" * 40
        if change in {"build", "version", "revision"}:
            field, value = {
                "build": ("build_id", "0" * 64),
                "version": ("package_version", "2.4.1"),
                "revision": ("source_revision", "0" * 40),
            }[change]
            document[field] = value
            result.stdout = canonical_json_bytes(document).decode()
        if change == "failed-process":
            result.returncode = 1
        return result

    with pytest.raises((ValueError, RuntimeError)):
        module().prepare_caller_native_execution(configuration, native, runner=changed_runner)


def test_policy_arguments_require_correct_conditional_background(tmp_path):
    configuration, native = native_case(tmp_path)
    execution = module().prepare_caller_native_execution(configuration, native, runner=runner)
    with pytest.raises(ValueError, match="background"):
        module().caller_native_policy_argv(configuration, execution, None)


def test_sqlite_wal_header_is_refused_before_native_probe(tmp_path):
    from pathlib import Path

    configuration, native = native_case(tmp_path)
    Path(native.model_snapshot).write_bytes(b"SQLite format 3\x00" + b"\x00\x00\x02\x02")
    with pytest.raises(ValueError, match="rollback"):
        module().prepare_caller_native_execution(configuration, native, runner=runner)


def test_missing_native_command_and_changed_execution_are_refused(tmp_path):
    configuration, native = native_case(tmp_path)
    with pytest.raises(ValueError, match="executable"):
        module().prepare_caller_native_execution(configuration, replace(native, tools={}), runner=runner)
    execution = module().prepare_caller_native_execution(configuration, native, runner=runner)
    for changed in (replace(execution, model_sha256="0" * 64), replace(execution, background_bytes=b"changed")):
        with pytest.raises(ValueError):
            module().caller_native_policy_argv(configuration, changed, tmp_path / "background.json")
    with pytest.raises(ValueError, match="producer"):
        module().prepare_caller_native_execution(
            replace(configuration, bundle=replace(configuration.bundle, advntr_policy=None)), native, runner=runner
        )
    with pytest.raises(ValueError, match="capture"):
        module().prepare_caller_native_execution(
            replace(configuration, context=replace(configuration.context, advntr_capture_policy=None)),
            native,
            runner=runner,
        )


def test_changed_background_snapshot_cannot_drive_native_call(tmp_path):
    configuration, native = native_case(tmp_path)
    execution = module().prepare_caller_native_execution(configuration, native, runner=runner)
    path = tmp_path / "background.json"
    path.write_bytes(b"changed")
    with pytest.raises(ValueError, match="background"):
        module().caller_native_policy_argv(configuration, execution, path)
