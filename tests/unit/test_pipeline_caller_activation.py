"""Approved caller snapshots and output protection precede native activation."""

from importlib import import_module

import pytest

from tests.unit.test_pipeline_caller_configuration import configuration, dual_configuration
from vntyper.scripts.calibration_caller_bundle import load_caller_model_bundle
from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def module():
    return import_module("vntyper.scripts.pipeline_caller_activation")


@pytest.mark.parametrize("mode", ["kestrel", "legacy", "exact"])
def test_snapshot_uses_frozen_bytes_and_exposes_only_closed_portable_inventory(tmp_path, mode):
    source = tmp_path / "source"
    source.mkdir()
    paths = configuration(source) if mode == "kestrel" else dual_configuration(source, mode)
    config = resolve_run_configuration(calibration_bundle=paths[0], calibration_context=paths[1]).caller_calibration
    expected = {p.name: p.read_bytes() for p in paths[0].iterdir()}
    for path in paths[0].iterdir():
        path.unlink()
    output = tmp_path / "run"
    output.mkdir()
    background = module().snapshot_caller_calibration(config, output)
    assert {p.name: p.read_bytes() for p in (output / "caller_calibration").iterdir()} == expected
    assert load_caller_model_bundle(output / "caller_calibration").sha256 == config.bundle.sha256
    assert (background is not None) is (mode == "exact")
    if background is not None:
        assert background.read_bytes() == config.bundle.background_bytes
    assert (
        load_strict_json_object((output / "caller_calibration_context.json").read_bytes())["schema_version"]
        == "caller-runtime-context-v1"
    )
    assert module().snapshot_caller_calibration(config, output) == background
    assert module().caller_calibration_identities(config) == {
        "caller_calibration_bundle_sha256": config.bundle.sha256,
        "caller_calibration_context_sha256": config.context.sha256,
    }


@pytest.mark.parametrize("change", ["context", "bundle", "context-symlink", "bundle-symlink"])
def test_existing_snapshot_requires_same_verified_bundle_and_context(tmp_path, change):
    source = tmp_path / "source"
    source.mkdir()
    paths = configuration(source)
    config = resolve_run_configuration(calibration_bundle=paths[0], calibration_context=paths[1]).caller_calibration
    output = tmp_path / "run"
    output.mkdir()
    module().snapshot_caller_calibration(config, output)
    if change == "context":
        (output / "caller_calibration_context.json").write_bytes(b"changed")
    elif change == "bundle":
        (output / "caller_calibration" / "candidate.json").write_bytes(b"changed")
    else:
        name = "caller_calibration" if change == "bundle-symlink" else "caller_calibration_context.json"
        original = output / name
        moved = output / (name + ".moved")
        original.rename(moved)
        original.symlink_to(moved, target_is_directory=change == "bundle-symlink")
    with pytest.raises(ValueError):
        module().snapshot_caller_calibration(config, output)


@pytest.mark.parametrize("target", ["bundle-child", "bundle-parent", "context", "context-alias"])
def test_outputs_never_overlap_operator_calibration_inputs(tmp_path, target):
    paths = configuration(tmp_path)
    config = resolve_run_configuration(calibration_bundle=paths[0], calibration_context=paths[1]).caller_calibration
    if target == "bundle-child":
        output = paths[0] / "out"
    elif target == "bundle-parent":
        output = tmp_path
    elif target == "context":
        output = paths[1]
    else:
        output = tmp_path / "alias"
        output.symlink_to(paths[1])
    with pytest.raises(ValueError, match="operator"):
        module().validate_caller_output_destination(config, output)
    module().validate_caller_output_destination(config, tmp_path / "separate-run")


def test_pipeline_preflight_checks_actual_options_before_read_or_model_io(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    paths = dual_configuration(source, assembly="hg19")
    run = resolve_run_configuration(calibration_bundle=paths[0], calibration_context=paths[1])
    kwargs = {
        "assembly": "hg19",
        "extra_modules": ["advntr"],
        "threads": 2,
        "additional_commands": None,
        "output": tmp_path / "run",
    }
    assert module().validate_caller_pipeline_request(run, **kwargs) == run.caller_calibration.operator_paths
    for change in (
        {"extra_modules": []},
        {"threads": 1},
        {"additional_commands": "--prune-reverse"},
        {"assembly": "hg38"},
    ):
        with pytest.raises(ValueError):
            module().validate_caller_pipeline_request(run, **(kwargs | change))
    assert not kwargs["output"].exists()
    assert module().validate_caller_pipeline_request(resolve_run_configuration(), **kwargs) == ()
