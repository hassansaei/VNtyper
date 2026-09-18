"""Simulation command produces inputs only and requires an external manifest commitment."""

import json
from importlib import import_module
from pathlib import Path

import pytest

from tests.unit.test_calibration_sim_protocol import simulation_protocol

pytestmark = pytest.mark.unit


def test_generate_and_verify_commands_preserve_external_manifest_binding(tmp_path: Path, capsys):
    module = import_module("calibration_simulation")
    declaration = tmp_path / "declaration.json"
    declaration.write_text(json.dumps(simulation_protocol()))
    output = tmp_path / "inputs"
    assert module.main(["generate", "--protocol", str(declaration), "--output", str(output)]) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["evidence_status"] == "generated-inputs-only"
    digest = generated["manifest_sha256"]
    assert module.main(["verify-inputs", "--input", str(output), "--manifest-sha256", digest]) == 0
    verified = json.loads(capsys.readouterr().out)
    assert verified["manifest_sha256"] == digest
    assert verified["protocol_sha256"] == generated["protocol_sha256"]
    assert verified["evidence_status"] == "verified-inputs-only"
    assert "passed" not in verified and "promotion_eligible" not in verified
    assert module.main(["verify-inputs", "--input", str(output), "--manifest-sha256", "0" * 64]) == 1
    assert not capsys.readouterr().out


@pytest.mark.parametrize("arguments", [[], ["run"], ["generate"], ["verify-inputs", "--input", "missing"]])
def test_usage_errors_remain_argparse_two(arguments):
    with pytest.raises(SystemExit) as caught:
        import_module("calibration_simulation").main(arguments)
    assert caught.value.code == 2


def test_missing_or_linked_protocol_fails_without_output(tmp_path: Path, capsys):
    module = import_module("calibration_simulation")
    missing = tmp_path / "missing.json"
    output = tmp_path / "output"
    assert module.main(["generate", "--protocol", str(missing), "--output", str(output)]) == 1
    source = tmp_path / "source.json"
    source.write_text(json.dumps(simulation_protocol()))
    missing.symlink_to(source)
    assert module.main(["generate", "--protocol", str(missing), "--output", str(output)]) == 1
    assert not output.exists()
    assert capsys.readouterr().out == ""


def test_existing_output_is_untouched_and_reports_failure(tmp_path: Path, capsys):
    module = import_module("calibration_simulation")
    declaration = tmp_path / "declaration.json"
    declaration.write_text(json.dumps(simulation_protocol()))
    output = tmp_path / "inputs"
    output.mkdir()
    marker = output / "keep"
    marker.write_text("original")
    assert module.main(["generate", "--protocol", str(declaration), "--output", str(output)]) == 1
    assert marker.read_text() == "original"
    assert not capsys.readouterr().out
