"""CLI admits generated callers only through complete approved portable bundles."""

import pytest

from tests.unit.test_pipeline_caller_configuration import configuration
from vntyper import cli

pytestmark = pytest.mark.unit


def test_cli_resolves_approved_bundle_and_context_before_dispatch(tmp_path, monkeypatch):
    path, context = configuration(tmp_path)
    observed = []
    monkeypatch.setattr(cli, "load_config", lambda _path=None: {"cli_defaults": {"log_file": None}})
    monkeypatch.setattr(cli, "setup_logging", lambda **_kwargs: None)
    monkeypatch.setitem(cli.HANDLERS, "pipeline", lambda args, **_kwargs: observed.append(args.run_configuration))
    cli.main(
        ["pipeline", "--bam", "input.bam", "--calibration-bundle", str(path), "--calibration-context", str(context)]
    )
    assert len(observed) == 1
    assert observed[0].caller_calibration is not None
    assert observed[0].decision_profile == observed[0].caller_calibration.bundle.profile


@pytest.mark.parametrize(
    "case", ["missing-context", "missing-bundle", "invalid-bundle", "extra-profile", "log-alias", "log-inside-bundle"]
)
def test_invalid_bundle_or_output_alias_is_refused_before_output_or_log_creation(tmp_path, case):
    path, context = configuration(tmp_path)
    original = context.read_bytes()
    output = tmp_path / "must-not-exist"
    args = ["pipeline", "--bam", "input.bam", "--output-dir", str(output)]
    if case != "missing-bundle":
        args += ["--calibration-bundle", str(path if case != "invalid-bundle" else tmp_path / "absent")]
    if case != "missing-context":
        args += ["--calibration-context", str(context)]
    if case == "extra-profile":
        args += ["--decision-profile", str(path / "decision-profile.json")]
    if case == "log-alias":
        args = ["--log-file", str(context)] + args
    elif case == "log-inside-bundle":
        args = ["--log-file", str(path / "new.log")] + args
    with pytest.raises(SystemExit) as error:
        cli.main(args)
    assert error.value.code == 1
    assert not output.exists()
    assert context.read_bytes() == original
    assert not (path / "new.log").exists()
