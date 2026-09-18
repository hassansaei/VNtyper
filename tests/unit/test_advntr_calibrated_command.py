"""Approved native arguments are consumed exactly once without legacy overrides."""

import shlex

import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def invoke(tmp_path, monkeypatch, *, calibrated=True, arguments=None, additional="", tool="advntr"):
    db = tmp_path / "model.db"
    bam = tmp_path / "input.bam"
    db.touch()
    bam.touch()
    component = dict(resolve_run_configuration().advntr)
    if calibrated:
        component["calibrated_calling"] = {"mode": "legacy", "cutoff": 0.002, "minimum_read_support": 4}
    calls = []
    monkeypatch.setattr(advntr, "run_command", lambda *args, **kwargs: calls.append((args, kwargs)) or True)
    result = advntr.run_advntr(
        str(db),
        str(bam),
        str(tmp_path),
        "output",
        {"tools": {"advntr": tool}},
        resolved_component=component,
        runtime_component={"settings": {"threads": 2, "additional_commands": additional}},
        calibrated_policy_arguments=arguments,
    )
    return result, calls


def test_calibrated_command_consumes_full_explicit_policy_without_duplicate_threads(tmp_path, monkeypatch):
    arguments = (
        "-t",
        "2",
        "--frameshift-pvalue-cutoff",
        "0.002",
        "--min-frameshift-read-support",
        "4",
        "--min-read-match-ratio",
        "0.6",
    )
    result, calls = invoke(tmp_path, monkeypatch, arguments=arguments, tool="mamba run -n envadvntr advntr")
    assert result == 0
    words = shlex.split(calls[0][0][0])
    assert tuple(words[-len(arguments) :]) == arguments
    assert words.count("-t") == 1
    assert calls[0][1]["critical"] is True


@pytest.mark.parametrize("change", ["missing", "unapproved", "extra", "wrong-threads", "non-tuple"])
def test_calibrated_native_command_refuses_unbound_or_conflicting_arguments(tmp_path, monkeypatch, change):
    arguments = ("-t", "2", "--frameshift-pvalue-cutoff", "0.002")
    kwargs = {"arguments": arguments}
    if change == "missing":
        kwargs["arguments"] = None
    elif change == "unapproved":
        kwargs["calibrated"] = False
    elif change == "extra":
        kwargs["additional"] = "--prune-reverse"
    elif change == "wrong-threads":
        kwargs["arguments"] = ("-t", "3")
    else:
        kwargs["arguments"] = list(arguments)
    with pytest.raises(ValueError):
        invoke(tmp_path, monkeypatch, **kwargs)


def test_calibrated_prefix_and_paths_are_shell_quoted_as_observed_tokens(tmp_path, monkeypatch):
    arguments = (
        "-t",
        "2",
        "--frameshift-pvalue-cutoff",
        "0.002",
        "--frameshift-background",
        "background with spaces.json",
    )
    _, calls = invoke(tmp_path, monkeypatch, arguments=arguments, tool="'tool with spaces' --literal '$VALUE'")
    words = shlex.split(calls[0][0][0])
    assert words[:3] == ["tool with spaces", "--literal", "$VALUE"]
    assert words[-1] == "background with spaces.json"
