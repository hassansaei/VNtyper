"""Strict answer-line grammar for the adVNTR version probe."""

from __future__ import annotations

import subprocess
from unittest.mock import patch

import pytest

from vntyper.modules.advntr.model_provenance import (
    AdvntrProbeStatus,
    AdvntrVersionProbe,
    detect_advntr_version,
    require_verified_advntr_version,
)

pytestmark = pytest.mark.unit


COMMAND = "mamba run -n envadvntr advntr"
CONFIG = {"tools": {"advntr": COMMAND}}
ARGV = ["mamba", "run", "-n", "envadvntr", "advntr", "--version"]


def _result(*, stdout: str = "", stderr: str = "") -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(ARGV, 0, stdout=stdout, stderr=stderr)


def test_a_diagnostic_semver_before_the_answer_cannot_authorize_the_run() -> None:
    """Only the bare adVNTR answer line is a candidate, never a diagnostic token."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="warning libmamba 2.1.0\n2.0.3\n"),
    ) as runner:
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert require_verified_advntr_version(outcome) == (2, 0, 3)
    runner.assert_called_once_with(ARGV, capture_output=True, text=True, check=False)


def test_usage_text_and_python_semver_cannot_authorize_the_run() -> None:
    """A semver in help or runtime diagnostics is not an adVNTR version answer."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="usage: advntr [options]\nPython 3.12.1 required\n"),
    ) as runner:
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert outcome.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
    runner.assert_called_once_with(ARGV, capture_output=True, text=True, check=False)
    with pytest.raises(RuntimeError, match="unparseable or ambiguous"):
        require_verified_advntr_version(outcome)


def test_conflicting_strict_answers_fail_closed() -> None:
    """Two different adVNTR answers are ambiguous even when each line is valid alone."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stdout="2.0.4\n", stderr="adVNTR 2.0.3: legacy help\n"),
    ) as runner:
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert outcome.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
    runner.assert_called_once_with(ARGV, capture_output=True, text=True, check=False)
    with pytest.raises(RuntimeError, match="unparseable or ambiguous"):
        require_verified_advntr_version(outcome)


def test_benign_stdout_does_not_hide_a_bare_stderr_answer() -> None:
    """Both streams contribute candidates; stdout's mere presence has no priority."""
    probe = AdvntrVersionProbe()
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stdout="launcher note\n", stderr="2.0.4\n"),
    ) as runner:
        first = detect_advntr_version(CONFIG, probe=probe)
        second = detect_advntr_version(CONFIG, probe=probe)

    assert require_verified_advntr_version(first) == (2, 0, 4)
    assert second == first
    runner.assert_called_once_with(ARGV, capture_output=True, text=True, check=False)


@pytest.mark.parametrize(
    ("stdout", "stderr", "expected"),
    [
        ("2.0.4", "", (2, 0, 4)),
        ("", "2.0.3\n", (2, 0, 3)),
        ("", "adVNTR 2.1.0: legacy help", (2, 1, 0)),
        ("", "adVNTR 0.9.0: historical help", (0, 9, 0)),
        ("usage: advntr [options]", "", None),
        ("Python 3.12.1 required", "", None),
        ("warning libmamba 2.1.0", "", None),
        ("advntr 2.1.0", "", None),
        ("", "", None),
    ],
)
def test_actual_probe_preserves_the_strict_answer_line_contract(
    stdout: str, stderr: str, expected: tuple[int, int, int] | None
) -> None:
    """The production parser accepts only bare or explicitly tagged adVNTR answers."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stdout=stdout, stderr=stderr),
    ):
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    if expected is None:
        assert outcome.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
        assert outcome.version is None
    else:
        assert outcome.status is AdvntrProbeStatus.VERIFIED
        assert require_verified_advntr_version(outcome) == expected


@pytest.mark.parametrize(
    "answer",
    [
        "02.00.004\n",
        "adVNTR 02.0.4: legacy help\n",
        "٢.٠.٤\n",
        "adVNTR ٢.0.4: legacy help\n",
    ],
)
def test_non_ascii_or_leading_zero_identifiers_are_rejected_and_never_cached(answer: str) -> None:
    """Relaxing any numeric identifier would authorize a non-SemVer answer line."""
    probe = AdvntrVersionProbe()
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr=answer),
    ) as runner:
        first = detect_advntr_version(CONFIG, probe=probe)
        second = detect_advntr_version(CONFIG, probe=probe)

    assert first.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
    assert first.version is None
    assert second == first
    assert runner.call_count == 2
