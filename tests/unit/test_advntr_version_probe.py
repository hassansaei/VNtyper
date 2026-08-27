"""Retry and cache contract for the run-scoped adVNTR version probe."""

from __future__ import annotations

import logging
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import cast
from unittest.mock import call, patch

import pytest

from vntyper.modules.advntr.model_provenance import (
    AdvntrModelError,
    AdvntrProbeStatus,
    AdvntrVersionOutcome,
    AdvntrVersionProbe,
    detect_advntr_version,
    require_compatible_advntr,
    require_verified_advntr_version,
)

pytestmark = pytest.mark.unit


COMMAND = "mamba run -n envadvntr advntr"
CONFIG = {"tools": {"advntr": COMMAND}}
ARGV = ["mamba", "run", "-n", "envadvntr", "advntr", "--version"]
RESOURCE_LOCK_OUTPUT = "error    libmamba Could not set lock (Resource temporarily unavailable)\n"
PROCESS_LOCK_OUTPUT = (
    "warning  libmamba Cannot lock '/home/test/.cache/mamba/proc'\n    Waiting for other mamba process to finish\n"
)
TRANSIENT_LOCK_OUTPUT = RESOURCE_LOCK_OUTPUT + PROCESS_LOCK_OUTPUT
MEASURED_LOCK_WITH_VERSION = f"{TRANSIENT_LOCK_OUTPUT}2.0.4\n{TRANSIENT_LOCK_OUTPUT}"
V2_MODEL = {"schema_version": "v2", "path": "model.db", "window_bp": 3525}


def _result(*, status: int = 0, stdout: str = "", stderr: str = "") -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(ARGV, status, stdout=stdout, stderr=stderr)


def _assert_verified(outcome: AdvntrVersionOutcome, version: tuple[int, int, int]) -> None:
    """Assert the typed successful-probe contract and its exact version."""
    assert outcome.status is AdvntrProbeStatus.VERIFIED
    assert outcome.version == version
    assert require_verified_advntr_version(outcome) == version


@pytest.mark.parametrize(
    ("status", "version", "message", "expected_error"),
    [
        ("verified", (2, 0, 4), "", "status must be an AdvntrProbeStatus"),
        (AdvntrProbeStatus.VERIFIED, None, "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, (99,), "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, (2, 0, 4, 1), "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, [2, 0, 4], "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, (True, 0, 4), "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, (2, -1, 4), "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, (2, "0", 4), "", "version must be exactly three non-negative integers"),
        (AdvntrProbeStatus.VERIFIED, (2, 0, 4), "launch failed", "failure message must be empty"),
        (AdvntrProbeStatus.LAUNCH_FAILURE, (2, 0, 4), "adVNTR version launch failed", "version must be None"),
        (AdvntrProbeStatus.UNPARSEABLE_SUCCESS, (2, 0, 4), "response was ambiguous", "version must be None"),
        (AdvntrProbeStatus.TRANSIENT_EXHAUSTED, (2, 0, 4), "retries exhausted", "version must be None"),
        (AdvntrProbeStatus.LAUNCH_FAILURE, None, "", "message must start with 'adVNTR version launch failed'"),
        (
            AdvntrProbeStatus.UNPARSEABLE_SUCCESS,
            None,
            "",
            "message must start with 'adVNTR version command succeeded'",
        ),
        (
            AdvntrProbeStatus.TRANSIENT_EXHAUSTED,
            None,
            "",
            "message must start with 'adVNTR version detection exhausted'",
        ),
        (
            AdvntrProbeStatus.LAUNCH_FAILURE,
            None,
            "adVNTR version command succeeded but its response was ambiguous.",
            "message must start with 'adVNTR version launch failed'",
        ),
    ],
)
def test_outcome_constructor_rejects_every_inconsistent_public_state(
    status: object,
    version: object,
    message: object,
    expected_error: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The public typed result must make contradictory discriminant states impossible."""
    caplog.set_level(logging.ERROR, logger="vntyper.modules.advntr.model_provenance")

    with pytest.raises(ValueError, match=expected_error) as exc_info:
        AdvntrVersionOutcome(
            cast(AdvntrProbeStatus, status),
            version=cast(tuple[int, int, int] | None, version),
            message=cast(str, message),
        )

    assert any(record.getMessage() == str(exc_info.value) for record in caplog.records)


def test_optimized_interpreter_guard_rejects_a_forged_invalid_verified_outcome() -> None:
    """Removing assertions cannot let a malformed verified payload reach compatibility."""
    script = "\n".join(
        [
            "from vntyper.modules.advntr.model_provenance import (",
            "    AdvntrProbeStatus, AdvntrVersionOutcome, require_verified_advntr_version,",
            ")",
            "outcome = object.__new__(AdvntrVersionOutcome)",
            "object.__setattr__(outcome, 'status', AdvntrProbeStatus.VERIFIED)",
            "object.__setattr__(outcome, 'version', (99,))",
            "object.__setattr__(outcome, 'message', '')",
            "try:",
            "    require_verified_advntr_version(outcome)",
            "except ValueError:",
            "    print('REJECTED_BEFORE_COMPATIBILITY')",
            "else:",
            "    print('COMPATIBILITY_REACHED')",
        ]
    )

    result = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True, check=False)

    assert result.returncode == 0
    assert result.stdout == "REJECTED_BEFORE_COMPATIBILITY\n"
    assert "COMPATIBILITY_REACHED" not in result.stdout


@pytest.mark.parametrize(
    ("result", "expected_status", "expected_message"),
    [
        (
            _result(status=127, stderr="not found\n"),
            AdvntrProbeStatus.LAUNCH_FAILURE,
            "adVNTR version launch failed: command exited with status 127.",
        ),
        (
            _result(stdout="usage: advntr\n", stderr="Python 3.12.1 required\n"),
            AdvntrProbeStatus.UNPARSEABLE_SUCCESS,
            "adVNTR version command succeeded but its response was unparseable or ambiguous.",
        ),
    ],
)
def test_terminal_process_results_have_explicit_distinct_outcomes(
    result: subprocess.CompletedProcess[str],
    expected_status: AdvntrProbeStatus,
    expected_message: str,
) -> None:
    """A permanent exit and a malformed success cannot collapse to the same None."""
    with patch("vntyper.modules.advntr.model_provenance.subprocess.run", return_value=result):
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert outcome == AdvntrVersionOutcome(expected_status, message=expected_message)
    with pytest.raises(RuntimeError, match=expected_message.replace(".", r"\.")):
        require_verified_advntr_version(outcome)


def test_a_transient_mamba_lock_failure_is_retried_then_cached() -> None:
    """Removing the lock-output classification would return unknown on attempt one."""
    probe = AdvntrVersionProbe()
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            side_effect=[_result(stderr=RESOURCE_LOCK_OUTPUT), _result(stderr="2.0.4\n")],
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))

    assert runner.call_args_list == [
        call(ARGV, capture_output=True, text=True, check=False),
        call(ARGV, capture_output=True, text=True, check=False),
    ]
    sleep.assert_called_once()


def test_the_measured_process_lock_failure_is_independently_retried() -> None:
    """The path-specific process-lock diagnostic is a positive retry marker on its own."""
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            side_effect=[_result(stderr=PROCESS_LOCK_OUTPUT), _result(stderr="2.0.4\n")],
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        _assert_verified(detect_advntr_version(CONFIG, probe=AdvntrVersionProbe()), (2, 0, 4))

    assert runner.call_count == 2
    sleep.assert_called_once()


def test_a_valid_version_surrounded_by_mamba_lock_noise_is_accepted_without_retry() -> None:
    """Rejecting the measured success shape would turn lock warnings into exhaustion."""
    probe = AdvntrVersionProbe()
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=_result(stderr=MEASURED_LOCK_WITH_VERSION),
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))

    assert runner.call_count == 1
    sleep.assert_not_called()


def test_transient_mamba_lock_failure_exhaustion_aborts_and_is_not_cached() -> None:
    """Reducing exhaustion to None would misreport launcher contention as incompatibility."""
    probe = AdvntrVersionProbe()
    responses = [_result(stderr=TRANSIENT_LOCK_OUTPUT) for _ in range(3)] + [_result(stderr="2.0.4\n")]
    with (
        patch("vntyper.modules.advntr.model_provenance.subprocess.run", side_effect=responses) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep"),
    ):
        exhausted = detect_advntr_version(CONFIG, probe=probe)
        assert exhausted.status is AdvntrProbeStatus.TRANSIENT_EXHAUSTED
        with pytest.raises(RuntimeError, match="transient mamba launch failures"):
            require_verified_advntr_version(exhausted)
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))

    assert runner.call_count == 4


def test_status_zero_malformed_output_is_not_retried_cached_or_accepted() -> None:
    """Broadly retrying every unknown answer would hide genuine unsupported output."""
    probe = AdvntrVersionProbe()
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=_result(stderr="usage: advntr [options]\n"),
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        first = detect_advntr_version(CONFIG, probe=probe)
        second = detect_advntr_version(CONFIG, probe=probe)

    assert runner.call_count == 2
    sleep.assert_not_called()
    assert first.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
    assert second.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
    with pytest.raises(RuntimeError, match="unparseable or ambiguous"):
        require_verified_advntr_version(first)


def test_a_diagnostic_semver_before_the_answer_cannot_authorize_the_run() -> None:
    """Only the bare adVNTR answer line is a candidate, never a diagnostic token."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="warning libmamba 2.1.0\n2.0.3\n"),
    ) as runner:
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    version = require_verified_advntr_version(outcome)
    assert version == (2, 0, 3)
    assert runner.call_count == 1
    with pytest.raises(AdvntrModelError, match="adVNTR 2.0.3"):
        require_compatible_advntr(V2_MODEL, version)


def test_usage_text_and_python_semver_cannot_authorize_the_run() -> None:
    """A semver in help or runtime diagnostics is not an adVNTR version answer."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="usage: advntr [options]\nPython 3.12.1 required\n"),
    ) as runner:
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert outcome.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS
    assert runner.call_count == 1
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
    assert runner.call_count == 1
    with pytest.raises(RuntimeError, match="unparseable or ambiguous"):
        require_verified_advntr_version(outcome)


def test_benign_stdout_does_not_hide_a_bare_stderr_answer() -> None:
    """Both streams contribute candidates; stdout's mere presence has no priority."""
    probe = AdvntrVersionProbe()
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stdout="launcher note\n", stderr="2.0.4\n"),
    ) as runner:
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))
        _assert_verified(detect_advntr_version(CONFIG, probe=probe), (2, 0, 4))

    assert runner.call_count == 1


def test_an_explicit_legacy_advntr_answer_is_still_parsed() -> None:
    """The historical help-banner answer remains supported when explicitly tagged."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="usage: advntr\noptions:\nadVNTR 1.4.0: legacy help\n"),
    ):
        _assert_verified(detect_advntr_version(CONFIG, probe=AdvntrVersionProbe()), (1, 4, 0))


def test_a_parsed_incompatible_version_is_not_retried_or_accepted() -> None:
    """A syntactically valid old version remains a compatibility failure, not a launch failure."""
    probe = AdvntrVersionProbe()
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=_result(stderr="2.0.3\n"),
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        outcome = detect_advntr_version(CONFIG, probe=probe)

    version = require_verified_advntr_version(outcome)
    assert version == (2, 0, 3)
    assert runner.call_count == 1
    sleep.assert_not_called()
    with pytest.raises(AdvntrModelError, match="adVNTR 2.0.3"):
        require_compatible_advntr(V2_MODEL, version)


def test_concurrent_callers_share_one_successful_probe() -> None:
    """Removing the probe lock would let every simultaneous caller launch mamba."""
    callers = 16
    barrier = threading.Barrier(callers)
    started = threading.Event()
    release = threading.Event()
    probe = AdvntrVersionProbe()

    def run_version(*_args: object, **_kwargs: object) -> subprocess.CompletedProcess[str]:
        started.set()
        assert release.wait(timeout=5)
        return _result(stderr="2.0.4\n")

    def detect() -> AdvntrVersionOutcome:
        barrier.wait(timeout=5)
        return detect_advntr_version(CONFIG, probe=probe)

    with (
        patch("vntyper.modules.advntr.model_provenance.subprocess.run", side_effect=run_version) as runner,
        ThreadPoolExecutor(max_workers=callers) as executor,
    ):
        futures = [executor.submit(detect) for _ in range(callers)]
        assert started.wait(timeout=5)
        release.set()
        expected = AdvntrVersionOutcome(AdvntrProbeStatus.VERIFIED, version=(2, 0, 4))
        assert [future.result(timeout=5) for future in futures] == [expected] * callers

    assert runner.call_count == 1


def test_a_new_run_scoped_probe_starts_with_an_empty_cache() -> None:
    """Making the cache process-global would reuse stale executable state across runs."""
    first_run = AdvntrVersionProbe()
    second_run = AdvntrVersionProbe()
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        side_effect=[_result(stderr="2.0.4\n"), _result(stderr="2.1.0\n")],
    ) as runner:
        _assert_verified(detect_advntr_version(CONFIG, probe=first_run), (2, 0, 4))
        _assert_verified(detect_advntr_version(CONFIG, probe=first_run), (2, 0, 4))
        _assert_verified(detect_advntr_version(CONFIG, probe=second_run), (2, 1, 0))

    assert runner.call_count == 2


def test_lock_text_from_a_non_mamba_command_is_not_classified_as_transient() -> None:
    """The retry classification belongs to the measured mamba launcher failure only."""
    config = {"tools": {"advntr": "/opt/advntr"}}
    probe = AdvntrVersionProbe()
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=_result(stderr=TRANSIENT_LOCK_OUTPUT),
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        outcome = detect_advntr_version(config, probe=probe)

    assert runner.call_count == 1
    sleep.assert_not_called()
    assert outcome.status is AdvntrProbeStatus.UNPARSEABLE_SUCCESS


def test_a_package_cache_permission_error_is_not_process_lock_contention() -> None:
    """A generic Cannot-lock error must return unknown without retrying or sleeping."""
    output = "error libmamba Cannot lock '/opt/conda/pkgs/cache' (Permission denied)\n"
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=_result(status=1, stderr=output),
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert runner.call_count == 1
    sleep.assert_not_called()
    assert outcome.status is AdvntrProbeStatus.LAUNCH_FAILURE


def test_nonzero_clean_output_is_not_accepted_as_a_version() -> None:
    """A failed process cannot become compatible merely by printing a version token."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(status=2, stderr="2.0.4\n"),
    ):
        outcome = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())
    assert outcome.status is AdvntrProbeStatus.LAUNCH_FAILURE


@pytest.mark.parametrize(
    ("failure", "message_fragment"),
    [
        (FileNotFoundError("missing"), "command not found"),
        (PermissionError("denied"), "permission denied"),
        (OSError("exec failed"), "version launch failed for"),
    ],
)
def test_permanent_launch_exceptions_return_unknown_without_being_cached(
    failure: OSError,
    message_fragment: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Removing a permanent-launch handler would bypass the compatibility refusal."""
    probe = AdvntrVersionProbe()
    caplog.set_level(logging.ERROR, logger="vntyper.modules.advntr.model_provenance")
    with patch("vntyper.modules.advntr.model_provenance.subprocess.run", side_effect=failure) as runner:
        first = detect_advntr_version(CONFIG, probe=probe)
        second = detect_advntr_version(CONFIG, probe=probe)
        with pytest.raises(RuntimeError):
            require_verified_advntr_version(first)

    assert runner.call_count == 2
    assert first.status is AdvntrProbeStatus.LAUNCH_FAILURE
    assert second.status is AdvntrProbeStatus.LAUNCH_FAILURE
    assert message_fragment in first.message
    assert second.message == first.message
    assert any(record.getMessage() == first.message for record in caplog.records)


def test_concurrent_failures_are_never_cached() -> None:
    """The shared lock serializes callers but only a verified tuple suppresses launches."""
    callers = 4
    barrier = threading.Barrier(callers)
    probe = AdvntrVersionProbe()

    def detect() -> AdvntrVersionOutcome:
        barrier.wait(timeout=5)
        return detect_advntr_version(CONFIG, probe=probe)

    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=_result(stderr="usage: advntr [options]\n"),
        ) as runner,
        ThreadPoolExecutor(max_workers=callers) as executor,
    ):
        outcomes = [future.result(timeout=5) for future in [executor.submit(detect) for _ in range(callers)]]

    assert [outcome.status for outcome in outcomes] == [AdvntrProbeStatus.UNPARSEABLE_SUCCESS] * callers
    assert runner.call_count == callers


def test_missing_command_returns_unknown_without_launching() -> None:
    """A partial configuration is not allowed to invent a bare executable name."""
    with patch("vntyper.modules.advntr.model_provenance.subprocess.run") as runner:
        outcome = detect_advntr_version({"tools": {}}, probe=AdvntrVersionProbe())
    runner.assert_not_called()
    assert outcome.status is AdvntrProbeStatus.LAUNCH_FAILURE
    assert outcome.message == "adVNTR version launch failed: no command is configured."
