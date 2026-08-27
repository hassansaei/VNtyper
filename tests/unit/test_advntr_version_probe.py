"""Retry and cache contract for the run-scoped adVNTR version probe."""

from __future__ import annotations

import logging
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import call, patch

import pytest

from vntyper.modules.advntr.model_provenance import (
    AdvntrModelError,
    AdvntrVersionProbe,
    detect_advntr_version,
    require_compatible_advntr,
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
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)

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
        assert detect_advntr_version(CONFIG, probe=AdvntrVersionProbe()) == (2, 0, 4)

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
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)

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
        with pytest.raises(RuntimeError, match="transient mamba launch failure"):
            detect_advntr_version(CONFIG, probe=probe)
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)

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
        assert detect_advntr_version(CONFIG, probe=probe) is None
        assert detect_advntr_version(CONFIG, probe=probe) is None

    assert runner.call_count == 2
    sleep.assert_not_called()
    with pytest.raises(AdvntrModelError, match="adVNTR unknown"):
        require_compatible_advntr(V2_MODEL, None)


def test_a_diagnostic_semver_before_the_answer_cannot_authorize_the_run() -> None:
    """Only the bare adVNTR answer line is a candidate, never a diagnostic token."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="warning libmamba 2.1.0\n2.0.3\n"),
    ) as runner:
        version = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

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
        version = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert version is None
    assert runner.call_count == 1
    with pytest.raises(AdvntrModelError, match="adVNTR unknown"):
        require_compatible_advntr(V2_MODEL, version)


def test_conflicting_strict_answers_fail_closed() -> None:
    """Two different adVNTR answers are ambiguous even when each line is valid alone."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stdout="2.0.4\n", stderr="adVNTR 2.0.3: legacy help\n"),
    ) as runner:
        version = detect_advntr_version(CONFIG, probe=AdvntrVersionProbe())

    assert version is None
    assert runner.call_count == 1
    with pytest.raises(AdvntrModelError, match="adVNTR unknown"):
        require_compatible_advntr(V2_MODEL, version)


def test_benign_stdout_does_not_hide_a_bare_stderr_answer() -> None:
    """Both streams contribute candidates; stdout's mere presence has no priority."""
    probe = AdvntrVersionProbe()
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stdout="launcher note\n", stderr="2.0.4\n"),
    ) as runner:
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)
        assert detect_advntr_version(CONFIG, probe=probe) == (2, 0, 4)

    assert runner.call_count == 1


def test_an_explicit_legacy_advntr_answer_is_still_parsed() -> None:
    """The historical help-banner answer remains supported when explicitly tagged."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(stderr="usage: advntr\noptions:\nadVNTR 1.4.0: legacy help\n"),
    ):
        assert detect_advntr_version(CONFIG, probe=AdvntrVersionProbe()) == (1, 4, 0)


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
        version = detect_advntr_version(CONFIG, probe=probe)

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

    def detect() -> tuple[int, int, int] | None:
        barrier.wait(timeout=5)
        return detect_advntr_version(CONFIG, probe=probe)

    with (
        patch("vntyper.modules.advntr.model_provenance.subprocess.run", side_effect=run_version) as runner,
        ThreadPoolExecutor(max_workers=callers) as executor,
    ):
        futures = [executor.submit(detect) for _ in range(callers)]
        assert started.wait(timeout=5)
        release.set()
        assert [future.result(timeout=5) for future in futures] == [(2, 0, 4)] * callers

    assert runner.call_count == 1


def test_a_new_run_scoped_probe_starts_with_an_empty_cache() -> None:
    """Making the cache process-global would reuse stale executable state across runs."""
    first_run = AdvntrVersionProbe()
    second_run = AdvntrVersionProbe()
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        side_effect=[_result(stderr="2.0.4\n"), _result(stderr="2.1.0\n")],
    ) as runner:
        assert detect_advntr_version(CONFIG, probe=first_run) == (2, 0, 4)
        assert detect_advntr_version(CONFIG, probe=first_run) == (2, 0, 4)
        assert detect_advntr_version(CONFIG, probe=second_run) == (2, 1, 0)

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
        assert detect_advntr_version(config, probe=probe) is None

    assert runner.call_count == 1
    sleep.assert_not_called()


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
        assert detect_advntr_version(CONFIG, probe=AdvntrVersionProbe()) is None

    assert runner.call_count == 1
    sleep.assert_not_called()


def test_nonzero_clean_output_is_not_accepted_as_a_version() -> None:
    """A failed process cannot become compatible merely by printing a version token."""
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=_result(status=2, stderr="2.0.4\n"),
    ):
        assert detect_advntr_version(CONFIG, probe=AdvntrVersionProbe()) is None


@pytest.mark.parametrize(
    ("failure", "logged"),
    [
        (FileNotFoundError("missing"), "Command not found"),
        (PermissionError("denied"), "Permission denied"),
        (OSError("exec failed"), "Failed to get adVNTR version"),
    ],
)
def test_permanent_launch_exceptions_return_unknown_without_being_cached(
    failure: OSError,
    logged: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Removing a permanent-launch handler would bypass the compatibility refusal."""
    probe = AdvntrVersionProbe()
    caplog.set_level(logging.ERROR, logger="vntyper.modules.advntr.model_provenance")
    with patch("vntyper.modules.advntr.model_provenance.subprocess.run", side_effect=failure) as runner:
        assert detect_advntr_version(CONFIG, probe=probe) is None
        assert detect_advntr_version(CONFIG, probe=probe) is None

    assert runner.call_count == 2
    assert logged in caplog.text


def test_missing_command_returns_unknown_without_launching() -> None:
    """A partial configuration is not allowed to invent a bare executable name."""
    with patch("vntyper.modules.advntr.model_provenance.subprocess.run") as runner:
        assert detect_advntr_version({"tools": {}}, probe=AdvntrVersionProbe()) is None
    runner.assert_not_called()
