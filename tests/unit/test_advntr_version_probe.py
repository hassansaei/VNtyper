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
TRANSIENT_LOCK_OUTPUT = (
    "error    libmamba Could not set lock (Resource temporarily unavailable)\n"
    "warning  libmamba Cannot lock '/home/test/.cache/mamba/proc'\n"
    "    Waiting for other mamba process to finish\n"
)
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
            side_effect=[_result(stderr=TRANSIENT_LOCK_OUTPUT), _result(stderr="2.0.4\n")],
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
