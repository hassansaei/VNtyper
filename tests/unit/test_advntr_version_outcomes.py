"""Public invariant contract for typed adVNTR version-probe outcomes."""

from __future__ import annotations

import logging
import subprocess
import sys
from typing import cast

import pytest

from vntyper.modules.advntr.model_provenance import (
    AdvntrProbeStatus,
    AdvntrVersionOutcome,
)

pytestmark = pytest.mark.unit


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
        (
            AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE,
            None,
            "adVNTR version launch failed: command exited with status 2.",
            "version must be exactly three non-negative integers",
        ),
        (
            AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE,
            (2, 0, 3),
            "",
            "message must start with 'adVNTR version launch failed'",
        ),
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
