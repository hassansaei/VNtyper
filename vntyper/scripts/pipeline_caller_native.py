"""Observed native assets and explicit argv for approved caller calibration."""

from __future__ import annotations

import hashlib
import logging
import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path

from vntyper.modules.advntr.advntr_calibration_policy import (
    AdvntrCapabilities,
    AdvntrToolPin,
    ProcessRunner,
    probe_advntr_capabilities,
    require_advntr_capabilities,
)
from vntyper.modules.advntr.advntr_capture import calibrated_policy_argv
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.pipeline_advntr_run_context import AdvntrRunContext
from vntyper.scripts.pipeline_caller_configuration import CallerPipelineConfiguration

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CallerNativeExecution:
    """Observed producer and model bytes, with immutable portable background bytes."""

    argv_prefix: tuple[str, ...]
    capabilities: AdvntrCapabilities
    model_sha256: str
    background_bytes: bytes | None


def _pin(configuration: CallerPipelineConfiguration) -> AdvntrToolPin:
    native = configuration.bundle.advntr_policy
    versions = configuration.bundle.candidate.producer.tool_versions
    if native is None or "advntr" not in versions or "advntr_build_id" not in versions:
        raise ValueError("calibrated native execution requires the complete frozen producer identity")
    return AdvntrToolPin(versions["advntr"], versions["advntr_build_id"], native.advntr_revision)


def _model_bytes(path: Path) -> bytes:
    if any(os.path.lexists(str(path) + suffix) for suffix in ("-wal", "-shm", "-journal")):
        raise ValueError("calibrated native model snapshot has unbound SQLite sidecars")
    raw = read_regular_path(path)
    if raw.startswith(b"SQLite format 3\x00") and (len(raw) < 20 or raw[18:20] != b"\x01\x01"):
        raise ValueError("calibrated native model snapshot requires SQLite rollback-journal format")
    return raw


def prepare_caller_native_execution(
    configuration: CallerPipelineConfiguration,
    native_context: AdvntrRunContext,
    *,
    runner: ProcessRunner = subprocess.run,
) -> CallerNativeExecution:
    """Verify actual selected model bytes and installed native build before read work.

    Args:
        configuration: Approved bundle and resolved capture context.
        native_context: Pipeline-owned model snapshot and actual executable prefix.
        runner: Injectable shell-free process boundary used to probe capabilities.

    Returns:
        Observed immutable native execution identity.

    Raises:
        ValueError: If model bytes, model claim, producer, or executable differ.
        RuntimeError: If the installed capability command fails.
    """
    pin = _pin(configuration)
    native = configuration.bundle.advntr_policy
    assert native is not None
    model_hash = hashlib.sha256(_model_bytes(Path(native_context.model_snapshot))).hexdigest()
    if model_hash != native.model_sha256 or native_context.model.get("sha256") != model_hash:
        raise ValueError("calibrated native model bytes differ from the frozen model identity")
    command = native_context.tools.get("advntr")
    if not isinstance(command, str) or not command.strip():
        raise ValueError("calibrated native execution requires an explicit executable prefix")
    prefix = tuple(shlex.split(command))
    capabilities = probe_advntr_capabilities(prefix, pin, runner=runner)
    # Reject unsupported CLI construction during preflight, before any read work.
    policy = configuration.context.advntr_capture_policy
    if policy is None:
        raise ValueError("calibrated native execution requires its full capture policy")
    calibrated_policy_argv(
        policy, configuration.bundle.caller_policy, Path("background.json") if native.mode == "exact" else None
    )
    return CallerNativeExecution(prefix, capabilities, model_hash, configuration.bundle.background_bytes)


def caller_native_policy_argv(
    configuration: CallerPipelineConfiguration,
    execution: CallerNativeExecution,
    background_path: Path | None,
) -> tuple[str, ...]:
    """Render exactly the selected full native policy using the shared capture adapter.

    Args:
        configuration: Approved bundle and complete capture policy.
        execution: Previously observed native execution identity.
        background_path: Run-owned exact background snapshot, or None for legacy mode.

    Returns:
        Explicit native policy arguments including the frozen thread count.

    Raises:
        ValueError: If identity, policy, or conditional background presence differs.
    """
    require_advntr_capabilities(execution.capabilities, _pin(configuration))
    native = configuration.bundle.advntr_policy
    policy = configuration.context.advntr_capture_policy
    if native is None or policy is None or execution.model_sha256 != native.model_sha256:
        raise ValueError("calibrated native execution differs from the frozen policy")
    if execution.background_bytes != configuration.bundle.background_bytes:
        raise ValueError("calibrated native execution background differs from the frozen bytes")
    if background_path is not None and read_regular_path(background_path) != execution.background_bytes:
        raise ValueError("calibrated native background snapshot differs from the frozen bytes")
    return tuple(calibrated_policy_argv(policy, configuration.bundle.caller_policy, background_path))
