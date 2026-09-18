"""Snapshot and protect approved caller runtime inputs at the pipeline boundary."""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from vntyper.modules.advntr.advntr_decision_config import project_advntr_settings
from vntyper.modules.advntr.advntr_genotyping import resolve_advntr_threads
from vntyper.scripts.calibration_atomic_io import atomic_output
from vntyper.scripts.calibration_caller_bundle import caller_model_bundle_files, load_caller_model_bundle
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.pipeline_caller_configuration import (
    CallerPipelineConfiguration,
    caller_runtime_context_document,
    validate_caller_run_options,
)
from vntyper.scripts.run_configuration import RunConfiguration

logger = logging.getLogger(__name__)


def caller_calibration_identities(configuration: CallerPipelineConfiguration) -> dict[str, str]:
    """Return the two exact approval/context commitments used by report and resume.

    Args:
        configuration: Resolved approved caller calibration configuration.

    Returns:
        Path-free bundle and context digests for analysis settings.

    Raises:
        ValueError: If the typed runtime context no longer matches its commitment.
    """
    caller_runtime_context_document(configuration.context)
    return {
        "caller_calibration_bundle_sha256": configuration.bundle.sha256,
        "caller_calibration_context_sha256": configuration.context.sha256,
    }


def validate_caller_output_destination(configuration: CallerPipelineConfiguration, destination: Path) -> None:
    """Refuse output or log destinations overlapping operator calibration inputs.

    Args:
        configuration: Resolved caller inputs whose paths remain operator-owned.
        destination: Pipeline output directory or proposed log file.

    Raises:
        ValueError: If output contains, equals, or falls within a protected input.
    """
    output = destination.resolve(strict=False)
    for path in configuration.operator_paths:
        protected = path.resolve(strict=False)
        if output == protected or output in protected.parents or protected in output.parents:
            raise ValueError("caller runtime output overlaps operator-owned calibration inputs")
        if output.exists() and protected.exists() and os.path.samefile(output, protected):
            raise ValueError("caller runtime output aliases operator-owned calibration inputs")


def snapshot_caller_calibration(configuration: CallerPipelineConfiguration, output: Path) -> Path | None:
    """Install the exact frozen portable bundle and separate context without clobber.

    Args:
        configuration: Approved and resolved caller configuration.
        output: Existing run-owned output directory.

    Returns:
        Exact-mode background snapshot path, otherwise None.

    Raises:
        ValueError: If ownership or existing snapshot identity differs.
        OSError: If the separate no-clobber context file cannot be written.
    """
    validate_caller_output_destination(configuration, output)
    files = caller_model_bundle_files(configuration.bundle)
    context = canonical_json_bytes(caller_runtime_context_document(configuration.context))
    destination = output / "caller_calibration"
    context_path = output / "caller_calibration_context.json"
    if os.path.lexists(context_path) and read_regular_path(context_path) != context:
        raise ValueError("caller runtime context snapshot differs from the resolved context")
    if os.path.lexists(destination):
        if load_caller_model_bundle(destination).sha256 != configuration.bundle.sha256:
            raise ValueError("caller runtime bundle snapshot differs from the approved bundle")
    else:

        def write_bundle(staging: Path) -> bool:
            for name, raw in files.items():
                (staging / name).write_bytes(raw)
            if load_caller_model_bundle(staging).sha256 != configuration.bundle.sha256:
                raise ValueError("caller runtime staged bundle identity differs")
            return True

        atomic_output(destination, write_bundle)
    if not os.path.lexists(context_path):
        descriptor = os.open(context_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(context)
            stream.flush()
            os.fsync(stream.fileno())
    return destination / "background.json" if configuration.bundle.background_bytes is not None else None


def validate_caller_pipeline_request(
    configuration: RunConfiguration,
    *,
    assembly: str,
    extra_modules: Sequence[str],
    threads: int,
    additional_commands: str | None,
    output: Path,
) -> tuple[Path, ...]:
    """Check effective caller options and ownership before pipeline asset/read I/O.

    Args:
        configuration: One complete resolved pipeline configuration.
        assembly: Actual selected reference assembly.
        extra_modules: Explicit optional-stage selection.
        threads: Pipeline thread budget.
        additional_commands: Optional native operator override.
        output: Proposed run output directory.

    Returns:
        Operator-owned calibration paths to protect throughout the pipeline.

    Raises:
        ValueError: If options or output ownership differ from the approved contract.
    """
    calibration = configuration.caller_calibration
    if calibration is None:
        return ()
    settings = project_advntr_settings(configuration.advntr, configuration.advntr_runtime)
    effective_additional = settings.additional_commands if additional_commands is None else additional_commands
    validate_caller_run_options(
        calibration,
        assembly=assembly,
        extra_modules=extra_modules,
        threads=resolve_advntr_threads(settings.command_mapping(), threads),
        additional_commands=effective_additional,
    )
    validate_caller_output_destination(calibration, output)
    return calibration.operator_paths
