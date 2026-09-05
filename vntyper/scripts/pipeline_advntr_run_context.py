"""Prepare immutable, run-owned adVNTR inputs before pipeline stage work."""

from __future__ import annotations

import logging
import os
import shutil
import stat
import tempfile
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

from vntyper.modules.advntr import model_provenance
from vntyper.scripts.archive_safety import revoke_public_archive
from vntyper.scripts.pipeline_advntr_cleanup import AdvntrCleanupPlan, plan_advntr_cleanup

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AdvntrRunContext:
    """Immutable adVNTR state shared by startup checks and later execution."""

    model_source: str
    model_snapshot: str
    model: Mapping[str, Any]
    version: tuple[int, int, int]
    tools: Mapping[str, Any]


def _revoke_prior_outputs(
    cleanup_plan: AdvntrCleanupPlan,
    *,
    protected_paths: Iterable[str | Path],
) -> None:
    """Withdraw only the public names that can advertise a previous adVNTR attempt."""
    from vntyper.modules.advntr.advntr_result_io import invalidate_advntr_artifact

    for destination in cleanup_plan.public_outputs:
        invalidate_advntr_artifact(destination)
    if cleanup_plan.archive is not None:
        if not os.path.lexists(cleanup_plan.archive.destination.parent):
            return
        revoke_public_archive(
            cleanup_plan.archive.base_name,
            cleanup_plan.archive.shutil_format,
            protected_paths=protected_paths,
        )


def _revoke_published_outputs(
    cleanup_plan: AdvntrCleanupPlan,
    *,
    protected_paths: Iterable[str | Path],
) -> None:
    """Withdraw published reports and archives without mutating stage checkpoints."""
    from vntyper.modules.advntr.advntr_result_io import invalidate_advntr_artifact

    for destination in cleanup_plan.published_reports:
        invalidate_advntr_artifact(destination)
    if cleanup_plan.archive is not None:
        if not os.path.lexists(cleanup_plan.archive.destination.parent):
            return
        revoke_public_archive(
            cleanup_plan.archive.base_name,
            cleanup_plan.archive.shutil_format,
            protected_paths=protected_paths,
        )


def _snapshot_model(source: str | Path, destination: str | Path) -> dict[str, Any]:
    """Copy, validate, hash, and atomically install one run-owned model snapshot."""
    source_path = Path(source)
    destination_path = Path(destination)
    temporary_path: Path | None = None
    installed = False
    try:
        if not stat.S_ISREG(source_path.stat().st_mode):
            raise ValueError(f"adVNTR model source is not a regular file: {source_path}")
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination_path.parent,
            prefix=f".{destination_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as candidate:
            temporary_path = Path(candidate.name)
            os.fchmod(candidate.fileno(), 0o644)
            with source_path.open("rb") as model_source:
                shutil.copyfileobj(model_source, candidate)
            candidate.flush()
            os.fsync(candidate.fileno())

        provenance = model_provenance.describe_model(temporary_path)
        os.replace(temporary_path, destination_path)
        installed = True
    except model_provenance.AdvntrModelError:
        raise
    except Exception as error:
        message = f"Failed to install run-owned adVNTR model snapshot {destination_path}: {error}"
        logger.error(message)
        raise RuntimeError(message) from error
    finally:
        if temporary_path is not None and not installed:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError as cleanup_error:
                logger.error(f"Failed to remove incomplete adVNTR model snapshot {temporary_path}: {cleanup_error}")

    provenance["path"] = str(destination_path)
    provenance["source_path"] = str(source_path)
    return provenance


def prepare_advntr_run_context(
    output_dir: str | Path,
    model_source: str | Path,
    config: dict[str, Any],
    *,
    archive_results: bool,
    archive_format: str,
    protected_paths: Iterable[str | Path] = (),
    revoke_outputs: bool = True,
    revoke_published: bool = False,
) -> AdvntrRunContext:
    """Prepare the exact adVNTR model, command, and version used by this run.

    Input and destination ownership must already be validated. Cleanup then happens
    before every model or version refusal, and the returned mappings cannot observe
    later caller mutations.

    Args:
        output_dir: Pipeline result root.
        model_source: Operator-selected adVNTR database.
        config: Complete pipeline configuration whose tools mapping is copied.
        archive_results: Whether this run selected an archive destination.
        archive_format: Selected public archive format, ``zip`` or ``tar.gz``.
        protected_paths: Operator-owned paths the selected archive cannot revoke.
        revoke_outputs: Whether to revoke prior outputs before starting.
        revoke_published: Whether to revoke only published reports and archives
            (used during resume to avoid invalidating stage checkpoints).

    Returns:
        Frozen run context binding the validated snapshot and verified command.

    Raises:
        RuntimeError: If cleanup, snapshot I/O, or version probing fails.
        AdvntrModelError: If the snapshot or binary/model combination is incompatible.
        ValueError: If the version outcome violates its typed invariant.
    """
    cleanup_plan = plan_advntr_cleanup(
        output_dir,
        archive_results=archive_results,
        archive_format=archive_format,
    )
    if revoke_outputs:
        _revoke_prior_outputs(cleanup_plan, protected_paths=protected_paths)
    elif revoke_published:
        _revoke_published_outputs(cleanup_plan, protected_paths=protected_paths)

    snapshot_path = cleanup_plan.model_snapshot
    model = _snapshot_model(model_source, snapshot_path)

    tools = MappingProxyType(dict(config.get("tools", {})))
    probe_config = {**config, "tools": dict(tools)}
    outcome = model_provenance.detect_advntr_version(probe_config, probe=model_provenance.AdvntrVersionProbe())
    if outcome.status not in {
        model_provenance.AdvntrProbeStatus.VERIFIED,
        model_provenance.AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE,
    }:
        model_provenance.require_verified_advntr_version(outcome)
    version = model_provenance.require_compatible_advntr_outcome(model, outcome)

    return AdvntrRunContext(
        model_source=str(Path(model_source)),
        model_snapshot=str(snapshot_path),
        model=MappingProxyType(dict(model)),
        version=version,
        tools=tools,
    )
