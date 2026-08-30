"""Plan adVNTR preflight cleanup and protect the active pipeline log."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.pipeline_inputs import archive_base_name

logger = logging.getLogger(__name__)

ADVNTR_MODEL_SNAPSHOT = Path("advntr/advntr_model.db")

_PUBLIC_OUTPUTS = (
    Path("advntr/output_adVNTR_result.tsv"),
    Path("advntr/output_adVNTR.tsv"),
    Path("advntr/output_adVNTR.vcf"),
    Path("advntr/cross_match_results.tsv"),
    Path("advntr/output_advntr.log"),
    Path("pipeline_summary.json"),
    Path("pipeline_summary.csv"),
    Path("pipeline_summary.tsv"),
    Path("summary_report.html"),
)

_ARCHIVE_FORMATS = {
    "zip": (".zip", "zip"),
    "tar.gz": (".tar.gz", "gztar"),
}


@dataclass(frozen=True)
class AdvntrArchiveCleanup:
    """The one selected sibling archive destination and revocation arguments."""

    destination: Path
    base_name: str
    shutil_format: str


@dataclass(frozen=True)
class AdvntrCleanupPlan:
    """Exact adVNTR cleanup names plus its separately destructive snapshot."""

    public_outputs: tuple[Path, ...]
    archive: AdvntrArchiveCleanup | None
    model_snapshot: Path

    @property
    def cleanup_destinations(self) -> tuple[Path, ...]:
        """Return only names stale-result cleanup may revoke."""
        if self.archive is None:
            return self.public_outputs
        return (*self.public_outputs, self.archive.destination)

    @property
    def destructive_destinations(self) -> tuple[Path, ...]:
        """Return every preflight name that may be unlinked or replaced."""
        return (*self.cleanup_destinations, self.model_snapshot)


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _path_variants(path: str | Path) -> tuple[Path, Path]:
    absolute = _absolute(path)
    return absolute, absolute.resolve(strict=False)


def _paths_alias(left: str | Path, right: str | Path) -> bool:
    return any(
        left_variant == right_variant
        for left_variant in _path_variants(left)
        for right_variant in _path_variants(right)
    ) or _same_file(Path(left), Path(right))


def plan_advntr_cleanup(
    output_dir: str | Path,
    *,
    archive_results: bool,
    archive_format: str,
) -> AdvntrCleanupPlan:
    """Build the single authoritative adVNTR destructive-preflight plan.

    Args:
        output_dir: Pipeline result root.
        archive_results: Whether this run selected a sibling result archive.
        archive_format: Selected public archive format, ``zip`` or ``tar.gz``.

    Returns:
        Exact stale-cleanup outputs and the separately destructive model snapshot.

    Raises:
        ValueError: If a selected archive format is unsupported.
    """
    output_root = Path(output_dir)
    public_outputs = tuple(output_root / relative_path for relative_path in _PUBLIC_OUTPUTS)
    model_snapshot = output_root / ADVNTR_MODEL_SNAPSHOT
    if not archive_results:
        return AdvntrCleanupPlan(public_outputs=public_outputs, archive=None, model_snapshot=model_snapshot)

    try:
        suffix, shutil_format = _ARCHIVE_FORMATS[archive_format]
    except KeyError:
        message = f"Unsupported archive format: {archive_format}"
        logger.error(message)
        raise ValueError(message) from None
    base_name = archive_base_name(output_root)
    archive = AdvntrArchiveCleanup(
        destination=Path(f"{base_name}{suffix}"),
        base_name=base_name,
        shutil_format=shutil_format,
    )
    return AdvntrCleanupPlan(public_outputs=public_outputs, archive=archive, model_snapshot=model_snapshot)


def validate_pipeline_log_outside_advntr_preflight(
    log_file: str | Path | None,
    cleanup_plan: AdvntrCleanupPlan,
) -> None:
    """Refuse an active log that preflight could unlink or replace by any alias.

    Args:
        log_file: Active application-log destination, or None when logging to a file
            was not requested.
        cleanup_plan: Exact cleanup and snapshot destinations selected for this run.

    Raises:
        ValueError: If the log names or aliases a destructive destination.
    """
    if log_file is None:
        return
    log_path = Path(log_file)
    for destination in cleanup_plan.destructive_destinations:
        if _paths_alias(log_path, destination):
            message = f"Pipeline log file aliases an adVNTR destructive preflight destination: {log_path}"
            logger.error(message)
            raise ValueError(message)


def validate_pipeline_log_outside_selected_advntr_model(
    log_file: str | Path | None,
    model_source: str | Path | None,
) -> None:
    """Refuse a selected-model log alias without emitting into that model.

    This direct-call guard runs before the pipeline's first log record. Its refusal
    deliberately raises without logging because an attached handler may already point
    at the operator-owned model being protected.

    Args:
        log_file: Active application-log destination, or None when absent.
        model_source: Selected operator adVNTR model, or None when disabled.

    Raises:
        ValueError: If the active log names or aliases the selected model.
    """
    if log_file is None or model_source is None:
        return
    log_path = Path(log_file)
    if _paths_alias(log_path, model_source):
        raise ValueError(f"Pipeline log file aliases selected operator adVNTR model: {log_path}")
