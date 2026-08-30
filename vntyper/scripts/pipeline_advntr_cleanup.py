"""Plan adVNTR preflight cleanup and protect the active pipeline log."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path

from vntyper.scripts.pipeline_inputs import archive_base_name

logger = logging.getLogger(__name__)

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
    """Every exact pathname an adVNTR preflight cleanup may revoke."""

    public_outputs: tuple[Path, ...]
    archive: AdvntrArchiveCleanup | None

    @property
    def destinations(self) -> tuple[Path, ...]:
        """Return the complete cleanup-owned pathname set."""
        if self.archive is None:
            return self.public_outputs
        return (*self.public_outputs, self.archive.destination)


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


def plan_advntr_cleanup(
    output_dir: str | Path,
    *,
    archive_results: bool,
    archive_format: str,
) -> AdvntrCleanupPlan:
    """Build the single authoritative adVNTR preflight-cleanup destination set.

    Args:
        output_dir: Pipeline result root.
        archive_results: Whether this run selected a sibling result archive.
        archive_format: Selected public archive format, ``zip`` or ``tar.gz``.

    Returns:
        The exact public outputs and optional selected archive cleanup owns.

    Raises:
        ValueError: If a selected archive format is unsupported.
    """
    output_root = Path(output_dir)
    public_outputs = tuple(output_root / relative_path for relative_path in _PUBLIC_OUTPUTS)
    if not archive_results:
        return AdvntrCleanupPlan(public_outputs=public_outputs, archive=None)

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
    return AdvntrCleanupPlan(public_outputs=public_outputs, archive=archive)


def validate_pipeline_log_outside_advntr_cleanup(
    log_file: str | Path | None,
    cleanup_plan: AdvntrCleanupPlan,
) -> None:
    """Refuse an active log that cleanup could unlink through any path alias.

    Args:
        log_file: Active application-log destination, or None when logging to a file
            was not requested.
        cleanup_plan: Exact destinations selected for this run's preflight cleanup.

    Raises:
        ValueError: If the log names or aliases a cleanup-owned destination.
    """
    if log_file is None:
        return
    log_path = Path(log_file)
    log_variants = _path_variants(log_path)
    for destination in cleanup_plan.destinations:
        destination_variants = _path_variants(destination)
        if any(
            log_variant == destination_variant
            for log_variant in log_variants
            for destination_variant in destination_variants
        ) or _same_file(log_path, destination):
            message = f"Pipeline log file aliases an adVNTR preflight cleanup destination: {log_path}"
            logger.error(message)
            raise ValueError(message)
