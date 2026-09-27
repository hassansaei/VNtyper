"""Why a report has no alignment view, decided once and stated the same way in the log and the report.

The report's alignment view is drawn around the position of a Kestrel call:
``kestrel/output.bed`` names that position, and the Kestrel stage writes it only when a
variant survives the final filter. So a sample with no call has no view, and that is the
normal state, not a failure.

For years that state logged ``BED file does not exist or not provided. Skipping IGV
report generation.`` at WARNING, and the report said only that "this run produced no
alignment session". An external user screening a 3,500-exome cohort read both as a
possible installation fault and asked before trusting the cohort. The two strings could
not tell the normal case apart from the two cases that are not normal:

* Kestrel called a variant and the region file is gone, so a view that should exist was
  not built;
* the region file exists but Kestrel called nothing, so the file was left in the output
  directory by an earlier run and would draw a position this run never called.

:func:`decide_igv_absence` separates the five states. The log line and the report
paragraph both follow from its ``reason``, so they cannot disagree.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from vntyper.scripts import report_assets
from vntyper.scripts.artifact_names import PIPELINE_BASENAME
from vntyper.scripts.summary_steps import STEP_READ

logger = logging.getLogger(__name__)

#: ``--report-igv off`` was in effect. The operator chose this.
IGV_SWITCHED_OFF: Final[str] = "switched_off"
#: Kestrel ran, read cleanly and called no variant: no position to draw. Normal.
IGV_NO_KESTREL_CALL: Final[str] = "no_kestrel_call"
#: The region file exists but Kestrel called nothing, so an earlier run left it behind.
IGV_STALE_REGION_FILE: Final[str] = "stale_region_file"
#: Kestrel called a variant but its region file is absent or was never given.
IGV_REGION_FILE_MISSING: Final[str] = "region_file_missing"
#: The summary holds no readable Kestrel result and no region file was found.
IGV_NO_KESTREL_RESULT: Final[str] = "no_kestrel_result"

#: A normal absence is INFO: a cohort's logs must not carry one WARNING per negative sample.
LOG_LEVELS: Final[dict[str, int]] = {
    IGV_SWITCHED_OFF: logging.INFO,
    IGV_NO_KESTREL_CALL: logging.INFO,
    IGV_STALE_REGION_FILE: logging.WARNING,
    IGV_REGION_FILE_MISSING: logging.WARNING,
    IGV_NO_KESTREL_RESULT: logging.WARNING,
}


@dataclass(frozen=True)
class IgvAbsence:
    """Whether to build the alignment view and, if not, why.

    Attributes:
        reason: One of the ``IGV_*`` constants, or None when the view is built.
        level: The logging level of ``message``.
        message: The log line explaining the absence; empty when the view is built.
        region_file_given: Whether a region file path was handed to the report.
        region_file_from_stage: Whether that path is the Kestrel stage's own file. The
            report names ``kestrel/output.bed`` and suggests a re-run only then.
    """

    reason: str | None
    level: int = logging.INFO
    message: str = ""
    region_file_given: bool = False
    region_file_from_stage: bool = False

    @property
    def build_view(self) -> bool:
        """Whether the report should run the IGV generator."""
        return self.reason is None


def decide_igv_absence(
    *,
    report_igv: str,
    bed_file: str | os.PathLike[str] | None,
    kestrel_state: str,
    kestrel_has_call: bool,
    kestrel_no_call_recorded: bool,
    bed_from_kestrel_stage: bool | None = None,
) -> IgvAbsence:
    """Decide whether a report gets an alignment view, and word the reason when it does not.

    Args:
        report_igv: The ``--report-igv`` mode in effect.
        bed_file: The region file the view would be anchored on, or None.
        kestrel_state: :func:`~vntyper.scripts.summary_steps.get_step_state` for Kestrel.
        kestrel_has_call: Whether the Kestrel step holds at least one row that is not
            the empty-result placeholder.
        kestrel_no_call_recorded: Whether the Kestrel step holds the empty-result
            placeholder a run that genotyped and called nothing writes. Only that is
            evidence of a no-call sample: a step with no rows at all recorded nothing,
            and must not be told "this is expected".
        bed_from_kestrel_stage: Whether ``bed_file`` is the Kestrel stage's own region
            file (the pipeline's, or one ``vntyper report`` discovered) rather than one
            the operator named. Only the stage's file can be left over by an earlier
            run. None falls back to the path: :func:`is_kestrel_region_file`.

    Returns:
        IgvAbsence: ``reason`` None when the view should be built.
    """
    bed_exists = bed_file is not None and str(bed_file) != "" and os.path.exists(bed_file)
    kestrel_read = kestrel_state == STEP_READ
    stage_owned = is_kestrel_region_file(bed_file) if bed_from_kestrel_stage is None else bed_from_kestrel_stage

    if report_igv == report_assets.REPORT_IGV_OFF:
        return _absence(
            IGV_SWITCHED_OFF,
            "--report-igv off: no alignment browser is produced for this run.",
        )
    if kestrel_read and kestrel_no_call_recorded and not kestrel_has_call:
        if bed_exists and not stage_owned:
            # `vntyper report --bed-file` naming the operator's own regions: their choice.
            return IgvAbsence(reason=None)
        if bed_exists:
            return _absence(
                IGV_STALE_REGION_FILE,
                f"No alignment view: Kestrel called no variant in this run, but {bed_file} exists. It was left "
                "by an earlier run into the same output directory and names a position this run did not call, "
                "so it was not used. Delete it, or use a fresh --output-dir per sample.",
            )
        return _absence(
            IGV_NO_KESTREL_CALL,
            "No alignment view for this sample: Kestrel called no variant, and the view is drawn around the "
            "position of a Kestrel call. This is expected for a sample without a call; it is not a sign of a "
            "failed step or an incomplete installation.",
        )
    if bed_exists:
        return IgvAbsence(reason=None)
    if kestrel_read and kestrel_has_call:
        if not bed_file:
            where, fix = "no region file was given", "pass one to `vntyper report --bed-file`"
        elif stage_owned:
            where, fix = (
                f"its region file {bed_file} does not exist",
                "re-run the pipeline for this sample, or pass a region file to `vntyper report --bed-file`",
            )
        else:
            where, fix = f"the region file given, {bed_file}, does not exist", "check the --bed-file path"
        return IgvAbsence(
            reason=IGV_REGION_FILE_MISSING,
            level=LOG_LEVELS[IGV_REGION_FILE_MISSING],
            message=f"No alignment view although Kestrel called a variant: {where}. The variant tables are "
            f"unaffected; to restore the view, {fix}.",
            region_file_given=bool(bed_file),
            region_file_from_stage=stage_owned,
        )
    return _absence(
        IGV_NO_KESTREL_RESULT,
        f"No alignment view: this run's summary holds no readable Kestrel result (state: {kestrel_state}, "
        f"{'no rows recorded' if kestrel_read else 'nothing read'}) and no region file was found, so there is "
        "no position to draw. See the report's Kestrel section.",
    )


def is_kestrel_region_file(bed_file: str | os.PathLike[str] | None) -> bool:
    """Whether ``bed_file`` is the Kestrel stage's own ``kestrel/output.bed``.

    Only that file can be stale: the stage writes it, and a later run into the same
    output directory that calls nothing used to leave it in place. A BED the operator
    names with ``vntyper report --bed-file`` is their choice and is drawn as given.

    Args:
        bed_file: The region file handed to the report, or None.

    Returns:
        bool: True for a path ending in ``kestrel/output.bed``.
    """
    if not bed_file:
        return False
    path = Path(bed_file)
    return path.name == f"{PIPELINE_BASENAME}.bed" and path.parent.name == "kestrel"


def _absence(reason: str, message: str) -> IgvAbsence:
    return IgvAbsence(reason=reason, level=LOG_LEVELS[reason], message=message)
