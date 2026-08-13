"""Build post-preflight samtools commands from a proven alignment plan."""

from __future__ import annotations

from pathlib import Path

from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.command_builders import (
    build_cram_unmapped_filter_command,
    build_cram_unmapped_indexed_command,
    build_samtools_slice_command,
)


def build_plan_slice_command(
    *,
    samtools_path: str,
    plan: AlignmentPlan,
    output_bam: str | Path,
    region: str | None,
    bed_file: str | Path | None,
    threads: int,
    fast_mode: bool,
    needs_advntr: bool = False,
) -> str:
    """Build a target slice that selects the plan's retained index.

    Args:
        samtools_path: Configured samtools invocation.
        plan: Proven alignment and retained-index plan.
        output_bam: Slice destination.
        region: Region string when no BED file is supplied.
        bed_file: Optional target BED file.
        threads: Samtools thread count.
        fast_mode: Whether this slice is the run's final alignment. In fast mode there
            is no unmapped recovery and no merge, so the slice survives as
            ``<name>_sliced.bam`` rather than being replaced.
        needs_advntr: Whether adVNTR will read that alignment. The slice index is only
            written in fast mode, and even there its only consumers are ``run_advntr``
            and ``downsample_bam_if_needed``; coverage reads the plan's own view. So
            indexing is gated on both.

    Returns:
        The complete samtools slice command.
    """
    return build_samtools_slice_command(
        samtools_path=samtools_path,
        in_bam=plan.view_path,
        output_bam=output_bam,
        region=region,
        bed_file=bed_file,
        reference_path=plan.reference_path,
        index_path=plan.stable_index_path,
        threads=threads,
        index_output=fast_mode and needs_advntr,
        exclude_unmapped=not fast_mode,
        # Not unconditionally True. In fast mode there is no unmapped recovery and no
        # merge, so this slice is never replaced: it *survives* as `<name>_sliced.bam`,
        # is an enumerated artifact and is shipped to users inside the archive. Writing
        # that at BGZF level 0 would roughly triple it.
        uncompressed=not fast_mode,
    )


def build_plan_unmapped_command(
    *,
    samtools_path: str,
    plan: AlignmentPlan,
    output_bam: str | Path,
    threads: int,
    uncompressed: bool = False,
) -> str:
    """Build the plan-selected indexed or streaming unmapped-read command.

    Args:
        samtools_path: Configured samtools invocation.
        plan: Proven alignment, reference, index, and scan decision.
        output_bam: Unmapped-read BAM destination.
        threads: Samtools thread count.
        uncompressed: Write the unmapped BAM at BGZF level 0. Forwarded to whichever
            scan the plan proved, so the decision cannot depend on that choice.

    Returns:
        The complete samtools view command.
    """
    if plan.unmapped_scan == "indexed":
        return build_cram_unmapped_indexed_command(
            samtools_path=samtools_path,
            in_bam=plan.view_path,
            unmapped_bam=output_bam,
            threads=threads,
            reference_path=plan.reference_path,
            index_path=plan.stable_index_path,
            uncompressed=uncompressed,
        )
    return build_cram_unmapped_filter_command(
        samtools_path=samtools_path,
        in_bam=plan.view_path,
        unmapped_bam=output_bam,
        threads=threads,
        reference_path=plan.reference_path,
        uncompressed=uncompressed,
    )
