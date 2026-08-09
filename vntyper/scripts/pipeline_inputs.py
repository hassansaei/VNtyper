"""Pipeline input selection and pre-write ownership dispatch."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import cast

from vntyper.scripts.alignment_target_io import protect_alignment_inputs, validate_fastq_pipeline_destinations

logger = logging.getLogger(__name__)


def resolve_pipeline_input(
    fastq1: str | Path | None,
    fastq2: str | Path | None,
    bam: str | Path | None,
    cram: str | Path | None,
    bwa_reference: str | Path | None,
    extra_modules: list[str],
) -> tuple[str, dict[str, str]]:
    """Validate the selected input family and describe it for the run summary.

    Args:
        fastq1: First FASTQ, when FASTQ mode is selected.
        fastq2: Optional second FASTQ.
        bam: BAM input, when selected.
        cram: CRAM input, when selected.
        bwa_reference: Required alignment reference for FASTQ mode.
        extra_modules: Optional pipeline modules requested for the run.

    Returns:
        The input type and basename-only summary mapping.

    Raises:
        ValueError: If the input family is missing, conflicting, or unsupported.
    """
    supplied_groups = [
        input_type
        for input_type, supplied in (
            ("FASTQ", fastq1 is not None or fastq2 is not None),
            ("BAM", bam is not None),
            ("CRAM", cram is not None),
        )
        if supplied
    ]
    if len(supplied_groups) > 1:
        message = "Provide either BAM, CRAM, or FASTQ files, not multiples."
        logger.error(message)
        raise ValueError(message)
    if not supplied_groups:
        message = "No input files provided."
        logger.error(message)
        raise ValueError(message)

    input_type = supplied_groups[0]
    if input_type == "FASTQ" and not fastq1:
        message = "When not providing BAM/CRAM, --fastq1 must be specified; --fastq2 is optional."
        logger.error(message)
        raise ValueError(message)
    if input_type == "FASTQ" and not bwa_reference:
        message = "BWA reference not provided or determined from configuration."
        logger.error(message)
        raise ValueError(message)
    if input_type == "FASTQ" and not fastq2 and "shark" in extra_modules:
        message = "SHARK requires paired-end FASTQ input; provide --fastq2 or remove the shark module."
        logger.error(message)
        raise ValueError(message)

    if input_type == "FASTQ":
        input_files = {"fastq1": os.path.basename(cast(str | Path, fastq1))}
        if fastq2:
            input_files["fastq2"] = os.path.basename(fastq2)
    elif input_type == "BAM":
        input_files = {"bam": os.path.basename(cast(str | Path, bam))}
    else:
        input_files = {"cram": os.path.basename(cast(str | Path, cram))}
    return input_type, input_files


def protect_pipeline_input_ownership(
    output_dir: str | Path,
    input_type: str,
    fastq1: str | Path | None,
    fastq2: str | Path | None,
    bam: str | Path | None,
    cram: str | Path | None,
    bed_file: str | Path | None,
    reference_fasta: str | Path | None,
    bwa_reference: str | Path | None,
    config: dict,
    reference_assembly: str,
) -> None:
    """Reject unsafe existing output entries before any pipeline mutation.

    Args:
        output_dir: Pipeline output root.
        input_type: Validated input type, ``FASTQ``, ``BAM``, or ``CRAM``.
        fastq1: First FASTQ when selected.
        fastq2: Optional second FASTQ.
        bam: BAM input when selected.
        cram: CRAM input when selected.
        bed_file: Optional operator-provided BED.
        reference_fasta: Optional explicit CRAM reference.
        bwa_reference: BWA reference selected for FASTQ mode.
        config: Pipeline configuration.
        reference_assembly: Declared reference assembly.

    Raises:
        ValueError: If an input or existing output-tree entry is unsafe.
    """
    if input_type in {"BAM", "CRAM"}:
        input_alignment = cast(str | Path, bam if input_type == "BAM" else cram)
        protect_alignment_inputs(
            output_dir,
            input_alignment,
            input_type.lower(),
            bed_file,
            reference_fasta,
            config,
            reference_assembly,
        )
        return
    validate_fastq_pipeline_destinations(
        output_dir,
        cast(str | Path, fastq1),
        fastq2,
        bed_file,
        cast(str | Path, bwa_reference),
        config,
    )
