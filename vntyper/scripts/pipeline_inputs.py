"""Pipeline input selection and pre-write ownership dispatch."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import cast

from vntyper.scripts.alignment_target_io import (
    alignment_operator_paths,
    fastq_operator_paths,
    protect_alignment_inputs,
    validate_fastq_pipeline_destinations,
)

logger = logging.getLogger(__name__)

_ARCHIVE_SUFFIXES = {"zip": ".zip", "tar.gz": ".tar.gz"}


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def archive_base_name(output_dir: str | Path) -> str:
    """Return the normalized suffix-free base name for a result archive.

    Args:
        output_dir: Pipeline output root, possibly ending in path separators.

    Returns:
        The output root without syntactically redundant trailing separators.
    """
    return str(Path(output_dir))


def _validate_archive_destination(
    output_dir: str | Path,
    archive_format: str,
    operator_paths: tuple[str | Path, ...],
    alignment_input: str | Path | None,
) -> None:
    """Reject an optional archive destination that aliases operator-owned state.

    This is a pre-existing-state check at the pipeline boundary. The archive writer
    remains responsible for atomic installation and for rechecking its source tree;
    this function does not attempt to redesign concurrent hostile mutation handling.

    Args:
        output_dir: Pipeline output root used as the suffix-free archive name.
        archive_format: CLI format name, ``zip`` or ``tar.gz``.
        operator_paths: Every selected input, index, reference, and sidecar.
        alignment_input: BAM or CRAM whose containing patient tree is protected.

    Raises:
        ValueError: If the format is unsupported or the destination aliases protected state.
    """
    try:
        suffix = _ARCHIVE_SUFFIXES[archive_format]
    except KeyError:
        message = f"Unsupported archive format: {archive_format}"
        logger.error(message)
        raise ValueError(message) from None

    destination = Path(f"{archive_base_name(output_dir)}{suffix}")
    destination_absolute = _absolute(destination)
    destination_variants = (destination_absolute, destination_absolute.resolve(strict=False))
    for operator_path in operator_paths:
        protected = _absolute(operator_path)
        protected_variants = (protected, protected.resolve(strict=False))
        if any(
            destination_variant == protected_variant
            for destination_variant in destination_variants
            for protected_variant in protected_variants
        ) or _same_file(destination, protected):
            message = f"Archive destination aliases protected operator input: {destination}"
            logger.error(message)
            raise ValueError(message)

    if alignment_input is None:
        return
    alignment = _absolute(alignment_input)
    patient_trees = (alignment.parent, alignment.resolve(strict=False).parent)
    if any(
        _is_within(destination_variant, patient_tree)
        for destination_variant in destination_variants
        for patient_tree in patient_trees
    ):
        message = f"Archive destination must stay outside the patient input tree: {destination}"
        logger.error(message)
        raise ValueError(message)


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
    archive_results: bool,
    archive_format: str,
    additional_operator_paths: tuple[str | Path, ...] = (),
) -> None:
    """Reject unsafe existing output entries before any pipeline mutation.

    This boundary validates the filesystem snapshot that exists before the run.
    Individual writers still revalidate their destinations immediately before
    replacement; concurrent hostile tree mutation is outside this snapshot contract.

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
        archive_results: Whether a sibling result archive will be installed.
        archive_format: Requested archive format, ``zip`` or ``tar.gz``.
        additional_operator_paths: Other operator-owned paths selected before writes.

    Raises:
        ValueError: If an input or existing output-tree entry is unsafe.
    """
    if input_type in {"BAM", "CRAM"}:
        input_alignment = cast(str | Path, bam if input_type == "BAM" else cram)
        operator_paths = (
            *alignment_operator_paths(
                input_alignment,
                input_type.lower(),
                bed_file,
                reference_fasta,
                config,
                reference_assembly,
            ),
            *additional_operator_paths,
        )
        if archive_results:
            _validate_archive_destination(output_dir, archive_format, operator_paths, input_alignment)
        protect_alignment_inputs(
            output_dir,
            input_alignment,
            input_type.lower(),
            bed_file,
            reference_fasta,
            config,
            reference_assembly,
            additional_operator_paths=additional_operator_paths,
        )
        return
    operator_paths = (
        *fastq_operator_paths(
            cast(str | Path, fastq1),
            fastq2,
            bed_file,
            cast(str | Path, bwa_reference),
            config,
        ),
        *additional_operator_paths,
    )
    if archive_results:
        _validate_archive_destination(output_dir, archive_format, operator_paths, None)
    validate_fastq_pipeline_destinations(
        output_dir,
        cast(str | Path, fastq1),
        fastq2,
        bed_file,
        cast(str | Path, bwa_reference),
        config,
        additional_operator_paths=additional_operator_paths,
    )
