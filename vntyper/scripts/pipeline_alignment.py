"""Prepare pipeline alignment targets and exact preflight requests."""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

from vntyper.scripts import preflight_error_io as error_io
from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import AlignmentPlan, preflight_error_payload
from vntyper.scripts.alignment_preflight import run_preflight
from vntyper.scripts.alignment_target_io import install_generated_bed
from vntyper.scripts.fastq_bam_processing import parse_contigs_from_header
from vntyper.scripts.pipeline_guards import enforce_declared_assembly, read_alignment_header
from vntyper.scripts.preflight_input_io import configured_preflight_text_limit, read_bounded_regular_text
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution, restore_reference_resolution
from vntyper.scripts.reference_uri_policy import (
    LocalHeaderReference,
    allow_ambient_reference_resolution,
    enforce_header_reference_policy,
    first_remote_header_reference,
    local_header_reference_paths,
    local_header_references,
)
from vntyper.scripts.region_utils import get_region_string_with_fallback

logger = logging.getLogger(__name__)


class AlignmentPreflightKwargs(TypedDict):
    """Keyword arguments shared by input and post-alignment preflight calls."""

    in_path: str
    output_dir: str
    output_name: str
    file_format: str
    config: dict
    threads: int
    region: None
    bed_file: Path
    coverage_region: str | None
    reference_assembly: str
    reference_fasta: str | None
    header_reference_paths: tuple[str, ...]
    header_references: tuple[LocalHeaderReference, ...]
    has_remote_header_reference: bool
    header_contigs: tuple[str, ...]
    m5: str | None
    header_m5s: tuple[tuple[str, str], ...]
    fast_mode: bool
    error_output_dir: str


@dataclass(frozen=True)
class PreparedAlignmentPreflight:
    """Header, exact target and proven plan produced by the input seam."""

    alignment_header: str | None
    bed_file: Path
    coverage_region: str
    plan: AlignmentPlan
    previous_ref_path: str | None


@dataclass(frozen=True)
class SummaryReferenceProvenance:
    """The reference identity a run's summary and report may honestly claim.

    Attributes:
        key_used: The `reference_data` config key that supplied the reference, or
            None when the run has no comparable "configured key" to report - which is
            every case but FASTQ (see `resolve_summary_reference_provenance`).
        path: The reference path the run actually used, or None when it read no
            reference file at all (a BAM run) or when the winning candidate was an
            ambient htslib resolution with no run-local path (a CRAM run).
        source_effective: For FASTQ, the reference source ("ucsc", "ncbi" or
            "ensembl") the resolved BWA reference belongs to. For BAM and CRAM, the
            alignment plan's own source label (e.g. "not-required", "cli",
            "config_bwa_reference", "htslib_resolved", "header_ur", "private M5
            cache") - a different vocabulary, but the same field, because it already
            carries exactly the meaning "where this run's reference identity came
            from".
    """

    key_used: str | None
    path: str | None
    source_effective: str | None


def resolve_summary_reference_provenance(
    *,
    input_type: str,
    bwa_reference_key: str | None,
    bwa_reference_path: str | None,
    bwa_reference_source: str | None,
    alignment_plan: AlignmentPlan | None,
) -> SummaryReferenceProvenance:
    """Describe the reference a run actually used, not the one BWA would have.

    `pipeline.py` used to record the resolved BWA reference in `pipeline_summary.json`
    for every input type, including BAM and CRAM runs that never open it at all (MAJOR
    5, milestone-5 PR-2 review): a BAM run reads no reference whatsoever, and a CRAM run
    decodes against whatever `alignment_plan` resolved - an explicit `--reference-fasta`,
    a different configured CRAM candidate, an ambient htslib resolution, or a private M5
    cache - which can differ entirely from the configured BWA path. Naming the BWA
    selection in either case is a false provenance claim in a machine-readable artefact,
    in a milestone whose thesis is that the recorded reference is the one actually used.

    A FASTQ run is the one case where the BWA reference *is* what was used to align, so
    its three values pass through unchanged - `key_used`'s only remaining reason to be a
    field distinct from `path` at all, since neither BAM nor CRAM has a comparable
    "configured key" to report.

    For BAM and CRAM, `alignment_plan.reference_path` and `.reference_source` already
    describe exactly what was proven at preflight - `None`/`"not-required"` for BAM,
    and whichever candidate `alignment_preflight.resolve_reference` actually proved for
    CRAM - so they are read through directly rather than re-derived. `key_used` is
    always None for both: there is no `reference_data` key equivalent to report, and
    inventing one would misrepresent an alignment-plan source label as a config key.

    Args:
        input_type: Pipeline input kind (``"FASTQ"``, ``"BAM"`` or ``"CRAM"``).
        bwa_reference_key: The `reference_data` key that resolved the BWA reference, as
            returned by `cli_handlers._resolve_bwa_reference`.
        bwa_reference_path: The resolved BWA reference path.
        bwa_reference_source: The effective reference source for the resolved BWA
            reference.
        alignment_plan: The proven `AlignmentPlan` for a BAM/CRAM run. Only FASTQ is
            expected to call this with None - `run_pipeline` always has a plan for
            BAM/CRAM by the time it builds the summary - and that case degrades to
            "nothing to claim" rather than raising, since a missing plan is itself
            already a sign nothing was proven yet.

    Returns:
        SummaryReferenceProvenance: What this run actually used, not what BWA was
        configured with.
    """
    if input_type == "FASTQ":
        return SummaryReferenceProvenance(bwa_reference_key, bwa_reference_path, bwa_reference_source)
    if alignment_plan is None:
        return SummaryReferenceProvenance(None, None, None)
    return SummaryReferenceProvenance(None, alignment_plan.reference_path, alignment_plan.reference_source)


def format_regions_as_bed(regions: str, *, convert_from_one_based: bool) -> str:
    """Validate comma-separated regions and format them as BED rows.

    User-supplied inline regions are 1-based and fully closed, so their starts are
    decremented for BED. Shipped default windows retain their literal starts because
    established real-data routing contracts include that historical interpretation.

    Args:
        regions: Comma-separated target regions.
        convert_from_one_based: Whether to decrement each validated start for BED.

    Returns:
        BED-formatted text with one target per line.

    Raises:
        ValueError: If a region has an invalid chromosome or coordinates.
    """
    rows: list[str] = []
    for region in regions.split(","):
        try:
            chrom, positions = region.strip().split(":")
            start_text, end_text = positions.strip().split("-")
            if (
                not chrom
                or any(character.isspace() for character in chrom)
                or not start_text.isascii()
                or not start_text.isdigit()
                or not end_text.isascii()
                or not end_text.isdigit()
            ):
                raise ValueError("region fields must be safe ASCII tokens")
        except ValueError as error:
            message = f"Invalid region format: {region}. Expected format 'chr:start-end'."
            logger.error(message)
            raise ValueError(message) from error
        start, end = int(start_text), int(end_text)
        if start < 1 or end < start:
            message = f"Invalid region coordinates: {region}. Expected 1 <= start <= end."
            logger.error(message)
            raise ValueError(message)
        bed_start = start - 1 if convert_from_one_based else start
        rows.append(f"{chrom}\t{bed_start}\t{end}\n")
    return "".join(rows)


def _first_active_bed_contig(bed_text: str) -> str | None:
    for line in bed_text.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return stripped.split(maxsplit=1)[0]
    return None


def _header_m5s(header: str) -> tuple[tuple[str, str], ...]:
    result: list[tuple[str, str]] = []
    for line in header.splitlines():
        fields = line.split("\t")
        if not fields or fields[0] != "@SQ":
            continue
        tags: dict[str, str] = {}
        for field in fields[1:]:
            key, separator, value = field.partition(":")
            if separator:
                tags[key] = value
        contig = tags.get("SN")
        checksum = tags.get("M5")
        if contig is not None and checksum is not None:
            result.append((contig, checksum))
    return tuple(result)


def prepare_alignment_target(
    *,
    input_type: str,
    bam: str | None,
    cram: str | None,
    output_dir: str | Path,
    reference_assembly: str,
    config: dict,
    bed_file: str | Path | None,
    custom_regions: str | None,
) -> Path:
    """Resolve and materialize the exact BED target later given to preflight.

    Args:
        input_type: Pipeline input kind (``BAM``, ``CRAM``, or ``FASTQ``).
        bam: BAM input path when selected.
        cram: CRAM input path when selected.
        output_dir: Pipeline output root for generated BED files.
        reference_assembly: Declared input assembly.
        config: Pipeline configuration.
        bed_file: Operator-provided BED path, when present.
        custom_regions: Operator-provided comma-separated regions, when present.

    Returns:
        The existing or generated target BED path.

    Raises:
        FileNotFoundError: If an operator-provided BED does not exist.
        ValueError: If a generated target is absent or malformed.
    """
    input_path = bam if input_type == "BAM" else cram if input_type == "CRAM" else None
    file_format = input_type.lower() if input_path is not None else None
    if bed_file:
        bed_file_path = Path(bed_file)
        if not bed_file_path.exists():
            logger.error(f"Provided BED file does not exist: {bed_file_path}")
            raise FileNotFoundError(f"BED file not found: {bed_file_path}")
        bed_text = read_bounded_regular_text(
            bed_file_path,
            max_bytes=configured_preflight_text_limit(config),
            description="alignment target BED",
        )
        owned_bed_path = Path(output_dir) / "operator_regions.bed"
        install_generated_bed(
            owned_bed_path,
            bed_text,
            input_path=input_path,
            file_format=file_format,
        )
        logger.info(f"Materialized provided BED file for this run: {owned_bed_path}")
        return owned_bed_path

    if custom_regions:
        bed_file_path = Path(output_dir) / "custom_regions.bed"
        regions = custom_regions
        description = "Custom regions"
        convert_from_one_based = True
    else:
        if input_type in {"BAM", "CRAM"}:
            input_file = bam if input_type == "BAM" else cram
            if input_file is None:
                message = f"{input_type} input path is required for target resolution."
                logger.error(message)
                raise ValueError(message)
            regions = get_region_string_with_fallback(
                bam_file=input_file,
                reference_assembly=reference_assembly,
                region_type="bam_region",
                config=config,
            )
        else:
            region_key = f"bam_region_{reference_assembly}"
            regions = config.get("bam_processing", {}).get(region_key)
            if not regions:
                message = f"Missing configuration for region: {region_key}"
                logger.error(f"Region key '{region_key}' not found in config.json under 'bam_processing'.")
                raise ValueError(message)
        bed_file_path = Path(output_dir) / f"predefined_regions_{reference_assembly}.bed"
        description = "Predefined regions"
        convert_from_one_based = False

    install_generated_bed(
        bed_file_path,
        format_regions_as_bed(regions, convert_from_one_based=convert_from_one_based),
        input_path=input_path,
        file_format=file_format,
    )
    logger.info(f"{description} converted to BED file: {bed_file_path}")
    return bed_file_path


def build_alignment_preflight_kwargs(
    *,
    in_path: str | Path,
    output_dir: str | Path,
    output_name: str,
    file_format: str,
    config: dict,
    threads: int,
    bed_file: Path,
    coverage_region: str | None = None,
    reference_assembly: str,
    fast_mode: bool,
    alignment_header: str | None = None,
    reference_fasta: str | Path | None = None,
    error_output_dir: str | Path | None = None,
) -> AlignmentPreflightKwargs:
    """Build the exact, BED-shaped request used by a real alignment preflight.

    Args:
        in_path: Alignment path to prove.
        output_dir: Run-local alignment working directory.
        output_name: Basename for the view and preflight logs.
        file_format: ``bam`` or ``cram``.
        config: Pipeline configuration.
        threads: Samtools thread count.
        bed_file: Exact BED later used for the real slice.
        coverage_region: Independently resolved region used by the later depth consumer.
        reference_assembly: Declared input assembly.
        fast_mode: Whether downstream processing permits any htslib index.
        alignment_header: Header already read by the assembly guard, if available.
        reference_fasta: Explicit CRAM reference candidate, when supplied.
        error_output_dir: Run output root for the curated failure artifact.

    Returns:
        Typed keyword arguments for :func:`alignment_preflight.run_preflight`.
    """
    header_contigs = (
        tuple(str(contig["name"]) for contig in parse_contigs_from_header(alignment_header) if "name" in contig)
        if alignment_header is not None
        else ()
    )
    bed_text = read_bounded_regular_text(
        bed_file,
        max_bytes=configured_preflight_text_limit(config),
        description="alignment target BED",
    )
    target_contig = _first_active_bed_contig(bed_text)
    header_m5s = _header_m5s(alignment_header) if alignment_header is not None else ()
    target_m5 = dict(header_m5s).get(target_contig) if target_contig is not None else None
    header_reference_paths = (
        local_header_reference_paths(alignment_header, in_path)
        if file_format == "cram" and alignment_header is not None
        else ()
    )
    header_references = (
        local_header_references(alignment_header, in_path)
        if file_format == "cram" and alignment_header is not None
        else ()
    )
    has_remote_header_reference = (
        first_remote_header_reference(alignment_header) is not None
        if file_format == "cram" and alignment_header is not None
        else False
    )
    return {
        "in_path": str(in_path),
        "output_dir": str(output_dir),
        "output_name": output_name,
        "file_format": file_format,
        "config": config,
        "threads": threads,
        "region": None,
        "bed_file": bed_file,
        "coverage_region": coverage_region,
        "reference_assembly": reference_assembly,
        "reference_fasta": str(reference_fasta) if reference_fasta is not None else None,
        "header_reference_paths": header_reference_paths,
        "header_references": header_references,
        "has_remote_header_reference": has_remote_header_reference,
        "header_contigs": header_contigs,
        "m5": target_m5,
        "header_m5s": header_m5s,
        "fast_mode": fast_mode,
        "error_output_dir": str(error_output_dir) if error_output_dir is not None else str(output_dir),
    }


def prepare_input_alignment_preflight(
    *,
    in_path: str | Path,
    input_type: str,
    output_dir: str | Path,
    config: dict,
    threads: int,
    reference_assembly: str,
    bed_file: str | Path | None,
    custom_regions: str | None,
    reference_fasta: str | Path | None,
    fast_mode: bool,
    alignment_validator: Callable[..., object] | None = None,
    validation_cwd: str | None = None,
) -> PreparedAlignmentPreflight:
    """Bind the input, then own validation, header, target and proof preparation.

    Args:
        in_path: Validated BAM or CRAM input path.
        input_type: Alignment kind, ``BAM`` or ``CRAM``.
        output_dir: Pipeline output root.
        config: Pipeline configuration.
        threads: Samtools thread count.
        reference_assembly: Declared input assembly.
        bed_file: Operator-provided BED path, when present.
        custom_regions: Operator-provided comma-separated regions, when present.
        reference_fasta: Explicit CRAM reference candidate, when supplied.
        fast_mode: Whether downstream processing permits any htslib index.
        alignment_validator: Quickcheck-compatible alignment validation callable.
        validation_cwd: Working directory for the validation subprocess.

    Returns:
        The reusable header, exact BED, proven plan, and REF_PATH restoration token.

    Raises:
        ValueError: If ``input_type`` is not an alignment format or preparation fails.
        RuntimeError: If alignment preflight cannot prove a required input.
    """
    if input_type not in {"BAM", "CRAM"}:
        raise ValueError(f"Alignment preflight requires BAM or CRAM input, got: {input_type}")
    input_path = str(in_path)
    is_cram = input_type == "CRAM"
    failure_context = error_io.PreflightErrorContext(output_dir)
    with error_io.persist_preflight_failure(failure_context):
        previous_ref_path = None
        binding: AlignmentBinding | None = None
        if is_cram:
            failure_context.phase = error_io.REFERENCE_POLICY_FAILURE
            previous_ref_path = pin_reference_resolution(config)
        try:
            failure_context.phase = error_io.VIEW_INDEX_FAILURE
            file_format = input_type.lower()
            stage_dir = Path(output_dir) / "fastq_bam_processing"
            stage_dir.mkdir(parents=True, exist_ok=True)
            bound_view = stage_dir / f"input.{file_format}"
            binding = AlignmentBinding(input_path)
            binding.install_view(bound_view)
            failure_context.phase = error_io.HEADER_PREPARATION_FAILURE
            if alignment_validator is not None:
                alignment_validator(str(bound_view), cwd=validation_cwd, log_dir=str(output_dir))
            alignment_header = read_alignment_header(str(bound_view), config)
            if input_type == "CRAM":
                if alignment_header is None:
                    message = (
                        "CRAM header could not be read; reference URI policy and declared assembly cannot be proven."
                    )
                    logger.error(message)
                    raise ValueError(message)
                failure_context.phase = error_io.REFERENCE_POLICY_FAILURE
                try:
                    enforce_header_reference_policy(
                        alignment_header,
                        allow_ambient=allow_ambient_reference_resolution(config),
                    )
                except ValueError as error:
                    failure_context.payload = preflight_error_payload("reference_policy_invalid", str(error), ())
                    raise
            failure_context.phase = error_io.HEADER_PREPARATION_FAILURE
            enforce_declared_assembly(reference_assembly, alignment_header)
            failure_context.phase = error_io.TARGET_PREPARATION_FAILURE
            exact_bed = prepare_alignment_target(
                input_type=input_type,
                bam=str(bound_view) if input_type == "BAM" else None,
                cram=str(bound_view) if input_type == "CRAM" else None,
                output_dir=output_dir,
                reference_assembly=reference_assembly,
                config=config,
                bed_file=bed_file,
                custom_regions=custom_regions,
            )
            coverage_region = get_region_string_with_fallback(
                bam_file=str(bound_view),
                reference_assembly=reference_assembly,
                region_type="vntr_region",
                config=config,
            )
            plan = run_preflight(
                **build_alignment_preflight_kwargs(
                    in_path=input_path,
                    output_dir=Path(output_dir) / "fastq_bam_processing",
                    output_name="input",
                    file_format=input_type.lower(),
                    config=config,
                    threads=threads,
                    bed_file=exact_bed,
                    coverage_region=coverage_region,
                    reference_assembly=reference_assembly,
                    fast_mode=fast_mode,
                    alignment_header=alignment_header,
                    reference_fasta=reference_fasta,
                    error_output_dir=output_dir,
                ),
                failure_context=failure_context,
                binding=binding,
                bound_view_path=str(bound_view),
            )
        except BaseException:
            if binding is not None and binding.is_open:
                try:
                    binding.close()
                except Exception:
                    logger.exception("Alignment binding cleanup failed while preserving preflight failure.")
            if is_cram:
                restore_reference_resolution(previous_ref_path)
            raise
    return PreparedAlignmentPreflight(alignment_header, exact_bed, coverage_region, plan, previous_ref_path)
