# vntyper/scripts/fastq_bam_processing.py

from __future__ import annotations

import json
import logging
import os
import subprocess
from pathlib import Path

from vntyper.scripts.alignment_consumer_commands import build_plan_slice_command, build_plan_unmapped_command
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.alignment_fastq_conversion import (
    plan_alignment_fastq_conversion,
    run_alignment_fastq_conversion,
)
from vntyper.scripts.alignment_target_io import (
    remove_validated_slice_indexes,
    validate_alignment_conversion_destinations,
    validate_fastq_processing_destinations,
)
from vntyper.scripts.artifact_publish import discard_partial, partial_path, publish_partial
from vntyper.scripts.command_builders import (
    build_fastp_command,
    build_samtools_depth_command,
    build_samtools_index_command,
    build_samtools_merge_command,
    build_threaded_samtools_index_argv,
)
from vntyper.scripts.coverage_qc import evaluate_coverage_qc
from vntyper.scripts.coverage_stats import (
    format_coverage_summary,
    parse_region_length,
    read_depth_positions,
    read_depth_values,
    summarise_coverage,
    vntr_geometry,
)
from vntyper.scripts.region_utils import get_region_string_with_fallback
from vntyper.scripts.utils import run_command

logger = logging.getLogger(__name__)


def process_fastq(fastq_1, fastq_2, threads, output, output_name, config):
    """
    Process one or two FASTQ files using fastp for quality control.

    Args:
        fastq_1 (str or Path): Path to the first FASTQ file.
        fastq_2 (str or Path or None): Optional path to the second FASTQ file.
        threads (int): Requested thread count. Dedup-enabled fastp is serialized
            stage-locally; dedup-disabled fastp uses this value unchanged.
        output (str or Path): Output directory.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary containing tool paths and parameters.

    Raises:
        ValueError: If a derived FASTQ destination is unsafe or aliases an input.
        RuntimeError: If FASTQ quality control fails.
    """
    validate_fastq_processing_destinations(output, output_name, fastq_1, fastq_2)
    qc_command = build_fastp_command(
        fastp_path=config["tools"]["fastp"],
        threads=threads,
        fastq_1=fastq_1,
        fastq_2=fastq_2,
        output=output,
        output_name=output_name,
        compression_level=config["bam_processing"]["compression_level"],
        qualified_quality_phred=config["bam_processing"]["qualified_quality_phred"],
        dup_calc_accuracy=config["bam_processing"]["dup_calc_accuracy"],
        length_required=config["bam_processing"]["length_required"],
        disable_adapter_trimming=config["bam_processing"]["disable_adapter_trimming"],
        deduplication=config["bam_processing"]["deduplication"],
    )

    log_file = Path(output) / f"{output_name}_fastp.log"
    logger.info(f"Executing FASTQ quality control with command: {qc_command}")

    success = run_command(str(qc_command), str(log_file), critical=True)
    if not success:
        logger.error("FASTQ quality control failed.")
        raise RuntimeError("FASTQ quality control failed.")

    logger.info("Quality control passed for FASTQ files.")


def process_bam_to_fastq(
    output,
    output_name,
    threads,
    config,
    plan: AlignmentPlan,
    reference_assembly="hg19",
    fast_mode=False,
    delete_intermediates=True,
    bed_file=None,
    needs_advntr: bool = False,
):
    """
    Process alignment files by slicing, filtering, and converting to FASTQ.

    Args:
        output (str or Path): Output directory.
        output_name (str): Base name for the output files.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary containing tool paths and parameters.
        plan: Proven alignment paths, index, reference, and unmapped-scan decision.
        reference_assembly (str, optional): Reference assembly used
            ("hg19", "hg38", "GRCh37", or "GRCh38"). Defaults to "hg19".
        fast_mode (bool, optional): If True, skips filtering of unmapped and partially
            mapped reads. Defaults to False.
        delete_intermediates (bool, optional): If True, deletes intermediate files after
            processing. Defaults to True.
        bed_file (Path, optional): Path to a BED file specifying regions for MUC1 analysis.
        needs_advntr (bool, optional): Whether adVNTR will read ``<name>_sliced.bam``.
            Its index has exactly one consumer -- ``run_advntr`` and
            ``downsample_bam_if_needed``, both behind ``--extra-modules advntr`` -- so
            writing it in Kestrel-only mode is dead work. Coverage reads the alignment
            plan's own view, not this file. A validated boolean rather than the module
            list, because this stage has no business knowing about module state beyond
            that one question. Defaults to False.
    Returns:
        tuple: Paths to the generated FASTQ files (R1, R2, other, single).

    Raises:
        ValueError: If any derived output path is unsafe to overwrite.
        RuntimeError: If any step in the processing fails.
    """
    protected_inputs = (bed_file,) if bed_file is not None else ()
    validate_alignment_conversion_destinations(output, output_name, plan, protected_inputs=protected_inputs)
    remove_validated_slice_indexes(output, output_name)
    samtools_path = config["tools"]["samtools"]

    if bed_file:
        if not bed_file.exists():
            logger.error(f"Provided BED file does not exist: {bed_file}")
            raise FileNotFoundError(f"BED file not found: {bed_file}")
        bam_region = f"-L {bed_file}"
        logger.debug(f"BAM regions set using BED file: {bam_region}")
    else:
        # Use dynamic region resolution with fallback to legacy format
        bam_region = get_region_string_with_fallback(
            bam_file=plan.view_path, reference_assembly=reference_assembly, region_type="bam_region", config=config
        )
        logger.debug(f"BAM region set to: {bam_region}")

    final_bam = Path(output) / f"{output_name}_sliced.bam"
    partial_slice = partial_path(final_bam)

    command_slice = build_plan_slice_command(
        samtools_path=samtools_path,
        plan=plan,
        output_bam=partial_slice,
        region=None if bed_file else bam_region,
        bed_file=bed_file,
        threads=threads,
        fast_mode=fast_mode,
        needs_advntr=needs_advntr,
    )
    log_file_slice = Path(output) / f"{output_name}_slice.log"
    logger.info(f"Executing region slicing with command: {command_slice}")

    discard_partial(partial_slice)
    try:
        success = run_command(str(command_slice), str(log_file_slice), critical=True)
        if not success:
            logger.error(f"{plan.file_format.upper()} region slicing failed.")
            discard_partial(partial_slice)
            raise RuntimeError(f"{plan.file_format.upper()} region slicing failed.")
    except BaseException:
        discard_partial(partial_slice)
        raise

    if fast_mode:
        publish_partial(partial_slice, final_bam)
    logger.info("BAM/CRAM region slicing completed.")

    # Extract & merge unmapped reads if not in fast_mode
    if not fast_mode:
        unmapped_bam = Path(output) / f"{output_name}_unmapped.bam"
        partial_unmapped = partial_path(unmapped_bam)

        command_filter = build_plan_unmapped_command(
            samtools_path=samtools_path,
            plan=plan,
            output_bam=partial_unmapped,
            threads=threads,
            uncompressed=delete_intermediates,
        )
        log_file_filter = Path(output) / f"{output_name}_filter.log"
        logger.info(f"Executing filtering with command: {command_filter}")

        discard_partial(partial_unmapped)
        try:
            success = run_command(str(command_filter), str(log_file_filter), critical=True)
            if not success:
                logger.error("BAM/CRAM filtering failed.")
                discard_partial(partial_unmapped)
                raise RuntimeError("BAM/CRAM filtering failed.")
        except BaseException:
            discard_partial(partial_unmapped)
            raise

        publish_partial(partial_unmapped, unmapped_bam)

        # Merge sliced + unmapped
        merged_bam = Path(output) / f"{output_name}_sliced_unmapped.bam"
        partial_merged = partial_path(merged_bam)
        command_merge = build_samtools_merge_command(
            samtools_path=samtools_path,
            merged_bam=partial_merged,
            sliced_bam=partial_slice,
            unmapped_bam=unmapped_bam,
            threads=threads,
        )
        log_file_merge = Path(output) / f"{output_name}_merge.log"
        logger.info(f"Executing BAM merging with command: {command_merge}")

        discard_partial(partial_merged)
        try:
            success = run_command(str(command_merge), str(log_file_merge), critical=True)
            if not success:
                logger.error("BAM merging failed.")
                discard_partial(partial_merged)
                raise RuntimeError("BAM merging failed.")
        except BaseException:
            discard_partial(partial_merged)
            raise
        finally:
            discard_partial(partial_slice)

        publish_partial(partial_merged, merged_bam)
        final_bam = merged_bam
        logger.info("BAM/CRAM filtering and merging completed.")

        # Rename merged BAM for adVNTR consistency
        final_bam_renamed = Path(output) / f"{output_name}_sliced.bam"
        os.replace(final_bam, final_bam_renamed)
        final_bam = final_bam_renamed
        logger.info(f"Renamed merged BAM file to {final_bam}")

    # Sliced BAM indexing: required only when adVNTR will read the alignment.
    # Converged across fast_mode and non-fast mode.
    if needs_advntr:
        final_bai = Path(f"{final_bam}.bai")
        partial_bai = partial_path(final_bai)
        command_index = build_samtools_index_command(
            samtools_path=samtools_path,
            bam_file=final_bam,
            output_bai=partial_bai,
            threads=threads,
        )
        log_file_index = Path(output) / f"{output_name}_index.log"
        logger.info(f"Re-indexing BAM file with command: {command_index}")

        discard_partial(partial_bai)
        try:
            if not run_command(command_index, str(log_file_index), critical=True):
                logger.error("Re-indexing BAM file failed.")
                discard_partial(partial_bai)
                raise RuntimeError("Re-indexing BAM file failed.")
        except BaseException:
            discard_partial(partial_bai)
            raise

        if not partial_bai.exists() or partial_bai.stat().st_size == 0:
            logger.error(f"BAM index file {partial_bai} not created or empty.")
            discard_partial(partial_bai)
            raise RuntimeError(f"BAM index file {partial_bai} not created or empty.")

        publish_partial(partial_bai, final_bai)
        logger.info(f"Published BAM index to {final_bai}")

    paths = plan_alignment_fastq_conversion(output=output, output_name=output_name)
    conversion_result = run_alignment_fastq_conversion(
        paths=paths,
        final_bam=final_bam,
        samtools_path=samtools_path,
        threads=threads,
        command_runner=run_command,
    )

    # Clean up intermediates if requested
    if delete_intermediates:
        logger.info("Removing intermediate BAM files...")
        intermediate_files = [
            Path(output) / f"{output_name}_unmapped.bam",
        ]
        for file in intermediate_files:
            if file.exists():
                file.unlink()
                logger.debug(f"Removed intermediate file: {file}")
        logger.info("Intermediate BAM files removed.")

    return conversion_result


def _region_is_window(region: str, window: tuple[int, int]) -> bool:
    """Is this region string exactly the configured VNTR window?

    Compared by coordinates rather than by string, because the contig prefix varies with
    the reference's naming convention (``chr1``, ``1``, ``NC_000001.11``) while the window
    does not. A mismatch means an operator-supplied region, not a naming difference.
    """
    _, _, coords = region.rpartition(":")
    start, _, end = coords.partition("-")
    try:
        return (int(start), int(end)) == window
    except ValueError:
        return False


def calculate_vntr_coverage(
    bam_file,
    region,
    threads,
    config,
    output_dir,
    output_name,
    summary_filename=None,
    reference_path=None,
    index_path=None,
    assembly_config=None,
):
    """
    Calculate the coverage over the VNTR region using samtools depth and write a TSV summary.

    Args:
        bam_file (str or Path): Path to the BAM file.
        region (str): Genomic region in 'chr:start-end' format.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary containing tool paths.
        output_dir (str or Path): Directory to store the coverage output.
        output_name (str): Base name for the coverage output file.
        summary_filename (str or Path, optional): File name for the TSV coverage summary.
            Defaults to "<output_name>_summary.tsv" in output_dir.
        reference_path (str or Path, optional): Proven reference FASTA for CRAM decoding.
        index_path (str or Path, optional): Exact retained BAI or CRAI for custom-index depth.
        assembly_config (dict, optional): The ``bam_processing.assemblies`` entry for this
            run. Supplies ``vntr_array_coords``, without which the build-comparable columns
            are recorded as not-measured rather than computed (#222).

    Returns:
        dict: Exactly the keys in :data:`~vntyper.scripts.coverage_stats.COVERAGE_COLUMNS`
            - the region-wide statistics from
            :func:`~vntyper.scripts.coverage_stats.summarise_coverage` plus the
            ``coverage_qc`` verdict this function merges in.

    Raises:
        RuntimeError: If coverage calculation fails.

    Note:
        The QC verdict is evaluated on the **published** figures - ``mean`` and
        ``percent_uncovered`` rounded to the two decimal places the TSV carries - so the
        emitted column can never contradict the report, which re-reads those rounded
        strings and recomputes the same verdict (#172).
    """
    samtools_path = config["tools"]["samtools"]
    # Parse and bound the region BEFORE samtools runs, not after. `depth -a` writes one row
    # per declared base, so a reversed or oversized region exhausts disk before Python ever
    # allocates anything -- validating it ten lines further down, after `run_command`, is too
    # late to prevent the failure it is meant to catch (#224).
    #
    # This raises outside the `try` below, so an invalid region propagates as the ValueError
    # `parse_region_length` raises rather than being re-wrapped as "Error calculating coverage
    # summary". That is the more accurate report: a misconfigured region is a configuration
    # error, not a failure to summarise. `pipeline_coverage` closes the alignment plan in a
    # `finally`, so cleanup is unaffected by the change of type.
    total_region_length = parse_region_length(region)
    if total_region_length <= 0:
        # `parse_region_length` degrades an unparseable region to 0 and only warns, which is
        # its documented contract and is right for a pure function. It is NOT sufficient here:
        # `samtools -r` accepts a bare contig, so passing `chr1` through unchecked emits
        # chromosome-wide depth and recreates the exact disk exhaustion the span bound exists
        # to prevent. The bound has to be enforced where the command is launched (#224).
        msg = (
            f"Region {region!r} could not be parsed into coordinates, so its span cannot be bounded. "
            "Refusing to run `samtools depth` on it: a region without an end is unbounded, and "
            "`depth -a` writes one row per base in whatever it is given."
        )
        logger.error(msg)
        raise ValueError(msg)
    coverage_output = Path(output_dir) / f"{output_name}_vntr_coverage.txt"
    depth_command = build_samtools_depth_command(
        samtools_path=samtools_path,
        threads=threads,
        region=region,
        bam_file=bam_file,
        coverage_output=coverage_output,
        reference_path=reference_path,
        index_path=index_path,
    )
    logger.info(f"Calculating VNTR coverage with command: {depth_command}")
    success = run_command(
        str(depth_command),
        str(coverage_output.with_suffix(".depth.log")),
        critical=True,
    )
    if not success:
        logger.error("VNTR coverage calculation failed.")
        raise RuntimeError("VNTR coverage calculation failed.")

    try:
        # The build-comparable columns are slices of this same depth table, never a second
        # `samtools depth` call: two calls over one BAM could disagree, and a wider one would
        # read outside the region CRAM preflight proved a reference against.
        #
        # Geometry is used only when the region actually is the configured VNTR window. A
        # `--custom-regions` or `--bed-file` run resolves a different region, and the array
        # coordinates do not describe it - so those runs record the columns as not-measured
        # rather than reporting figures for an interval nobody asked about (#222).
        geometry = vntr_geometry(assembly_config) if assembly_config else None
        if geometry is not None and not _region_is_window(region, geometry.window):
            logger.info(
                "Coverage region %s is not the configured VNTR window %s-%s, so the "
                "build-comparable columns are recorded as not measured.",
                region,
                geometry.window[0],
                geometry.window[1],
            )
            geometry = None

        if geometry is None:
            coverage_values = read_depth_values(coverage_output)
            stats = summarise_coverage(coverage_values, total_region_length)
        else:
            pairs = read_depth_positions(coverage_output)
            coverage_values = [depth for _, depth in pairs]
            depth_at = dict(pairs)
            # `.get(..., 0)` mirrors `summarise_coverage`'s own zero-padding: a depth table
            # truncated at a contig end is short, and the missing positions are uncovered.
            array_depths = [depth_at.get(position, 0) for position in range(geometry.array[0], geometry.array[1] + 1)]
            flank_depths = [
                depth_at.get(position, 0) for start, end in geometry.flank for position in range(start, end + 1)
            ]
            stats = summarise_coverage(
                coverage_values,
                total_region_length,
                array_depths=array_depths,
                flank_depths=flank_depths,
                flank_bases=geometry.flank_bases,
            )

        thresholds = config.get("thresholds", {})
        # `.get` with the shipped defaults rather than `[...]`: `--config-path` replaces
        # the whole config instead of merging (AGENTS.md trap 2), so a caller-supplied
        # config legitimately lacks these keys and a KeyError here would abort a run over
        # a display threshold.
        #
        # The `round` calls are the point: the verdict is evaluated on the same figures
        # `format_coverage_summary` is about to write, so the emitted column and the
        # report's recomputed screening axis cannot disagree at a boundary (#172).
        qc = evaluate_coverage_qc(
            round(stats["mean"], 2),
            round(stats["percent_uncovered"], 2),
            thresholds.get("mean_vntr_coverage", 100),
            thresholds.get("percent_vntr_uncovered", 50.0),
        )
        stats["coverage_qc"] = qc.status

        logger.info(f"Mean VNTR coverage: {stats['mean']:.2f}")
        logger.info(f"Median VNTR coverage: {stats['median']:.2f}")
        logger.info(f"Standard deviation: {stats['stdev']:.2f}")
        logger.info(f"Min coverage: {stats['min']}")
        logger.info(f"Max coverage: {stats['max']}")
        logger.info(f"VNTR region total length: {stats['region_length']} bp")
        logger.info(f"VNTR region uncovered bases: {stats['uncovered_bases']} bp")
        logger.info(f"Percentage of VNTR region with zero coverage: {stats['percent_uncovered']:.2f}%")
        logger.info(f"Coverage QC: {stats['coverage_qc']}")

        if summary_filename is None:
            summary_filename = Path(output_dir) / f"{output_name}_summary.tsv"
        else:
            summary_filename = Path(summary_filename)

        with open(summary_filename, "w") as out_f:
            out_f.write(format_coverage_summary(stats))
        logger.info(f"Coverage summary written to: {summary_filename}")

        return stats
    except Exception as e:
        logger.error(f"Error calculating coverage summary: {e}")
        raise RuntimeError(f"Error calculating coverage summary: {e}") from e


def downsample_bam_if_needed(
    bam_path,
    max_coverage,
    reference_assembly,
    threads,
    config,
    coverage_dir,
    coverage_prefix,
):
    """
    Check the current coverage of 'bam_path' in the VNTR region and
    downsample if coverage exceeds 'max_coverage'.

    Args:
        bam_path (str or Path): Path to the input BAM file.
        max_coverage (int): The maximum coverage threshold (e.g., 300).
        reference_assembly (str): Reference assembly ("hg19", "hg38", "GRCh37", or "GRCh38")
            to pick correct region.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary with samtools paths, etc.
        coverage_dir (Path or str): Directory where coverage logs can be written.
        coverage_prefix (str): Prefix for coverage logs (e.g., 'advntr_precheck').

    Returns:
        Path: The path to the (optionally) downsampled BAM.
    """
    from pathlib import Path

    bam_path = Path(bam_path)  # ensure it's a Path object

    # Use dynamic region resolution with fallback to legacy format
    region = get_region_string_with_fallback(
        bam_file=str(bam_path), reference_assembly=reference_assembly, region_type="vntr_region", config=config
    )

    current_coverage = calculate_vntr_coverage(
        bam_file=str(bam_path),
        region=region,
        threads=threads,
        config=config,
        output_dir=coverage_dir,
        output_name=coverage_prefix,
    )["mean"]

    if current_coverage <= max_coverage:
        logger.info(
            f"Current coverage ({current_coverage:.2f}) <= max_coverage ({max_coverage}). No downsampling needed."
        )
        return bam_path

    fraction = max_coverage / current_coverage
    logger.info(
        f"Current coverage: {current_coverage:.2f}, max coverage: {max_coverage}, downsampling fraction: {fraction:.4f}"
    )

    samtools_path = config["tools"]["samtools"]
    downsampled_bam = bam_path.parent / (bam_path.stem + "_downsampled.bam")
    partial_downsampled_bam = partial_path(downsampled_bam)
    seed = 42
    subsample_param = f"{seed}.{int(fraction * 1000):03d}"

    cmd_view = [
        samtools_path,
        "view",
        "-s",
        subsample_param,
        "-@",
        str(threads),
        "-b",
        "-o",
        str(partial_downsampled_bam),
        str(bam_path),
    ]
    logger.info(f"Downsampling BAM with command: {' '.join(cmd_view)}")
    discard_partial(partial_downsampled_bam)
    try:
        subprocess.run(cmd_view, check=True)
        publish_partial(partial_downsampled_bam, downsampled_bam)
    except (subprocess.CalledProcessError, OSError) as err:
        discard_partial(partial_downsampled_bam)
        logger.error(f"Downsampling failed: {err}")
        return bam_path

    sorted_down_bam = downsampled_bam.with_suffix(".sorted.bam")
    partial_sorted_down_bam = partial_path(sorted_down_bam)
    cmd_sort = [
        samtools_path,
        "sort",
        "-@",
        str(threads),
        "-o",
        str(partial_sorted_down_bam),
        str(downsampled_bam),
    ]
    discard_partial(partial_sorted_down_bam)
    final_down_bai = Path(f"{sorted_down_bam}.bai")
    partial_down_bai = partial_path(final_down_bai)
    discard_partial(partial_down_bai)
    try:
        subprocess.run(cmd_sort, check=True)
        publish_partial(partial_sorted_down_bam, sorted_down_bam)
        downsampled_bam.unlink(missing_ok=True)
        cmd_index = build_threaded_samtools_index_argv(
            samtools_path=samtools_path,
            bam_file=sorted_down_bam,
            threads=threads,
            output_bai=partial_down_bai,
        )
        subprocess.run(cmd_index, check=True)
        if not partial_down_bai.exists() or partial_down_bai.stat().st_size == 0:
            raise OSError(f"Downsampled index file {partial_down_bai} not created or empty.")
        publish_partial(partial_down_bai, final_down_bai)
    except (subprocess.CalledProcessError, OSError) as err:
        discard_partial(partial_sorted_down_bam)
        discard_partial(partial_down_bai)
        logger.error(f"Sorting/indexing failed after downsampling: {err}")
        return bam_path

    logger.info(f"Downsampling complete. Using BAM: {sorted_down_bam}")
    return sorted_down_bam


def parse_contigs_from_header(header: str) -> list:
    """
    Parses the BAM header to extract contig information from lines starting with '@SQ'.
    Returns a list of dictionaries with keys 'name' and 'length'.
    """
    contigs = []
    for line in header.splitlines():
        if line.startswith("@SQ"):
            parts = line.split("\t")
            contig_info: dict[str, str | int | None] = {}
            for part in parts:
                if part.startswith("SN:"):
                    contig_info["name"] = part.replace("SN:", "")
                elif part.startswith("LN:"):
                    try:
                        contig_info["length"] = int(part.replace("LN:", ""))
                    except ValueError:
                        contig_info["length"] = None
            if "name" in contig_info and contig_info.get("length") is not None:
                contigs.append(contig_info)
    return contigs


def detect_assembly_from_contigs(header: str, config: dict, threshold: float | None = None) -> str:
    """
    Detects the reference genome assembly by comparing contig information from the BAM header
    against known assemblies from config. Returns the detected assembly name if the match percentage
    is above the threshold, otherwise returns 'Not detected'.

    Args:
        header (str): BAM header string
        config (dict): Configuration dictionary containing assembly_detection section
        threshold (float, optional): Match percentage threshold (0.0-1.0).
            If not provided, uses value from config (default 0.9)

    Returns:
        str: Detected assembly name or "Not detected"
    """
    # Load known assemblies and threshold from config
    assembly_config = config.get("assembly_detection", {})
    known_assemblies = assembly_config.get("known_assemblies", {})

    if threshold is None:
        threshold = assembly_config.get("detection_threshold", 0.9)

    if not known_assemblies:
        logger.warning("No known_assemblies found in config. Assembly detection disabled.")
        return "Not detected"

    bam_contigs = parse_contigs_from_header(header)
    for assembly_data in known_assemblies.values():
        expected_contigs = assembly_data["contigs"]
        match_count = 0
        for expected in expected_contigs:
            for contig in bam_contigs:
                if contig["name"] == expected["name"] and contig["length"] == expected["length"]:
                    match_count += 1
                    break
        match_percentage = match_count / len(expected_contigs)
        if match_percentage >= threshold:
            return assembly_data["name"]
    return "Not detected"


def parse_header_pipeline_info(
    header: str, output_dir: Path, config: dict, output_name: str = "pipeline_info.json"
) -> None:
    """
    Parses the BAM header to extract assembly and alignment pipeline information.
    Uses both text matching and contig matching to detect the assembly.
    Warns if the Dragen or CLC pipeline is detected or if the alignment pipeline cannot be detected,
    recommending the use of BWA aligner.
    Writes the extracted information as a JSON file to the specified output directory.

    Parameters:
        header (str): The BAM header.
        output_dir (Path): Directory where the output file will be written.
        config (dict): Configuration dictionary containing assembly_detection section.
        output_name (str): Base name for the output file. Defaults to 'pipeline_info.json'.
    """
    lower_header = header.lower()

    # Text matching for assembly detection
    if "hg19" in lower_header or "hs37" in lower_header or "grch37" in lower_header:
        assembly_text = "hg19"
    elif "hg38" in lower_header or "hs38" in lower_header or "grch38" in lower_header or "hs38dh" in lower_header:
        assembly_text = "hg38"
    else:
        assembly_text = "Not detected"  # Contig matching
    assembly_contig = detect_assembly_from_contigs(header, config)

    # Determine the alignment pipeline.
    if "dragen" in lower_header:
        pipeline = "Dragen"
    elif "clc" in lower_header or "clcbio" in lower_header:
        pipeline = "CLC"
    elif "bwa" in lower_header:
        pipeline = "BWA"
    else:
        pipeline = "Unknown"

    warning_message = ""
    if pipeline.lower() == "dragen":
        warning_message = (
            "WARNING: The Dragen pipeline has known issues aligning reads in the VNTR region. "
            "It is recommended to use normal mode."
        )
        logger.warning(warning_message)
    elif pipeline.lower() == "clc":
        warning_message = (
            "WARNING: The CLC pipeline may have issues aligning reads in the VNTR region (Not tested). "
            "It is recommended to use normal mode and verify results carefully."
        )
        logger.warning(warning_message)
    elif pipeline.lower() == "unknown":
        warning_message = (
            "WARNING: Alignment pipeline could not be detected from the header. "
            "It is recommended to use the BWA aligner."
        )
        logger.warning(warning_message)

    result = {
        "assembly_text": assembly_text,
        "assembly_contig": assembly_contig,
        "alignment_pipeline": pipeline,
    }
    if warning_message:
        result["warning"] = warning_message

    output_path = output_dir / output_name
    try:
        with open(output_path, "w", encoding="utf-8") as out_f:
            json.dump(result, out_f, indent=4)
        logger.info(f"Pipeline info written to {output_path}")
    except Exception as e:
        logger.error(f"Failed to write pipeline info file: {e}")
        raise


def extract_bam_header(bam_file: str, config: dict) -> str:
    """
    Extracts the header from a BAM file using samtools view -H.

    Args:
        bam_file (str): Path to the BAM file.
        config (dict): Configuration dictionary containing tool paths.

    Returns:
        str: The BAM header.

    Raises:
        subprocess.CalledProcessError: If samtools fails.
    """
    samtools_path = config["tools"]["samtools"]
    cmd = [samtools_path, "view", "-H", bam_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout
