# vntyper/scripts/fastq_bam_processing.py

from __future__ import annotations

import json
import logging
import os
import subprocess
from pathlib import Path

from vntyper.scripts.alignment_consumer_commands import build_plan_slice_command, build_plan_unmapped_command
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.alignment_target_io import (
    remove_validated_slice_indexes,
    validate_alignment_conversion_destinations,
    validate_fastq_processing_destinations,
)
from vntyper.scripts.command_builders import (
    build_bam_to_fastq_command,
    build_fastp_command,
    build_samtools_depth_command,
    build_samtools_index_command,
    build_samtools_merge_command,
)
from vntyper.scripts.coverage_qc import evaluate_coverage_qc
from vntyper.scripts.coverage_stats import (
    format_coverage_summary,
    parse_region_length,
    read_depth_values,
    summarise_coverage,
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
        threads (int): Number of threads to use.
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
    keep_intermediates=False,
    bed_file=None,
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
            processing, overriding ``keep_intermediates``. Defaults to True.
        keep_intermediates (bool, optional): If True, keeps intermediate files for later
            use unless ``delete_intermediates`` is True. Defaults to False.
        bed_file (Path, optional): Path to a BED file specifying regions for MUC1 analysis.
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

    command_slice = build_plan_slice_command(
        samtools_path=samtools_path,
        plan=plan,
        output_bam=final_bam,
        region=None if bed_file else bam_region,
        bed_file=bed_file,
        threads=threads,
        fast_mode=fast_mode,
    )
    log_file_slice = Path(output) / f"{output_name}_slice.log"
    logger.info(f"Executing region slicing with command: {command_slice}")

    success = run_command(str(command_slice), str(log_file_slice), critical=True)
    if not success:
        logger.error(f"{plan.file_format.upper()} region slicing failed.")
        raise RuntimeError(f"{plan.file_format.upper()} region slicing failed.")
    logger.info("BAM/CRAM region slicing completed.")

    # Extract & merge unmapped reads if not in fast_mode
    if not fast_mode:
        unmapped_bam = Path(output) / f"{output_name}_unmapped.bam"

        command_filter = build_plan_unmapped_command(
            samtools_path=samtools_path,
            plan=plan,
            output_bam=unmapped_bam,
            threads=threads,
        )
        log_file_filter = Path(output) / f"{output_name}_filter.log"
        logger.info(f"Executing filtering with command: {command_filter}")

        success = run_command(str(command_filter), str(log_file_filter), critical=True)
        if not success:
            logger.error("BAM/CRAM filtering failed.")
            raise RuntimeError("BAM/CRAM filtering failed.")

        # Merge sliced + unmapped
        merged_bam = Path(output) / f"{output_name}_sliced_unmapped.bam"
        command_merge = build_samtools_merge_command(
            samtools_path=samtools_path,
            merged_bam=merged_bam,
            sliced_bam=final_bam,
            unmapped_bam=unmapped_bam,
            threads=threads,
        )
        log_file_merge = Path(output) / f"{output_name}_merge.log"
        logger.info(f"Executing BAM merging with command: {command_merge}")

        success = run_command(str(command_merge), str(log_file_merge), critical=True)
        if not success:
            logger.error("BAM merging failed.")
            raise RuntimeError("BAM merging failed.")

        final_bam = merged_bam
        logger.info("BAM/CRAM filtering and merging completed.")

        # Rename merged BAM for adVNTR consistency and re-index
        final_bam_renamed = Path(output) / f"{output_name}_sliced.bam"
        os.replace(final_bam, final_bam_renamed)
        final_bam = final_bam_renamed
        logger.info(f"Renamed merged BAM file to {final_bam}")

        # No output_bai here, deliberately: final_bam is the merged BAM this stage
        # just wrote inside `output`, so samtools' default destination beside it is
        # already inside the output directory (#162).
        command_index = build_samtools_index_command(
            samtools_path=samtools_path,
            bam_file=final_bam,
            threads=threads,
        )
        log_file_index = Path(output) / f"{output_name}_index.log"
        logger.info(f"Re-indexing BAM file with command: {command_index}")
        if not run_command(command_index, str(log_file_index), critical=True):
            logger.error("Re-indexing BAM file failed.")
            raise RuntimeError("Re-indexing BAM file failed.")

    # Convert final BAM to FASTQ
    final_fastq_1 = Path(output) / f"{output_name}_R1.fastq.gz"
    final_fastq_2 = Path(output) / f"{output_name}_R2.fastq.gz"
    final_fastq_other = Path(output) / f"{output_name}_other.fastq.gz"
    final_fastq_single = Path(output) / f"{output_name}_single.fastq.gz"

    command_sort_fastq = build_bam_to_fastq_command(
        samtools_path=samtools_path,
        in_bam=final_bam,
        threads=threads,
        fastq_r1=final_fastq_1,
        fastq_r2=final_fastq_2,
        fastq_other=final_fastq_other,
        fastq_single=final_fastq_single,
    )
    log_file_sort_fastq = Path(output) / f"{output_name}_sort_fastq.log"
    logger.info(f"Executing BAM to FASTQ conversion with command: {command_sort_fastq}")

    success = run_command(str(command_sort_fastq), str(log_file_sort_fastq), critical=True)
    if not success:
        logger.error("BAM to FASTQ conversion failed.")
        raise RuntimeError("BAM to FASTQ conversion failed.")
    logger.info("BAM to FASTQ conversion completed.")

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

    return (
        str(final_fastq_1),
        str(final_fastq_2),
        str(final_fastq_other),
        str(final_fastq_single),
    )


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
        total_region_length = parse_region_length(region)
        coverage_values = read_depth_values(coverage_output)
        stats = summarise_coverage(coverage_values, total_region_length)

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
        str(downsampled_bam),
        str(bam_path),
    ]
    logger.info(f"Downsampling BAM with command: {' '.join(cmd_view)}")
    try:
        subprocess.run(cmd_view, check=True)
    except subprocess.CalledProcessError as err:
        logger.error(f"Downsampling failed: {err}")
        return bam_path

    sorted_down_bam = downsampled_bam.with_suffix(".sorted.bam")
    cmd_sort = [
        samtools_path,
        "sort",
        "-@",
        str(threads),
        "-o",
        str(sorted_down_bam),
        str(downsampled_bam),
    ]
    try:
        subprocess.run(cmd_sort, check=True)
        downsampled_bam.unlink()
        cmd_index = [samtools_path, "index", str(sorted_down_bam)]
        subprocess.run(cmd_index, check=True)
    except subprocess.CalledProcessError as err:
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
