#!/usr/bin/env python3
# vntyper/scripts/fastq_bam_processing.py

import logging
import os
from pathlib import Path

from vntyper.scripts.utils import run_command


def process_fastq(
    fastq_1, fastq_2, threads, output, output_name, config
):
    """
    Process FASTQ files using fastp for quality control.

    Args:
        fastq_1 (str or Path): Path to the first FASTQ file.
        fastq_2 (str or Path): Path to the second FASTQ file.
        threads (int): Number of threads to use.
        output (str or Path): Output directory.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary containing tool paths and parameters.

    Raises:
        RuntimeError: If FASTQ quality control fails.
    """
    fastp_path = config["tools"]["fastp"]
    compression_level = config["bam_processing"]["compression_level"]
    disable_adapter_trimming = config["bam_processing"]["disable_adapter_trimming"]
    deduplication = config["bam_processing"]["deduplication"]
    dup_calc_accuracy = config["bam_processing"]["dup_calc_accuracy"]
    length_required = config["bam_processing"]["length_required"]
    qualified_quality_phred = config["bam_processing"]["qualified_quality_phred"]

    qc_command = (
        f"{fastp_path} --thread {threads} --in1 {fastq_1} --in2 {fastq_2} "
        f"--out1 {output}/{output_name}_R1.fastq.gz --out2 {output}/{output_name}_R2.fastq.gz "
        f"--compression {compression_level} "
        f"--qualified_quality_phred {qualified_quality_phred} "
        f"--dup_calc_accuracy {dup_calc_accuracy} "
        f"--length_required {length_required} "
        f"--html {output}/{output_name}.html "
        f"--json {output}/{output_name}.json "
    )

    if disable_adapter_trimming:
        qc_command += " --disable_adapter_trimming"
    if deduplication:
        qc_command += " --dedup"

    log_file = Path(output) / f"{output_name}_fastp.log"
    logging.info(
        f"Executing FASTQ quality control with command: {qc_command}"
    )

    success = run_command(str(qc_command), str(log_file), critical=True)
    if not success:
        logging.error("FASTQ quality control failed.")
        raise RuntimeError("FASTQ quality control failed.")

    logging.info("Quality control passed for FASTQ files.")


def process_bam_to_fastq(
    in_bam,
    output,
    output_name,
    threads,
    config,
    reference_assembly="hg19",
    fast_mode=False,
    delete_intermediates=True,
    keep_intermediates=False,
    bed_file=None,
    file_format="bam"
):
    """
    Process alignment files by slicing, filtering, and converting to FASTQ.

    Args:
        in_bam (str or Path): Path to the input BAM/CRAM file.
        output (str or Path): Output directory.
        output_name (str): Base name for the output files.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary containing tool paths and parameters.
        reference_assembly (str, optional): Reference assembly used ("hg19" or "hg38").
            Defaults to "hg19".
        fast_mode (bool, optional): If True, skips filtering of unmapped and partially
            mapped reads. Defaults to False.
        delete_intermediates (bool, optional): If True, deletes intermediate files after
            processing. Defaults to True.
        keep_intermediates (bool, optional): If True, keeps intermediate files for later
            use. Defaults to False.
        bed_file (Path, optional): Path to a BED file specifying regions for MUC1 analysis.
        file_format (str, optional): "bam" or "cram". Default is "bam". This parameter
            enables CRAM support.

    Returns:
        tuple: Paths to the generated FASTQ files (R1, R2, other, single).

    Raises:
        RuntimeError: If any step in the processing fails.
    """
    samtools_path = config["tools"]["samtools"]

    if bed_file:
        if not bed_file.exists():
            logging.error(f"Provided BED file does not exist: {bed_file}")
            raise FileNotFoundError(f"BED file not found: {bed_file}")
        bam_region = f"-L {bed_file}"
        logging.debug(f"BAM regions set using BED file: {bam_region}")
    else:
        bam_region = (
            config["bam_processing"]["bam_region_hg38"]
            if reference_assembly == "hg38"
            else config["bam_processing"]["bam_region_hg19"]
        )
        logging.debug(f"BAM region set to: {bam_region}")

    cram_ref_option = ""
    final_bam = Path(output) / f"{output_name}_sliced.bam"

    if keep_intermediates and final_bam.exists():
        logging.info(f"Reusing existing BAM slice: {final_bam}")
    else:
        if bed_file:
            command_slice = (
                f"{samtools_path} view -P -b {cram_ref_option} {in_bam} -L {bed_file} -o {final_bam} && "
                f"{samtools_path} index {final_bam}"
            )
        else:
            command_slice = (
                f"{samtools_path} view -P -b {cram_ref_option} {in_bam} {bam_region} -o {final_bam} && "
                f"{samtools_path} index {final_bam}"
            )
        log_file_slice = Path(output) / f"{output_name}_slice.log"
        logging.info(f"Executing region slicing with command: {command_slice}")

        success = run_command(
            str(command_slice), str(log_file_slice), critical=True
        )
        if not success:
            logging.error(f"{file_format.upper()} region slicing failed.")
            raise RuntimeError(f"{file_format.upper()} region slicing failed.")
        logging.info("BAM/CRAM region slicing completed.")

    if not fast_mode:
        command_filter = (
            f"{samtools_path} view {cram_ref_option} -@ {threads} -h {in_bam} | tee "
            f" >(samtools view -b -f 4 -F 264 -@ {threads} - -o {output}/"
            f"{output_name}_unmapped1.bam) "
            f" >(samtools view -b -f 8 -F 260 -@ {threads} - -o {output}/"
            f"{output_name}_unmapped2.bam) "
            f" >(samtools view -b -f 12 -F 265 -@ {threads} - -o {output}/"
            f"{output_name}_unmapped3.bam) "
            f"> /dev/null"
        )
        log_file_filter = Path(output) / f"{output_name}_filter.log"
        logging.info(f"Executing filtering with command: {command_filter}")

        success = run_command(
            str(command_filter), str(log_file_filter), critical=True
        )
        if not success:
            logging.error("BAM/CRAM filtering failed.")
            raise RuntimeError("BAM/CRAM filtering failed.")

        merged_bam = Path(output) / f"{output_name}_sliced_unmapped.bam"
        command_merge = (
            f"{samtools_path} merge -f -@ {threads} {merged_bam} "
            f"{final_bam} {output}/{output_name}_unmapped1.bam "
            f"{output}/{output_name}_unmapped2.bam "
            f"{output}/{output_name}_unmapped3.bam"
        )
        log_file_merge = Path(output) / f"{output_name}_merge.log"
        logging.info(f"Executing BAM merging with command: {command_merge}")

        success = run_command(
            str(command_merge), str(log_file_merge), critical=True
        )
        if not success:
            logging.error("BAM merging failed.")
            raise RuntimeError("BAM merging failed.")

        final_bam = merged_bam
        logging.info("BAM/CRAM filtering and merging completed.")

    final_fastq_1 = Path(output) / f"{output_name}_R1.fastq.gz"
    final_fastq_2 = Path(output) / f"{output_name}_R2.fastq.gz"
    final_fastq_other = Path(output) / f"{output_name}_other.fastq.gz"
    final_fastq_single = Path(output) / f"{output_name}_single.fastq.gz"

    if keep_intermediates and all(
        p.exists()
        for p in [
            final_fastq_1,
            final_fastq_2,
            final_fastq_other,
            final_fastq_single,
        ]
    ):
        logging.info(
            f"Reusing existing FASTQ files: {final_fastq_1}, {final_fastq_2}, "
            f"{final_fastq_other}, and {final_fastq_single}"
        )
    else:
        command_sort_fastq = (
            f"{samtools_path} sort -n -@ {threads} {final_bam} | "
            f"{samtools_path} fastq -@ {threads} - -1 {final_fastq_1} "
            f"-2 {final_fastq_2} -0 {final_fastq_other} "
            f"-s {final_fastq_single}"
        )
        log_file_sort_fastq = Path(output) / f"{output_name}_sort_fastq.log"
        logging.info(
            f"Executing BAM to FASTQ conversion with command: {command_sort_fastq}"
        )

        success = run_command(
            str(command_sort_fastq), str(log_file_sort_fastq), critical=True
        )
        if not success:
            logging.error("BAM to FASTQ conversion failed.")
            raise RuntimeError("BAM to FASTQ conversion failed.")
        logging.info("BAM to FASTQ conversion completed.")

    if delete_intermediates and not keep_intermediates:
        logging.info("Removing intermediate BAM files...")
        intermediate_files = [
            Path(output) / f"{output_name}_unmapped1.bam",
            Path(output) / f"{output_name}_unmapped2.bam",
            Path(output) / f"{output_name}_unmapped3.bam",
        ]
        for file in intermediate_files:
            if file.exists():
                file.unlink()
                logging.debug(f"Removed intermediate file: {file}")
        logging.info("Intermediate BAM files removed.")

    return (
        str(final_fastq_1),
        str(final_fastq_2),
        str(final_fastq_other),
        str(final_fastq_single),
    )


def calculate_vntr_coverage(bam_file, region, threads, config, output_dir, output_name):
    """
    Calculate the mean coverage over the VNTR region using samtools depth.

    Args:
        bam_file (str or Path): Path to the BAM file.
        region (str): Genomic region in 'chr:start-end' format.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary containing tool paths.
        output_dir (str or Path): Directory to store the coverage log.
        output_name (str): Base name for the coverage log file.

    Returns:
        float: Mean coverage over the VNTR region.

    Raises:
        RuntimeError: If coverage calculation fails.
    """
    samtools_path = config["tools"]["samtools"]
    coverage_output = Path(output_dir) / f"{output_name}_vntr_coverage.txt"
    depth_command = (
        f"{samtools_path} depth -@ {threads} -r {region} {bam_file} > {coverage_output}"
    )
    logging.info(f"Calculating VNTR coverage with command: {depth_command}")
    success = run_command(
        str(depth_command), str(coverage_output.with_suffix('.depth.log')), critical=True
    )
    if not success:
        logging.error("VNTR coverage calculation failed.")
        raise RuntimeError("VNTR coverage calculation failed.")

    try:
        with open(coverage_output, 'r') as f:
            coverage_values = [int(line.strip().split('\t')[2]) for line in f if line.strip()]
        mean_coverage = sum(coverage_values) / len(coverage_values) if coverage_values else 0
        logging.info(f"Mean VNTR coverage: {mean_coverage:.2f}")
        return mean_coverage
    except Exception as e:
        logging.error(f"Error calculating mean coverage: {e}")
        raise RuntimeError(f"Error calculating mean coverage: {e}")
