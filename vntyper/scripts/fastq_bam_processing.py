import os
import logging
from vntyper.scripts.utils import run_command


def process_fastq(fastq_1, fastq_2, threads, output, output_name, config):
    """
    Process FASTQ files using fastp for quality control.

    Args:
        fastq_1 (str): Path to the first FASTQ file.
        fastq_2 (str): Path to the second FASTQ file.
        threads (int): Number of threads to use.
        output (str): Output directory.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary containing tool paths and parameters.
    """
    fastp_path = config["tools"]["fastp"]
    compression_level = config["bam_processing"]["compression_level"]
    disable_adapter_trimming = config["bam_processing"]["disable_adapter_trimming"]
    deduplication = config["bam_processing"]["deduplication"]
    dup_calc_accuracy = config["bam_processing"]["dup_calc_accuracy"]
    length_required = config["bam_processing"]["length_required"]

    QC_command = (
        f"{fastp_path} --thread {threads} --in1 {fastq_1} --in2 {fastq_2} "
        f"--out1 {output}/{output_name}_R1.fastq.gz --out2 {output}/{output_name}_R2.fastq.gz "
        f"--compression {compression_level} --dup_calc_accuracy {dup_calc_accuracy} --length_required {length_required} "
        f"--html {output}/{output_name}.html"
    )

    if disable_adapter_trimming:
        QC_command += " --disable_adapter_trimming"
    if deduplication:
        QC_command += " --dedup"

    log_file = os.path.join(output, f"{output_name}_fastp.log")
    if not run_command(QC_command, log_file, critical=True):
        logging.error("FASTQ quality control failed.")
        raise RuntimeError("FASTQ quality control failed.")

    logging.info('Quality control passed for FASTQ files.')


def process_bam_to_fastq(in_bam, output, output_name, threads, config, reference_assembly="hg19", fast_mode=False, delete_intermediates=True, keep_intermediates=False):
    """
    Process BAM files by slicing, filtering, and converting to FASTQ.

    Args:
        in_bam (str): Path to the input BAM file.
        output (str): Output directory.
        output_name (str): Base name for the output files.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary containing tool paths and parameters.
        reference_assembly (str): Reference assembly used ("hg19" or "hg38").
        fast_mode (bool): If True, skips filtering of unmapped and partially mapped reads.
        delete_intermediates (bool): If True, deletes intermediate files after processing.
        keep_intermediates (bool): If True, keeps intermediate files for later use.

    Returns:
        Tuple[str, str, str]: Paths to the generated FASTQ files (R1, R2, and other).
    """
    samtools_path = config["tools"]["samtools"]

    if reference_assembly == "hg38":
        bam_region = config["bam_processing"]["bam_region_hg38"]
    else:
        bam_region = config["bam_processing"]["bam_region_hg19"]

    # Slicing BAM region
    final_bam = f"{output}/{output_name}_sliced.bam"
    if keep_intermediates and os.path.exists(final_bam):
        logging.info(f"Reusing existing BAM slice: {final_bam}")
    else:
        command_slice = f"{samtools_path} view -P -b {in_bam} {bam_region} -o {final_bam}"
        log_file = os.path.join(output, f"{output_name}_slice.log")
        if not run_command(command_slice, log_file, critical=True):
            logging.error("BAM region slicing failed.")
            return None, None, None
        logging.info("BAM region slicing completed.")

    if fast_mode:
        logging.info("Fast mode enabled: Skipping filtering of unmapped and partially mapped reads.")
    else:
        # Filtering step
        command_filter = (
            f"{samtools_path} view -@ {threads} -h {in_bam} | tee "
            f">(samtools view -b -f 4 -F 264 -@ {threads} - -o {output}/{output_name}_unmapped1.bam) "
            f">(samtools view -b -f 8 -F 260 -@ {threads} - -o {output}/{output_name}_unmapped2.bam) "
            f">(samtools view -b -f 12 -F 256 -@ {threads} - -o {output}/{output_name}_unmapped3.bam) "
            f"> /dev/null"
        )
        log_file = os.path.join(output, f"{output_name}_filter.log")
        if not run_command(command_filter, log_file, critical=True):
            logging.error("BAM filtering failed.")
            return None, None, None

        command_merge = (
            f"{samtools_path} merge -@ {threads} {output}/{output_name}_sliced_unmapped.bam "
            f"{final_bam} {output}/{output_name}_unmapped1.bam {output}/{output_name}_unmapped2.bam {output}/{output_name}_unmapped3.bam"
        )
        log_file = os.path.join(output, f"{output_name}_merge.log")
        if not run_command(command_merge, log_file, critical=True):
            logging.error("BAM merging failed.")
            return None, None, None

        final_bam = f"{output}/{output_name}_sliced_unmapped.bam"
        logging.info("BAM filtering and merging completed.")

    # Sorting and converting BAM to FASTQ using pipes
    final_fastq_1 = f"{output}/{output_name}_R1.fastq.gz"
    final_fastq_2 = f"{output}/{output_name}_R2.fastq.gz"
    final_fastq_other = f"{output}/{output_name}_other.fastq.gz"
    final_fastq_single = f"{output}/{output_name}_single.fastq.gz"

    if keep_intermediates and os.path.exists(final_fastq_1) and os.path.exists(final_fastq_2) and os.path.exists(final_fastq_other):
        logging.info(f"Reusing existing FASTQ files: {final_fastq_1}, {final_fastq_2}, {final_fastq_other} and {final_fastq_single}")
    else:
        # FUTURE: Check if the single and other FASTQ files are required
        command_sort_fastq = (
            f"{samtools_path} sort -n -@ {threads} {final_bam} | "
            f"{samtools_path} fastq -@ {threads} - -1 {final_fastq_1} -2 {final_fastq_2} -0 {final_fastq_other} -s {final_fastq_single}"
        )
        log_file = os.path.join(output, f"{output_name}_sort_fastq.log")
        if not run_command(command_sort_fastq, log_file, critical=True):
            logging.error("BAM to FASTQ conversion failed.")
            return None, None, None
        logging.info("BAM to FASTQ conversion completed.")

    # Remove intermediate BAM files if required
    if delete_intermediates and not keep_intermediates:
        logging.info("Removing intermediate BAM files...")
        intermediate_files = [
            f"{output}/{output_name}_unmapped1.bam",
            f"{output}/{output_name}_unmapped2.bam",
            f"{output}/{output_name}_unmapped3.bam",
        ]
        for file in intermediate_files:
            if os.path.exists(file):
                os.remove(file)
        logging.info("Intermediate files removed.")

    # Return the paths to the generated FASTQ files
    return final_fastq_1, final_fastq_2, final_fastq_other, final_fastq_single
