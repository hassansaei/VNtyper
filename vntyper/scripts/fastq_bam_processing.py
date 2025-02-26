#!/usr/bin/env python3
# vntyper/scripts/fastq_bam_processing.py

import logging
import os
from pathlib import Path
import subprocess
import statistics  # for median and stdev calculations
import json

from vntyper.scripts.utils import run_command


def process_fastq(fastq_1, fastq_2, threads, output, output_name, config):
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
    logging.info(f"Executing FASTQ quality control with command: {qc_command}")

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
    file_format="bam",
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

        success = run_command(str(command_slice), str(log_file_slice), critical=True)
        if not success:
            logging.error(f"{file_format.upper()} region slicing failed.")
            raise RuntimeError(f"{file_format.upper()} region slicing failed.")
        logging.info("BAM/CRAM region slicing completed.")

    if not fast_mode:
        command_filter = (
            f"{samtools_path} view {cram_ref_option} -@ {threads} -h {in_bam} | tee "
            f" >(samtools view -b -f 12 -@ {threads} - -o {output}/{output_name}_unmapped.bam) "
            f"> /dev/null"
        )
        log_file_filter = Path(output) / f"{output_name}_filter.log"
        logging.info(f"Executing filtering with command: {command_filter}")

        success = run_command(str(command_filter), str(log_file_filter), critical=True)
        if not success:
            logging.error("BAM/CRAM filtering failed.")
            raise RuntimeError("BAM/CRAM filtering failed.")

        merged_bam = Path(output) / f"{output_name}_sliced_unmapped.bam"
        command_merge = (
            f"{samtools_path} merge -f -@ {threads} {merged_bam} "
            f"{output}/{output_name}_sliced.bam "
            f"{output}/{output_name}_unmapped.bam"
        )
        log_file_merge = Path(output) / f"{output_name}_merge.log"
        logging.info(f"Executing BAM merging with command: {command_merge}")

        success = run_command(str(command_merge), str(log_file_merge), critical=True)
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
            Path(output) / f"{output_name}_unmapped.bam",
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


def calculate_vntr_coverage(
    bam_file, region, threads, config, output_dir, output_name, summary_filename=None
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

    Returns:
        dict: A dictionary containing mean, median, standard deviation, min, and max coverage.

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
        str(depth_command),
        str(coverage_output.with_suffix(".depth.log")),
        critical=True,
    )
    if not success:
        logging.error("VNTR coverage calculation failed.")
        raise RuntimeError("VNTR coverage calculation failed.")

    try:
        with open(coverage_output, "r") as f:
            coverage_values = [
                int(line.strip().split("\t")[2]) for line in f if line.strip()
            ]
        if not coverage_values:
            raise RuntimeError("No coverage data found.")

        mean_coverage = sum(coverage_values) / len(coverage_values)
        median_coverage = statistics.median(coverage_values)
        stdev_coverage = (
            statistics.stdev(coverage_values) if len(coverage_values) > 1 else 0
        )
        min_coverage = min(coverage_values)
        max_coverage = max(coverage_values)

        logging.info(f"Mean VNTR coverage: {mean_coverage:.2f}")
        logging.info(f"Median VNTR coverage: {median_coverage:.2f}")
        logging.info(f"Standard deviation: {stdev_coverage:.2f}")
        logging.info(f"Min coverage: {min_coverage}")
        logging.info(f"Max coverage: {max_coverage}")

        if summary_filename is None:
            summary_filename = Path(output_dir) / f"{output_name}_summary.tsv"
        else:
            summary_filename = Path(summary_filename)

        with open(summary_filename, "w") as out_f:
            out_f.write("mean\tmedian\tstdev\tmin\tmax\n")
            out_f.write(
                f"{mean_coverage:.2f}\t{median_coverage:.2f}\t{stdev_coverage:.2f}\t{min_coverage}\t{max_coverage}\n"
            )
        logging.info(f"Coverage summary written to: {summary_filename}")

        return {
            "mean": mean_coverage,
            "median": median_coverage,
            "stdev": stdev_coverage,
            "min": min_coverage,
            "max": max_coverage,
        }
    except Exception as e:
        logging.error(f"Error calculating coverage summary: {e}")
        raise RuntimeError(f"Error calculating coverage summary: {e}")


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
        reference_assembly (str): "hg19" or "hg38" to pick correct region.
        threads (int): Number of threads to use.
        config (dict): Configuration dictionary with samtools paths, etc.
        coverage_dir (Path or str): Directory where coverage logs can be written.
        coverage_prefix (str): Prefix for coverage logs (e.g., 'advntr_precheck').

    Returns:
        Path: The path to the (optionally) downsampled BAM.
    """
    from pathlib import Path

    bam_path = Path(bam_path)  # ensure it's a Path object

    if reference_assembly == "hg38":
        region = config["bam_processing"]["vntr_region_hg38"]
    else:
        region = config["bam_processing"]["vntr_region_hg19"]

    current_coverage = calculate_vntr_coverage(
        bam_file=str(bam_path),
        region=region,
        threads=threads,
        config=config,
        output_dir=coverage_dir,
        output_name=coverage_prefix,
    )["mean"]

    if current_coverage <= max_coverage:
        logging.info(
            f"Current coverage ({current_coverage:.2f}) <= max_coverage ({max_coverage}). No downsampling needed."
        )
        return bam_path

    fraction = max_coverage / current_coverage
    logging.info(
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
    logging.info(f"Downsampling BAM with command: {' '.join(cmd_view)}")
    try:
        subprocess.run(cmd_view, check=True)
    except subprocess.CalledProcessError as err:
        logging.error(f"Downsampling failed: {err}")
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
        logging.error(f"Sorting/indexing failed after downsampling: {err}")
        return bam_path

    logging.info(f"Downsampling complete. Using BAM: {sorted_down_bam}")
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
            contig_info = {}
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


def detect_assembly_from_contigs(header: str, threshold: float = 0.9) -> str:
    """
    Detects the reference genome assembly by comparing contig information from the BAM header
    against known assemblies. Returns the detected assembly name if the match percentage
    is above the threshold, otherwise returns 'Not detected'.
    """
    known_assemblies = {
        "hg19": {
            "name": "hg19",
            "contigs": [
                {"name": "chr1", "length": 249250621},
                {"name": "chr2", "length": 243199373},
                {"name": "chr3", "length": 198022430},
                {"name": "chr4", "length": 191154276},
                {"name": "chr5", "length": 180915260},
                {"name": "chr6", "length": 171115067},
                {"name": "chr7", "length": 159138663},
                {"name": "chr8", "length": 146364022},
                {"name": "chr9", "length": 141213431},
                {"name": "chr10", "length": 135534747},
                {"name": "chr11", "length": 135006516},
                {"name": "chr12", "length": 133851895},
                {"name": "chr13", "length": 115169878},
                {"name": "chr14", "length": 107349540},
                {"name": "chr15", "length": 102531392},
                {"name": "chr16", "length": 90354753},
                {"name": "chr17", "length": 81195210},
                {"name": "chr18", "length": 78077248},
                {"name": "chr19", "length": 59128983},
                {"name": "chr20", "length": 63025520},
                {"name": "chr21", "length": 48129895},
                {"name": "chr22", "length": 51304566},
                {"name": "chrX", "length": 155270560},
                {"name": "chrY", "length": 59373566},
            ],
        },
        "hg38": {
            "name": "hg38",
            "contigs": [
                {"name": "chr1", "length": 248956422},
                {"name": "chr2", "length": 242193529},
                {"name": "chr3", "length": 198295559},
                {"name": "chr4", "length": 190214555},
                {"name": "chr5", "length": 181538259},
                {"name": "chr6", "length": 170805979},
                {"name": "chr7", "length": 159345973},
                {"name": "chr8", "length": 145138636},
                {"name": "chr9", "length": 138394717},
                {"name": "chr10", "length": 133797422},
                {"name": "chr11", "length": 135086622},
                {"name": "chr12", "length": 133275309},
                {"name": "chr13", "length": 114364328},
                {"name": "chr14", "length": 107043718},
                {"name": "chr15", "length": 101991189},
                {"name": "chr16", "length": 90338345},
                {"name": "chr17", "length": 83257441},
                {"name": "chr18", "length": 80373285},
                {"name": "chr19", "length": 58617616},
                {"name": "chr20", "length": 64444167},
                {"name": "chr21", "length": 46709983},
                {"name": "chr22", "length": 50818468},
                {"name": "chrX", "length": 156040895},
                {"name": "chrY", "length": 57227415},
            ],
        },
    }

    bam_contigs = parse_contigs_from_header(header)
    for assembly_key, assembly_data in known_assemblies.items():
        expected_contigs = assembly_data["contigs"]
        match_count = 0
        for expected in expected_contigs:
            for contig in bam_contigs:
                if (
                    contig["name"] == expected["name"]
                    and contig["length"] == expected["length"]
                ):
                    match_count += 1
                    break
        match_percentage = match_count / len(expected_contigs)
        if match_percentage >= threshold:
            return assembly_data["name"]
    return "Not detected"


def parse_header_pipeline_info(
    header: str, output_dir: Path, output_name: str = "pipeline_info.json"
) -> None:
    """
    Parses the BAM header to extract assembly and alignment pipeline information.
    Uses both text matching and contig matching to detect the assembly.
    Warns if the Dragen pipeline is detected or if the alignment pipeline cannot be detected,
    recommending the use of BWA aligner.
    Writes the extracted information as a JSON file to the specified output directory.

    Parameters:
        header (str): The BAM header.
        output_dir (Path): Directory where the output file will be written.
        output_name (str): Base name for the output file. Defaults to 'pipeline_info.json'.
    """
    lower_header = header.lower()

    # Text matching for assembly detection
    if "hg19" in lower_header or "hs37" in lower_header or "grch37" in lower_header:
        assembly_text = "hg19"
    elif (
        "hg38" in lower_header
        or "hs38" in lower_header
        or "grch38" in lower_header
        or "hs38dh" in lower_header
    ):
        assembly_text = "hg38"
    else:
        assembly_text = "Not detected"

    # Contig matching for assembly (assumes detect_assembly_from_contigs is defined/imported)
    assembly_contig = detect_assembly_from_contigs(header)

    # Determine the alignment pipeline.
    if "dragen" in lower_header:
        pipeline = "Dragen"
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
        logging.warning(warning_message)
    elif pipeline.lower() == "unknown":
        warning_message = (
            "WARNING: Alignment pipeline could not be detected from the header. "
            "It is recommended to use the BWA aligner."
        )
        logging.warning(warning_message)

    # Compose the JSON result.
    result = {
        "assembly_text": assembly_text,
        "assembly_contig": assembly_contig,
        "alignment_pipeline": pipeline,
    }
    if warning_message:
        result["warning"] = warning_message

    # Write the result as JSON.
    output_path = output_dir / output_name
    try:
        with open(output_path, "w", encoding="utf-8") as out_f:
            json.dump(result, out_f, indent=4)
        logging.info(f"Pipeline info written to {output_path}")
    except Exception as e:
        logging.error(f"Failed to write pipeline info file: {e}")
        raise e


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
