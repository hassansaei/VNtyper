#!/usr/bin/env python3
# vntyper/scripts/pipeline.py

import logging
import os
import shutil
import sys
import timeit
from pathlib import Path
from datetime import datetime

from vntyper.scripts.alignment_processing import align_and_sort_fastq
from vntyper.scripts.fastq_bam_processing import (
    process_bam_to_fastq,
    process_fastq,
    calculate_vntr_coverage,
    downsample_bam_if_needed,
    extract_bam_header,
    parse_header_pipeline_info,
)
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.utils import (
    create_output_directories,
    get_tool_versions,
    validate_bam_file,
    validate_fastq_file,
)
from vntyper.version import __version__ as VERSION

# Import our new summary functions (including end_summary and CSV/TSV conversion functions)
from vntyper.scripts.summary import (
    start_summary,
    record_step,
    write_summary,
    end_summary,
    convert_summary_to_csv,
    convert_summary_to_tsv,
)

# Import cross-match functions from cross_match.py
from vntyper.scripts.cross_match import (
    cross_match_variants,
    write_results_tsv,
    extract_results_from_pipeline_summary,
)


def write_bed_file(regions, bed_file_path):
    """
    Writes regions to a BED file in the correct format.

    Parameters:
      - regions (str): Comma-separated regions in 'chr:start-end' format.
      - bed_file_path (Path): Path to the BED file to be written.
    """
    with open(bed_file_path, "w", encoding="utf-8") as bed_fh:
        for region in regions.split(","):
            try:
                chrom, positions = region.strip().split(":")
                start, end = positions.strip().split("-")
                bed_fh.write(f"{chrom}\t{start}\t{end}\n")
            except ValueError:
                logging.error(
                    f"Invalid region format: {region}. Expected format 'chr:start-end'."
                )
                raise ValueError(
                    f"Invalid region format: {region}. Expected format 'chr:start-end'."
                )


def run_pipeline(
    bwa_reference,
    output_dir,
    extra_modules,
    module_args,
    config,
    fastq1=None,
    fastq2=None,
    bam=None,
    cram=None,
    threads=4,
    reference_assembly="hg19",
    fast_mode=False,
    keep_intermediates=False,
    delete_intermediates=False,
    archive_results=False,
    archive_format="zip",
    custom_regions=None,
    bed_file=None,
    log_level=logging.INFO,
    sample_name=None,
    log_file=None,
    summary_formats=None,  # New parameter: list of additional summary output formats (e.g., ['csv', 'tsv'])
):
    """
    Main pipeline function that orchestrates the genotyping process.

    Args:
        bwa_reference (str): Path to the genome reference FASTA file for BWA.
        output_dir (Path): Path to the output directory.
        extra_modules (list): Optional modules to include (e.g., ['advntr', 'shark']).
        module_args (dict): Dictionary containing module-specific arguments.
        config (dict): Configuration dictionary.
        fastq1 (str, optional): Path to the first FASTQ file.
        fastq2 (str, optional): Path to the second FASTQ file.
        bam (str, optional): Path to the BAM file.
        cram (str, optional): Path to the CRAM file.
        threads (int, optional): Number of threads to use. Default is 4.
        reference_assembly (str, optional): Reference assembly ("hg19" or "hg38").
        fast_mode (bool, optional): Skip filtering steps if True.
        keep_intermediates (bool, optional): Keep intermediate files.
        delete_intermediates (bool, optional): Delete intermediate files after processing.
        archive_results (bool, optional): Archive results after completion.
        archive_format (str, optional): Format for archiving (zip or tar.gz).
        custom_regions (str, optional): Comma-separated custom regions.
        bed_file (Path, optional): BED file for MUC1 analysis.
        log_level (int, optional): Logging level.
        sample_name (str, optional): Sample name for labeling results.
        log_file (str, optional): Path to the log file.
        summary_formats (list, optional): Additional summary output formats to generate.
            Supported formats are 'csv' and 'tsv'. JSON summary is always generated.

    Raises:
        ValueError: Various input validation errors.
        FileNotFoundError: If the specified BED file is not found.
        RuntimeError: If alignment fails due to missing indexes.
    """
    if not bwa_reference:
        logging.error("BWA reference not provided or determined from configuration.")
        raise ValueError("BWA reference not provided or determined from configuration.")

    logging.debug(f"BWA reference set to: {bwa_reference}")
    logging.debug(f"Output directory set to: {output_dir}")

    dirs = create_output_directories(output_dir)
    logging.info(f"Created output directories in: {output_dir}")

    tool_versions = get_tool_versions(config)
    logging.info(
        f"VNtyper pipeline {VERSION} started with tool versions: {tool_versions}"
    )

    overall_start = timeit.default_timer()
    logging.info("Pipeline execution started.")

    # Initialize summary to record pipeline steps, including vntyper version and empty input_files.
    summary = start_summary(version=VERSION, input_files={})
    # Define the summary file path to be used for writing after each step.
    summary_file_path = os.path.join(output_dir, "pipeline_summary.json")

    input_type = None
    if fastq1 and fastq2:
        input_type = "FASTQ"
    elif bam:
        input_type = "BAM"
    elif cram:
        input_type = "CRAM"
    else:
        logging.error("No input files provided.")
        raise ValueError("No input files provided.")

    input_files = {}
    if input_type == "FASTQ":
        input_files["fastq1"] = os.path.basename(fastq1)
        input_files["fastq2"] = os.path.basename(fastq2)
    elif input_type == "BAM":
        input_files["bam"] = os.path.basename(bam)
    elif input_type == "CRAM":
        input_files["cram"] = os.path.basename(cram)
    # Update the summary with input file information
    summary["input_files"] = input_files

    try:
        input_count = sum(
            [
                1 if input_type == "FASTQ" else 0,
                1 if input_type == "BAM" else 0,
                1 if input_type == "CRAM" else 0,
            ]
        )
        if input_count > 1:
            logging.error(
                "Multiple input types provided. Provide only one: FASTQ, BAM, or CRAM."
            )
            raise ValueError("Provide either BAM, CRAM, or FASTQ files, not multiples.")

        if not bam and not cram and (not fastq1 or not fastq2):
            logging.error(
                "When not providing BAM/CRAM, both --fastq1 and --fastq2 must be specified for paired-end sequencing."
            )
            raise ValueError(
                "When not providing BAM/CRAM, both --fastq1 and --fastq2 must be specified for paired-end sequencing."
            )

        # Validate input files
        if input_type == "BAM":
            validate_bam_file(bam)
        elif input_type == "CRAM":
            validate_bam_file(cram)
        elif input_type == "FASTQ":
            validate_fastq_file(fastq1)
            validate_fastq_file(fastq2)
        else:
            logging.error("Incomplete FASTQ inputs provided.")
            raise ValueError(
                "Both FASTQ files must be provided for paired-end sequencing."
            )  # BED file logic
        if bed_file:
            bed_file_path = Path(bed_file)
            if not bed_file_path.exists():
                logging.error(f"Provided BED file does not exist: {bed_file_path}")
                raise FileNotFoundError(f"BED file not found: {bed_file_path}")
            logging.info(f"Using provided BED file: {bed_file_path}")
        elif custom_regions:
            bed_file_path = Path(output_dir) / "custom_regions.bed"
            write_bed_file(custom_regions, bed_file_path)
            logging.info(f"Custom regions converted to BED file: {bed_file_path}")
        else:
            region_key = f"bam_region_{reference_assembly}"
            predefined_regions = config.get("bam_processing", {}).get(region_key)
            if not predefined_regions:
                logging.error(
                    f"Region key '{region_key}' not found in config.json under 'bam_processing'."
                )
                raise ValueError(f"Missing configuration for region: {region_key}")
            bed_file_path = (
                Path(output_dir) / f"predefined_regions_{reference_assembly}.bed"
            )
            write_bed_file(predefined_regions, bed_file_path)
            logging.info(f"Predefined regions converted to BED file: {bed_file_path}")

        logging.debug(f"Final bed_file_path => {bed_file_path}")
        logging.debug(f"bed_file_path exists? {bed_file_path.exists()}")

        # --- Input Conversion ---
        if input_type in ["BAM", "CRAM"]:
            logging.info(
                f"Starting {input_type} to FASTQ conversion with specified regions."
            )
            conversion_start = datetime.utcnow()
            if input_type == "BAM":
                if bam is None or str(bam).strip().lower() == "none":
                    logging.error("Invalid BAM input (None).")
                    raise ValueError("Invalid BAM file input.")

                fastq1, fastq2, _, _ = process_bam_to_fastq(
                    in_bam=bam,
                    output=dirs["fastq_bam_processing"],
                    output_name="output",
                    threads=threads,
                    config=config,
                    reference_assembly=reference_assembly,
                    fast_mode=fast_mode,
                    delete_intermediates=delete_intermediates,
                    keep_intermediates=keep_intermediates,
                    bed_file=bed_file_path,
                )
                conversion_command = f"process_bam_to_fastq(in_bam={bam}, ...)"
                header_parse_start = datetime.utcnow()
                header = extract_bam_header(bam, config)
                parse_header_pipeline_info(header, Path(dirs["fastq_bam_processing"]))
                header_parse_end = datetime.utcnow()
                record_step(
                    summary,
                    "BAM Header Parsing",
                    str(Path(dirs["fastq_bam_processing"]) / "pipeline_info.json"),
                    "json",
                    "parse_header_pipeline_info(extracted header)",
                    header_parse_start,
                    header_parse_end,
                    write_summary_path=summary_file_path,
                )
            else:  # CRAM branch
                if cram is None or str(cram).strip().lower() == "none":
                    logging.error("Invalid CRAM input (None).")
                    raise ValueError("Invalid CRAM file input.")

                fastq1, fastq2, _, _ = process_bam_to_fastq(
                    in_bam=cram,
                    output=dirs["fastq_bam_processing"],
                    output_name="output",
                    threads=threads,
                    config=config,
                    reference_assembly=reference_assembly,
                    fast_mode=fast_mode,
                    delete_intermediates=delete_intermediates,
                    keep_intermediates=keep_intermediates,
                    bed_file=bed_file_path,
                    file_format="cram",
                )
                conversion_command = (
                    f"process_bam_to_fastq(in_bam={cram}, file_format='cram', ...)"
                )
            conversion_end = datetime.utcnow()
            record_step(
                summary,
                f"{input_type} to FASTQ Conversion",
                str(Path(dirs["fastq_bam_processing"]) / "output_R1.fastq.gz"),
                "fastq",
                conversion_command,
                conversion_start,
                conversion_end,
                write_summary_path=summary_file_path,
            )
            if not fastq1 or not fastq2:
                logging.error("Failed to generate FASTQ files from input. Exiting.")
                raise ValueError("Failed to generate FASTQ files from input.")

        elif input_type == "FASTQ":
            # --- SHARK Filtering Module ---
            if "shark" in extra_modules:
                from vntyper.modules.shark.shark_filtering import (
                    run_shark_filter,
                    load_shark_config,
                )

                shark_config = load_shark_config()
                logging.info("SHARK module included. Running SHARK filtering first.")
                run_sample_name = sample_name if sample_name else "sample"
                shark_start = datetime.utcnow()
                fastq1, fastq2 = run_shark_filter(
                    fastq_1=fastq1,
                    fastq_2=fastq2,
                    output_dir=dirs["fastq_bam_processing"],
                    config=shark_config,
                    main_config=config,
                    sample_name=run_sample_name,
                    reference_assembly=reference_assembly,
                    threads=threads,
                )
                shark_end = datetime.utcnow()
                record_step(
                    summary,
                    "SHARK Filtering",
                    os.path.join(dirs["fastq_bam_processing"], "filtered_R1.fastq.gz"),
                    "fastq",
                    "run_shark_filter(...)",
                    shark_start,
                    shark_end,
                    write_summary_path=summary_file_path,
                )
            logging.info("Starting FASTQ quality control.")
            qc_start = datetime.utcnow()
            process_fastq(
                fastq1,
                fastq2,
                threads,
                dirs["fastq_bam_processing"],
                "output",
                config,
            )
            qc_end = datetime.utcnow()
            record_step(
                summary,
                "FASTQ Quality Control",
                os.path.join(dirs["fastq_bam_processing"], "output.json"),
                "json",
                "process_fastq(...)",
                qc_start,
                qc_end,
                write_summary_path=summary_file_path,
            )
            logging.info("FASTQ quality control completed.")
            fastq1 = os.path.join(dirs["fastq_bam_processing"], "output_R1.fastq.gz")
            fastq2 = os.path.join(dirs["fastq_bam_processing"], "output_R2.fastq.gz")
            logging.info("Starting FASTQ alignment.")
            align_start = datetime.utcnow()
            sorted_bam = align_and_sort_fastq(
                fastq1,
                fastq2,
                bwa_reference,
                dirs["alignment_processing"],
                "output",
                threads,
                config,
            )
            align_end = datetime.utcnow()
            record_step(
                summary,
                "FASTQ Alignment",
                sorted_bam,
                "bam",
                "align_and_sort_fastq(...)",
                align_start,
                align_end,
                write_summary_path=summary_file_path,
            )
            if not sorted_bam:
                logging.error(
                    "Alignment failed: BWA index files for the provided reference "
                    "are missing or incomplete. Please run 'bwa index <reference.fa>' "
                    "to generate them."
                )
                raise RuntimeError(
                    "Alignment failed due to missing or incomplete BWA reference indices."
                )
            logging.info("FASTQ alignment completed.")
            logging.info("Starting BAM to FASTQ conversion (Post-alignment).")
            conv2_start = datetime.utcnow()
            fastq1, fastq2, _, _ = process_bam_to_fastq(
                in_bam=sorted_bam,
                output=dirs["fastq_bam_processing"],
                output_name="output",
                threads=threads,
                config=config,
                reference_assembly=reference_assembly,
                fast_mode=fast_mode,
                delete_intermediates=delete_intermediates,
                keep_intermediates=keep_intermediates,
                bed_file=bed_file_path,
            )
            conv2_end = datetime.utcnow()
            record_step(
                summary,
                "BAM to FASTQ Conversion (Post-alignment)",
                os.path.join(dirs["fastq_bam_processing"], "output_R1.fastq.gz"),
                "fastq",
                "process_bam_to_fastq(sorted_bam, ...)",
                conv2_start,
                conv2_end,
                write_summary_path=summary_file_path,
            )
            if not fastq1 or not fastq2:
                logging.error("Failed to generate FASTQ files from BAM. Exiting.")
                raise ValueError(
                    "Failed to generate FASTQ files from BAM."
                )  # --- Coverage Calculation ---
        logging.info("Calculating mean coverage over the VNTR region.")
        if input_type == "BAM":
            input_bam = Path(bam)
        elif input_type == "CRAM":
            input_bam = Path(cram)
        else:
            input_bam = Path(dirs["alignment_processing"]) / "output_sorted.bam"

        vntr_region_key = f"vntr_region_{reference_assembly}"
        vntr_region = config.get("bam_processing", {}).get(vntr_region_key)
        if not vntr_region:
            logging.error(
                f"Region key '{vntr_region_key}' not found in config.json under 'bam_processing'."
            )
            raise ValueError(
                f"Missing configuration for VNTR region: {vntr_region_key}"
            )

        cov_start = datetime.utcnow()
        coverage_stats = calculate_vntr_coverage(
            bam_file=str(input_bam),
            region=vntr_region,
            threads=threads,
            config=config,
            output_dir=dirs["coverage"],
            output_name="coverage",
        )
        cov_end = datetime.utcnow()
        mean_coverage = coverage_stats["mean"]
        record_step(
            summary,
            "Coverage Calculation",
            str(Path(dirs["coverage"]) / "coverage_summary.tsv"),
            "tsv",
            "calculate_vntr_coverage(...)",
            cov_start,
            cov_end,
            write_summary_path=summary_file_path,
        )

        # --- Kestrel Genotyping ---
        logging.info("Starting Kestrel genotyping.")
        vcf_out = os.path.join(dirs["kestrel"], "output.vcf")
        kestrel_path = config["tools"]["kestrel"]
        reference_vntr = config["reference_data"]["muc1_reference_vntr"]

        kestrel_start = datetime.utcnow()
        if fastq1 and fastq2:
            run_kestrel(
                vcf_path=Path(vcf_out),
                output_dir=Path(dirs["kestrel"]),
                fastq_1=fastq1,
                fastq_2=fastq2,
                reference_vntr=reference_vntr,
                kestrel_path=kestrel_path,
                config=config,
                sample_name=sample_name,
                log_level=log_level,
            )
        else:
            logging.error("FASTQ files required for Kestrel genotyping not provided.")
            raise ValueError(
                "FASTQ files required for Kestrel genotyping not provided."
            )
        kestrel_end = datetime.utcnow()
        record_step(
            summary,
            "Kestrel Genotyping",
            os.path.join(dirs["kestrel"], "kestrel_result.tsv"),
            "tsv",
            "run_kestrel(...)",
            kestrel_start,
            kestrel_end,
            write_summary_path=summary_file_path,
        )
        logging.info(
            "Kestrel genotyping completed."
        )  # --- adVNTR Genotyping and Cross-Match (only if advntr requested and performed) ---
        if "advntr" in extra_modules:
            logging.info("adVNTR module included. Starting adVNTR genotyping.")
            try:
                from vntyper.modules.advntr.advntr_genotyping import (
                    load_advntr_config,
                    process_advntr_output,
                    run_advntr,
                )
            except ImportError as exc:
                logging.error(f"adVNTR module import failed: {exc}")
                sys.exit(1)

            advntr_config = load_advntr_config()
            advntr_settings = advntr_config.get("advntr_settings", {})
            advntr_reference = module_args.get("advntr", {}).get("advntr_reference")

            if not advntr_reference:
                ref_map = {
                    "hg19": "hg19",
                    "GRCh37": "hg19",
                    "hg38": "hg38",
                    "GRCh38": "hg38",
                }
                ucsc_style_ref = ref_map.get(reference_assembly, "hg19")
                advntr_key = f"advntr_reference_vntr_{ucsc_style_ref}"
                advntr_reference = config.get("reference_data", {}).get(advntr_key)
            else:
                if advntr_reference == "hg19":
                    advntr_reference = config.get("reference_data", {}).get(
                        "advntr_reference_vntr_hg19"
                    )
                elif advntr_reference == "hg38":
                    advntr_reference = config.get("reference_data", {}).get(
                        "advntr_reference_vntr_hg38"
                    )
                else:
                    logging.error(f"Invalid advntr_reference: {advntr_reference}")
                    raise ValueError(f"Invalid advntr_reference: {advntr_reference}")

            if not advntr_reference:
                logging.error("adVNTR reference path not found in configuration.")
                raise ValueError("adVNTR reference path not found in configuration.")

            logging.debug(f"adVNTR reference set to: {advntr_reference}")

            max_cov = module_args.get("advntr", {}).get("max_coverage")
            sorted_bam = Path(dirs["fastq_bam_processing"]) / "output_sliced.bam"
            if sorted_bam and sorted_bam.exists():
                if max_cov:
                    logging.info(
                        f"Using quick adVNTR mode with max coverage = {max_cov}"
                    )
                    sorted_bam = downsample_bam_if_needed(
                        bam_path=sorted_bam,
                        max_coverage=max_cov,
                        reference_assembly=reference_assembly,
                        threads=threads,
                        config=config,
                        coverage_dir=dirs["coverage"],
                        coverage_prefix="advntr_precheck",
                    )
                advntr_start = datetime.utcnow()
                run_advntr(
                    advntr_reference,
                    sorted_bam,
                    dirs["advntr"],
                    "output",
                    config=config,
                )
                output_format = advntr_settings.get("output_format", "tsv")
                output_ext = ".vcf" if output_format == "vcf" else ".tsv"
                output_path = os.path.join(dirs["advntr"], f"output_adVNTR{output_ext}")
                process_advntr_output(
                    output_path, dirs["advntr"], "output", config=config
                )
                advntr_end = datetime.utcnow()
                record_step(
                    summary,
                    "adVNTR Genotyping",
                    os.path.join(dirs["advntr"], "output_adVNTR_result.tsv"),
                    "tsv",
                    "run_advntr(...), process_advntr_output(...)",
                    advntr_start,
                    advntr_end,
                    write_summary_path=summary_file_path,
                )
                logging.info("adVNTR genotyping completed.")

                # --- Cross-Match Variant Comparison ---
                logging.info(
                    "Starting cross-match of Kestrel and adVNTR variant calls."
                )
                cross_start = datetime.utcnow()
                # Extract parsed results from the summary
                kestrel_records, advntr_records = extract_results_from_pipeline_summary(
                    summary
                )
                if not kestrel_records:
                    logging.error(
                        "Kestrel genotyping results not found for cross-match."
                    )
                    raise ValueError(
                        "Kestrel genotyping results not found for cross-match."
                    )
                if not advntr_records:
                    logging.error(
                        "adVNTR genotyping results not found for cross-match."
                    )
                    raise ValueError(
                        "adVNTR genotyping results not found for cross-match."
                    )
                crossmatch_summary = cross_match_variants(
                    kestrel_records, advntr_records, config=config
                )
                cross_match_output = os.path.join(
                    dirs["advntr"], "cross_match_results.tsv"
                )
                write_results_tsv(crossmatch_summary["matches"], cross_match_output)
                cross_end = datetime.utcnow()
                record_step(
                    summary,
                    "Cross-Match Variant Comparison",
                    cross_match_output,
                    "tsv",
                    "cross_match_variants(kestrel_results, advntr_results)",
                    cross_start,
                    cross_end,
                    write_summary_path=summary_file_path,
                )
                logging.info(
                    f"Cross-match variant comparison completed. Overall match: {crossmatch_summary['overall_match']}"
                )
            else:
                logging.error("Sorted BAM required for adVNTR not provided.")
                raise ValueError("Sorted BAM required for adVNTR not provided.")
        else:
            logging.info("adVNTR module not included. Skipping adVNTR genotyping.")

        # --- Generate Summary Report and Archiving ---
        logging.info("Generating summary report.")
        report_file = "summary_report.html"
        template_dir = config.get("paths", {}).get("template_dir", "vntyper/templates")

        sorted_vcf = os.path.join(dirs["kestrel"], "output_indel.vcf.gz")
        bam_out = os.path.join(dirs["kestrel"], "output.bam")
        bed_out = os.path.join(dirs["kestrel"], "output.bed")
        fasta_reference = config["reference_data"]["muc1_reference_vntr"]

        generate_summary_report(
            output_dir,
            template_dir,
            report_file,
            log_file,
            bed_file=bed_out,
            bam_file=bam_out,
            fasta_file=fasta_reference,
            flanking=config.get("default_values", {}).get("flanking", 50),
            vcf_file=sorted_vcf,
            config=config,
        )
        logging.info(f"Summary report generated: {report_file}")

        if archive_results:
            logging.info("Archiving the results folder.")
            if archive_format == "zip":
                fmt = "zip"
            elif archive_format == "tar.gz":
                fmt = "gztar"
            else:
                logging.error(f"Unsupported archive format: {archive_format}")
                raise ValueError(f"Unsupported archive format: {archive_format}")

            archive_name = f"{output_dir}"
            try:
                archive_path = shutil.make_archive(
                    base_name=archive_name,
                    format=fmt,
                    root_dir=output_dir,
                    base_dir=".",
                )
                logging.info(f"Results folder archived at: {archive_path}")
            except Exception as exc:
                logging.error(f"Failed to archive results folder: {exc}")

        logging.info("Pipeline finished successfully.")

        # Mark pipeline end in summary
        end_summary(summary)

        # Write out the complete pipeline summary
        summary_file_path = os.path.join(output_dir, "pipeline_summary.json")
        write_summary(summary, summary_file_path)
        logging.info(f"Pipeline summary written to: {summary_file_path}")

        # Generate additional summary output formats if specified
        if summary_formats:
            if "csv" in summary_formats:
                csv_path = os.path.join(output_dir, "pipeline_summary.csv")
                convert_summary_to_csv(summary, csv_path)
                logging.info(f"Pipeline summary CSV written to: {csv_path}")
            if "tsv" in summary_formats:
                tsv_path = os.path.join(output_dir, "pipeline_summary.tsv")
                convert_summary_to_tsv(summary, tsv_path)
                logging.info(f"Pipeline summary TSV written to: {tsv_path}")

    except Exception as exc:
        logging.error(f"An error occurred: {exc}", exc_info=True)
        sys.exit(1)

    overall_stop = timeit.default_timer()
    elapsed_time = (overall_stop - overall_start) / 60
    logging.info(f"Pipeline completed in {elapsed_time:.2f} minutes.")
