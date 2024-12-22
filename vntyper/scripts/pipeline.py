#!/usr/bin/env python3
# vntyper/scripts/pipeline.py

import logging
import os
import shutil
import sys
import timeit
from pathlib import Path

from vntyper.scripts.alignment_processing import align_and_sort_fastq
from vntyper.scripts.fastq_bam_processing import (
    process_bam_to_fastq,
    process_fastq,
    calculate_vntr_coverage
)
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.utils import (
    create_output_directories,
    get_tool_versions,
    setup_logging,
    validate_bam_file,
    validate_fastq_file
)
from vntyper.version import __version__ as VERSION


def write_bed_file(regions, bed_file_path):
    """
    Writes regions to a BED file in the correct format.

    Parameters:
    - regions (str): Comma-separated regions in 'chr:start-end' format.
    - bed_file_path (Path): Path to the BED file to be written.
    """
    with open(bed_file_path, 'w', encoding='utf-8') as bed_fh:
        for region in regions.split(','):
            try:
                chrom, positions = region.strip().split(':')
                start, end = positions.strip().split('-')
                bed_fh.write(f"{chrom}\t{start}\t{end}\n")
            except ValueError:
                logging.error(
                    f"Invalid region format: {region}. "
                    "Expected format 'chr:start-end'."
                )
                raise ValueError(
                    f"Invalid region format: {region}. "
                    "Expected format 'chr:start-end'."
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
    archive_format='zip',
    custom_regions=None,
    bed_file=None,
    log_level=logging.INFO,
    sample_name=None
):
    """
    Main pipeline function that orchestrates the genotyping process.

    Args:
        bwa_reference (str): Path to the genome reference FASTA file for BWA.
        output_dir (Path): Path to the output directory.
        extra_modules (list): Optional modules to include (e.g., ['advntr']).
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

    Raises:
        ValueError: Various input validation errors.
        FileNotFoundError: If the specified BED file is not found.
        RuntimeError: If alignment fails due to missing indexes.
    """
    if not bwa_reference:
        logging.error(
            "BWA reference not provided or determined from configuration."
        )
        raise ValueError(
            "BWA reference not provided or determined from configuration."
        )

    logging.debug(f"BWA reference set to: {bwa_reference}")
    logging.debug(f"Output directory set to: {output_dir}")

    dirs = create_output_directories(output_dir)
    logging.info(f"Created output directories in: {output_dir}")

    log_file = os.path.join(output_dir, "pipeline.log")
    setup_logging(log_level=log_level, log_file=log_file)
    logging.info(f"Logging to file: {log_file}")

    tool_versions = get_tool_versions(config)
    logging.info(
        f"VNtyper pipeline {VERSION} started with tool versions: {tool_versions}"
    )

    start_time = timeit.default_timer()
    logging.info("Pipeline execution started.")

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
        input_files['fastq1'] = os.path.basename(fastq1)
        input_files['fastq2'] = os.path.basename(fastq2)
    elif input_type == "BAM":
        input_files['bam'] = os.path.basename(bam)
    elif input_type == "CRAM":
        input_files['cram'] = os.path.basename(cram)

    try:
        input_count = sum([
            1 if input_type == "FASTQ" else 0,
            1 if input_type == "BAM" else 0,
            1 if input_type == "CRAM" else 0
        ])
        if input_count > 1:
            logging.error(
                "Multiple input types provided. Provide only one: FASTQ, BAM, or CRAM."
            )
            raise ValueError(
                "Provide either BAM, CRAM, or FASTQ files, not multiples."
            )

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
            )

        if bed_file:
            bed_file_path = Path(bed_file)
            if not bed_file_path.exists():
                logging.error(
                    f"Provided BED file does not exist: {bed_file_path}"
                )
                raise FileNotFoundError(
                    f"BED file not found: {bed_file_path}"
                )
            logging.info(f"Using provided BED file: {bed_file_path}")
        elif custom_regions:
            bed_file_path = Path(output_dir) / "custom_regions.bed"
            write_bed_file(custom_regions, bed_file_path)
            logging.info(f"Custom regions converted to BED file: {bed_file_path}")
        else:
            if reference_assembly == "hg38":
                predefined_regions = config["bam_processing"]["bam_region_hg38"]
            else:
                predefined_regions = config["bam_processing"]["bam_region_hg19"]

            bed_file_path = Path(output_dir) / f"predefined_regions_{reference_assembly}.bed"
            write_bed_file(predefined_regions, bed_file_path)
            logging.info(f"Predefined regions converted to BED file: {bed_file_path}")

        if input_type == "BAM":
            logging.info("Starting BAM to FASTQ conversion with specified regions.")
            if bam is None or str(bam).strip().lower() == "none":
                logging.error("Invalid BAM input (None).")
                raise ValueError("Invalid BAM file input.")

            fastq1, fastq2, _, _ = process_bam_to_fastq(
                in_bam=bam,
                output=dirs['fastq_bam_processing'],
                output_name="output",
                threads=threads,
                config=config,
                reference_assembly=reference_assembly,
                fast_mode=fast_mode,
                delete_intermediates=delete_intermediates,
                keep_intermediates=keep_intermediates,
                bed_file=bed_file_path
            )
            if not fastq1 or not fastq2:
                logging.error("Failed to generate FASTQ files from BAM. Exiting.")
                raise ValueError("Failed to generate FASTQ files from BAM.")

        elif input_type == "CRAM":
            logging.info("Starting CRAM to FASTQ conversion with specified regions.")
            if cram is None or str(cram).strip().lower() == "none":
                logging.error("Invalid CRAM input (None).")
                raise ValueError("Invalid CRAM file input.")

            fastq1, fastq2, _, _ = process_bam_to_fastq(
                in_bam=cram,
                output=dirs['fastq_bam_processing'],
                output_name="output",
                threads=threads,
                config=config,
                reference_assembly=reference_assembly,
                fast_mode=fast_mode,
                delete_intermediates=delete_intermediates,
                keep_intermediates=keep_intermediates,
                bed_file=bed_file_path,
                file_format="cram"
            )
            if not fastq1 or not fastq2:
                logging.error("Failed to generate FASTQ files from CRAM. Exiting.")
                raise ValueError("Failed to generate FASTQ files from CRAM.")

        elif input_type == "FASTQ":
            if 'shark' in extra_modules:
                from vntyper.modules.shark.shark_filtering import (
                    run_shark_filter, load_shark_config
                )
                shark_config = load_shark_config()
                logging.info("SHARK module included. Running SHARK filtering first.")
                run_sample_name = sample_name if sample_name else "sample"
                fastq1, fastq2 = run_shark_filter(
                    fastq_1=fastq1,
                    fastq_2=fastq2,
                    output_dir=dirs['fastq_bam_processing'],
                    config=shark_config,
                    main_config=config,
                    sample_name=run_sample_name,
                    reference_assembly=reference_assembly,
                    threads=threads
                )

            logging.info("Starting FASTQ quality control.")
            process_fastq(
                fastq1,
                fastq2,
                threads,
                dirs['fastq_bam_processing'],
                "output",
                config,
            )
            logging.info("FASTQ quality control completed.")

            fastq1 = os.path.join(dirs['fastq_bam_processing'], "output_R1.fastq.gz")
            fastq2 = os.path.join(dirs['fastq_bam_processing'], "output_R2.fastq.gz")

            logging.info("Starting FASTQ alignment.")
            sorted_bam = align_and_sort_fastq(
                fastq1,
                fastq2,
                bwa_reference,
                dirs['alignment_processing'],
                "output",
                threads,
                config,
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

            logging.info("Starting BAM to FASTQ conversion with specified regions.")
            if sorted_bam is None or str(sorted_bam).lower() == "none":
                logging.error("Alignment produced a None value for sorted_bam.")
                raise ValueError("No valid sorted_bam file provided.")

            fastq1, fastq2, _, _ = process_bam_to_fastq(
                in_bam=sorted_bam,
                output=dirs['fastq_bam_processing'],
                output_name="output",
                threads=threads,
                config=config,
                reference_assembly=reference_assembly,
                fast_mode=fast_mode,
                delete_intermediates=delete_intermediates,
                keep_intermediates=keep_intermediates,
                bed_file=bed_file_path
            )
            if not fastq1 or not fastq2:
                logging.error("Failed to generate FASTQ files from BAM. Exiting.")
                raise ValueError("Failed to generate FASTQ files from BAM.")

        logging.info("Calculating mean coverage over the VNTR region.")
        if input_type == "BAM":
            input_bam = Path(bam)
        elif input_type == "CRAM":
            input_bam = Path(cram)
        else:
            input_bam = Path(dirs['alignment_processing']) / "output_sorted.bam"

        if reference_assembly == "hg38":
            vntr_region = config["bam_processing"]["vntr_region_hg38"]
        else:
            vntr_region = config["bam_processing"]["vntr_region_hg19"]

        mean_coverage = calculate_vntr_coverage(
            bam_file=str(input_bam),
            region=vntr_region,
            threads=threads,
            config=config,
            output_dir=dirs['coverage'],
            output_name="coverage"
        )

        logging.info("Starting Kestrel genotyping.")
        vcf_out = os.path.join(dirs['kestrel'], "output.vcf")
        kestrel_path = config["tools"]["kestrel"]
        reference_vntr = config["reference_data"]["muc1_reference_vntr"]

        if fastq1 and fastq2:
            run_kestrel(
                vcf_path=Path(vcf_out),
                output_dir=Path(dirs['kestrel']),
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
            raise ValueError("FASTQ files required for Kestrel genotyping not provided.")

        logging.info("Kestrel genotyping completed.")

        if 'advntr' in extra_modules:
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
            advntr_reference = module_args.get('advntr', {}).get('advntr_reference')

            if not advntr_reference:
                if reference_assembly == "hg19":
                    advntr_reference = config.get("reference_data", {}).get(
                        "advntr_reference_vntr_hg19"
                    )
                else:
                    advntr_reference = config.get("reference_data", {}).get(
                        "advntr_reference_vntr_hg38"
                    )
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
                    raise ValueError(
                        f"Invalid advntr_reference: {advntr_reference}"
                    )

            if not advntr_reference:
                logging.error("adVNTR reference path not found in configuration.")
                raise ValueError("adVNTR reference path not found in configuration.")

            logging.debug(f"adVNTR reference set to: {advntr_reference}")
            sorted_bam = Path(dirs['fastq_bam_processing']) / "output_sliced.bam"
            if sorted_bam and sorted_bam.exists():
                run_advntr(
                    advntr_reference,
                    sorted_bam,
                    dirs['advntr'],
                    "output",
                    config,
                )

                output_format = advntr_settings.get("output_format", "tsv")
                output_ext = ".vcf" if output_format == "vcf" else ".tsv"
                output_path = os.path.join(
                    dirs['advntr'], f"output_adVNTR{output_ext}"
                )
                process_advntr_output(output_path, dirs['advntr'], "output")
                logging.info("adVNTR genotyping completed.")
            else:
                logging.error("Sorted BAM required for adVNTR not provided.")
                raise ValueError("Sorted BAM required for adVNTR not provided.")
        else:
            logging.info("adVNTR module not included. Skipping adVNTR genotyping.")

        logging.info("Generating summary report.")
        report_file = "summary_report.html"
        template_dir = config.get(
            'paths', {}
        ).get('template_dir', 'vntyper/templates')

        sorted_vcf = os.path.join(dirs['kestrel'], "output_indel.vcf.gz")
        bam_out = os.path.join(dirs['kestrel'], "output.bam")
        bed_out = os.path.join(dirs['kestrel'], "output.bed")
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
            input_files=input_files,
            pipeline_version=VERSION,
            mean_vntr_coverage=mean_coverage,
            vcf_file=sorted_vcf,
            config=config
        )
        logging.info(f"Summary report generated: {report_file}")

        if archive_results:
            logging.info("Archiving the results folder.")
            if archive_format == 'zip':
                fmt = 'zip'
            elif archive_format == 'tar.gz':
                fmt = 'gztar'
            else:
                logging.error(f"Unsupported archive format: {archive_format}")
                raise ValueError(f"Unsupported archive format: {archive_format}")

            archive_name = f"{output_dir}"
            try:
                archive_path = shutil.make_archive(
                    base_name=archive_name,
                    format=fmt,
                    root_dir=output_dir,
                    base_dir='.',
                )
                logging.info(f"Results folder archived at: {archive_path}")
            except Exception as exc:
                logging.error(f"Failed to archive results folder: {exc}")

        logging.info("Pipeline finished successfully.")

    except Exception as exc:
        logging.error(f"An error occurred: {exc}", exc_info=True)
        sys.exit(1)

    stop_time = timeit.default_timer()
    elapsed_time = (stop_time - start_time) / 60
    logging.info(f"Pipeline completed in {elapsed_time:.2f} minutes.")
