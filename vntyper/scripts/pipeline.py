import logging
import os
import sys
import timeit
from datetime import datetime, timezone
from pathlib import Path

from vntyper.scripts.alignment_preflight import run_preflight
from vntyper.scripts.alignment_processing import align_and_sort_fastq
from vntyper.scripts.archive_safety import create_safe_archive
from vntyper.scripts.artifact_names import select_best_vcf_file

# Import cross-match functions from cross_match.py
from vntyper.scripts.cross_match import (
    cross_match_variants,
    extract_results_from_pipeline_summary,
    write_results_tsv,
)
from vntyper.scripts.fastq_bam_processing import (
    calculate_vntr_coverage,
    downsample_bam_if_needed,
    extract_bam_header,
    parse_header_pipeline_info,
    process_bam_to_fastq,
    process_fastq,
)
from vntyper.scripts.generate_report import generate_summary_report
from vntyper.scripts.kestrel_genotyping import run_kestrel
from vntyper.scripts.pipeline_alignment import (
    build_alignment_preflight_kwargs,
    prepare_alignment_target,
    prepare_input_alignment_preflight,
    resolve_summary_reference_provenance,
)
from vntyper.scripts.pipeline_cleanup import close_alignment_plan
from vntyper.scripts.pipeline_coverage import calculate_alignment_coverage
from vntyper.scripts.pipeline_inputs import archive_base_name, protect_pipeline_input_ownership, resolve_pipeline_input
from vntyper.scripts.pipeline_kestrel import run_kestrel_stage
from vntyper.scripts.pipeline_read_routing import route_converted_fastqs
from vntyper.scripts.reference_resolution import resolve_from_mapping
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution as pin_reference_resolution
from vntyper.scripts.reference_resolution_environment import restore_reference_resolution
from vntyper.scripts.region_utils import get_region_string_with_fallback

# Import our new summary functions (including end_summary and CSV/TSV conversion functions)
from vntyper.scripts.summary import (
    convert_summary_to_csv,
    convert_summary_to_tsv,
    end_summary,
    record_step,
    start_summary,
    write_summary,
)

# The five step names are matched by exact string comparison in generate_report.py,
# cohort_summary.py and cross_match.py. A typo does not fail - it silently drops a
# report section (AGENTS.md trap 5), so they are named, never spelled out.
from vntyper.scripts.summary_steps import (
    STEP_ADVNTR,
    STEP_BAM_HEADER,
    STEP_COVERAGE,
    STEP_CROSS_MATCH,
)
from vntyper.scripts.utils import (
    create_output_directories,
    get_tool_versions,
    validate_bam_file,
    validate_fastq_file,
)
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


def select_advntr_reference(config: dict, reference_assembly: str) -> str | None:
    """The adVNTR database for an assembly, by coordinate system.

    Args:
        config: Pipeline configuration.
        reference_assembly: Supported assembly label.

    Returns:
        str | None: Database path, or None when no key is configured.
    """
    resolved = resolve_from_mapping("advntr", reference_assembly, config.get("reference_data", {}))
    return resolved.value if resolved is not None else None


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
    reference_fasta=None,
    threads=4,
    reference_assembly="hg19",
    reference_key_used=None,
    reference_source_effective=None,
    fast_mode=False,
    keep_intermediates=False,
    delete_intermediates=False,
    archive_results=False,
    archive_format="zip",
    custom_regions=None,
    bed_file=None,
    log_level=logging.INFO,
    sample_name=None,
    sample_name_is_explicit=False,
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
        reference_fasta (Path, optional): Explicit reference FASTA for CRAM decoding.
        threads (int, optional): Number of threads to use. Default is 4.
        reference_assembly (str, optional): Reference assembly ("hg19" or "hg38").
        reference_key_used (str, optional): The `reference_data` config key that
            supplied `bwa_reference`, as resolved by
            :func:`vntyper.scripts.cli_handlers.select_bwa_reference`. Recorded in the
            run summary so a UCSC-family fallback is visible in the report.
        reference_source_effective (str, optional): The reference source
            ("ucsc", "ncbi" or "ensembl") the run actually used, which can differ from
            `reference_assembly`'s own source when a fallback was taken. Recorded in
            the run summary alongside `reference_key_used`.
        fast_mode (bool, optional): Skip filtering steps if True.
        keep_intermediates (bool, optional): Keep intermediate files.
        delete_intermediates (bool, optional): Delete intermediate files after processing.
        archive_results (bool, optional): Archive results after completion.
        archive_format (str, optional): Format for archiving (zip or tar.gz).
        custom_regions (str, optional): Comma-separated custom regions.
        bed_file (Path, optional): BED file for MUC1 analysis.
        log_level (int, optional): Logging level.
        sample_name (str, optional): Sample name for labeling results.
        sample_name_is_explicit (bool, optional): Whether `sample_name` is the
            operator's own `--sample-name` rather than a name `handle_pipeline`
            derived from an input path. Recorded in the run summary beside the name
            itself, because the report has to finish deriving a derived one
            (`S1_R1.fastq` -> `S1`) and must leave an explicit one alone - and the
            string on its own cannot say which it is (#242).
        log_file (str, optional): Path to the log file.
        summary_formats (list, optional): Additional summary output formats to generate.
            Supported formats are 'csv' and 'tsv'. JSON summary is always generated.

    Raises:
        ValueError: Various input validation errors.
        FileNotFoundError: If the specified BED file is not found.
        RuntimeError: If alignment fails due to missing indexes.
    """
    # Capture the working directory at the start to ensure it remains valid
    # throughout the pipeline execution, especially for tools like Java that need it
    try:
        project_root = os.getcwd()
        logger.debug(f"Captured project root directory: {project_root}")
    except (OSError, FileNotFoundError) as e:
        # If we can't get the current directory, try to use the absolute path of the script
        project_root = str(Path(__file__).parent.parent.parent)
        logger.warning(f"Could not determine current working directory ({e}), using fallback: {project_root}")

    logger.debug(f"BWA reference set to: {bwa_reference}")
    logger.debug(f"Output directory set to: {output_dir}")

    overall_start = timeit.default_timer()
    logger.info("Pipeline execution started.")

    input_type, input_files = resolve_pipeline_input(fastq1, fastq2, bam, cram, bwa_reference, extra_modules)
    archive_protected_paths = tuple(
        path for path in (bam, cram, fastq1, fastq2, reference_fasta, bed_file, bwa_reference) if path
    )
    previous_ref_path = None
    reference_resolution_pinned = False
    alignment_plan = None
    primary_outcome_is_active = False
    try:
        protect_pipeline_input_ownership(
            output_dir,
            input_type,
            fastq1,
            fastq2,
            bam,
            cram,
            bed_file,
            reference_fasta,
            bwa_reference,
            config,
            reference_assembly,
            archive_results,
            archive_format,
        )
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        if input_type == "FASTQ":
            validate_fastq_file(fastq1)
            if fastq2:
                validate_fastq_file(fastq2)
        elif input_type not in {"BAM", "CRAM"}:
            logger.error("No supported input was provided.")
            raise ValueError("No supported input was provided.")

        alignment_header = None
        vntr_region = None
        if input_type in ["BAM", "CRAM"]:
            input_alignment = bam if input_type == "BAM" else cram
            prepared = prepare_input_alignment_preflight(
                in_path=str(input_alignment),
                input_type=input_type,
                output_dir=output_dir,
                config=config,
                threads=threads,
                reference_assembly=reference_assembly,
                bed_file=bed_file,
                custom_regions=custom_regions,
                reference_fasta=reference_fasta,
                fast_mode=fast_mode,
                alignment_validator=validate_bam_file,
                validation_cwd=project_root,
            )
            alignment_header = prepared.alignment_header
            bed_file_path = prepared.bed_file
            vntr_region = prepared.coverage_region
            alignment_plan = prepared.plan
            previous_ref_path = prepared.previous_ref_path
            reference_resolution_pinned = input_type == "CRAM"
        else:
            bed_file_path = prepare_alignment_target(
                input_type=input_type,
                bam=bam,
                cram=cram,
                output_dir=output_dir,
                reference_assembly=reference_assembly,
                config=config,
                bed_file=bed_file,
                custom_regions=custom_regions,
            )

        dirs = create_output_directories(output_dir)
        logger.info(f"Created output directories in: {output_dir}")

        tool_versions = get_tool_versions(config)
        logger.info(f"VNtyper pipeline {VERSION} started with tool versions: {tool_versions}")

        # What the run actually used, not what BWA was configured with (MAJOR 5,
        # milestone-5 PR-2 review): a BAM run never reads a reference, and a CRAM run
        # decodes against whatever `alignment_plan` resolved above, which can differ
        # entirely from the configured BWA path. Only FASTQ's own BWA resolution is
        # correct as recorded, so it passes through unchanged.
        reference_provenance = resolve_summary_reference_provenance(
            input_type=input_type,
            bwa_reference_key=reference_key_used,
            bwa_reference_path=bwa_reference,
            bwa_reference_source=reference_source_effective,
            alignment_plan=alignment_plan,
        )
        summary = start_summary(
            version=VERSION,
            input_files=input_files,
            # The same string Kestrel embeds below, so the report and the VCF cannot
            # name the same run differently (#242). `start_summary` runs before any
            # step, so this is on disk from the first `record_step` onwards. The
            # provenance flag travels with it because the string alone cannot say
            # whether it is a name or a `Path.stem` still to be derived.
            sample_name=sample_name,
            sample_name_is_explicit=sample_name_is_explicit,
            reference_assembly_requested=reference_assembly,
            reference_key_used=reference_provenance.key_used,
            reference_path=reference_provenance.path,
            reference_source_effective=reference_provenance.source_effective,
        )
        summary_file_path = os.path.join(output_dir, "pipeline_summary.json")
        if input_type in ["BAM", "CRAM"]:
            logger.info(f"Starting {input_type} to FASTQ conversion with specified regions.")
            conversion_start = datetime.now(timezone.utc).replace(tzinfo=None)
            if input_type == "BAM":
                if bam is None or str(bam).strip().lower() == "none":
                    logger.error("Invalid BAM input (None).")
                    raise ValueError("Invalid BAM file input.")
            else:  # CRAM branch
                if cram is None or str(cram).strip().lower() == "none":
                    logger.error("Invalid CRAM input (None).")
                    raise ValueError("Invalid CRAM file input.")

            produced_fastqs = process_bam_to_fastq(
                output=dirs["fastq_bam_processing"],
                output_name="output",
                threads=threads,
                config=config,
                plan=alignment_plan,
                reference_assembly=reference_assembly,
                fast_mode=fast_mode,
                delete_intermediates=delete_intermediates,
                keep_intermediates=keep_intermediates,
                bed_file=bed_file_path,
            )
            kestrel_fastq_files = route_converted_fastqs(produced_fastqs, config)
            if not kestrel_fastq_files:
                raise ValueError("FASTQ routing produced no inputs for Kestrel.")
            conversion_command = "process_bam_to_fastq(plan=alignment_plan, ...)"

            if input_type == "BAM":
                header_parse_start = datetime.now(timezone.utc).replace(tzinfo=None)
                # Re-read a failed guard from the same proven alignment view.
                header = alignment_header or extract_bam_header(alignment_plan.view_path, config)
                parse_header_pipeline_info(header, Path(dirs["fastq_bam_processing"]), config)
                header_parse_end = datetime.now(timezone.utc).replace(tzinfo=None)
                record_step(
                    summary,
                    STEP_BAM_HEADER,
                    str(Path(dirs["fastq_bam_processing"]) / "pipeline_info.json"),
                    "json",
                    "parse_header_pipeline_info(extracted header)",
                    header_parse_start,
                    header_parse_end,
                    write_summary_path=summary_file_path,
                )
            conversion_end = datetime.now(timezone.utc).replace(tzinfo=None)
            record_step(
                summary,
                f"{input_type} to FASTQ Conversion",
                str(kestrel_fastq_files[0]),
                "fastq",
                conversion_command,
                conversion_start,
                conversion_end,
                write_summary_path=summary_file_path,
            )

        elif input_type == "FASTQ":
            # --- SHARK Filtering Module ---
            if "shark" in extra_modules:
                from vntyper.modules.shark.shark_filtering import (
                    load_shark_config,
                    run_shark_filter,
                )

                shark_config = load_shark_config()
                logger.info("SHARK module included. Running SHARK filtering first.")
                run_sample_name = sample_name or "sample"
                shark_start = datetime.now(timezone.utc).replace(tzinfo=None)
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
                shark_end = datetime.now(timezone.utc).replace(tzinfo=None)
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
            logger.info("Starting FASTQ quality control.")
            qc_start = datetime.now(timezone.utc).replace(tzinfo=None)
            paired_fastq_input = fastq2 is not None
            process_fastq(
                fastq1,
                fastq2,
                threads,
                dirs["fastq_bam_processing"],
                "output",
                config,
            )
            qc_end = datetime.now(timezone.utc).replace(tzinfo=None)
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
            logger.info("FASTQ quality control completed.")
            fastq1 = os.path.join(dirs["fastq_bam_processing"], "output_R1.fastq.gz")
            fastq2 = os.path.join(dirs["fastq_bam_processing"], "output_R2.fastq.gz") if paired_fastq_input else None
            logger.info("Starting FASTQ alignment.")
            align_start = datetime.now(timezone.utc).replace(tzinfo=None)
            sorted_bam = align_and_sort_fastq(
                fastq1,
                fastq2,
                bwa_reference,
                dirs["alignment_processing"],
                "output",
                threads,
                config,
            )
            align_end = datetime.now(timezone.utc).replace(tzinfo=None)
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
                logger.error(
                    "Alignment failed: BWA index files for the provided reference "
                    "are missing or incomplete. Please run 'bwa index <reference.fa>' "
                    "to generate them."
                )
                raise RuntimeError("Alignment failed due to missing or incomplete BWA reference indices.")
            logger.info("FASTQ alignment completed.")
            alignment_plan = run_preflight(
                **build_alignment_preflight_kwargs(
                    in_path=str(sorted_bam),
                    output_dir=Path(output_dir) / "fastq_bam_processing",
                    output_name="post_alignment",
                    file_format="bam",
                    config=config,
                    threads=threads,
                    bed_file=bed_file_path,
                    reference_assembly=reference_assembly,
                    fast_mode=fast_mode,
                )
            )
            logger.info("Starting BAM to FASTQ conversion (Post-alignment).")
            conv2_start = datetime.now(timezone.utc).replace(tzinfo=None)
            produced_fastqs = process_bam_to_fastq(
                output=dirs["fastq_bam_processing"],
                output_name="output",
                threads=threads,
                config=config,
                plan=alignment_plan,
                reference_assembly=reference_assembly,
                fast_mode=fast_mode,
                delete_intermediates=delete_intermediates,
                keep_intermediates=keep_intermediates,
                bed_file=bed_file_path,
            )
            kestrel_fastq_files = route_converted_fastqs(produced_fastqs, config)
            if not kestrel_fastq_files:
                raise ValueError("FASTQ routing produced no inputs for Kestrel.")
            conv2_end = datetime.now(timezone.utc).replace(tzinfo=None)
            record_step(
                summary,
                "BAM to FASTQ Conversion (Post-alignment)",
                str(kestrel_fastq_files[0]),
                "fastq",
                "process_bam_to_fastq(plan=post_alignment_plan, ...)",
                conv2_start,
                conv2_end,
                write_summary_path=summary_file_path,
            )
            # --- Coverage Calculation ---
        logger.info("Calculating mean coverage over the VNTR region.")
        if alignment_plan is None:
            raise RuntimeError("Alignment preflight did not produce a plan for coverage.")
        cov_start = datetime.now(timezone.utc).replace(tzinfo=None)
        vntr_region = calculate_alignment_coverage(
            plan=alignment_plan,
            region=vntr_region,
            reference_assembly=reference_assembly,
            threads=threads,
            config=config,
            output_dir=dirs["coverage"],
            coverage_calculator=calculate_vntr_coverage,
            region_resolver=get_region_string_with_fallback,
        )
        # The exact span the coverage stage consumed - resolved here and, until
        # #242, thrown away. The report could not otherwise state it: reading
        # `config["default_values"]["reference_assembly"]` back would mislabel any
        # `--reference-assembly` override and cannot reconstruct `--custom-regions`
        # at all. `record_step` below serialises the whole summary dict, not just
        # the step, so setting it here puts it on disk immediately.
        summary["region_resolved"] = vntr_region
        cov_end = datetime.now(timezone.utc).replace(tzinfo=None)
        record_step(
            summary,
            STEP_COVERAGE,
            str(Path(dirs["coverage"]) / "coverage_summary.tsv"),
            "tsv",
            "calculate_vntr_coverage(...)",
            cov_start,
            cov_end,
            write_summary_path=summary_file_path,
        )

        logger.info("Starting Kestrel genotyping.")
        run_kestrel_stage(
            fastq_files=kestrel_fastq_files,
            dirs=dirs,
            config=config,
            sample_name=sample_name or "sample",
            log_level=log_level,
            cwd=project_root,
            summary=summary,
            summary_file_path=summary_file_path,
            runner=run_kestrel,
        )
        logger.info(
            "Kestrel genotyping completed."
        )  # --- adVNTR Genotyping and Cross-Match (only if advntr requested and performed) ---
        if "advntr" in extra_modules:
            logger.info("adVNTR module included. Starting adVNTR genotyping.")
            try:
                from vntyper.modules.advntr.advntr_genotyping import (
                    advntr_output_extension,
                    load_advntr_config,
                    process_advntr_output,
                    run_advntr,
                )
            except ImportError as exc:
                logger.error(f"adVNTR module import failed: {exc}")
                sys.exit(1)

            advntr_config = load_advntr_config()
            advntr_settings = advntr_config.get("advntr_settings", {})
            advntr_reference = module_args.get("advntr", {}).get("advntr_reference")

            if not advntr_reference:
                advntr_reference = select_advntr_reference(config, reference_assembly)
            else:
                if advntr_reference == "hg19":
                    advntr_reference = config.get("reference_data", {}).get("advntr_reference_vntr_hg19")
                elif advntr_reference == "hg38":
                    advntr_reference = config.get("reference_data", {}).get("advntr_reference_vntr_hg38")
                else:
                    logger.error(f"Invalid advntr_reference: {advntr_reference}")
                    raise ValueError(f"Invalid advntr_reference: {advntr_reference}")

            if not advntr_reference:
                logger.error("adVNTR reference path not found in configuration.")
                raise ValueError("adVNTR reference path not found in configuration.")

            logger.debug(f"adVNTR reference set to: {advntr_reference}")

            max_cov = module_args.get("advntr", {}).get("max_coverage")
            sorted_bam = Path(dirs["fastq_bam_processing"]) / "output_sliced.bam"
            if sorted_bam and sorted_bam.exists():
                if max_cov:
                    logger.info(f"Using quick adVNTR mode with max coverage = {max_cov}")
                    sorted_bam = downsample_bam_if_needed(
                        bam_path=sorted_bam,
                        max_coverage=max_cov,
                        reference_assembly=reference_assembly,
                        threads=threads,
                        config=config,
                        coverage_dir=dirs["coverage"],
                        coverage_prefix="advntr_precheck",
                    )
                advntr_start = datetime.now(timezone.utc).replace(tzinfo=None)
                run_advntr(
                    advntr_reference,
                    sorted_bam,
                    dirs["advntr"],
                    "output",
                    config=config,
                    cwd=project_root,
                    pipeline_threads=threads,
                )
                # Shared with run_advntr, which built the path adVNTR actually wrote. This
                # used to repeat `.get("output_format", "tsv")` here, so the producer and the
                # consumer of that path each carried their own fallback (#247).
                output_ext = advntr_output_extension(advntr_settings)
                output_path = os.path.join(dirs["advntr"], f"output_adVNTR{output_ext}")
                process_advntr_output(output_path, dirs["advntr"], "output", config=config)
                advntr_end = datetime.now(timezone.utc).replace(tzinfo=None)
                record_step(
                    summary,
                    STEP_ADVNTR,
                    os.path.join(dirs["advntr"], "output_adVNTR_result.tsv"),
                    "tsv",
                    "run_advntr(...), process_advntr_output(...)",
                    advntr_start,
                    advntr_end,
                    write_summary_path=summary_file_path,
                )
                logger.info("adVNTR genotyping completed.")

                # --- Cross-Match Variant Comparison ---
                logger.info("Starting cross-match of Kestrel and adVNTR variant calls.")
                cross_start = datetime.now(timezone.utc).replace(tzinfo=None)
                # Extract parsed results from the summary
                kestrel_records, advntr_records = extract_results_from_pipeline_summary(summary)
                if not kestrel_records:
                    logger.error("Kestrel genotyping results not found for cross-match.")
                    raise ValueError("Kestrel genotyping results not found for cross-match.")
                if not advntr_records:
                    logger.error("adVNTR genotyping results not found for cross-match.")
                    raise ValueError("adVNTR genotyping results not found for cross-match.")
                crossmatch_summary = cross_match_variants(kestrel_records, advntr_records, config=config)
                cross_match_output = os.path.join(dirs["advntr"], "cross_match_results.tsv")
                write_results_tsv(crossmatch_summary["matches"], cross_match_output)
                cross_end = datetime.now(timezone.utc).replace(tzinfo=None)
                record_step(
                    summary,
                    STEP_CROSS_MATCH,
                    cross_match_output,
                    "tsv",
                    "cross_match_variants(kestrel_results, advntr_results)",
                    cross_start,
                    cross_end,
                    write_summary_path=summary_file_path,
                )
                logger.info(
                    f"Cross-match variant comparison completed. Overall match: {crossmatch_summary['overall_match']}"
                )
            else:
                logger.error("Sorted BAM required for adVNTR not provided.")
                raise ValueError("Sorted BAM required for adVNTR not provided.")
        else:
            logger.info("adVNTR module not included. Skipping adVNTR genotyping.")

        # --- Generate Summary Report and Archiving ---
        logger.info("Generating summary report.")
        report_file = "summary_report.html"
        template_dir = config.get("paths", {}).get("template_dir", "vntyper/templates")

        # Select best available VCF file (compressed preferred, uncompressed fallback)
        vcf_file = select_best_vcf_file(dirs["kestrel"])
        bam_out = os.path.join(dirs["kestrel"], "output.bam")
        bed_out = os.path.join(dirs["kestrel"], "output.bed")
        fasta_reference = config["reference_data"]["muc1_reference_vntr"]

        # `generate_summary_report` reads `pipeline_summary.json` back **from
        # disk**, and the final `write_summary` below runs after it. Every
        # top-level key set since the last `record_step` would therefore be
        # missing from the report even though the finished file carries it. This
        # makes the file match the summary in hand at the moment it is read,
        # rather than relying on a later `record_step` happening to fire (#242).
        write_summary(summary, summary_file_path)

        generate_summary_report(
            output_dir,
            template_dir,
            report_file,
            log_file,
            bed_file=bed_out,
            bam_file=bam_out,
            fasta_file=fasta_reference,
            flanking=config.get("default_values", {}).get("flanking", 50),
            vcf_file=vcf_file,  # Can be .gz, .vcf, or None - handled gracefully
            config=config,
        )
        logger.info(f"Summary report generated: {report_file}")

        # Mark pipeline end in summary
        end_summary(summary)

        # Write out the complete pipeline summary
        summary_file_path = os.path.join(output_dir, "pipeline_summary.json")
        write_summary(summary, summary_file_path)
        logger.info(f"Pipeline summary written to: {summary_file_path}")

        # Generate additional summary output formats if specified
        if summary_formats:
            if "csv" in summary_formats:
                csv_path = os.path.join(output_dir, "pipeline_summary.csv")
                convert_summary_to_csv(summary, csv_path)
                logger.info(f"Pipeline summary CSV written to: {csv_path}")
            if "tsv" in summary_formats:
                tsv_path = os.path.join(output_dir, "pipeline_summary.tsv")
                convert_summary_to_tsv(summary, tsv_path)
                logger.info(f"Pipeline summary TSV written to: {tsv_path}")

        if archive_results:
            logger.info("Archiving the results folder.")
            formats = {"zip": "zip", "tar.gz": "gztar"}
            if archive_format not in formats:
                raise ValueError(f"Unsupported archive format: {archive_format}")
            archive_path = create_safe_archive(
                archive_base_name(output_dir),
                formats[archive_format],
                output_dir,
                protected_paths=archive_protected_paths,
            )
            logger.info(f"Results folder archived at: {archive_path}")

        logger.info("Pipeline finished successfully.")

    except Exception:
        primary_outcome_is_active = True
        logger.exception("An error occurred")
        sys.exit(1)
    except BaseException:
        primary_outcome_is_active = True
        raise
    finally:
        try:
            close_alignment_plan(alignment_plan, preserve_primary=primary_outcome_is_active)
        finally:
            if reference_resolution_pinned:
                restore_reference_resolution(previous_ref_path)

    overall_stop = timeit.default_timer()
    elapsed_time = (overall_stop - overall_start) / 60
    logger.info(f"Pipeline completed in {elapsed_time:.2f} minutes.")
