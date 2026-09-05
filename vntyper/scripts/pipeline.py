import contextlib
import logging
import os
import shutil
import sys
import tempfile
import timeit
from collections.abc import Mapping
from datetime import datetime, timezone
from functools import partial
from pathlib import Path
from typing import Any

from vntyper.scripts.alignment_preflight import run_preflight
from vntyper.scripts.alignment_processing import align_and_sort_fastq
from vntyper.scripts.archive_safety import create_safe_archive
from vntyper.scripts.artifact_names import (
    ADVNTR_EVIDENCE_SNAPSHOT_RELATIVE,
    DECISION_PROFILE_SNAPSHOT_RELATIVE,
    select_best_vcf_file,
)
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

# Import cross-match functions from cross_match.py
from vntyper.scripts.nomenclature_annotate import DominanceSeamOutcome, reconcile_caller_outputs
from vntyper.scripts.pipeline_advntr_cleanup import (
    plan_advntr_cleanup,
    validate_pipeline_log_outside_advntr_preflight,
    validate_pipeline_log_outside_selected_advntr_model,
)
from vntyper.scripts.pipeline_advntr_preflight import plan_advntr_preflight, plan_valid_advntr_preflight
from vntyper.scripts.pipeline_advntr_preflight import select_advntr_reference as select_advntr_reference
from vntyper.scripts.pipeline_advntr_run_context import AdvntrRunContext, prepare_advntr_run_context
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
from vntyper.scripts.pipeline_resume_planning import (
    build_analysis_settings,
    build_canonical_inputs_and_fingerprints,
    evaluate_resume_compatibility,
    initial_stage_carry_forward,
    record_reused_stage,
    resolve_effective_advntr_runtime,
    resolve_effective_kestrel_runtime,
    resolve_effective_shark_runtime,
)
from vntyper.scripts.profile_provenance import snapshot_decision_profile
from vntyper.scripts.reference_resolution_environment import pin_reference_resolution as pin_reference_resolution
from vntyper.scripts.reference_resolution_environment import restore_reference_resolution
from vntyper.scripts.region_utils import get_region_string_with_fallback
from vntyper.scripts.report_assets import DEFAULT_REPORT_IGV
from vntyper.scripts.resume import (
    fingerprint_file,
    load_prior_summary,
    resolve_reused_artifact_path,
    resume_refusals,
    step_is_reusable,
)
from vntyper.scripts.run_configuration import RunConfiguration, resolve_run_configuration

# Import our new summary functions (including end_summary and CSV/TSV conversion functions)
from vntyper.scripts.summary import (
    convert_summary_rows_to_delimited,
    convert_summary_to_csv,
    convert_summary_to_tsv,
    end_summary,
    record_step,
    refresh_step,
    start_summary,
    write_summary,
)

# The five step names are matched by exact string comparison in generate_report.py,
# cohort_summary.py and cross_match.py. A typo does not fail - it silently drops a
# report section (AGENTS.md trap 5), so they are named, never spelled out.
from vntyper.scripts.summary_steps import (
    STEP_ADVNTR,
    STEP_BAM_HEADER,
    STEP_BAM_TO_FASTQ,
    STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
    STEP_COVERAGE,
    STEP_CRAM_TO_FASTQ,
    STEP_CROSS_MATCH,
    STEP_FASTQ_ALIGNMENT,
    STEP_FASTQ_QC,
    STEP_KESTREL,
    STEP_SHARK,
)
from vntyper.scripts.utils import (
    create_output_directories,
    get_tool_versions,
    run_command,
    validate_bam_file,
    validate_fastq_file,
)
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


def _record_reused_step(summary: dict[str, Any], record: dict[str, Any]) -> None:
    st_name = record.get("step")
    idx = next((i for i, s in enumerate(summary["steps"]) if s.get("step") == st_name), None)
    if idx is not None:
        summary["steps"][idx] = record
    else:
        summary["steps"].append(record)


def _run_and_record_fastq_qc(
    fastq1: str | Path,
    fastq2: str | Path | None,
    threads: int,
    dirs: Mapping[str, str],
    config: dict[str, Any],
    summary: dict[str, Any],
    summary_file_path: str,
    *,
    qc_only: bool = False,
) -> None:
    logger.info("Starting FASTQ quality control.")
    qc_start = datetime.now(timezone.utc).replace(tzinfo=None)
    if qc_only:
        with tempfile.TemporaryDirectory(dir=dirs["fastq_bam_processing"], prefix=".vntyper_qc_") as tmpdir:
            process_fastq(
                fastq1,
                fastq2,
                threads,
                tmpdir,
                "output",
                config,
            )
            for artifact in ("output.json", "output_fastp.html", "output_fastp.log"):
                src = Path(tmpdir) / artifact
                if src.is_file():
                    shutil.copy2(src, Path(dirs["fastq_bam_processing"]) / artifact)
    else:
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
        STEP_FASTQ_QC,
        os.path.join(dirs["fastq_bam_processing"], "output.json"),
        "json",
        "process_fastq(...)",
        qc_start,
        qc_end,
        write_summary_path=summary_file_path,
    )
    logger.info("FASTQ quality control completed.")


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
    report_igv=DEFAULT_REPORT_IGV,
    run_configuration=None,
    resume=False,
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
        report_igv (str, optional): How the report carries its alignment browser -
            one of `report_assets.REPORT_IGV_MODES`. Threaded through from
            `--report-igv` on the `pipeline` subcommand, because an ordinary run
            reaches `generate_summary_report` here rather than through
            `vntyper report` (#242).
        run_configuration: Immutable decision profile and stage components resolved
            before the run. Direct compatibility callers load the packaged profile.

    Raises:
        ValueError: Various input validation errors.
        FileNotFoundError: If the specified BED file is not found.
        RuntimeError: If alignment fails due to missing indexes.
    """
    if run_configuration is None:
        run_configuration = resolve_run_configuration()
    elif not isinstance(run_configuration, RunConfiguration):
        raise ValueError("pipeline run_configuration must be a resolved RunConfiguration")

    if log_file is not None:
        early_advntr_preflight = plan_valid_advntr_preflight(
            config,
            extra_modules,
            module_args,
            reference_assembly,
        )
        if early_advntr_preflight is not None:
            validate_pipeline_log_outside_selected_advntr_model(log_file, early_advntr_preflight.reference)

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
        advntr_preflight = plan_advntr_preflight(config, extra_modules, module_args, reference_assembly)
        needs_advntr = advntr_preflight.enabled
        advntr_reference = advntr_preflight.reference
        additional_operator_paths = (advntr_reference,) if advntr_reference is not None else ()
        canonical_input_files, input_fingerprints = build_canonical_inputs_and_fingerprints(
            input_type, fastq1, fastq2, bam, cram, bed_file
        )

        analysis_settings = build_analysis_settings(
            reference_assembly=reference_assembly,
            fast_mode=fast_mode,
            custom_regions=custom_regions,
            bed_file=bed_file,
            advntr_reference=advntr_reference,
            module_args=module_args,
            config=config,
            extra_modules=extra_modules,
            input_type=input_type,
        )

        effective_reference_path = None
        if input_type == "FASTQ" and bwa_reference:
            effective_reference_path = str(Path(os.path.join(project_root, bwa_reference)).resolve())
        elif input_type == "CRAM" and reference_fasta:
            effective_reference_path = str(Path(os.path.join(project_root, reference_fasta)).resolve())

        effective_reference_fingerprint = (
            fingerprint_file(Path(effective_reference_path))
            if effective_reference_path and Path(effective_reference_path).is_file()
            else None
        )

        shark_reference_path = None
        shark_reference_fingerprint = None
        if input_type == "FASTQ" and "shark" in extra_modules:
            try:
                from vntyper.modules.shark.shark_filtering import select_muc1_region_fasta

                raw_shark = select_muc1_region_fasta(
                    dict(run_configuration.shark_runtime),
                    config,
                    reference_assembly,
                )
                if raw_shark:
                    shark_p = Path(os.path.join(project_root, raw_shark)).resolve()
                    shark_reference_path = str(shark_p)
                    if shark_p.is_file():
                        shark_reference_fingerprint = fingerprint_file(shark_p)
            except (ValueError, KeyError, OSError) as err:
                logger.debug("Could not resolve SHARK reference for fingerprinting: %s", err)

        raw_muc1 = config.get("reference_data", {}).get("muc1_reference_vntr")
        kestrel_reference_path = str(Path(os.path.join(project_root, raw_muc1)).resolve()) if raw_muc1 else None
        kestrel_reference_fingerprint = (
            fingerprint_file(Path(kestrel_reference_path))
            if kestrel_reference_path and Path(kestrel_reference_path).is_file()
            else None
        )

        raw_motifs = config.get("reference_data", {}).get("muc1_motifs_rev_com")
        kestrel_motifs_path = str(Path(os.path.join(project_root, raw_motifs)).resolve()) if raw_motifs else None
        kestrel_motifs_fingerprint = (
            fingerprint_file(Path(kestrel_motifs_path))
            if kestrel_motifs_path and Path(kestrel_motifs_path).is_file()
            else None
        )

        raw_advntr_rus = config.get("reference_data", {}).get("code_adVNTR_RUs")
        advntr_rus_path = str(Path(os.path.join(project_root, raw_advntr_rus)).resolve()) if raw_advntr_rus else None
        advntr_rus_fingerprint = (
            fingerprint_file(Path(advntr_rus_path)) if advntr_rus_path and Path(advntr_rus_path).is_file() else None
        )
        (
            kestrel_counting_mode,
            effective_kestrel_runtime,
            kestrel_runtime_fingerprint,
        ) = resolve_effective_kestrel_runtime(run_configuration, config, project_root)
        (
            effective_advntr_runtime,
            advntr_runtime_fingerprint,
        ) = resolve_effective_advntr_runtime(run_configuration, config) if "advntr" in extra_modules else ({}, None)
        (
            effective_shark_runtime,
            shark_runtime_fingerprint,
        ) = (
            resolve_effective_shark_runtime(run_configuration, config)
            if input_type == "FASTQ" and "shark" in extra_modules
            else ({}, None)
        )

        summary_file_path = os.path.join(output_dir, "pipeline_summary.json")
        prior_summary = None
        if resume:
            prior_summary = load_prior_summary(summary_file_path)
            if prior_summary is None:
                msg = f"Cannot resume: prior summary not found or invalid at {summary_file_path}"
                logger.error(msg)
                raise ValueError(msg)

            refusals = resume_refusals(
                prior_summary,
                version=VERSION,
                input_files=input_files,
                sample_name=sample_name,
                reference_key_used=reference_key_used if input_type == "FASTQ" else None,
                decision_profile_sha256=run_configuration.decision_profile.digest,
                canonical_input_files=canonical_input_files,
                reference_assembly=reference_assembly,
                analysis_settings=analysis_settings,
                reference_path=effective_reference_path,
                input_fingerprints=input_fingerprints,
                reference_fingerprint=effective_reference_fingerprint,
                shark_reference_path=shark_reference_path,
                shark_reference_fingerprint=shark_reference_fingerprint,
            )

            if refusals:
                for refusal in refusals:
                    logger.error("Resume refused: %s", refusal)
                raise ValueError(f"Cannot resume pipeline: {'; '.join(refusals)}")
            logger.info("Resuming execution from prior summary (started %s)", prior_summary.get("pipeline_start"))

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
            additional_operator_paths,
        )

        if resume:
            # On resume, revoke published reports and export tables from prior run so
            # an early failure during re-execution does not leave stale public outputs.
            # This runs strictly after ownership validation to protect input files.
            for stale_name in (
                "summary_report.html",
                "pipeline_summary.csv",
                "pipeline_summary.tsv",
                "pipeline_summary_rows.csv",
                "pipeline_summary_rows.tsv",
            ):
                stale_path = Path(output_dir) / stale_name
                if stale_path.is_file():
                    stale_path.unlink()

        # Refuse an unknown/incompatible adVNTR before alignment preparation,
        # conversion, coverage, or Kestrel. Keep the classified startup answer in this
        # run rather than probing again after those expensive stages: failed answers
        # deliberately remain absent from the probe's reusable success cache.
        advntr_context: AdvntrRunContext | None = None
        advntr_version_overrides = {}
        if needs_advntr:
            logger.debug(f"adVNTR reference set to: {advntr_reference}")
            advntr_cleanup = plan_advntr_cleanup(
                output_dir,
                archive_results=archive_results,
                archive_format=archive_format,
            )
            validate_pipeline_log_outside_advntr_preflight(log_file, advntr_cleanup)
            advntr_context = prepare_advntr_run_context(
                output_dir,
                advntr_reference,
                config,
                archive_results=archive_results,
                archive_format=archive_format,
                protected_paths=(*archive_protected_paths, *additional_operator_paths),
                revoke_outputs=not resume,
                revoke_published=resume,
            )
            advntr_version_overrides["advntr"] = ".".join(str(part) for part in advntr_context.version)
            (
                effective_advntr_runtime,
                advntr_runtime_fingerprint,
            ) = resolve_effective_advntr_runtime(
                run_configuration,
                config,
                advntr_version=advntr_context.version,
            )
        out_path = Path(output_dir)
        if not resume and out_path.exists() and any(out_path.iterdir()):
            logger.warning("Output directory %s is non-empty; prior results may be overwritten.", output_dir)

        out_path.mkdir(parents=True, exist_ok=True)
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
                alignment_validator=partial(validate_bam_file, samtools_path=config["tools"]["samtools"]),
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

        snapshot_decision_profile(
            run_configuration.decision_profile,
            Path(output_dir) / DECISION_PROFILE_SNAPSHOT_RELATIVE,
        )

        advntr_evidence = None
        if needs_advntr:
            from vntyper.modules.advntr.artifact_evidence import (
                load_packaged_artifact_evidence,
                snapshot_artifact_evidence,
            )

            advntr_evidence = load_packaged_artifact_evidence()
            snapshot_artifact_evidence(
                advntr_evidence,
                Path(output_dir) / ADVNTR_EVIDENCE_SNAPSHOT_RELATIVE,
            )

        # Probing every configured tool shelled out to adVNTR (315 ms) and SHARK (36 ms)
        # on every Kestrel-only run for a value that is only logged. The set is derived
        # from the input type as well as the modules: fastp and BWA belong to the FASTQ
        # path and are never invoked for BAM or CRAM, which `extra_modules` cannot say.
        # A requested module is named only when the config declares a tool of that name,
        # so this cannot assert the existence of an entry it has not read (trap 2).
        # Whether anything will read `<name>_sliced.bam.bai`. The conversion stage is
        # given the answer rather than `extra_modules`, so it cannot grow further
        # dependencies on module state it has no business knowing.

        tools_in_use = {"samtools", "kestrel", "java_path"}
        if input_type == "FASTQ":
            tools_in_use |= {"fastp", "bwa"}
        tools_in_use |= {module for module in extra_modules if module in config.get("tools", {})}
        tool_versions = get_tool_versions(
            config,
            tools_in_use=tools_in_use,
            version_overrides=advntr_version_overrides,
        )
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
        if input_type == "CRAM" and reference_provenance.path:
            effective_reference_path = reference_provenance.path
            effective_reference_fingerprint = (
                fingerprint_file(Path(effective_reference_path)) if Path(effective_reference_path).is_file() else None
            )

        summary = start_summary(
            version=VERSION,
            input_files=input_files,
            canonical_input_files=canonical_input_files,
            input_fingerprints=input_fingerprints,
            analysis_settings=analysis_settings,
            kestrel_reference_path=kestrel_reference_path,
            kestrel_reference_fingerprint=kestrel_reference_fingerprint,
            kestrel_motifs_path=kestrel_motifs_path,
            kestrel_motifs_fingerprint=kestrel_motifs_fingerprint,
            advntr_rus_path=advntr_rus_path,
            advntr_rus_fingerprint=advntr_rus_fingerprint,
            kestrel_runtime_fingerprint=kestrel_runtime_fingerprint,
            advntr_runtime_fingerprint=advntr_runtime_fingerprint,
            reference_fingerprint=effective_reference_fingerprint,
            shark_reference_path=shark_reference_path,
            shark_reference_fingerprint=shark_reference_fingerprint,
            shark_runtime_fingerprint=shark_runtime_fingerprint,
            sample_name=sample_name,
            sample_name_is_explicit=sample_name_is_explicit,
            reference_assembly_requested=reference_assembly,
            reference_key_used=reference_provenance.key_used,
            reference_path=reference_provenance.path,
            persistent_reference_path=reference_provenance.path,
            reference_consumer_path=alignment_plan.reference_path if alignment_plan else None,
            reference_source_effective=reference_provenance.source_effective,
            advntr_evidence_digest=advntr_evidence.digest if advntr_evidence is not None else None,
            decision_profile=run_configuration.decision_profile,
        )

        compatibility = evaluate_resume_compatibility(
            prior_summary,
            input_type=input_type,
            kestrel_reference_path=kestrel_reference_path,
            kestrel_reference_fingerprint=kestrel_reference_fingerprint,
            kestrel_motifs_path=kestrel_motifs_path,
            kestrel_motifs_fingerprint=kestrel_motifs_fingerprint,
            kestrel_runtime_fingerprint=kestrel_runtime_fingerprint,
            kestrel_counting_mode=kestrel_counting_mode,
            advntr_model_sha=advntr_context.model.get("sha256") if advntr_context else None,
            advntr_rus_path=advntr_rus_path,
            advntr_rus_fingerprint=advntr_rus_fingerprint,
            advntr_runtime_fingerprint=advntr_runtime_fingerprint,
            shark_reference_path=shark_reference_path,
            shark_reference_fingerprint=shark_reference_fingerprint,
            shark_runtime_fingerprint=shark_runtime_fingerprint,
            effective_reference_path=effective_reference_path,
            effective_reference_fingerprint=effective_reference_fingerprint,
            advntr_version=advntr_version_overrides.get("advntr"),
            current_preprocessing_tools=analysis_settings.get("preprocessing_tools"),
        )

        if resume and prior_summary:
            donor_summary_path = Path(output_dir) / "pipeline_summary.donor.json"
            if not donor_summary_path.is_file() and os.path.isfile(summary_file_path):
                with contextlib.suppress(OSError):
                    shutil.copy2(summary_file_path, donor_summary_path)
            initial_stage_carry_forward(
                summary,
                prior_summary,
                output_dir,
                compatibility=compatibility,
                needs_advntr=needs_advntr,
                extra_modules=extra_modules,
            )
        else:
            donor_summary_path = Path(output_dir) / "pipeline_summary.donor.json"
            if donor_summary_path.is_file():
                with contextlib.suppress(OSError):
                    donor_summary_path.unlink()

        if input_type in ["BAM", "CRAM"]:
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

            conversion_step = STEP_BAM_TO_FASTQ if input_type == "BAM" else STEP_CRAM_TO_FASTQ
            can_reuse_conversion = (
                resume
                and prior_summary is not None
                and not (input_type == "CRAM" and compatibility.inval_cram)
                and not compatibility.inval_align
                and step_is_reusable(
                    prior_summary,
                    conversion_step,
                    output_dir,
                    needs_advntr=needs_advntr,
                )
            )
            fq_dir = Path(dirs["fastq_bam_processing"])
            candidate_fastqs = (
                str(fq_dir / "output_R1.fastq.gz"),
                str(fq_dir / "output_R2.fastq.gz"),
                str(fq_dir / "output_other.fastq.gz"),
                str(fq_dir / "output_single.fastq.gz"),
            )
            kestrel_fastq_files = None
            if can_reuse_conversion and all(os.path.isfile(p) for p in candidate_fastqs):
                try:
                    kestrel_fastq_files = route_converted_fastqs(candidate_fastqs, config)
                except (ValueError, OSError) as err:
                    logger.warning("Existing converted FASTQs failed routing on resume: %s", err)
                    kestrel_fastq_files = None

            if kestrel_fastq_files is not None:
                logger.info("Reusing previous %s step results.", conversion_step)
                record_reused_stage(summary, prior_summary, conversion_step, output_root=output_dir)
                write_summary(summary, summary_file_path)

            else:
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
                    bed_file=bed_file_path,
                    needs_advntr=needs_advntr,
                )
                kestrel_fastq_files = route_converted_fastqs(produced_fastqs, config)
                if not kestrel_fastq_files:
                    raise ValueError("FASTQ routing produced no inputs for Kestrel.")
                conversion_command = "process_bam_to_fastq(plan=alignment_plan, ...)"
                conversion_end = datetime.now(timezone.utc).replace(tzinfo=None)
                record_step(
                    summary,
                    STEP_BAM_TO_FASTQ if input_type == "BAM" else STEP_CRAM_TO_FASTQ,
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
                    run_shark_filter,
                    write_shark_step_summary,
                )

                logger.info("SHARK module included. Running SHARK filtering first.")
                run_sample_name = sample_name or "sample"
                shark_start = datetime.now(timezone.utc).replace(tzinfo=None)
                fastq1, fastq2 = run_shark_filter(
                    fastq_1=fastq1,
                    fastq_2=fastq2,
                    output_dir=dirs["fastq_bam_processing"],
                    config=run_configuration.shark_runtime,
                    main_config=config,
                    sample_name=run_sample_name,
                    reference_assembly=reference_assembly,
                    threads=threads,
                    resolved_component=run_configuration.shark,
                    custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
                )
                shark_step_file = os.path.join(dirs["fastq_bam_processing"], f"{run_sample_name}_shark_step.json")
                write_shark_step_summary(
                    fastq1,
                    fastq2,
                    shark_step_file,
                    config=run_configuration.shark_runtime,
                    shark_version=tool_versions.get("shark", "unknown"),
                )
                # Counting and sidecar creation are part of this stage. Capture its end
                # only after both succeed, so a recorded duration cannot exclude work
                # required to make the step readable.
                shark_end = datetime.now(timezone.utc).replace(tzinfo=None)
                record_step(
                    summary,
                    STEP_SHARK,
                    shark_step_file,
                    "json",
                    "run_shark_filter(...), write_shark_step_summary(...)",
                    shark_start,
                    shark_end,
                    write_summary_path=summary_file_path,
                )
            paired_fastq_input = fastq2 is not None
            can_reuse_alignment_and_conv = (
                resume
                and prior_summary is not None
                and not compatibility.inval_align
                and step_is_reusable(prior_summary, STEP_FASTQ_ALIGNMENT, output_dir)
                and step_is_reusable(
                    prior_summary,
                    STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                    output_dir,
                    needs_advntr=needs_advntr,
                )
            )
            fq_dir = Path(dirs["fastq_bam_processing"])
            candidate_post_fastqs = (
                str(fq_dir / "output_R1.fastq.gz"),
                str(fq_dir / "output_R2.fastq.gz"),
                str(fq_dir / "output_other.fastq.gz"),
                str(fq_dir / "output_single.fastq.gz"),
            )
            kestrel_fastq_files = None
            if can_reuse_alignment_and_conv and all(os.path.isfile(p) for p in candidate_post_fastqs):
                try:
                    kestrel_fastq_files = route_converted_fastqs(candidate_post_fastqs, config)
                except (ValueError, OSError) as err:
                    logger.warning("Existing post-alignment FASTQs failed routing: %s", err)
                    kestrel_fastq_files = None

            if kestrel_fastq_files is not None:
                logger.info("Reusing previous alignment and post-alignment conversion.")
                record_reused_stage(summary, prior_summary, STEP_FASTQ_ALIGNMENT, output_root=output_dir)
                record_reused_stage(summary, prior_summary, STEP_BAM_TO_FASTQ_POST_ALIGNMENT, output_root=output_dir)
                can_reuse_qc = (
                    resume
                    and prior_summary is not None
                    and not compatibility.inval_qc
                    and step_is_reusable(prior_summary, STEP_FASTQ_QC, output_dir)
                )
                if can_reuse_qc:
                    logger.info("Reusing previous %s step results.", STEP_FASTQ_QC)
                    record_reused_stage(summary, prior_summary, STEP_FASTQ_QC, output_root=output_dir)
                else:
                    _run_and_record_fastq_qc(
                        fastq1, fastq2, threads, dirs, config, summary, summary_file_path, qc_only=True
                    )
                write_summary(summary, summary_file_path)
                prior_align_st = next(
                    (s for s in prior_summary.get("steps", []) if s.get("step") == STEP_FASTQ_ALIGNMENT), None
                )
                candidate_bam = None
                if prior_align_st and prior_align_st.get("result_file"):
                    candidate_bam = resolve_reused_artifact_path(prior_align_st["result_file"], output_dir)
                sorted_bam = candidate_bam or (Path(dirs["alignment_processing"]) / "output_sorted.bam")
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
            else:
                can_reuse_alignment = (
                    resume
                    and prior_summary is not None
                    and not compatibility.inval_align
                    and step_is_reusable(prior_summary, STEP_FASTQ_ALIGNMENT, output_dir)
                )
                if can_reuse_alignment:
                    logger.info("Reusing previous %s step results.", STEP_FASTQ_ALIGNMENT)
                    record_reused_stage(summary, prior_summary, STEP_FASTQ_ALIGNMENT, output_root=output_dir)
                    can_reuse_qc = (
                        resume
                        and prior_summary is not None
                        and not compatibility.inval_qc
                        and step_is_reusable(prior_summary, STEP_FASTQ_QC, output_dir)
                    )
                    if can_reuse_qc:
                        logger.info("Reusing previous %s step results.", STEP_FASTQ_QC)
                        record_reused_stage(summary, prior_summary, STEP_FASTQ_QC, output_root=output_dir)
                    else:
                        _run_and_record_fastq_qc(
                            fastq1, fastq2, threads, dirs, config, summary, summary_file_path, qc_only=True
                        )
                    write_summary(summary, summary_file_path)
                    prior_align_st = next(
                        (s for s in prior_summary.get("steps", []) if s.get("step") == STEP_FASTQ_ALIGNMENT), None
                    )
                    candidate_bam = None
                    if prior_align_st and prior_align_st.get("result_file"):
                        candidate_bam = resolve_reused_artifact_path(prior_align_st["result_file"], output_dir)
                    sorted_bam = candidate_bam or (Path(dirs["alignment_processing"]) / "output_sorted.bam")
                else:
                    _run_and_record_fastq_qc(fastq1, fastq2, threads, dirs, config, summary, summary_file_path)

                    fastq1 = os.path.join(dirs["fastq_bam_processing"], "output_R1.fastq.gz")
                    fastq2 = (
                        os.path.join(dirs["fastq_bam_processing"], "output_R2.fastq.gz") if paired_fastq_input else None
                    )
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
                        STEP_FASTQ_ALIGNMENT,
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
                    bed_file=bed_file_path,
                    needs_advntr=needs_advntr,
                )
                kestrel_fastq_files = route_converted_fastqs(produced_fastqs, config)
                if not kestrel_fastq_files:
                    raise ValueError("FASTQ routing produced no inputs for Kestrel.")
                conv2_end = datetime.now(timezone.utc).replace(tzinfo=None)
                record_step(
                    summary,
                    STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
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

        if (
            resume
            and prior_summary
            and compatibility.kestrel_ref_matches
            and step_is_reusable(prior_summary, STEP_KESTREL, output_dir)
        ):
            logger.info("Reusing previous %s step results.", STEP_KESTREL)
            record_reused_stage(summary, prior_summary, STEP_KESTREL, output_root=output_dir)
            if "kestrel_counting_mode" in prior_summary:
                summary["kestrel_counting_mode"] = prior_summary["kestrel_counting_mode"]
            write_summary(summary, summary_file_path)

        else:
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
                threads=threads,
                resolved_component=run_configuration.kestrel,
                nomenclature_component=run_configuration.nomenclature,
                dominance_component=run_configuration.dominance,
                runtime_component=run_configuration.kestrel_runtime,
                custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
            )
            logger.info(
                "Kestrel genotyping completed."
            )  # --- adVNTR Genotyping and Cross-Match (only if advntr requested and performed) ---
        if "advntr" not in extra_modules and run_configuration.dominance.get("enabled") is True:
            dominance_outcome = reconcile_caller_outputs(
                os.path.join(dirs["kestrel"], "kestrel_result.tsv"),
                None,
                resolved_component=run_configuration.nomenclature,
                dominance_component=run_configuration.dominance,
                custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
            )
            if not isinstance(dominance_outcome, DominanceSeamOutcome) or not dominance_outcome.evaluated:
                raise RuntimeError("enabled dominance seam did not evaluate the selected run policy")
            if dominance_outcome.rewritten:
                refresh_step(summary, STEP_KESTREL, write_summary_path=summary_file_path)
        if "advntr" in extra_modules:
            logger.info("adVNTR module included. Starting adVNTR genotyping.")
            if advntr_context is None:
                raise RuntimeError("adVNTR execution reached without a verified run context.")
            # Which model a run resolved decides which reads adVNTR can ever see: the
            # fetch window comes from the model's own content. Validate before running,
            # because the failure mode is a confident result over a truncated locus
            # rather than an error (#268).
            # Top-level, not inside record_step: this is run state, and a step record is
            # not where state belongs.
            advntr_model = dict(advntr_context.model)
            summary["advntr_model"] = advntr_model
            logger.info(
                "adVNTR model %s (%s), window %s bp over %s",
                advntr_model["sha256"][:12],
                advntr_model["schema_version"],
                advntr_model["window_bp"],
                advntr_model["genomic_interval"],
            )

            if (
                resume
                and prior_summary
                and compatibility.advntr_model_matches
                and step_is_reusable(prior_summary, STEP_ADVNTR, output_dir)
            ):
                logger.info("Reusing previous %s step results.", STEP_ADVNTR)
                record_reused_stage(summary, prior_summary, STEP_ADVNTR, output_root=output_dir)
                write_summary(summary, summary_file_path)

            else:
                try:
                    from vntyper.modules.advntr.advntr_genotyping import (
                        advntr_output_extension,
                        process_advntr_output,
                        run_advntr,
                    )
                    from vntyper.modules.advntr.advntr_result_io import invalidate_advntr_artifact
                except ImportError as exc:
                    logger.error(f"adVNTR module import failed: {exc}")
                    sys.exit(1)

                for stage_rel in (
                    "output_adVNTR_result.tsv",
                    "output_adVNTR.tsv",
                    "output_adVNTR.vcf",
                    "cross_match_results.tsv",
                    "output_advntr.log",
                ):
                    invalidate_advntr_artifact(Path(dirs["advntr"]) / stage_rel)

                advntr_settings = run_configuration.advntr.get("settings")
                if not isinstance(advntr_settings, Mapping):
                    raise ValueError("resolved adVNTR settings must be a mapping")
                # Shared with run_advntr, which builds the path adVNTR writes. Resolve it for
                # this run's parser. Both supported producer names and the derived result
                # were already invalidated immediately after input-ownership validation.
                output_ext = advntr_output_extension(advntr_settings)

                max_cov = module_args.get("advntr", {}).get("max_coverage")
                sorted_bam = Path(dirs["fastq_bam_processing"]) / "output_sliced.bam"
                if sorted_bam and sorted_bam.exists():
                    sorted_bai = Path(f"{sorted_bam}.bai")
                    if sorted_bam.stat().st_size > 0 and (not sorted_bai.exists() or sorted_bai.stat().st_size == 0):
                        logger.info("Rebuilding missing or empty sliced BAM index for adVNTR: %s", sorted_bai)
                        index_cmd = f"{config['tools']['samtools']} index -@ {threads} {sorted_bam}"
                        run_command(
                            index_cmd,
                            str(Path(dirs["fastq_bam_processing"]) / "reindex_sliced_bam.log"),
                            critical=False,
                        )
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
                    advntr_execution_config = {**config, "tools": dict(advntr_context.tools)}
                    advntr_status = run_advntr(
                        advntr_context.model_snapshot,
                        sorted_bam,
                        dirs["advntr"],
                        "output",
                        config=advntr_execution_config,
                        cwd=project_root,
                        pipeline_threads=threads,
                        resolved_component=run_configuration.advntr,
                        runtime_component=run_configuration.advntr_runtime,
                        custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
                    )
                    if advntr_status != 0:
                        msg = f"adVNTR genotyping returned non-zero status {advntr_status}; result parsing was not attempted."
                        logger.error(msg)
                        raise RuntimeError(msg)
                    output_path = os.path.join(dirs["advntr"], f"output_adVNTR{output_ext}")
                    process_advntr_output(
                        output_path,
                        dirs["advntr"],
                        "output",
                        config=config,
                        artifact_evidence=advntr_evidence,
                        resolved_component=run_configuration.advntr,
                        nomenclature_component=run_configuration.nomenclature,
                        custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
                    )
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
                else:
                    logger.error("Sorted BAM required for adVNTR not provided.")
                    raise ValueError("Sorted BAM required for adVNTR not provided.")

            # Tier A needs two independent callers agreeing, which no single
            # caller stage can see. Without this step production could never
            # emit a tier-A name however well the two agreed (#nomenclature).
            # Runs whenever adVNTR is requested, ensuring caller results are reconciled
            # even when adVNTR was reused (#20).
            kestrel_tsv = os.path.join(dirs["kestrel"], "kestrel_result.tsv")
            advntr_tsv = os.path.join(dirs["advntr"], "output_adVNTR_result.tsv")
            reconciliation_outcome = reconcile_caller_outputs(
                kestrel_tsv,
                advntr_tsv,
                artifact_evidence=advntr_evidence,
                resolved_component=run_configuration.nomenclature,
                dominance_component=run_configuration.dominance,
                custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
            )
            if run_configuration.dominance.get("enabled") is True and (
                not isinstance(reconciliation_outcome, DominanceSeamOutcome) or not reconciliation_outcome.evaluated
            ):
                raise RuntimeError("enabled dominance seam did not evaluate the selected run policy")
            if (
                isinstance(reconciliation_outcome, DominanceSeamOutcome) and reconciliation_outcome.rewritten
            ) or reconciliation_outcome is True:
                refresh_step(
                    summary,
                    STEP_KESTREL,
                    write_summary_path=summary_file_path,
                )
                refresh_step(
                    summary,
                    STEP_ADVNTR,
                    write_summary_path=summary_file_path,
                )

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
            crossmatch_summary = cross_match_variants(
                kestrel_records,
                advntr_records,
                resolved_component=run_configuration.cross_match,
                custom_context_active=run_configuration.decision_profile.source == "explicit-cli",
            )
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
            logger.info("adVNTR module not included. Skipping adVNTR genotyping.")

        # --- Generate Summary Report and Archiving ---
        logger.info("Generating summary report.")
        report_file = "summary_report.html"
        template_dir = config.get("paths", {}).get("template_dir")

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
            report_igv=report_igv,
        )
        logger.info(f"Summary report generated: {report_file}")

        # Mark pipeline end in summary
        end_summary(summary)

        # Write out the complete pipeline summary
        summary_file_path = os.path.join(output_dir, "pipeline_summary.json")
        write_summary(summary, summary_file_path)
        logger.info(f"Pipeline summary written to: {summary_file_path}")

        donor_summary_path = Path(output_dir) / "pipeline_summary.donor.json"
        if donor_summary_path.is_file():
            with contextlib.suppress(OSError):
                donor_summary_path.unlink()

        # Each requested summary format writes the one-row-per-step table and, beside
        # it, the long rows table (#119). Unknown format names are ignored, as before.
        summary_writers = (("csv", ",", convert_summary_to_csv), ("tsv", "\t", convert_summary_to_tsv))
        for fmt, delimiter, write_table in summary_writers:
            if fmt not in (summary_formats or ()):
                continue
            table_path = os.path.join(output_dir, f"pipeline_summary.{fmt}")
            write_table(summary, table_path)
            rows_path = os.path.join(output_dir, f"pipeline_summary_rows.{fmt}")
            convert_summary_rows_to_delimited(summary, rows_path, delimiter=delimiter)
            logger.info(f"Pipeline summary {fmt.upper()} written to: {table_path} and {rows_path}")

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
