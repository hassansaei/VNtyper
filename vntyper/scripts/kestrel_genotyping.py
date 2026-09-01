"""
kestrel_genotyping.py

This script orchestrates mapping-free genotyping using the Kestrel tool,
focusing on MUC1-VNTR analysis. It coordinates both Kestrel execution
and postprocessing steps (filtering, scoring, confidence assignment,
motif annotation).

--------------------------------------------------------------------------------
High-Level Flow:
  1) Construct & run the Kestrel command to genotype VNTR from provided FASTQs.
  2) Convert intermediate SAM→BAM, ensuring indexing for further steps.
  3) Postprocess the resulting VCF:
     - Filter to INDEL variants.
     - Split into insertion vs. deletion.
     - Merge with MUC1 motif references.
     - Apply frameshift logic & empirical coverage cutoffs from
       Saei et al., iScience 26, 107171 (2023).
  4) Assign each variant a confidence label (e.g., Low_Precision, High_Precision).
  5) Generate final output: `kestrel_result.tsv`.

--------------------------------------------------------------------------------
from Saei et al., iScience 26, 107171 (2023).
"""

import logging
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from vntyper.scripts.command_builders import build_sam_to_bam_command, build_samtools_index_command, quote_path
from vntyper.scripts.confidence_assignment import (
    calculate_depth_score_and_assign_confidence,
)
from vntyper.scripts.file_processing import filter_indel_vcf, filter_vcf
from vntyper.scripts.flagging import (
    KESTREL_FLAG_COLUMNS,
    CompiledFlagRules,
    add_artifact_gate,
    add_flags,
    compile_flag_rules,
    validate_duplicate_flagging_config,
)
from vntyper.scripts.identity_candidate_persistence import (
    IDENTITY_CAPTURE_COLUMNS,
    candidate_capture_cells,
    selected_candidate_cells,
)
from vntyper.scripts.identity_candidates import (
    IdentityTranslationComponent,
    capture_kestrel_observations,
    overlay_legacy_projection,
    translation_component_from_config,
    with_candidate_evidence,
)
from vntyper.scripts.kestrel_command import construct_kestrel_command as construct_kestrel_command  # noqa: F401
from vntyper.scripts.kestrel_counting import DEFAULT_KANALYZE_PATH, execute_attempt
from vntyper.scripts.kestrel_execution import KestrelCommandArguments, plan_kestrel_invocations
from vntyper.scripts.kestrel_vcf_contract import describe_unusable_vcf
from vntyper.scripts.molecular_identity_presentation import (
    IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS,
    identity_translation_diagnostic_cells,
)
from vntyper.scripts.motif_processing import (
    load_additional_motifs,
    load_muc1_reference,
    motif_correction_and_annotation,  # Moved from the main script
    preprocessing_deletion,
    preprocessing_insertion,
)
from vntyper.scripts.nomenclature import load_nomenclature_config
from vntyper.scripts.nomenclature_annotate import annotate_kestrel_frame
from vntyper.scripts.scoring import (
    extract_frameshifts,
    split_depth_and_calculate_frame_score,
    split_frame_score,
)
from vntyper.scripts.subthreshold import detect_from_file, format_note
from vntyper.scripts.utils import load_config, run_command

# Modularized functions for variant parsing/scoring/confidence
from vntyper.scripts.variant_parsing import (
    filter_by_alt_values_and_finalize,
    read_vcf_without_comments,
)
from vntyper.version import __version__ as VERSION

logger = logging.getLogger(__name__)


def load_kestrel_config(config_path=None):
    """
    Loads the Kestrel configuration file or defaults to the
    local `kestrel_config.json`.

    For example, thresholds for "Depth_Score" and coverage cutoffs
    are read here, referencing the empirical cutoffs from Saei et al.,
    iScience 26, 107171 (2023).

    Args:
        config_path (str, optional): Path to a custom kestrel_config.json file.
            If not provided, use the default located next to this script.

    Returns:
        dict: Configuration dictionary with Kestrel-specific settings.
    """
    if config_path is None:
        # Default path to kestrel_config.json
        config_path = os.path.join(os.path.dirname(__file__), "kestrel_config.json")
    return load_config(config_path)


# Load the Kestrel configuration globally so it can be used
# by run_kestrel() and subsequent steps
kestrel_config = load_kestrel_config()


#: The boolean gate columns :func:`filter_final_dataframe` requires and ANDs, in order.
#:
#: A module constant rather than a list inside that function because
#: :func:`vntyper.scripts.subthreshold.detect` is handed it: eligibility for the #266
#: below-reporting-floor note is "fails ``depth_confidence_pass`` and nothing else", and
#: restating the gates there would let a seventh gate be added here and silently widen it.
FILTER_COLUMNS: tuple[str, ...] = (
    "is_frameshift",
    "is_valid_frameshift",
    "depth_confidence_pass",
    "alt_filter_pass",
    "motif_filter_pass",
    "flag_filter_pass",
)


def generate_header(reference_vntr, version=VERSION):
    """
    Creates a list of header lines containing metadata about
    the Kestrel genotyping process. Intended to prepend to TSV outputs.

    Args:
        reference_vntr (str): Path to the MUC1-VNTR reference used by Kestrel.
        version (str): VNtyper's version (defaults to the package VERSION).

    Returns:
        list of str: Header lines with version info & analysis date.
    """
    header = [
        "## VNtyper Kestrel result",
        f"## VNtyper Version: {version}",
        f"## Analysis date: {datetime.now(timezone.utc).astimezone().strftime('%Y-%m-%d %H:%M:%S')}",
        f"## Reference file: {reference_vntr}",
    ]
    return header


def convert_sam_to_bam_and_index(sam_file, output_dir, samtools_path="samtools", threads=1):
    """
    Converts the Kestrel-generated SAM to a BAM, indexes it for
    potential downstream usage, then deletes the original SAM.

    Args:
        sam_file (str): Path to the "output.sam" created by Kestrel.
        output_dir (str): Directory for the resulting BAM.
        samtools_path (str): Configured samtools invocation. Defaults to the
            historical bare executable.
        threads (int): Samtools thread count. One preserves the historical
            command strings byte-for-byte.

    Returns:
        str: Path to the indexed BAM file.
    """
    bam_file = os.path.join(output_dir, "output.bam")
    bam_index = bam_file + ".bai"

    # Both results are checked. They used to be discarded, and success was then inferred
    # from `os.path.exists` on the two outputs -- which conflates "samtools succeeded" with
    # "a file is present". A failed `view` still leaves a truncated BAM behind, so the SAM
    # was deleted and the truncated BAM kept, and the only record of the failure was a log
    # file nobody reads. This BAM is what the report's IGV track shows (#255).
    logger.info(f"Converting SAM to BAM: {sam_file} -> {bam_file}")
    if not run_command(
        build_sam_to_bam_command(
            samtools_path=samtools_path,
            sam_file=sam_file,
            bam_file=bam_file,
            threads=threads,
        ),
        log_file=os.path.join(output_dir, "samtools_view.log"),
    ):
        msg = f"Converting Kestrel's SAM to BAM failed: {sam_file} -> {bam_file}."
        logger.error(msg)
        raise RuntimeError(msg)

    logger.info(f"Indexing BAM file: {bam_file}")
    if not run_command(
        build_samtools_index_command(samtools_path=samtools_path, bam_file=bam_file, threads=threads),
        log_file=os.path.join(output_dir, "samtools_index.log"),
    ):
        msg = f"Indexing Kestrel's BAM failed: {bam_file}."
        logger.error(msg)
        raise RuntimeError(msg)

    # Only now is deleting the SAM safe. The existence test is kept as a belt-and-braces
    # check on what samtools claims to have written, not as the success signal it used to be.
    if os.path.exists(bam_file) and os.path.exists(bam_index):
        os.remove(sam_file)
        logger.info(f"Deleted SAM file: {sam_file}")

    return bam_file


def _discard_attempt_artifacts(output_dir, kmer_size, reason, vcf_path=None):
    """Remove what a discarded k-mer attempt left behind, so the next one starts clean.

    Every configured k-mer size writes to the same ``output.vcf`` and the same
    ``output.sam``. Both discard paths need this: the unusable-VCF branch, and the branch
    where Kestrel exits 0 having written no VCF at all -- that one has no VCF to remove but
    may still have written a SAM before giving up.

    Isolation is the point rather than a demonstrated corruption. The pinned Kestrel opens
    its haplotype output through Java's truncating ``FileOutputStream``, so an attempt that
    reaches SAM initialisation overwrites the file anyway; the gap is an attempt that exits
    before reaching it, whose predecessor's SAM would then be converted into ``output.bam``
    by a later successful attempt -- the alignment the report's IGV track shows (#255).

    Args:
        output_dir (str): Directory the attempt wrote into.
        kmer_size (int): The k-mer size being discarded, for the message.
        reason (str): Why the attempt is being discarded, for the message.
        vcf_path (pathlib.Path, optional): The VCF to remove. ``None`` when there is none.

    Raises:
        RuntimeError: If an artefact exists and cannot be removed. This must not become a
            new way to leave the loop silently: carrying the file into the next attempt
            would report it against the wrong k-mer size.
    """
    stale_paths = [path for path in (vcf_path, Path(output_dir) / "output.sam") if path is not None]
    for stale in stale_paths:
        try:
            stale.unlink()
        except FileNotFoundError:
            # Kestrel need not have written it. Absent is the desired state, not an error.
            continue
        except OSError as exc:
            msg = (
                f"Kestrel's attempt at k-mer size {kmer_size} was discarded because {reason}, and "
                f"{stale} could not be removed ({exc}). Continuing would carry that file into the "
                "next k-mer size and report it against the wrong one."
            )
            logger.error(msg)
            raise RuntimeError(msg) from exc


def run_kestrel(
    vcf_path,
    output_dir,
    fastq_files,
    reference_vntr,
    kestrel_path,
    config,
    sample_name,
    log_level=logging.INFO,
    cwd=None,
    threads=4,
):
    """
    Main entry point to run the Kestrel jar for MUC1 genotyping, then
    postprocess the output. By default, iterates over a set of k-mer
    sizes from the config (often just a single size: 20).

    Steps:
      1) Remove any VCF left behind by an earlier run (#212).
      2) Construct command for each k-mer size.
      3) Run Kestrel, capturing logs.
      4) Convert Kestrel's SAM→BAM, index.
      5) Call `process_kestrel_output()` for final filtering.

    Args:
        vcf_path (Path): Where the VCF from Kestrel should be written.
        output_dir (str): Folder for intermediate & final results.
        fastq_files (Sequence[str | Path]): Ordered FASTQs for one sample.
        reference_vntr (str): MUC1 reference FASTA for Kestrel.
        kestrel_path (str): Path to the Kestrel jar.
        config (dict): Overall pipeline config (tools, references, etc.).
        sample_name (str): Name to embed in outputs/VCF header.
        log_level (int): Logging level (INFO, DEBUG, etc.).
        cwd (str, optional): Working directory to use when running Java/Kestrel.
            Important for Java initialization.
        threads (int, optional): The run's total thread budget, allocated across
            KAnalyze's three concurrent counting stages. Appended **after** ``cwd``
            rather than inserted before it: this signature permits positional calls,
            so a new parameter in front of the existing optional ones would silently
            rebind them for any external caller. Defaults to 4, matching
            ``config.json``'s ``default_values.threads``.

    Raises:
        RuntimeError: If Kestrel fails for a given k-mer size, or if no configured
            k-mer size produced a VCF (see #212 -- returning silently there let the
            pipeline report a manufactured negative).

    Returns:
        None
    """
    kestrel_settings = kestrel_config.get("kestrel_settings", {})
    java_path = config["tools"]["java_path"]
    java_memory = kestrel_settings.get("java_memory", "12g")
    kmer_sizes = kestrel_settings.get("kmer_sizes", [20])
    max_align_states = kestrel_settings.get("max_align_states", 30)
    max_hap_states = kestrel_settings.get("max_hap_states", 30)
    log_level_str = logging.getLevelName(log_level)

    additional_settings = kestrel_settings.get("additional_settings", "")
    # Every new key is read with `.get` and a shipped default: `--config-path` replaces
    # the whole config rather than merging, so a replacement config legitimately lacks
    # all of them (trap 2). `kestrel_config` itself is a module global read at import
    # time, so tests patch the global rather than passing a config (trap 1).
    java_opts_count = kestrel_settings.get("java_opts_count", "")
    java_opts_call = kestrel_settings.get("java_opts_call", "-XX:+UseSerialGC")
    split_counting = kestrel_settings.get("split_counting", True)
    keep_ikc = kestrel_settings.get("keep_ikc", False)
    kanalyze_path = config["tools"].get("kanalyze", DEFAULT_KANALYZE_PATH)

    # #212: a pre-existing output.vcf used to skip the whole stage and `return`, which
    # also skipped the two statements that turn a VCF into a result. Re-running into a
    # directory left by an interrupted run then either re-reported a stale result or
    # produced none at all -- and a step with no result file is rendered as a negative.
    # Deliberate reuse belongs behind the --resume flag proposed in #20; this ad-hoc,
    # unflagged version is unsafe. The stale file is removed rather than left for the JAR
    # to overwrite so the log says why, and the removal sits *before* the loop rather than
    # inside it: the only way to reach a second iteration is for the first to have written
    # no VCF, so a per-iteration unlink is a no-op at best and, if the loop ever grows a
    # `continue` after a successful iteration, deletes a good result.
    if vcf_path.is_file():
        logger.warning(f"Removing a VCF left by an earlier run before re-running Kestrel: {vcf_path}")
        vcf_path.unlink()

    command_arguments = KestrelCommandArguments(
        kestrel_path=kestrel_path,
        reference_vntr=reference_vntr,
        vcf_out=vcf_path,
        java_path=java_path,
        java_memory=java_memory,
        max_align_states=max_align_states,
        max_hap_states=max_hap_states,
        log_level=log_level_str,
        sample_name=sample_name,
        additional_settings=additional_settings,
        java_opts_call=java_opts_call,
        java_opts_count=java_opts_count,
        kanalyze_path=kanalyze_path,
        threads=threads,
        split_counting=split_counting,
    )
    invocations = plan_kestrel_invocations(
        fastq_files=fastq_files,
        kmer_sizes=kmer_sizes,
        output_dir=Path(output_dir),
        command_arguments=command_arguments,
    )

    # Whether some k-mer size both ran Kestrel and post-processed its VCF. `break` below
    # is reachable only inside `if vcf_path.is_file()`, so a Kestrel invocation that exits
    # 0 without writing a VCF used to fall out of the loop and return None.
    completed = False

    for invocation in invocations:
        kmer_size = invocation.kmer_size

        # Counting, calling and their cleanup live in `kestrel_counting` (AGENTS.md
        # rule 3): this file is well over the ~650-line guideline and the split is the
        # region that changed. `run_command` is passed so the module attribute the test
        # suite patches is the one that runs.
        execute_attempt(invocation, cwd=cwd, keep_ikc=keep_ikc, run_command=run_command)
        logger.info(f"Mapping-free genotyping of MUC1-VNTR with k-mer size {kmer_size} done!")

        if vcf_path.is_file():
            problem = describe_unusable_vcf(vcf_path)
            if problem is None:
                # Convert the intermediate SAM→BAM (for debugging or IGV)
                sam_file = os.path.join(output_dir, "output.sam")
                convert_sam_to_bam_and_index(
                    sam_file,
                    output_dir,
                    samtools_path=config["tools"]["samtools"],
                    threads=threads,
                )

                # Postprocess final output
                process_kestrel_output(output_dir, vcf_path, reference_vntr, kestrel_config, config)
                completed = True
                break  # Stop after the first successful k-mer size

            # A VCF nothing can parse is the same kind of event as no VCF at all: Kestrel
            # exited 0 and produced nothing this stage can read. Treat it the same way and let
            # the next configured k-mer size try. Aborting here would contradict the deliberate
            # fall-through below, and the terminal raise still catches the case where every
            # k-mer size fails -- so no path reaches a manufactured negative either way (#223).
            logger.warning(
                f"Kestrel exited successfully for k-mer size {kmer_size} and wrote {vcf_path}, but "
                f"{problem}. Removing it and trying the next configured k-mer size."
            )
            # Every k-mer size writes to the same `vcf_path` *and* the same `output.sam`, so a
            # discarded attempt must leave neither behind. A stale VCF would be re-examined by
            # the next iteration and reported against the wrong k-mer size.
            #
            # The SAM is narrower than it first looks, and the honest statement is worth having
            # here. The pinned Kestrel opens its haplotype output through Java's truncating
            # `FileOutputStream`, so an attempt that *does* initialise SAM output rewrites the
            # file rather than appending to it -- the common case cleans up after itself. What
            # remains is an attempt that exits 0 and never reaches that initialisation: its
            # predecessor's SAM then survives, and a later successful attempt converts whatever
            # occupies the path into `output.bam`, the alignment the report's IGV track shows.
            # The genotype comes from the VCF and is unaffected either way. Removing it is
            # attempt isolation -- the same reason the VCF is removed -- not a fix for a
            # demonstrated conversion of stale data (#255).
            #
            # Both removals are guarded because neither may become a *new* way to leave the
            # loop: an unhandled OSError would escape as a bare filesystem error, skipping the
            # terminal RuntimeError that is this function's whole contract.
            _discard_attempt_artifacts(output_dir, kmer_size, reason=problem, vcf_path=vcf_path)
            continue

        # Exit status 0 with no VCF is the silent path #212 is about. Say so, then let the
        # next configured k-mer size try.
        #
        # This branch discards an attempt too, so it needs the same artefact isolation as the
        # unusable-VCF branch above: there is no VCF to remove, but Kestrel may still have
        # written a SAM before giving up, and every attempt shares that path (#255).
        logger.warning(f"Kestrel exited successfully for k-mer size {kmer_size} but wrote no VCF to {vcf_path}.")
        _discard_attempt_artifacts(output_dir, kmer_size, reason="it wrote no VCF")

    if not completed:
        msg = (
            "Kestrel produced no usable VCF for any configured k-mer size, so no result file was "
            "written. Every configured k-mer size either wrote no VCF at all, or wrote one that could "
            "not be parsed into records. Reporting this as a negative would manufacture a confident "
            "negative genotype. See issues #212 and #223."
        )
        logger.error(msg)
        raise RuntimeError(msg)


def _try_compress_vcf_with_bcftools(input_vcf, output_vcf_gz, output_dir):
    """
    Attempt to compress and sort a VCF file using bcftools.

    This function follows the Single Responsibility Principle (SRP) by doing
    exactly one thing: attempting VCF compression. It gracefully handles the
    case where bcftools is not available.

    Args:
        input_vcf (str): Path to the input uncompressed VCF file.
        output_vcf_gz (str): Path to the desired compressed output file.
        output_dir (str): Directory for log files.

    Returns:
        bool: True if compression succeeded, False otherwise.

    Notes:
        - If bcftools is not in PATH, logs a WARNING and returns False
        - If bcftools command fails, logs an ERROR and returns False
        - This allows the pipeline to gracefully fall back to uncompressed VCF
    """
    # Check if bcftools is available (defensive programming)
    if shutil.which("bcftools") is None:
        logger.warning(
            "bcftools not found in PATH. VCF compression skipped. "
            "IGV report will use uncompressed VCF. "
            "For optimal performance, install bcftools: 'conda install bcftools'"
        )
        return False

    # Attempt compression using existing run_command infrastructure (DRY principle)
    log_file = os.path.join(output_dir, "bcftools_sort.log")
    success = run_command(
        f"bcftools sort {quote_path(input_vcf)} -o {quote_path(output_vcf_gz)} -W -O z",
        log_file=log_file,
    )

    if not success:
        logger.error(
            f"bcftools sort command failed. Check {log_file} for details. IGV report will use uncompressed VCF."
        )
        return False

    logger.info(f"VCF successfully compressed: {output_vcf_gz}")
    return True


def _subthreshold_note(output_dir, config):
    """The #266 banner line for a sample that called nothing, or None.

    Emitted only on the no-call path, and only from the one `output_empty_result` branch
    that has a scored frame behind it. On the 400-sample simulated benchmark the eligible
    set covers 22 of 22 false negatives and 1 of 200 true negatives -- but also 85 of 178
    *called* samples, where the eligible rows are weaker descriptions of the event that was
    called. Printing it there would be noise beside a call, which is the user confusion
    @hassansaei warned about on #266; which of several passing candidates is reported is
    #270's subject, not this one.

    The two earlier `output_empty_result` branches run before anything is scored, so no
    pre-result of this run exists there and reading one could only pick up a stale file from
    an earlier run into the same output directory.

    Args:
        output_dir (str): The Kestrel output directory, holding `kestrel_pre_result.tsv`.
        config (dict): The Kestrel configuration.

    Returns:
        str or None: The line, or None when the feature is off, nothing is eligible, or the
        evidence cannot be read. Never raises: the note is an annotation, and losing one
        must not cost a result.
    """
    settings = config.get("subthreshold_note", {})
    if not settings.get("enabled", False):
        return None
    template = settings.get("template")
    if not template:
        logger.warning("subthreshold_note is enabled but carries no template; no note will be written.")
        return None
    # A partial configuration is an operator error the pipeline refuses elsewhere
    # (`calculate_depth_score_and_assign_confidence` raises KeyError on one), but this is an
    # annotation: a config without the floor should lose the note, not the run.
    try:
        floor = config["confidence_assignment"]["depth_score_thresholds"]["low"]
    except (KeyError, TypeError):
        logger.warning(
            "No confidence_assignment.depth_score_thresholds.low in the Kestrel config; "
            "no subthreshold note will be written."
        )
        return None
    signal = detect_from_file(os.path.join(output_dir, "kestrel_pre_result.tsv"), FILTER_COLUMNS, floor)
    if signal is None:
        return None
    return format_note(signal, template)


def process_kestrel_output(output_dir, vcf_path, reference_vntr, kestrel_config, config):
    """
    Processes the Kestrel output VCF files after Kestrel finishes.

    Steps:
      1) Filter the VCF to extract only INDELs -> output_indel.vcf
      2) Fix the file format line from "VCF4.2" to "VCFv4.2"
      3) Sort & index (via bcftools) the resulting INDEL VCF
         (used for future expansions, e.g., IGV coverage).
      4) Split into insertion and deletion VCFs.
      5) Merge with reference motif data, apply scoring & annotation
         logic (frame shifts, coverage thresholds, confidence).
      6) Write out final `kestrel_result.tsv`.

    Args:
        output_dir (str): Where intermediate & final files live.
        vcf_path (Path): The original VCF from Kestrel.
        reference_vntr (str): MUC1 reference file for motif annotation.
        kestrel_config (dict): Contains thresholds for depth, alt coverage, etc.
        config (dict): Main pipeline config (paths, additional references).

    Returns:
        pd.DataFrame or None:
          The final processed DataFrame of variants, or None if no variants found.

    Raises:
        ValueError: If the configured flag rules are invalid for the Kestrel result schema.
    """
    flagging_rules = kestrel_config.get("flagging_rules", {})
    compiled_flag_rules = compile_flag_rules(flagging_rules, KESTREL_FLAG_COLUMNS)
    duplicates_config = kestrel_config.get("duplicate_flagging", {})
    validate_duplicate_flagging_config(duplicates_config, compiled_flag_rules)
    logger.info("Processing Kestrel VCF results...")

    # Step 1) Filter the original VCF to extract INDELs
    indel_vcf = os.path.join(output_dir, "output_indel.vcf")
    filter_vcf(vcf_path, indel_vcf)

    # Step 2) Fix the file format line from "VCF4.2" -> "VCFv4.2"
    fixed_indel_vcf = indel_vcf + ".fixed"
    with open(indel_vcf) as fin, open(fixed_indel_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##fileformat=VCF4.2"):
                fout.write("##fileformat=VCFv4.2\n")
            else:
                fout.write(line)
    os.replace(fixed_indel_vcf, indel_vcf)

    # Step 3) Compress & sort using bcftools (if available)
    # Uses modular helper function following SRP (Single Responsibility Principle)
    sorted_indel_vcf_gz = indel_vcf + ".gz"
    _try_compress_vcf_with_bcftools(indel_vcf, sorted_indel_vcf_gz, output_dir)
    # Note: Compression may fail if bcftools unavailable - this is gracefully handled
    # The pipeline will continue and use uncompressed VCF for IGV report

    # Step 4) Split into insertion and deletion sub-VCFs
    output_ins = os.path.join(output_dir, "output_insertion.vcf")
    output_del = os.path.join(output_dir, "output_deletion.vcf")
    filter_indel_vcf(indel_vcf, output_ins, output_del)

    # The derived files can fail independently of the raw VCF, and `read_vcf_without_comments`
    # converts any read failure into an empty frame -- which four lines below is the `Negative`
    # placeholder. This is the second half of #223: "Kestrel ran and found nothing" and "Kestrel
    # produced nothing readable" must not render identically. Kestrel has already succeeded by
    # this point, so there is nothing left to retry and this raises rather than falling through.
    #
    # `read_vcf_without_comments` itself is deliberately untouched. Its empty-frame fallback is
    # a reviewed disposition (`scripts/ble001_policy.json`, preserved-no-authorized-alternative)
    # with tests asserting it as correct, so the check belongs at the call site instead.
    #
    # A valid header with zero records passes: that is the legitimate empty result, and it must
    # still reach `output_empty_result` below.
    for derived_vcf in (output_ins, output_del):
        problem = describe_unusable_vcf(derived_vcf)
        if problem is not None:
            msg = (
                f"Kestrel's derived VCF {derived_vcf} cannot be parsed: {problem}. Reporting this as a "
                "negative would manufacture a confident negative genotype. See issue #223."
            )
            logger.error(msg)
            raise RuntimeError(msg)

    # Step 5) Read the insertion & deletion VCFs
    header = generate_header(reference_vntr)
    vcf_insertion = read_vcf_without_comments(output_ins)
    vcf_deletion = read_vcf_without_comments(output_del)

    # If both are empty, produce an empty result file
    if vcf_insertion.empty and vcf_deletion.empty:
        logger.warning("No insertion/deletion variants found. Skipping.")
        output_empty_result(output_dir, header)
        return None

    # Merge with MUC1 reference data
    muc1_ref = load_muc1_reference(reference_vntr)

    # Preprocess insertion/deletion
    insertion_df = preprocessing_insertion(vcf_insertion, muc1_ref) if not vcf_insertion.empty else pd.DataFrame()
    deletion_df = preprocessing_deletion(vcf_deletion, muc1_ref) if not vcf_deletion.empty else pd.DataFrame()

    combined_df = pd.concat([insertion_df, deletion_df], axis=0)
    # Sort deterministically
    sort_columns = list(combined_df.columns)
    combined_df = combined_df.sort_values(by=sort_columns).reset_index(drop=True)

    if combined_df.empty:
        logger.warning("Empty combined DataFrame of insertions+deletions.")
        output_empty_result(output_dir, header)
        return None

    # Load additional motif data from config
    merged_motifs = load_additional_motifs(config)

    # Perform frame scoring, depth scoring, confidence assignment, etc.
    identity_component = translation_component_from_config(load_nomenclature_config())
    processed_df = process_kmer_results(
        combined_df,
        merged_motifs,
        output_dir,
        kestrel_config,
        compiled_flag_rules=compiled_flag_rules,
        identity_component=identity_component,
    )

    if processed_df.empty:
        logger.warning("Final processed DataFrame is empty. Writing empty result.")
        # #266: the one empty-result branch with a scored frame behind it, and therefore the
        # only one that can say anything about what was suppressed.
        output_empty_result(output_dir, header, note=_subthreshold_note(output_dir, kestrel_config))
        return None

    # Flagging is now applied inside process_kmer_results() (step 6.5) so that
    # flags are available before variant selection (fixes #145). Only re-apply
    # here if the Flag column is missing (e.g., process_kmer_results was called
    # without flagging rules configured).
    if "Flag" not in processed_df.columns and (compiled_flag_rules.rules or duplicates_config.get("enabled", False)):
        processed_df = add_flags(processed_df, compiled_flag_rules, duplicates_config=duplicates_config)

    # Name the variants before writing. Doing it here rather than in a later stage is
    # what makes one edit reach every surface: the TSV below, the pipeline summary
    # built from this same frame, and the HTML report all inherit the columns.
    processed_df = annotate_kestrel_frame(processed_df, output_dir, identity_component=identity_component)

    # Write the final processed results
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")
    with open(final_output_path, "w") as f:
        f.write("\n".join(header) + "\n")
        processed_df.to_csv(f, sep="\t", index=False)

    logger.info("Kestrel VCF processing completed.")
    return processed_df


def output_empty_result(output_dir, header, note=None):
    """
    Creates an empty result file with correct headers and a
    placeholder 'Negative' variant row to indicate no variants
    passed the filtering heuristics.

    Args:
        output_dir (str): Path to where we place the resulting .tsv
        header (list of str): The header lines from generate_header().
        note (str, optional): One extra banner line, appended after `header` as a `##`
            comment. #266's below-reporting-floor note arrives this way. It is a comment
            and never a row, so the 10-column placeholder schema below is unchanged and
            no consumer that reads the table can mistake it for a call -- `parse_tsv`
            routes `#` lines into `comments`, and `data` never sees them.
    """
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")

    empty_result_data = {
        "Motif": ["None"],
        "Variant": ["None"],
        "POS": ["None"],
        "REF": ["None"],
        "ALT": ["None"],
        "Motif_sequence": ["None"],
        "Estimated_Depth_AlternateVariant": ["None"],
        "Estimated_Depth_Variant_ActiveRegion": ["None"],
        "Depth_Score": ["None"],
        "Confidence": ["Negative"],
    }
    empty_df = pd.DataFrame(empty_result_data)

    banner = list(header)
    if note:
        banner.append(f"## {note}")
        logger.info("Recording a below-reporting-floor note on an otherwise empty result.")

    with open(final_output_path, "w") as f:
        f.write("\n".join(banner) + "\n")
        empty_df.to_csv(f, sep="\t", index=False)

    logger.info(f"Empty result file with placeholder saved at {final_output_path}")


def process_kmer_results(
    combined_df,
    merged_motifs,
    output_dir,
    kestrel_config,
    compiled_flag_rules: CompiledFlagRules | None = None,
    identity_component: IdentityTranslationComponent | None = None,
):
    """
    Applies the main postprocessing heuristics:
      1) Split the depth from the 'Sample' column and compute frame score
      2) Split frame score into numeric parts, mark frameshifts vs. non-frameshifts
      3) Extract frameshift variants (3n+1 / 3n+2)
      4) Compute Depth_Score, assign confidence
      4.5) Sort by Depth_Score + add haplo_count
      5) ALT-based filtering logic (e.g., 'GG' threshold)
      6) Motif correction & annotation
      6.5) Apply flagging rules before selection (fixes #145)
      6.5b) Derive the artifact gate from the flags (fixes #174)
      7) Final filter + select single best variant
      8) Generate BED file for coverage

    References:
      - Saei et al., iScience 26, 107171 (2023) for empirical cutoffs

    Args:
        combined_df (pd.DataFrame):
            A combined DataFrame of insertions & deletions with some columns
            (POS, REF, ALT, Sample, etc.).
        merged_motifs (pd.DataFrame):
            Additional motif info loaded from config for final annotation.
        output_dir (str):
            Folder to store optional BED file / final outputs.
        kestrel_config (dict):
            Contains thresholds for depth score & confidence assignment.
        compiled_flag_rules: Rules already validated against the Kestrel postprocessing
            schema by the VCF-level caller. Direct callers may omit this argument.
        identity_component: Explicit complete-context translation component. The
            production VCF boundary supplies the checked-in nomenclature component;
            focused legacy callers may omit identity capture.

    Returns:
        pd.DataFrame: The final, fully annotated & filtered DataFrame. Could be empty.
    """
    if compiled_flag_rules is None:
        compiled_flag_rules = compile_flag_rules(kestrel_config.get("flagging_rules", {}), KESTREL_FLAG_COLUMNS)
    duplicates_config = kestrel_config.get("duplicate_flagging", {})
    validate_duplicate_flagging_config(duplicates_config, compiled_flag_rules)

    if combined_df.empty:
        return combined_df

    # (1) Split depth fields & calculate initial frame score
    df = split_depth_and_calculate_frame_score(combined_df)
    if df.empty:
        return df

    # (2) Split frame score into numeric 'direction'/'frameshift_amount'
    df = split_frame_score(df)
    if df.empty:
        return df

    # (3) Extract frameshifts by analyzing the pattern of 3n+1 or 3n+2
    df = extract_frameshifts(df)
    if df.empty:
        return df

    # (4) Assign confidence via coverage-based heuristics from config
    df = calculate_depth_score_and_assign_confidence(df, kestrel_config)
    if df.empty:
        return df

    identity_candidates = None
    if identity_component is not None:
        capture_records = df.to_dict("records")
        for record in capture_records:
            record["Motif_sequence"] = str(record["Motif_sequence"])
        identity_candidates = capture_kestrel_observations(capture_records, identity_component)
        capture_rows = [candidate_capture_cells(candidate) for candidate in identity_candidates.candidates]
        df = df.copy()
        for column in IDENTITY_CAPTURE_COLUMNS:
            df[column] = [cells[column] for cells in capture_rows]
        diagnostic_rows = [
            identity_translation_diagnostic_cells(candidate.observation.translation)
            for candidate in identity_candidates.candidates
        ]
        for column in IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS:
            df[column] = [cells[column] for cells in diagnostic_rows]

    # (4.5) Add haplo_count after confidence assignment
    df = add_haplo_count(df)

    # (5) Filter certain ALT values (e.g., discarding 'GG' if below threshold)
    df = filter_by_alt_values_and_finalize(df, kestrel_config)
    if df.empty:
        return df

    # (6) Motif correction & annotation
    df = motif_correction_and_annotation(df, merged_motifs, kestrel_config)
    if df.empty:
        return df

    # (6.5) Apply flagging BEFORE selection so flags inform variant choice.
    # All columns needed by flagging rules (Depth_Score, Motif, REF, ALT)
    # are available after step (6). Moving flagging here fixes #145: previously,
    # a flagged variant could be selected over an unflagged one because
    # select_single_best_variant ran before add_flags.
    if compiled_flag_rules.rules or duplicates_config.get("enabled", False):
        df = add_flags(df, compiled_flag_rules, duplicates_config=duplicates_config)

    # (6.5b) #174: derive the artifact gate. Unconditional, unlike add_flags above: a
    # frame that reached the final filter without `flag_filter_pass` would abort the run
    # on a missing required gate column (#185). Which flags are artifacts is
    # configuration; an absent or empty `artifact_flags` list excludes nothing, which is
    # exactly the pre-#174 behaviour.
    df = add_artifact_gate(df, kestrel_config.get("artifact_flags", []))

    evidenced_candidates = None
    passing_identity_ordinals: tuple[int, ...] = ()
    if identity_candidates is not None:
        evidenced_candidates = with_candidate_evidence(identity_candidates, df.to_dict("records"))
        passing_mask = df[list(FILTER_COLUMNS)].all(axis=1)
        passing_identity_ordinals = tuple(
            int(serialized) for serialized in df.loc[passing_mask, IDENTITY_CAPTURE_COLUMNS[5]]
        )

    # (7) Final Filter
    df = filter_final_dataframe(df, output_dir)
    if df.empty:
        logger.info("All rows failed one or more filter criteria. Returning empty.")
        return df

    if evidenced_candidates is not None:
        df = df.drop(columns=list(IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS))
        selected_ordinal = int(df.iloc[0][IDENTITY_CAPTURE_COLUMNS[5]])
        selected_candidates = overlay_legacy_projection(
            evidenced_candidates,
            passing_identity_ordinals,
            selected_ordinal,
        )
        for column, value in selected_candidate_cells(selected_candidates).items():
            df[column] = value

    # (8) Now generate the BED file from the fully filtered result
    bed_file_path = generate_bed_file(df, output_dir)
    if bed_file_path:
        logger.info(f"BED file created at: {bed_file_path}")

    return df


def generate_bed_file(df, output_dir):
    """
    Generates a BED file from the final Kestrel output DataFrame
    if 'Motif_fasta' and 'POS_fasta' columns exist. This can help
    with coverage visualizations in IGV or other genome browsers.

    Column 1 is 'Motif_fasta', the 120 bp pair record of
    All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa that the variant was called
    against -- not 'Motif', which is a half-motif name and not a record in that file.
    'POS_fasta' is the 1-based VCF position inside that record, and BED intervals are
    0-based and half-open, so it is written as [POS_fasta - 1, POS_fasta) (#203).

    Args:
        df (pd.DataFrame): The final processed DataFrame with motif info.
        output_dir (str): Where to save the resulting BED file.

    Returns:
        str or None: Path to the generated BED or None if data is missing/empty.
    """
    # We only generate a BED if columns 'Motif_fasta' & 'POS_fasta' exist
    if "Motif_fasta" not in df.columns or "POS_fasta" not in df.columns:
        logger.warning("Missing 'Motif_fasta' or 'POS_fasta' in DataFrame. Skipping BED file generation.")
        return None

    if df.empty:
        logger.warning("DataFrame is empty. No variants to generate a BED file.")
        return None

    bed_file_path = os.path.join(output_dir, "output.bed")

    # Each row: motif_fasta, start=pos_fasta-1, end=pos_fasta.
    #
    # `POS_fasta` is the 1-based VCF position within the 120 bp pair record named by
    # `Motif_fasta` (#203). BED intervals are 0-based and half-open, so position p is the
    # interval [p-1, p). This used to write [p, p+1), which named the base after the
    # variant - IGV highlighted the wrong position on every row.
    with open(bed_file_path, "w") as bed_file:
        for _, row in df.iterrows():
            motif_fasta = row["Motif_fasta"]
            pos = int(row["POS_fasta"])
            bed_file.write(f"{motif_fasta}\t{pos - 1}\t{pos}\n")

    logger.info(f"BED file generated at: {bed_file_path}")
    return bed_file_path


def add_haplo_count(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate haplo_count: number of variants sharing (POS, REF, ALT).

    Higher haplo_count = more supporting evidence across haplotypes.
    This represents how many times the exact same variant (position + alleles)
    appears across different haplotype calls.

    Example:
        POS=67, REF=G, ALT=GG appears 389 times → haplo_count=389
        POS=54, REF=C, ALT=GC appears 176 times → haplo_count=176

    Args:
        df: DataFrame with POS, REF, ALT columns

    Returns:
        DataFrame with haplo_count column added

    References:
        - Issue #136 fix implementation
        - Streamlined implementation plan (Day 2)
    """
    df = df.copy()

    if df.empty:
        df["haplo_count"] = 0
        return df

    # Group by exact variant (POS, REF, ALT) and count occurrences
    if all(col in df.columns for col in ["POS", "REF", "ALT"]):
        df["haplo_count"] = df.groupby(["POS", "REF", "ALT"])["ALT"].transform("size")
    else:
        logger.warning("Missing POS/REF/ALT columns for haplo_count calculation")
        df["haplo_count"] = 0

    return df


def select_single_best_variant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select the single best variant (Hassan's requirement: "one representative variant").

    Selection criteria (strict priority order):
        1. Highest Confidence level (High_Precision* > High_Precision > Low_Precision > Negative)
        2. Unflagged variants preferred over flagged (fixes #145)
        3. Highest Depth_Score (coverage tie-breaker)
        4. Highest haplo_count (most supporting evidence)
        5. Lowest POS (genomic position tie-breaker for reproducibility)

    This ensures deterministic, biologically-informed variant selection.

    Args:
        df: DataFrame with Confidence, haplo_count, Depth_Score, POS columns

    Returns:
        DataFrame with exactly 1 row (the best variant), or empty if input empty

    References:
        - Issue #136 fix implementation
        - PR #140 code review (Critical Issue #2: tie-breaking)
        - Streamlined implementation plan (Day 3)

    Examples:
        >>> df = pd.DataFrame(
        ...     {
        ...         "Confidence": ["High_Precision", "High_Precision", "Low_Precision"],
        ...         "haplo_count": [389, 120, 15],
        ...         "Depth_Score": [0.010, 0.009, 0.006],
        ...         "POS": [67, 67, 54],
        ...     }
        ... )
        >>> result = select_single_best_variant(df)
        >>> len(result)
        1
        >>> result.iloc[0]["haplo_count"]
        389
    """
    if df.empty:
        return df

    # Define Confidence priority mapping
    confidence_priority = {
        "High_Precision*": 3,
        "High_Precision": 2,
        "Low_Precision": 1,
        "Negative": 0,
    }

    df = df.copy()
    df["_priority"] = df["Confidence"].map(confidence_priority).fillna(0)

    # Ensure numeric types for sorting (create missing columns with default 0)
    if "haplo_count" not in df.columns:
        df["haplo_count"] = 0
    else:
        df["haplo_count"] = pd.to_numeric(df["haplo_count"], errors="coerce").fillna(0)

    if "Depth_Score" not in df.columns:
        df["Depth_Score"] = 0
    else:
        df["Depth_Score"] = pd.to_numeric(df["Depth_Score"], errors="coerce").fillna(0)

    if "POS" not in df.columns:
        df["POS"] = 0
    else:
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce").fillna(0)

    # Deprioritize flagged variants: unflagged (1) sorts before flagged (0).
    # If Flag column is absent, all variants are treated as unflagged.
    if "Flag" in df.columns:
        df["_is_unflagged"] = (df["Flag"] == "Not flagged").astype(int)
    else:
        df["_is_unflagged"] = 1

    # Multi-key sort (deterministic tie-breaking)
    # Priority: Confidence DESC, unflagged DESC, Depth_Score DESC, haplo_count DESC, POS ASC
    df = df.sort_values(
        ["_priority", "_is_unflagged", "Depth_Score", "haplo_count", "POS"],
        ascending=[False, False, False, False, True],
    )

    # Keep only the best variant
    result = df.head(1).drop(columns=["_priority", "_is_unflagged"])

    logger.info(
        "Selected best variant: Confidence=%s, haplo_count=%d, Depth_Score=%.5f, POS=%d",
        result.iloc[0]["Confidence"],
        int(result.iloc[0]["haplo_count"]),
        result.iloc[0]["Depth_Score"],
        int(result.iloc[0]["POS"]),
    )

    return result


def filter_final_dataframe(df: pd.DataFrame, output_dir: str) -> pd.DataFrame:
    """
    Final step: filter the DataFrame based on the boolean columns introduced
    by earlier steps ('is_frameshift', 'is_valid_frameshift',
    'depth_confidence_pass', 'alt_filter_pass', 'motif_filter_pass',
    'flag_filter_pass').

    Five of the six encode a *pathogenicity or quality* judgement. The sixth,
    'flag_filter_pass' (#174), encodes an *artifact* judgement: the row carries a
    flag that `kestrel_config.json` declares under `artifact_flags`, so it is not
    a candidate variant at all. Advisory flags such as
    'Low_Depth_Conserved_Motifs' leave that gate True and only deprioritise the
    row during selection.

    We keep rows where *all* filter columns are True. Every column is required:
    the earlier stages mark rather than filter, so a column missing from a
    non-empty frame would turn its safety gate into a permit. That is an error
    (issue #185), not something to skip.

    An empty frame is the one documented exception: it carries no candidate
    variant, so there is nothing for a missing gate to permit, and it is
    returned unchanged.

    Additionally, this function logs how many rows pass or fail each filter
    and writes the unfiltered DataFrame to 'kestrel_pre_result.tsv' directly
    in the 'output_dir' -- before any raise, so an aborted run keeps its
    evidence. The final filtered DataFrame is returned in memory.

    Args:
        df (pd.DataFrame): The postprocessed DataFrame, with all six
            boolean filter columns.
        output_dir (str): Path to the main output directory.

    Returns:
        pd.DataFrame: A copy of `df` containing only rows that pass
            all filter columns, or `df` itself when it is empty.

    Raises:
        ValueError: If a required filter column is missing from a non-empty
            DataFrame.
    """
    logger.info("Starting final filter of DataFrame with %d rows...", len(df))

    # Write the unfiltered DataFrame to 'kestrel_pre_result.tsv' in output_dir
    pre_result_path = os.path.join(output_dir, "kestrel_pre_result.tsv")
    df.to_csv(pre_result_path, sep="\t", index=False)
    logger.info("Wrote pre-filter DataFrame to %s", pre_result_path)

    # #185: the explicit, documented empty-result path. An empty frame holds no
    # candidate variant, so no gate can be silently turned into a permit here.
    # In production this is unreachable -- process_kmer_results guards every stage
    # ahead of the single call site with `if df.empty: return df` -- but stating it
    # keeps the carve-out explicit rather than resting on that structure.
    if df.empty:
        logger.info("Empty DataFrame reached the final filter; returning it unchanged.")
        return df

    # Columns every non-empty frame is required to carry
    filter_cols = FILTER_COLUMNS

    # Build a mask requiring all existing boolean filters == True
    final_mask = pd.Series(True, index=df.index)
    for col in filter_cols:
        if col in df.columns:
            before_count = final_mask.sum()
            final_mask &= df[col]
            after_count = final_mask.sum()
            logger.info(
                "Filter column '%s' exists; %d -> %d rows remain after requiring True.",
                col,
                before_count,
                after_count,
            )
        else:
            # #185: a missing gate is an error, not a permit. @hassansaei:
            # "a missing required gate column should raise (abort the run), not
            # be skipped [...] That is not acceptable for this pipeline."
            # Reachability: this function's only caller is process_kmer_results,
            # behind six `if df.empty: return df` guards, so any frame arriving
            # here is non-empty and has traversed every stage that adds a gate
            # column. An empty frame short-circuits above -- that is the explicit
            # empty-result path his decision carves out.
            msg = (
                f"Required filter column '{col}' is missing from a non-empty Kestrel result frame. "
                "An upstream stage stopped emitting it, so its safety gate would silently become a "
                "permit. Aborting rather than reporting unfiltered variants. See issue #185."
            )
            logger.error(msg)
            raise ValueError(msg)

    filtered_df = df[final_mask].copy()
    logger.info("Final DataFrame has %d rows after all filters.", len(filtered_df))

    # Select single best variant using multi-key priority sorting
    if len(filtered_df) > 1:
        filtered_df = select_single_best_variant(filtered_df)
        logger.info("Selected 1 best variant from %d candidates using priority sorting.", len(df[final_mask]))
    elif len(filtered_df) == 1:
        logger.info("Only 1 variant passed all filters (no selection needed).")
    else:
        logger.info("No variants passed all filters.")

    return filtered_df
