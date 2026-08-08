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

import pandas as pd

from vntyper.scripts.command_builders import quote_path
from vntyper.scripts.confidence_assignment import (
    calculate_depth_score_and_assign_confidence,
)
from vntyper.scripts.file_processing import filter_indel_vcf, filter_vcf
from vntyper.scripts.kestrel_command import construct_kestrel_command
from vntyper.scripts.motif_processing import (
    load_additional_motifs,
    load_muc1_reference,
    motif_correction_and_annotation,  # Moved from the main script
    preprocessing_deletion,
    preprocessing_insertion,
)
from vntyper.scripts.scoring import (
    extract_frameshifts,
    split_depth_and_calculate_frame_score,
    split_frame_score,
)
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


def convert_sam_to_bam_and_index(sam_file, output_dir):
    """
    Converts the Kestrel-generated SAM to a BAM, indexes it for
    potential downstream usage, then deletes the original SAM.

    Args:
        sam_file (str): Path to the "output.sam" created by Kestrel.
        output_dir (str): Directory for the resulting BAM.

    Returns:
        str: Path to the indexed BAM file.
    """
    bam_file = os.path.join(output_dir, "output.bam")
    bam_index = bam_file + ".bai"

    # Convert SAM to BAM using samtools
    logger.info(f"Converting SAM to BAM: {sam_file} -> {bam_file}")
    run_command(
        f"samtools view -Sb {quote_path(sam_file)} > {quote_path(bam_file)}",
        log_file=os.path.join(output_dir, "samtools_view.log"),
    )

    # Index the BAM file
    logger.info(f"Indexing BAM file: {bam_file}")
    run_command(
        f"samtools index {quote_path(bam_file)}",
        log_file=os.path.join(output_dir, "samtools_index.log"),
    )

    # Delete the SAM file if indexing succeeds
    if os.path.exists(bam_file) and os.path.exists(bam_index):
        os.remove(sam_file)
        logger.info(f"Deleted SAM file: {sam_file}")

    return bam_file


def run_kestrel(
    vcf_path,
    output_dir,
    fastq_1,
    fastq_2,
    reference_vntr,
    kestrel_path,
    config,
    sample_name,
    log_level=logging.INFO,
    cwd=None,
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
        fastq_1 (str): Path to FASTQ R1.
        fastq_2 (str | None): Optional path to FASTQ R2.
        reference_vntr (str): MUC1 reference FASTA for Kestrel.
        kestrel_path (str): Path to the Kestrel jar.
        config (dict): Overall pipeline config (tools, references, etc.).
        sample_name (str): Name to embed in outputs/VCF header.
        log_level (int): Logging level (INFO, DEBUG, etc.).
        cwd (str, optional): Working directory to use when running Java/Kestrel.
            Important for Java initialization.

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

    # Retrieve additional settings (defaults to empty)
    additional_settings = kestrel_settings.get("additional_settings", "")

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

    # Whether some k-mer size both ran Kestrel and post-processed its VCF. `break` below
    # is reachable only inside `if vcf_path.is_file()`, so a Kestrel invocation that exits
    # 0 without writing a VCF used to fall out of the loop and return None.
    completed = False

    # Try each k-mer size in sequence. Usually just [20], can be more.
    for kmer_size in kmer_sizes:
        kmer_command = construct_kestrel_command(
            kmer_size=kmer_size,
            kestrel_path=kestrel_path,
            reference_vntr=reference_vntr,
            output_dir=output_dir,
            fastq_1=fastq_1,
            fastq_2=fastq_2,
            vcf_out=vcf_path,
            java_path=java_path,
            java_memory=java_memory,
            max_align_states=max_align_states,
            max_hap_states=max_hap_states,
            log_level=log_level_str,
            sample_name=sample_name,
            additional_settings=additional_settings,
        )

        log_file = os.path.join(output_dir, f"kestrel_kmer_{kmer_size}.log")

        logger.info(f"Launching Kestrel with k-mer size {kmer_size}...")

        # Actually run the Kestrel command
        if not run_command(kmer_command, log_file, critical=True, cwd=cwd):
            logger.error(f"Kestrel failed for k-mer size {kmer_size}. Check {log_file} for details.")
            raise RuntimeError(f"Kestrel failed for kmer size {kmer_size}.")

        logger.info(f"Mapping-free genotyping of MUC1-VNTR with k-mer size {kmer_size} done!")

        # Now that Kestrel completed, confirm the VCF is present
        if vcf_path.is_file():
            # Convert the intermediate SAM→BAM (for debugging or IGV)
            sam_file = os.path.join(output_dir, "output.sam")
            convert_sam_to_bam_and_index(sam_file, output_dir)

            # Postprocess final output
            process_kestrel_output(output_dir, vcf_path, reference_vntr, kestrel_config, config)
            completed = True
            break  # Stop after the first successful k-mer size

        # Exit status 0 with no VCF is the silent path #212 is about. Say so, then let the
        # next configured k-mer size try.
        logger.warning(f"Kestrel exited successfully for k-mer size {kmer_size} but wrote no VCF to {vcf_path}.")

    if not completed:
        msg = (
            "Kestrel produced no VCF for any configured k-mer size, so no result file was written. "
            "Reporting this as a negative would manufacture a confident negative genotype. See issue #212."
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
    """
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
    processed_df = process_kmer_results(combined_df, merged_motifs, output_dir, kestrel_config)

    if processed_df.empty:
        logger.warning("Final processed DataFrame is empty. Writing empty result.")
        output_empty_result(output_dir, header)
        return None

    # Flagging is now applied inside process_kmer_results() (step 6.5) so that
    # flags are available before variant selection (fixes #145). Only re-apply
    # here if the Flag column is missing (e.g., process_kmer_results was called
    # without flagging rules configured).
    if "Flag" not in processed_df.columns:
        flagging_rules = kestrel_config.get("flagging_rules", {})
        duplicates_config = kestrel_config.get("duplicate_flagging", {})

        if flagging_rules or duplicates_config.get("enabled", False):
            from vntyper.scripts.flagging import add_flags

            processed_df = add_flags(processed_df, flagging_rules, duplicates_config=duplicates_config)

    # Write the final processed results
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")
    with open(final_output_path, "w") as f:
        f.write("\n".join(header) + "\n")
        processed_df.to_csv(f, sep="\t", index=False)

    logger.info("Kestrel VCF processing completed.")
    return processed_df


def output_empty_result(output_dir, header):
    """
    Creates an empty result file with correct headers and a
    placeholder 'Negative' variant row to indicate no variants
    passed the filtering heuristics.

    Args:
        output_dir (str): Path to where we place the resulting .tsv
        header (list of str): The header lines from generate_header().
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

    with open(final_output_path, "w") as f:
        f.write("\n".join(header) + "\n")
        empty_df.to_csv(f, sep="\t", index=False)

    logger.info(f"Empty result file with placeholder saved at {final_output_path}")


def process_kmer_results(combined_df, merged_motifs, output_dir, kestrel_config):
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

    Returns:
        pd.DataFrame: The final, fully annotated & filtered DataFrame. Could be empty.
    """
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
    flagging_rules = kestrel_config.get("flagging_rules", {})
    duplicates_config = kestrel_config.get("duplicate_flagging", {})
    if flagging_rules or duplicates_config.get("enabled", False):
        from vntyper.scripts.flagging import add_flags

        df = add_flags(df, flagging_rules, duplicates_config=duplicates_config)

    # (6.5b) #174: derive the artifact gate. Unconditional, unlike add_flags above: a
    # frame that reached the final filter without `flag_filter_pass` would abort the run
    # on a missing required gate column (#185). Which flags are artifacts is
    # configuration; an absent or empty `artifact_flags` list excludes nothing, which is
    # exactly the pre-#174 behaviour.
    from vntyper.scripts.flagging import add_artifact_gate

    df = add_artifact_gate(df, kestrel_config.get("artifact_flags", []))

    # (7) Final Filter
    df = filter_final_dataframe(df, output_dir)
    if df.empty:
        logger.info("All rows failed one or more filter criteria. Returning empty.")
        return df

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
    filter_cols = [
        "is_frameshift",
        "is_valid_frameshift",
        "depth_confidence_pass",
        "alt_filter_pass",
        "motif_filter_pass",
        "flag_filter_pass",
    ]

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
