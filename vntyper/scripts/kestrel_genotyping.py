#!/usr/bin/env python3

import logging
import os
from datetime import datetime
from pathlib import Path

import pandas as pd

from vntyper.scripts.file_processing import filter_indel_vcf, filter_vcf
from vntyper.scripts.motif_processing import (
    load_additional_motifs,
    load_muc1_reference,
    preprocessing_deletion,
    preprocessing_insertion,
)
from vntyper.scripts.utils import load_config, run_command
from vntyper.version import __version__ as VERSION


def load_kestrel_config(config_path=None):
    """
    Loads the Kestrel configuration file.

    Args:
        config_path (str, optional): Path to the kestrel_config.json file.
            Defaults to None.

    Returns:
        dict: Configuration dictionary loaded from kestrel_config.json.
    """
    if config_path is None:
        # Default path to kestrel_config.json
        config_path = os.path.join(
            os.path.dirname(__file__), 'kestrel_config.json'
        )
    return load_config(config_path)


# Load the Kestrel configuration
kestrel_config = load_kestrel_config()


def construct_kestrel_command(
    kmer_size,
    kestrel_path,
    reference_vntr,
    output_dir,
    fastq_1,
    fastq_2,
    vcf_out,
    java_path,
    java_memory,
    max_align_states,
    max_hap_states,
    log_level,
    sample_name,
):
    """
    Constructs the command for running Kestrel based on various settings.

    Args:
        kmer_size (int): Size of the kmer to use.
        kestrel_path (str): Path to the Kestrel jar file.
        reference_vntr (str): Path to the reference VNTR file.
        output_dir (str): Directory for the output files.
        fastq_1 (str): Path to the first FASTQ file.
        fastq_2 (str): Path to the second FASTQ file.
        vcf_out (str): Path to the VCF output file.
        java_path (str): Path to the Java executable.
        java_memory (str): Amount of memory to allocate to the JVM.
        max_align_states (int): Maximum alignment states.
        max_hap_states (int): Maximum haplotype states.
        log_level (str): Log level to use for Kestrel.
        sample_name (str): Sample name to be set using the -s/--sample option.

    Returns:
        str: The constructed Kestrel command.
    """
    if not fastq_1 or not fastq_2:
        raise ValueError("FASTQ input files are missing or invalid.")

    return (
        f"{java_path} -Xmx{java_memory} -jar {kestrel_path} -k {kmer_size} "
        f"--maxalignstates {max_align_states} --maxhapstates {max_hap_states} "
        f"-r {reference_vntr} -o {vcf_out} "
        f"-s{sample_name} {fastq_1} {fastq_2} "
        f"--hapfmt sam -p {output_dir}/output.sam --logstderr --logstdout "
        f"--loglevel {log_level.upper()} --temploc {output_dir}"
    )


def generate_header(reference_vntr, version=VERSION):
    """
    Generates a header for the output files with analysis information.

    Args:
        reference_vntr (str): Path to the reference VNTR file.
        version (str): Version of the tool being used.

    Returns:
        list: A list of header lines to be written at the top of the output file.
    """
    header = [
        "## VNtyper Kestrel result",
        f"## VNtyper Version: {version}",
        f"## Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"## Reference file: {reference_vntr}",
    ]
    return header


def convert_sam_to_bam_and_index(sam_file, output_dir):
    """
    Converts the SAM file to BAM format, indexes the BAM file, and deletes the SAM file.

    Args:
        sam_file (str): Path to the SAM file to be converted.
        output_dir (str): Directory where the BAM and BAM index files will be saved.

    Returns:
        str: Path to the indexed BAM file.
    """
    bam_file = os.path.join(output_dir, "output.bam")
    bam_index = bam_file + ".bai"

    # Convert SAM to BAM
    logging.info(f"Converting SAM to BAM: {sam_file} -> {bam_file}")
    run_command(
        f"samtools view -Sb {sam_file} > {bam_file}",
        log_file=os.path.join(output_dir, "samtools_view.log"),
    )

    # Index the BAM file
    logging.info(f"Indexing BAM file: {bam_file}")
    run_command(
        f"samtools index {bam_file}",
        log_file=os.path.join(output_dir, "samtools_index.log"),
    )

    # Delete the SAM file
    if os.path.exists(bam_file) and os.path.exists(bam_index):
        os.remove(sam_file)
        logging.info(f"Deleted SAM file: {sam_file}")

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
):
    """
    Orchestrates the Kestrel genotyping process by iterating through kmer sizes
    and processing the VCF output.

    Args:
        vcf_path (Path): Path to the output VCF file.
        output_dir (str): Directory where Kestrel outputs will be saved.
        fastq_1 (str): Path to the first FASTQ file.
        fastq_2 (str): Path to the second FASTQ file.
        reference_vntr (str): Path to the reference VNTR file.
        kestrel_path (str): Path to the Kestrel jar file.
        config (dict): Overall configuration dictionary.
        sample_name (str): The sample name to pass to Kestrel.
        log_level (int): Logging level.
    """
    global kestrel_config  # Ensure we use the global kestrel_config

    kestrel_settings = kestrel_config.get("kestrel_settings", {})
    java_path = config["tools"]["java_path"]
    java_memory = kestrel_settings.get("java_memory", "12g")
    kmer_sizes = kestrel_settings.get("kmer_sizes", [20])
    max_align_states = kestrel_settings.get("max_align_states", 30)
    max_hap_states = kestrel_settings.get("max_hap_states", 30)
    log_level_str = logging.getLevelName(log_level)  # Get the current log level as a string

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
            log_level=log_level_str,  # Pass the log level here
            sample_name=sample_name,
        )

        log_file = os.path.join(output_dir, f"kestrel_kmer_{kmer_size}.log")

        if vcf_path.is_file():
            logging.info("VCF file already exists, skipping Kestrel run...")
            return
        else:
            logging.info(f"Launching Kestrel with kmer size {kmer_size}...")

            # Run the command and log output to the specified log file
            if not run_command(kmer_command, log_file, critical=True):
                logging.error(
                    f"Kestrel failed for kmer size {kmer_size}. Check {log_file} for details."
                )
                raise RuntimeError(f"Kestrel failed for kmer size {kmer_size}.")

            logging.info(f"Mapping-free genotyping of MUC1-VNTR with kmer size {kmer_size} done!")

            if vcf_path.is_file():
                # Convert the SAM to BAM and index it
                sam_file = os.path.join(output_dir, "output.sam")
                bam_file = convert_sam_to_bam_and_index(sam_file, output_dir)

                process_kestrel_output(
                    output_dir, vcf_path, reference_vntr, kestrel_config, config
                )
                break


def process_kestrel_output(
    output_dir, vcf_path, reference_vntr, kestrel_config, config
):
    """
    Processes the Kestrel output VCF files, filters variants, fixes the VCF file format header,
    and writes the results to output files.

    This function performs the following steps:
    1. Filters the initial VCF to extract INDEL variants into `output_indel.vcf`.
    2. Fixes the incorrect file format line in `output_indel.vcf` from `##fileformat=VCF4.2` 
       to `##fileformat=VCFv4.2`.
    3. (New step) As a preparation for future integration, sorts and indexes the `output_indel.vcf` 
       using bcftools. However, this sorted/indexed VCF is currently not used to generate the 
       insertion/deletion VCFs for the rest of the pipeline. Instead, it is just produced and kept 
       for now, as requested.
    4. Filters the original (fixed-header) `output_indel.vcf` (not the sorted one) into insertion 
       and deletion VCFs, and proceeds with the pre-existing processing steps (motif annotation, 
       depth score calculation, etc.).
    5. Writes out intermediate and final TSV results.

    NOTE: The sorted and indexed VCF is created but not integrated further in the pipeline yet.
    Future steps might include using this sorted/indexed VCF for IGV.js session coverage improvements.

    Args:
        output_dir (str): Directory where Kestrel outputs are saved.
        vcf_path (Path): Path to the original Kestrel VCF output file.
        reference_vntr (str): Path to the reference VNTR file.
        kestrel_config (dict): Configuration dictionary for Kestrel settings.
        config (dict): Overall configuration dictionary.

    Returns:
        pd.DataFrame or None: The final processed DataFrame of variants, or None if no variants found.
    """
    logging.info("Processing Kestrel VCF results...")

    # Define file paths for intermediate and final files
    indel_vcf = os.path.join(output_dir, "output_indel.vcf")
    output_ins = os.path.join(output_dir, "output_insertion.vcf")
    output_del = os.path.join(output_dir, "output_deletion.vcf")

    # Filter the original VCF to extract indels into output_indel.vcf
    filter_vcf(vcf_path, indel_vcf)

    # Fix the fileformat line in the output_indel.vcf
    # Kestrel currently writes '##fileformat=VCF4.2', which should be '##fileformat=VCFv4.2'
    fixed_indel_vcf = indel_vcf + ".fixed"
    with open(indel_vcf, "r") as fin, open(fixed_indel_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##fileformat=VCF4.2"):
                fout.write("##fileformat=VCFv4.2\n")
            else:
                fout.write(line)
    # Replace the original indel VCF with the fixed version
    os.replace(fixed_indel_vcf, indel_vcf)

    # ------------------------------------------------------------------------
    # NEW STEP: Produce a sorted and indexed VCF for IGV integration
    # Note: We do NOT use this sorted/indexed VCF for the current pipeline steps, 
    #       just produce it and keep it in place for future integration.
    # ------------------------------------------------------------------------
    sorted_indel_vcf_gz = indel_vcf + ".gz"
    # Sort the VCF
    run_command(
        f"bcftools sort {indel_vcf} -o {sorted_indel_vcf_gz} -W -O z",
        log_file=os.path.join(output_dir, "bcftools_sort.log"),
    )
    # ------------------------------------------------------------------------

    # Filter the original (fixed) indel VCF to separate insertions and deletions
    # NOTE: We do NOT use the sorted and indexed file here to avoid breaking current logic.
    filter_indel_vcf(indel_vcf, output_ins, output_del)

    # Generate a custom header for result files
    header = generate_header(reference_vntr)

    def read_vcf_without_comments(vcf_file):
        """
        Utility function to read a VCF (possibly gzipped) without '##' comments, 
        returning a DataFrame of variants.
        """
        import gzip
        open_func = gzip.open if vcf_file.endswith(".gz") else open
        data = []
        header_line = None
        with open_func(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header_line = line.strip().split('\t')
                elif not line.startswith("##"):
                    data.append(line.strip().split('\t'))
        if data:
            return pd.DataFrame(data, columns=header_line)
        else:
            return pd.DataFrame()

    # Load insertion and deletion VCFs into DataFrames
    vcf_insertion = read_vcf_without_comments(output_ins)
    vcf_deletion = read_vcf_without_comments(output_del)

    if vcf_insertion.empty and vcf_deletion.empty:
        logging.warning(
            "No insertion or deletion variants found in the VCF. "
            "Skipping Kestrel processing."
        )
        output_empty_result(output_dir, header)
        return None

    # Load MUC1 reference VNTR motifs
    muc1_ref = load_muc1_reference(reference_vntr)

    # Preprocess insertion and deletion data
    insertion_df = (
        preprocessing_insertion(vcf_insertion, muc1_ref) 
        if not vcf_insertion.empty 
        else pd.DataFrame()
    )
    deletion_df = (
        preprocessing_deletion(vcf_deletion, muc1_ref) 
        if not vcf_deletion.empty 
        else pd.DataFrame()
    )

    # Combine insertion and deletion DataFrames
    combined_df = pd.concat([insertion_df, deletion_df], axis=0)
    # Sort combined results deterministically
    sort_columns = list(combined_df.columns)
    combined_df = combined_df.sort_values(by=sort_columns).reset_index(drop=True)

    if combined_df.empty:
        logging.warning(
            "Combined DataFrame of insertions and deletions is empty. "
            "Skipping further processing."
        )
        output_empty_result(output_dir, header)
        return None

    # Load additional motifs from the configuration
    merged_motifs = load_additional_motifs(config)

    # Process the combined DataFrame (calculate frame score, depth score, apply filters, etc.)
    processed_df = process_kmer_results(
        combined_df, merged_motifs, output_dir, kestrel_config
    )

    if processed_df.empty:
        logging.warning(
            "Final processed DataFrame is empty. Outputting empty results."
        )
        output_empty_result(output_dir, header)
        return None

    # Write the intermediate (pre-result) DataFrame to a TSV file
    pre_result_path = os.path.join(output_dir, "kestrel_pre_result.tsv")
    with open(pre_result_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        combined_df.to_csv(f, sep='\t', index=False)
    logging.info(f"Intermediate results saved as {pre_result_path}")

    # Write the final processed DataFrame to a TSV file
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")
    with open(final_output_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        processed_df.to_csv(f, sep='\t', index=False)

    logging.info("Kestrel VCF processing completed.")
    return processed_df


def output_empty_result(output_dir, header):
    """
    Creates an empty result file with the correct headers and a placeholder 'Negative' result row.

    Args:
        output_dir (str): Directory where the output file will be saved.
        header (list): List of header lines for the output file.
    """
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")

    # Create a DataFrame with one row containing "None" values and "Negative" in the Confidence column
    empty_result_data = {
        'Motif': ['None'],
        'Variant': ['None'],
        'POS': ['None'],
        'REF': ['None'],
        'ALT': ['None'],
        'Motif_sequence': ['None'],
        'Estimated_Depth_AlternateVariant': ['None'],
        'Estimated_Depth_Variant_ActiveRegion': ['None'],
        'Depth_Score': ['None'],
        'Confidence': ['Negative'],
    }
    empty_df = pd.DataFrame(empty_result_data)

    # Write the header and the empty DataFrame to the output file
    with open(final_output_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        empty_df.to_csv(f, sep='\t', index=False)

    logging.info(f"Empty result file with placeholder saved as {final_output_path}")


# Function 1: Split Depth and Calculate Frame Score
def split_depth_and_calculate_frame_score(df):
    """
    Splits the Depth column (Sample) into components and calculates the frame score.

    Args:
        df (pd.DataFrame): DataFrame containing variant information.

    Returns:
        pd.DataFrame: DataFrame with additional columns for depth and frame score.
    """
    if df.empty:
        return df

    # Remove the redundant logic of renaming 'Sample' to 'Depth'
    # Directly split the 'Sample' column
    df[['Del', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']] = df[
        'Sample'
    ].str.split(':', expand=True)

    # Select relevant columns
    df = df[
        [
            'Motifs',
            'Variant',
            'POS',
            'REF',
            'ALT',
            'Motif_sequence',
            'Estimated_Depth_AlternateVariant',
            'Estimated_Depth_Variant_ActiveRegion',
        ]
    ].copy()

    # Calculate reference and alternate allele lengths
    df["ref_len"] = df["REF"].str.len()
    df["alt_len"] = df["ALT"].str.len()

    # Calculate frame score
    df["Frame_Score"] = round((df["alt_len"] - df["ref_len"]) / 3, 2)
    df['Frame_Score'] = df['Frame_Score'].astype(str).apply(lambda x: x.replace('.0', 'C'))

    # Filter out non-frameshift variants
    df["TrueFalse"] = df['Frame_Score'].str.contains('C', regex=True)
    df = df[df["TrueFalse"] == False].copy()

    return df


# Function 2: Split Frame Score
def split_frame_score(df):
    """
    Splits the Frame_Score column into left and right parts and adjusts values.

    Args:
        df (pd.DataFrame): DataFrame with the 'Frame_Score' column.

    Returns:
        pd.DataFrame: DataFrame with 'left' and 'right' columns added.
    """
    if df.empty:
        return df

    # Split Frame_Score into left and right parts
    df[['left', 'right']] = df['Frame_Score'].str.split('.', expand=True)

    # Replace '-0' with '-1' in the 'left' column
    df['left'] = df['left'].replace('-0', '-1')

    # Drop unnecessary columns
    df.drop(['TrueFalse', 'ref_len', 'alt_len'], axis=1, inplace=True)

    return df


# Function 3: Extract Frameshifts
def extract_frameshifts(df):
    """
    Extracts insertion and deletion frameshifts based on 'left' and 'right' parts.

    Args:
        df (pd.DataFrame): DataFrame with 'left' and 'right' columns.

    Returns:
        pd.DataFrame: Filtered DataFrame with frameshift variants.
    """
    if df.empty:
        return df

    # Extract good frameshifts (3n+1 for insertion and 3n+2 for deletion)
    ins = df[
        df["left"].apply(lambda x: '-' not in x)
        & df["right"].apply(lambda y: '33' in y)
    ]
    del_ = df[
        df["left"].apply(lambda x: '-' in x)
        & df["right"].apply(lambda y: '67' in y)
    ]

    # Concatenate insertions and deletions
    frameshifts_df = pd.concat([ins, del_], axis=0)

    return frameshifts_df


# Function 4: Calculate Depth Score and Assign Confidence
def calculate_depth_score_and_assign_confidence(df, kestrel_config):
    """
    Calculates depth score and assigns confidence levels based on thresholds.

    Args:
        df (pd.DataFrame): DataFrame with depth information.
        kestrel_config (dict): Configuration dictionary for Kestrel settings.

    Returns:
        pd.DataFrame: DataFrame with 'Depth_Score' and 'Confidence' columns added.
    """
    if df.empty:
        return df

    # Convert depth-related columns to integers
    df['Estimated_Depth_AlternateVariant'] = df['Estimated_Depth_AlternateVariant'].astype(int)
    df['Estimated_Depth_Variant_ActiveRegion'] = df['Estimated_Depth_Variant_ActiveRegion'].astype(
        int
    )

    # Calculate depth score
    df['Depth_Score'] = (
        df['Estimated_Depth_AlternateVariant'] / df['Estimated_Depth_Variant_ActiveRegion']
    )

    # Extract thresholds and confidence levels from kestrel_config
    depth_score_thresholds = kestrel_config['confidence_assignment']['depth_score_thresholds']
    alt_depth_thresholds = kestrel_config['confidence_assignment']['alt_depth_thresholds']
    var_active_region_threshold = kestrel_config['confidence_assignment'][
        'var_active_region_threshold'
    ]
    confidence_levels = kestrel_config['confidence_assignment']['confidence_levels']

    # Define conditions for assigning confidence scores
    def assign_confidence(row):
        depth_score = row['Depth_Score']
        alt_depth = row['Estimated_Depth_AlternateVariant']
        var_active_region = row['Estimated_Depth_Variant_ActiveRegion']

        if (
            depth_score <= depth_score_thresholds['low']
            or var_active_region <= var_active_region_threshold
        ):
            return confidence_levels['low_precision']
        elif (
            alt_depth_thresholds['mid_low']
            <= alt_depth
            <= alt_depth_thresholds['mid_high']
            and depth_score_thresholds['low']
            <= depth_score
            <= depth_score_thresholds['high']
        ):
            return confidence_levels['low_precision']
        elif alt_depth > alt_depth_thresholds['mid_high']:
            return confidence_levels['high_precision']
        elif alt_depth <= alt_depth_thresholds['low']:
            return confidence_levels['low_precision']
        elif (
            alt_depth_thresholds['mid_low']
            <= alt_depth
            < alt_depth_thresholds['mid_high']
            and depth_score >= depth_score_thresholds['high']
        ):
            return confidence_levels['high_precision']
        elif (
            alt_depth >= alt_depth_thresholds['mid_high']
            and depth_score >= depth_score_thresholds['high']
        ):
            return confidence_levels['high_precision_star']
        else:
            return confidence_levels['low_precision']

    # Apply confidence score assignment
    df['Confidence'] = df.apply(assign_confidence, axis=1)

    return df


# Function 5: Filter by ALT Values and Finalize Data
def filter_by_alt_values_and_finalize(df, kestrel_config):
    """
    Filters DataFrame based on ALT values and finalizes data for output.

    Args:
        df (pd.DataFrame): DataFrame with variant information.
        kestrel_config (dict): Configuration dictionary for Kestrel settings.

    Returns:
        pd.DataFrame: Filtered and finalized DataFrame.
    """
    if df.empty:
        return df

    # Extract parameters from kestrel_config
    gg_alt_value = kestrel_config['alt_filtering']['gg_alt_value']
    gg_depth_score_threshold = kestrel_config['alt_filtering']['gg_depth_score_threshold']
    exclude_alts = kestrel_config['alt_filtering']['exclude_alts']

    # Filter based on specific ALT values (e.g., 'GG')
    if df['ALT'].str.contains(r'\b' + gg_alt_value + r'\b').any():
        gg_condition = df['ALT'] == gg_alt_value
        df = pd.concat(
            [
                df[~gg_condition],
                df[gg_condition & (df['Depth_Score'] >= gg_depth_score_threshold)],
            ]
        )

    # Exclude specified ALT values
    df = df[~df['ALT'].isin(exclude_alts)]

    # Drop unnecessary columns
    df.drop(['left', 'right'], axis=1, inplace=True)

    return df


# Function 6: Motif Correction and Annotation
def motif_correction_and_annotation(df, merged_motifs, kestrel_config):
    """
    Performs motif correction and annotates the DataFrame with motif information.

    Args:
        df (pd.DataFrame): DataFrame with variant information.
        merged_motifs (pd.DataFrame): DataFrame with merged motif information.
        kestrel_config (dict): Configuration dictionary for Kestrel settings.

    Returns:
        pd.DataFrame: Annotated DataFrame with motif corrections.
    """
    if df.empty:
        return df

    # Extract parameters from kestrel_config
    position_threshold = kestrel_config['motif_filtering']['position_threshold']
    exclude_motifs_right = kestrel_config['motif_filtering']['exclude_motifs_right']
    alt_for_motif_right_gg = kestrel_config['motif_filtering']['alt_for_motif_right_gg']
    motifs_for_alt_gg = kestrel_config['motif_filtering']['motifs_for_alt_gg']
    exclude_alts_combined = kestrel_config['motif_filtering']['exclude_alts_combined']
    exclude_motifs_combined = kestrel_config['motif_filtering']['exclude_motifs_combined']

    # Make a copy of the 'Motifs' column and rename it to 'Motif_fasta'
    if 'Motifs' in df.columns:
        df['Motif_fasta'] = df['Motifs']

    # Ensure that splitting will result in exactly two columns
    if df['Motifs'].str.count('-').max() == 1:
        df[['Motif_left', 'Motif_right']] = df['Motifs'].str.split('-', expand=True)
    else:
        logging.error("Unexpected format in 'Motifs' column during splitting.")
        return pd.DataFrame()  # Return an empty DataFrame or handle appropriately

    df['POS'] = df['POS'].astype(int)

    # Split into left and right motifs based on position
    motif_left = df[df['POS'] < position_threshold].copy()
    motif_right = df[df['POS'] >= position_threshold].copy()

    # Process the left motifs
    if not motif_left.empty:
        motif_left.rename(columns={'Motif_right': 'Motif'}, inplace=True)
        motif_left.drop(['Motif_sequence'], axis=1, inplace=True)
        motif_left = motif_left.merge(merged_motifs, on='Motif', how='left')
        motif_left = motif_left[
            [
                'Motif',
                'Motif_fasta',
                'Variant',
                'POS',
                'REF',
                'ALT',
                'Motif_sequence',
                'Estimated_Depth_AlternateVariant',
                'Estimated_Depth_Variant_ActiveRegion',
                'Depth_Score',
                'Confidence',
            ]
        ]
        motif_left = (
            motif_left.sort_values('Depth_Score', ascending=False)
            .drop_duplicates('ALT', keep='first')
            .sort_values('POS', ascending=False)
            .tail(1)
        )

    # Process the right motifs
    if not motif_right.empty:
        motif_right.rename(columns={'Motif_left': 'Motif'}, inplace=True)
        motif_right.drop(['Motif_sequence'], axis=1, inplace=True)
        motif_right = motif_right.merge(merged_motifs, on='Motif', how='left')
        motif_right = motif_right[
            [
                'Motif',
                'Motif_fasta',
                'Variant',
                'POS',
                'REF',
                'ALT',
                'Motif_sequence',
                'Estimated_Depth_AlternateVariant',
                'Estimated_Depth_Variant_ActiveRegion',
                'Depth_Score',
                'Confidence',
            ]
        ]

        if motif_right['ALT'].str.contains(r'\b' + alt_for_motif_right_gg + r'\b').any():
            motif_right = motif_right.loc[
                ~motif_right['Motif'].isin(exclude_motifs_right)
            ]
            motif_right = motif_right.loc[motif_right['ALT'] == alt_for_motif_right_gg]
            motif_right = (
                motif_right.sort_values('Depth_Score', ascending=False)
                .drop_duplicates('ALT', keep='first')
            )
            if motif_right['Motif'].isin(motifs_for_alt_gg).any():
                motif_right = motif_right[motif_right['Motif'].isin(motifs_for_alt_gg)]
        else:
            motif_right = (
                motif_right.sort_values('Depth_Score', ascending=False)
                .drop_duplicates('ALT', keep='first')
            )

        motif_right.drop_duplicates(subset=['REF', 'ALT'], inplace=True)

    # Combine the processed left and right motifs
    combined_df = pd.concat([motif_right, motif_left])

    # Additional filtering
    combined_df = combined_df.loc[~combined_df['ALT'].isin(exclude_alts_combined)]
    combined_df = combined_df[~combined_df['Motif'].isin(exclude_motifs_combined)]
    combined_df['POS'] = combined_df['POS'].astype(int)

    # Make a copy of the 'POS' column and rename it to 'POS_fasta'
    if 'POS' in combined_df.columns:
        combined_df['POS_fasta'] = combined_df['POS']

    # Adjust positions where necessary
    combined_df.update(
        combined_df['POS'].mask(
            combined_df['POS'] >= position_threshold, lambda x: x - position_threshold
        )
    )

    return combined_df


# Function 7: Generate BED File
def generate_bed_file(df, output_dir):
    """
    Generates a BED file from the processed Kestrel output DataFrame.

    Args:
        df (pd.DataFrame): Processed DataFrame with 'Motif_fasta' and 'POS_fasta' columns.
        output_dir (str): Directory to save the generated BED file.

    Returns:
        str or None: Path to the generated BED file, or None if the file wasn't created.
    """
    # Ensure the required columns are present
    if 'Motif_fasta' not in df.columns or 'POS_fasta' not in df.columns:
        logging.warning(
            "Missing 'Motif_fasta' or 'POS_fasta' columns in the DataFrame. Skipping BED file generation."
        )
        return None

    if df.empty:
        logging.warning("DataFrame is empty. No variants to generate a BED file.")
        return None

    bed_file_path = os.path.join(output_dir, "output.bed")

    with open(bed_file_path, 'w') as bed_file:
        for _, row in df.iterrows():
            motif_fasta = row['Motif_fasta']
            pos = row['POS_fasta']
            bed_file.write(f"{motif_fasta}\t{pos}\t{pos + 1}\n")

    logging.info(f"BED file generated at: {bed_file_path}")
    return bed_file_path


# Main Function: Process Kmer Results
def process_kmer_results(combined_df, merged_motifs, output_dir, kestrel_config):
    """
    Processes and filters Kestrel results by applying several steps including frame score calculation,
    depth score assignment, filtering based on ALT values and confidence scores, and final motif correction.

    Args:
        combined_df (pd.DataFrame): Combined DataFrame of insertion and deletion variants.
        merged_motifs (pd.DataFrame): DataFrame with merged motif information.
        output_dir (str): Directory where output files will be saved.
        kestrel_config (dict): Configuration dictionary for Kestrel settings.

    Returns:
        pd.DataFrame: Final processed DataFrame ready for output.
    """
    if combined_df.empty:
        return combined_df

    # Step 1: Split Depth and Calculate Frame Score
    df = split_depth_and_calculate_frame_score(combined_df)
    if df.empty:
        return df

    # Step 2: Split Frame Score
    df = split_frame_score(df)
    if df.empty:
        return df

    # Step 3: Extract Frameshifts
    df = extract_frameshifts(df)
    if df.empty:
        return df

    # Step 4: Calculate Depth Score and Assign Confidence
    df = calculate_depth_score_and_assign_confidence(df, kestrel_config)
    if df.empty:
        return df

    # Step 5: Filter by ALT Values and Finalize Data
    df = filter_by_alt_values_and_finalize(df, kestrel_config)
    if df.empty:
        return df

    # Step 6: Motif Correction and Annotation
    df = motif_correction_and_annotation(df, merged_motifs, kestrel_config)
    if df.empty:
        return df

    # Step 7: Generate BED File
    generate_bed_file(df, output_dir)

    return df
