import logging
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from datetime import datetime
from vntyper.scripts.utils import run_command
from vntyper.scripts.file_processing import filter_vcf, filter_indel_vcf
from vntyper.scripts.motif_processing import load_muc1_reference, load_additional_motifs, preprocessing_insertion, preprocessing_deletion
from vntyper.version import __version__ as VERSION


def construct_kestrel_command(kmer_size, kestrel_path, reference_vntr, output_dir, fastq_1, fastq_2, vcf_out, java_path, java_memory, max_align_states, max_hap_states):
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

    Returns:
        str: The constructed Kestrel command.
    """
    if not fastq_1 or not fastq_2:
        raise ValueError("FASTQ input files are missing or invalid.")

    return (
        f"{java_path} -Xmx{java_memory} -jar {kestrel_path} -k {kmer_size} "
        f"--maxalignstates {max_align_states} --maxhapstates {max_hap_states} "
        f"-r {reference_vntr} -o {vcf_out} "
        f"{fastq_1} {fastq_2} "
        f"--hapfmt sam -p {output_dir}/output.sam --logstderr --logstdout --temploc {output_dir}"
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
        f"## Reference file: {reference_vntr}"
    ]
    return header


def run_kestrel(vcf_path, output_dir, fastq_1, fastq_2, reference_vntr, kestrel_path, kestrel_settings, config):
    """
    Orchestrates the Kestrel genotyping process by iterating through kmer sizes and processing the VCF output.

    Args:
        vcf_path (Path): Path to the output VCF file.
        output_dir (str): Directory where Kestrel outputs will be saved.
        fastq_1 (str): Path to the first FASTQ file.
        fastq_2 (str): Path to the second FASTQ file.
        reference_vntr (str): Path to the reference VNTR file.
        kestrel_path (str): Path to the Kestrel jar file.
        kestrel_settings (dict): Settings for Kestrel.
        config (dict): Overall configuration dictionary.
    """
    java_path = kestrel_settings.get("java_path", "java")
    java_memory = kestrel_settings.get("java_memory", "15g")
    kmer_sizes = kestrel_settings.get("kmer_sizes", [20, 17, 25, 41])
    max_align_states = kestrel_settings.get("max_align_states", 30)
    max_hap_states = kestrel_settings.get("max_hap_states", 30)

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
            max_hap_states=max_hap_states
        )

        log_file = os.path.join(output_dir, f"kestrel_kmer_{kmer_size}.log")

        if vcf_path.is_file():
            logging.info("VCF file already exists, skipping Kestrel run...")
            return
        else:
            logging.info(f"Launching Kestrel with kmer size {kmer_size}...")

            # Run the command and log output to the specified log file
            if not run_command(kmer_command, log_file, critical=True):
                logging.error(f"Kestrel failed for kmer size {kmer_size}. Check {log_file} for details.")
                raise RuntimeError(f"Kestrel failed for kmer size {kmer_size}.")

            logging.info(f"Mapping-free genotyping of MUC1-VNTR with kmer size {kmer_size} done!")

            if vcf_path.is_file():
                process_kestrel_output(output_dir, vcf_path, reference_vntr, config)  # Pass config here
                break  # If successful, break out of the loop and stop trying other kmer sizes


def process_kestrel_output(output_dir, vcf_path, reference_vntr, config):
    logging.info("Processing Kestrel VCF results...")

    indel_vcf = os.path.join(output_dir, "output_indel.vcf")
    output_ins = os.path.join(output_dir, "output_insertion.vcf")
    output_del = os.path.join(output_dir, "output_deletion.vcf")

    # Filter the VCF to extract indels, insertions, and deletions
    filter_vcf(vcf_path, indel_vcf)
    filter_indel_vcf(indel_vcf, output_ins, output_del)

    # Generate the header
    header = generate_header(reference_vntr)

    # Manually read the VCF file to remove '##' comments and retain the header
    def read_vcf_without_comments(vcf_file):
        data = []
        header = None
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().split('\t')
                elif not line.startswith("##"):
                    data.append(line.strip().split('\t'))
        if data:  # Ensure that there is data before creating a DataFrame
            return pd.DataFrame(data, columns=header)
        else:
            return pd.DataFrame()  # Return an empty DataFrame if no data is found

    # Load the filtered VCF files into DataFrames
    vcf_insertion = read_vcf_without_comments(output_ins)
    vcf_deletion = read_vcf_without_comments(output_del)

    if vcf_insertion.empty and vcf_deletion.empty:
        logging.warning("No insertion or deletion variants found in the VCF. Skipping Kestrel processing.")
        output_empty_result(output_dir, header)
        return None  # Early exit if no data

    # Load MUC1 VNTR reference motifs
    MUC1_ref = load_muc1_reference(reference_vntr)

    # Preprocess insertion and deletion dataframes
    insertion_df = preprocessing_insertion(vcf_insertion, MUC1_ref) if not vcf_insertion.empty else pd.DataFrame()
    deletion_df = preprocessing_deletion(vcf_deletion, MUC1_ref) if not vcf_deletion.empty else pd.DataFrame()

    # Combine insertion and deletion dataframes
    combined_df = pd.concat([insertion_df, deletion_df], axis=0)

    if combined_df.empty:
        logging.warning("Combined DataFrame of insertions and deletions is empty. Skipping further processing.")
        output_empty_result(output_dir, header)
        return None  # Early exit if no data

    # Load additional motifs from the configuration
    merged_motifs = load_additional_motifs(config)

    # Process and filter results based on frameshifts and confidence scores, with motif correction
    processed_df = process_kmer_results(combined_df, merged_motifs)

    if processed_df.empty:
        logging.warning("Final processed DataFrame is empty. Outputting empty results.")
        output_empty_result(output_dir, header)
        return None  # Early exit if no data

    # Save the intermediate pre-result as `_pre_result.tsv`
    pre_result_path = os.path.join(output_dir, "kestrel_pre_result.tsv")
    with open(pre_result_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        combined_df.to_csv(f, sep='\t', index=False)
    logging.info(f"Intermediate results saved as {pre_result_path}")

    # Save the final processed dataframe as `kestrel_result.tsv`
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")
    with open(final_output_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        processed_df.to_csv(f, sep='\t', index=False)

    logging.info("Kestrel VCF processing completed.")
    return processed_df


def output_empty_result(output_dir, header):
    """
    Creates an empty result file with the correct headers.
    """
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")
    empty_df = pd.DataFrame(columns=[
        'Motifs', 'POS', 'REF', 'ALT', 'Variant', 'Motif_sequence',
        'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion',
        'Depth_Score', 'Confidence'
    ])
    with open(final_output_path, 'w') as f:
        f.write("\n".join(header) + "\n")
        empty_df.to_csv(f, sep='\t', index=False)
    logging.info(f"Empty result file saved as {final_output_path}")


# Function 1: Split Depth and Calculate Frame Score
def split_depth_and_calculate_frame_score(df):
    """
    Splits the Depth column (Sample) into Del, Estimated_Depth_AlternateVariant, and Estimated_Depth_Variant_ActiveRegion.
    Calculates the frame score and filters out non-frameshift variants.
    """
    if df.empty:
        return df

    # Rename 'Sample' column to 'Depth'
    df = df.rename(columns={'Sample': 'Depth'})

    # Split the 'Depth' column into components
    df[['Del', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']] = df['Depth'].str.split(':', expand=True)

    # Select relevant columns
    df = df[['Motifs', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']]

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
    Splits the Frame_Score column into left and right parts. Adjusts the left column values.
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
    Extracts insertion and deletion frameshifts based on the left and right parts of the Frame_Score.
    """
    if df.empty:
        return df

    # Extract good frameshifts (3n+1 for insertion and 3n+2 for deletion)
    ins = df[df["left"].apply(lambda x: '-' not in x) & df["right"].apply(lambda y: '33' in y)]
    del_ = df[df["left"].apply(lambda x: '-' in x) & df["right"].apply(lambda y: '67' in y)]

    # Concatenate insertions and deletions
    frameshifts_df = pd.concat([ins, del_], axis=0)

    return frameshifts_df


# Function 4: Calculate Depth Score and Assign Confidence
def calculate_depth_score_and_assign_confidence(df):
    """
    Calculates depth score and assigns confidence based on depth and variant region conditions.
    """
    if df.empty:
        return df

    # Convert depth-related columns to integers
    df['Estimated_Depth_AlternateVariant'] = df['Estimated_Depth_AlternateVariant'].astype(int)
    df['Estimated_Depth_Variant_ActiveRegion'] = df['Estimated_Depth_Variant_ActiveRegion'].astype(int)

    # Calculate depth score
    df['Depth_Score'] = df['Estimated_Depth_AlternateVariant'] / df['Estimated_Depth_Variant_ActiveRegion']

    # Define conditions for assigning confidence scores
    def assign_confidence(row):
        depth_score = row['Depth_Score']
        alt_depth = row['Estimated_Depth_AlternateVariant']
        var_active_region = row['Estimated_Depth_Variant_ActiveRegion']

        if depth_score <= 0.00469 or var_active_region <= 200:
            return 'Low_Precision'
        elif 21 <= alt_depth <= 100 and 0.00469 <= depth_score <= 0.00515:
            return 'Low_Precision'
        elif alt_depth > 100:
            return 'High_Precision'
        elif alt_depth <= 20:
            return 'Low_Precision'
        elif 21 <= alt_depth < 100 and depth_score >= 0.00515:
            return 'High_Precision'
        elif alt_depth >= 100 and depth_score >= 0.00515:
            return 'High_Precision*'
        else:
            return 'Low_Precision'

    # Apply confidence score assignment
    df['Confidence'] = df.apply(assign_confidence, axis=1)

    return df


# Function 5: Filter by ALT Values and Finalize Data
def filter_by_alt_values_and_finalize(df):
    """
    Filters based on specific ALT values and finalizes the dataframe.
    """
    if df.empty:
        return df

    # Filter based on specific ALT values (e.g., 'GG')
    if df['ALT'].str.contains(r'\bGG\b').any():
        gg_condition = df['ALT'] == 'GG'
        df = pd.concat([
            df[~gg_condition],
            df[gg_condition & (df['Depth_Score'] >= 0.00469)]
        ])

    # Filter out specific ALT sequences (CG, TG)
    df = df[~df['ALT'].isin(['CG', 'TG'])]

    # Drop unnecessary columns
    df.drop(['left', 'right'], axis=1, inplace=True)

    # Keep only high precision results not in the red zone
    if 'Confidence' in df.columns:
        df = df[df['Confidence'] != 'Red_Zone']

    return df


# Function 6: Motif Correction and Annotation
def motif_correction_and_annotation(df, merged_motifs):
    if df.empty:
        return df

    # Ensure that splitting will result in exactly two columns
    if df['Motifs'].str.count('-').max() == 1:
        df[['Motif_left', 'Motif_right']] = df['Motifs'].str.split('-', expand=True)
    else:
        logging.error("Unexpected format in 'Motifs' column during splitting.")
        return pd.DataFrame()  # Return an empty DataFrame or handle appropriately

    df['POS'] = df['POS'].astype(int)

    # Split into left and right motifs based on position
    motif_left = df[df['POS'] < 60].copy()  # Explicitly create a copy
    motif_right = df[df['POS'] >= 60].copy()  # Explicitly create a copy

    # Process the left motifs
    if not motif_left.empty:
        motif_left.rename(columns={'Motif_right': 'Motif'}, inplace=True)
        motif_left.drop(['Motif_sequence'], axis=1, inplace=True)
        motif_left = motif_left.merge(merged_motifs, on='Motif', how='left')
        motif_left = motif_left[['Motif', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence',
                                 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion',
                                 'Depth_Score', 'Confidence']]
        motif_left = motif_left.sort_values('Depth_Score', ascending=False).drop_duplicates('ALT', keep='first')
        motif_left = motif_left.sort_values('POS', ascending=False).tail(1)

    # Process the right motifs
    if not motif_right.empty:
        motif_right.rename(columns={'Motif_left': 'Motif'}, inplace=True)
        motif_right.drop(['Motif_sequence'], axis=1, inplace=True)
        motif_right = motif_right.merge(merged_motifs, on='Motif', how='left')
        motif_right = motif_right[['Motif', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence',
                                   'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion',
                                   'Depth_Score', 'Confidence']]

        if motif_right['ALT'].str.contains(r'\bGG\b').any():
            motif_right = motif_right.loc[(motif_right['Motif'] != 'Q') & (motif_right['Motif'] != '8') &
                                          (motif_right['Motif'] != '9') & (motif_right['Motif'] != '7') &
                                          (motif_right['Motif'] != '6p') & (motif_right['Motif'] != '6') &
                                          (motif_right['Motif'] != 'V') & (motif_right['Motif'] != 'J') &
                                          (motif_right['Motif'] != 'I') & (motif_right['Motif'] != 'G') &
                                          (motif_right['Motif'] != 'E') & (motif_right['Motif'] != 'A')]
            motif_right = motif_right.loc[motif_right['ALT'] == 'GG']
            motif_right = motif_right.sort_values('Depth_Score', ascending=False).drop_duplicates('ALT', keep='first')
            if motif_right['Motif'].str.contains('X').any():
                motif_right = motif_right[motif_right['Motif'] == 'X']
        else:
            motif_right = motif_right.sort_values('Depth_Score', ascending=False).drop_duplicates('ALT', keep='first')

        motif_right.drop_duplicates(subset=['REF', 'ALT'], inplace=True)

    # Combine the processed left and right motifs
    combined_df = pd.concat([motif_right, motif_left])

    # Additional filtering
    combined_df = combined_df.loc[(combined_df['ALT'] != 'CCGCC') & (combined_df['ALT'] != 'CGGCG') &
                                  (combined_df['ALT'] != 'CGGCC')]
    combined_df = combined_df[(combined_df['Motif'] != '6') & (combined_df['Motif'] != '6p') &
                              (combined_df['Motif'] != '7')]
    combined_df['POS'] = combined_df['POS'].astype(int)

    # Adjust positions where necessary
    combined_df.update(combined_df['POS'].mask(combined_df['POS'] >= 60, lambda x: x - 60))

    return combined_df


# Main Function: Process Kmer Results
def process_kmer_results(combined_df, merged_motifs):
    """
    Processes and filters Kestrel results by applying several steps including frame score calculation,
    depth score assignment, filtering based on ALT values and confidence scores, and final motif correction.
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
    df = calculate_depth_score_and_assign_confidence(df)
    if df.empty:
        return df

    # Step 5: Filter by ALT Values and Finalize Data
    df = filter_by_alt_values_and_finalize(df)
    if df.empty:
        return df

    # Step 6: Motif Correction and Annotation
    df = motif_correction_and_annotation(df, merged_motifs)

    return df