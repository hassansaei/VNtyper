import logging
import pandas as pd
import numpy as np
import os
import subprocess as sp
from pathlib import Path
from vntyper.scripts.utils import run_command


def run_advntr(db_file_hg19, sorted_bam, output, output_name, config, output_format="tsv"):
    """
    Run adVNTR genotyping using the specified database file and BAM file, fetching settings from the config.

    Args:
        db_file_hg19 (str): Path to the adVNTR VNTR database file.
        sorted_bam (str): Path to the sorted BAM file.
        output (str): Directory where the results will be saved.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary.
        output_format (str): The format of the output file (e.g., "vcf", "tsv").
    """
    advntr_path = config["tools"]["advntr"]
    advntr_settings = config["advntr_settings"]

    # Set the number of threads from the config or default to 8
    threads = advntr_settings.get("threads", 8)

    # Retrieve additional command parts from the config, if available
    additional_commands = advntr_settings.get("additional_commands", "-aln")

    # Determine the output file extension based on the specified format
    output_ext = ".vcf" if output_format == "vcf" else ".tsv"
    
    advntr_command = (
        f"{advntr_path} genotype -fs -vid {advntr_settings['vid']} "
        f"--alignment_file {sorted_bam} -o {output}/{output_name}_adVNTR{output_ext} "
        f"-m {db_file_hg19} --working_directory {output} -t {threads} {additional_commands}"
    )

    # Define log files for stdout and stderr
    log_file = os.path.join(output, f"{output_name}_advntr.log")

    logging.info("Launching adVNTR genotyping!")

    try:
        # Run the adVNTR command and log output to the specified log file
        if not run_command(advntr_command, log_file, critical=True):
            logging.error("adVNTR genotyping failed. Check the log for details.")
            # Do not raise an exception here, allowing the pipeline to continue
    except Exception as e:
        logging.error(f"adVNTR genotyping encountered an error: {e}")
        # Continue without raising an exception

    logging.info("adVNTR genotyping of MUC1-VNTR completed!")


def advntr_processing_del(df, config):
    """
    Process adVNTR deletions by calculating deletion length and filtering by frameshift.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.
        config (dict): Configuration dictionary containing advntr settings.

    Returns:
        pd.DataFrame: Filtered DataFrame containing deletions.
    """
    # Create a copy of the original DataFrame to avoid modifying it directly
    df1 = df.copy()

    # Rename the 'State' column to 'Variant' and 'Pvalue\n' to 'Pvalue' for clarity
    df1.rename(columns={'State': 'Variant', 'Pvalue\n': 'Pvalue'}, inplace=True)

    # Calculate the length of deletions in the 'Variant' column by counting occurrences of 'D'
    df1['Deletion_length'] = df1['Variant'].str.count('D')

    # Calculate the length of insertions in the 'Variant' column by counting occurrences of 'I'
    df1['Insertion'] = df1['Variant'].str.count('I')

    # Extract the insertion length from the 'Variant' column, handling cases with missing values
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1[['I', 'Insertion_len']] = df1['Insertion_len'].str.split('LEN', expand=True)
    df1['Insertion_len'] = df1['Insertion_len'].replace('', 0).astype(int)

    # Convert the 'Deletion_length' column to integer type
    df1['Deletion_length'] = df1['Deletion_length'].astype(int)

    # Calculate the frameshift value as the absolute difference between insertion and deletion lengths
    df1['frame'] = abs(df1['Insertion_len'] - df1['Deletion_length']).astype(str)

    # Define valid frameshift patterns based on settings from the config file
    max_frameshift = config["advntr_settings"]["max_frameshift"]
    frameshift_multiplier = config["advntr_settings"]["frameshift_multiplier"]
    del_frame = (np.arange(max_frameshift) * frameshift_multiplier + 2).astype(str)

    # Filter the DataFrame to keep only rows with valid deletions and frameshift values
    df1 = df1[(df1['Deletion_length'] >= 1) & df1['frame'].isin(del_frame)]

    # Return the filtered DataFrame containing deletions
    return df1


def advntr_processing_ins(df, config):
    """
    Process adVNTR insertions by calculating insertion length and filtering by frameshift.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.
        config (dict): Configuration dictionary containing advntr settings.

    Returns:
        pd.DataFrame: Filtered DataFrame containing insertions.
    """
    # Create a copy of the original DataFrame to avoid modifying it directly
    df1 = df.copy()

    # Rename the 'State' column to 'Variant' and 'Pvalue\n' to 'Pvalue' for clarity
    df1.rename(columns={'State': 'Variant', 'Pvalue\n': 'Pvalue'}, inplace=True)

    # Calculate the length of deletions in the 'Variant' column by counting occurrences of 'D'
    df1['Deletion_length'] = df1['Variant'].str.count('D')

    # Calculate the length of insertions in the 'Variant' column by counting occurrences of 'I'
    df1['Insertion'] = df1['Variant'].str.count('I')

    # Extract the insertion length from the 'Variant' column, handling cases with missing values
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1[['I', 'Insertion_len']] = df1['Insertion_len'].str.split('LEN', expand=True)
    df1['Insertion_len'] = df1['Insertion_len'].replace('', 0).astype(int)

    # Convert the 'Deletion_length' column to integer type
    df1['Deletion_length'] = df1['Deletion_length'].astype(int)

    # Calculate the frameshift value as the absolute difference between insertion and deletion lengths
    df1['frame'] = abs(df1['Insertion_len'] - df1['Deletion_length']).astype(str)

    # Define valid frameshift patterns based on settings from the config file
    max_frameshift = config["advntr_settings"]["max_frameshift"]
    frameshift_multiplier = config["advntr_settings"]["frameshift_multiplier"]
    ins_frame = (np.arange(max_frameshift) * frameshift_multiplier + 1).astype(str)

    # Filter the DataFrame to keep only rows with valid insertions and frameshift values
    df1 = df1[(df1['Insertion_len'] >= 1) & df1['frame'].isin(ins_frame)]

    # Return the filtered DataFrame containing insertions
    return df1


def process_advntr_output(output_path, output, output_name, config, file_format="tsv"):
    """
    Process the adVNTR output TSV to extract relevant information and generate final results.
    
    Args:
        output_path (str): Path to the adVNTR output file.
        output (str): Directory where the final results will be saved.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary.
        file_format (str): The format of the output file ("vcf" or "tsv").
    """
    if not os.path.exists(output_path):
        logging.error(f"adVNTR output file {output_path} not found!")
        return

    logging.info('Processing adVNTR result...')

    # Read the output file
    with open(output_path, 'r') as file:
        content = file.readlines()

    # Define default column names that match the expected structure of the adVNTR output
    default_columns_advntr_output = ['#VID', 'State', 'NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']

    # Check if the file is empty or only contains comments
    if not content or all(line.startswith('#') for line in content):
        logging.warning('No pathogenic variant was found with adVNTR!')

        # Write the header and the Negative line if the file is empty or only contains comments
        with open(output_path, 'w') as f:
            f.write(f'#Input File: {output_name} \n')
            f.write('#Reference file: None\n')
            f.write('#P-value cutoff: 0.001\n')
            f.write('\t'.join(default_columns_advntr_output) + '\n')
            f.write('Negative\t' + '\t'.join(['None'] * (len(default_columns_advntr_output) - 1)) + '\n')

        logging.info(f'Empty or comment-only output found, written header and negative line to {output_path}')
        return

    # Read the output file again to update the header
    with open(output_path, 'r') as file:
        content = file.readlines()

    # Replace #VID with VID in the header line to prevent it from being treated as a comment
    content = [line.replace('#VID', 'VID') if line.startswith('#VID') else line for line in content]

    # Write the modified content back to the file
    with open(output_path, 'w') as file:
        file.writelines(content)

    try:
        # Load the data using pandas, ensuring the correct header is used and ignoring comments
        logging.info('Loading data into DataFrame...')
        df = pd.read_csv(output_path, sep='\t', comment='#')
        logging.info(f"Data loaded successfully with shape: {df.shape}")

        # Log the first few rows of the DataFrame to verify the data
        logging.debug(f"First few rows of the DataFrame:\n{df.head()}")

    except Exception as e:
        logging.error(f"Error loading data into DataFrame: {e}")
        return

    try:
        # Process deletions and insertions using the respective functions
        logging.info('Processing deletions...')
        df_del = advntr_processing_del(df, config)

        logging.info('Processing insertions...')
        df_ins = advntr_processing_ins(df, config)

        # Concatenate the processed results for deletions and insertions
        logging.info('Concatenating deletions and insertions...')
        advntr_concat = pd.concat([df_del, df_ins], axis=0)

        # Ensure only the relevant columns are kept and processed for the final output
        default_columns_advntr_final = ['VID', 'Variant', 'NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']
        advntr_concat = advntr_concat[default_columns_advntr_final]

        # Remove duplicate records based on specific columns
        logging.info('Removing duplicates...')
        advntr_concat.drop_duplicates(subset=['VID', 'Variant', 'NumberOfSupportingReads'], inplace=True)

        # Save the processed adVNTR results to a TSV file
        output_result_path = os.path.join(output, f"{output_name}_adVNTR_result.tsv")
        advntr_concat.to_csv(output_result_path, sep='\t', index=False)
        logging.info(f"Processed adVNTR results saved to {output_result_path}")

    except Exception as e:
        logging.error(f"Error during processing of deletions and insertions: {e}")
        return

    # CLean up intermediate files as specified in the configuration
    # currently this is not used
    cleanup_files(output, output_name, config.get("advntr_cleanup_files", []))


def cleanup_files(output, output_name, files_to_remove):
    """
    Clean up intermediate files as specified in the configuration.

    Args:
        output (str): The output directory.
        output_name (str): The base name for the output files.
        files_to_remove (list): List of files to remove.
    """
    for db in files_to_remove:
        db_str = os.path.join(output, f"{output_name}{db}")
        if os.path.exists(db_str):
            rm_command = f"rm {db_str}"
            process = sp.Popen(rm_command, shell=True)
            process.wait()

    logging.info('Intermediate files cleaned up.')
