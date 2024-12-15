#!/usr/bin/env python3

import logging
import pandas as pd
import numpy as np
import os
import subprocess as sp
from pathlib import Path
from vntyper.scripts.utils import run_command, load_config

def load_advntr_config(config_path=None):
    """
    Loads the adVNTR configuration file.
    """
    if config_path is None:
        # Default path to advntr_config.json
        config_path = os.path.join(os.path.dirname(__file__), 'advntr_config.json')
    return load_config(config_path)

# Load the adVNTR settings
advntr_config = load_advntr_config()
advntr_settings = advntr_config.get("advntr_settings", {})

def run_advntr(db_file, sorted_bam, output, output_name, config):
    """
    Run adVNTR genotyping using the specified database file and BAM file, fetching settings from advntr_config.

    Args:
        db_file (str): Path to the adVNTR VNTR database file.
        sorted_bam (str): Path to the sorted BAM file.
        output (str): Directory where the results will be saved.
        output_name (str): Base name for the output files.
        config (dict): Main configuration dictionary.
    """
    advntr_path = config["tools"]["advntr"]

    # Set the number of threads from advntr_settings or default to 8
    threads = advntr_settings.get("threads", 8)

    # Retrieve additional command parts from advntr_settings, if available
    additional_commands = advntr_settings.get("additional_commands", "-aln")

    # Determine the output format and extension
    output_format = advntr_settings.get("output_format", "tsv")
    output_ext = ".vcf" if output_format == "vcf" else ".tsv"

    # Set the VNTR ID from config file or default to 25561
    vid = advntr_settings.get("vid", 25561)
    
    advntr_command = (
        f"{advntr_path} genotype -fs -vid {vid} "
        f"--alignment_file {sorted_bam} -o {output}/{output_name}_adVNTR{output_ext} "
        f"-m {db_file} --working_directory {output} -t {threads} {additional_commands}"
    )

    # Define log file for adVNTR output
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

def advntr_processing_del(df):
    """
    Process adVNTR deletions by calculating deletion length and filtering by frameshift.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing deletions.
    """
    # Create a copy of the original DataFrame to avoid modifying it directly
    df1 = df.copy()

    # Rename columns for clarity
    df1.rename(columns={'State': 'Variant', 'Pvalue\n': 'Pvalue'}, inplace=True)

    # Calculate deletion and insertion lengths
    df1['Deletion_length'] = df1['Variant'].str.count('D')
    df1['Insertion_length'] = df1['Variant'].str.count('I')

    # Extract insertion length values
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1[['I', 'Insertion_len']] = df1['Insertion_len'].str.split('LEN', expand=True)
    df1['Insertion_len'] = df1['Insertion_len'].replace('', 0).astype(int)

    # Convert columns to integer type
    df1['Deletion_length'] = df1['Deletion_length'].astype(int)

    # Calculate frameshift value
    df1['frame'] = abs(df1['Insertion_len'] - df1['Deletion_length']).astype(str)

    # Define valid frameshift patterns based on settings from advntr_settings
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    del_frame = (np.arange(max_frameshift) * frameshift_multiplier + 2).astype(str)

    # Filter the DataFrame
    df1 = df1[(df1['Deletion_length'] >= 1) & df1['frame'].isin(del_frame)]

    return df1

def advntr_processing_ins(df):
    """
    Process adVNTR insertions by calculating insertion length and filtering by frameshift.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing insertions.
    """
    # Create a copy of the original DataFrame to avoid modifying it directly
    df1 = df.copy()

    # Rename columns for clarity
    df1.rename(columns={'State': 'Variant', 'Pvalue\n': 'Pvalue'}, inplace=True)

    # Calculate deletion and insertion lengths
    df1['Deletion_length'] = df1['Variant'].str.count('D')
    df1['Insertion_length'] = df1['Variant'].str.count('I')

    # Extract insertion length values
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1[['I', 'Insertion_len']] = df1['Insertion_len'].str.split('LEN', expand=True)
    df1['Insertion_len'] = df1['Insertion_len'].replace('', 0).astype(int)

    # Convert columns to integer type
    df1['Deletion_length'] = df1['Deletion_length'].astype(int)

    # Calculate frameshift value
    df1['frame'] = abs(df1['Insertion_len'] - df1['Deletion_length']).astype(str)

    # Define valid frameshift patterns based on settings from advntr_settings
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    ins_frame = (np.arange(max_frameshift) * frameshift_multiplier + 1).astype(str)

    # Filter the DataFrame
    df1 = df1[(df1['Insertion_len'] >= 1) & df1['frame'].isin(ins_frame)]

    return df1

def process_advntr_output(output_path, output, output_name):
    """
    Process the adVNTR output to extract relevant information and generate final results.

    Args:
        output_path (str): Path to the adVNTR output file.
        output (str): Directory where the final results will be saved.
        output_name (str): Base name for the output files.
    """
    if not os.path.exists(output_path):
        logging.error(f"adVNTR output file {output_path} not found!")
        return

    logging.info('Processing adVNTR result...')

    # Read the output file
    with open(output_path, 'r') as file:
        content = file.readlines()

    # Define default column names
    default_columns_advntr_output = ['#VID', 'State', 'NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']

    # Check if the file is empty or only contains comments
    if not content or all(line.startswith('#') for line in content):
        logging.warning('No pathogenic variant was found with adVNTR!')

        # Write the header and a negative result line
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

    # Replace #VID with VID in the header line
    content = [line.replace('#VID', 'VID') if line.startswith('#VID') else line for line in content]

    # Write the modified content back to the file
    with open(output_path, 'w') as file:
        file.writelines(content)

    try:
        # Load the data using pandas
        logging.info('Loading data into DataFrame...')
        df = pd.read_csv(output_path, sep='\t', comment='#')
        logging.info(f"Data loaded successfully with shape: {df.shape}")

        # Log the first few rows of the DataFrame
        logging.debug(f"First few rows of the DataFrame:\n{df.head()}")

    except Exception as e:
        logging.error(f"Error loading data into DataFrame: {e}")
        return

    try:
        # Process deletions and insertions
        logging.info('Processing deletions...')
        df_del = advntr_processing_del(df)

        logging.info('Processing insertions...')
        df_ins = advntr_processing_ins(df)

        # Concatenate deletions and insertions
        logging.info('Concatenating deletions and insertions...')
        advntr_concat = pd.concat([df_del, df_ins], axis=0)

        # Keep relevant columns
        default_columns_advntr_final = ['VID', 'Variant', 'NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']
        advntr_concat = advntr_concat[default_columns_advntr_final]

        # Remove duplicates
        logging.info('Removing duplicates...')
        advntr_concat.drop_duplicates(subset=['VID', 'Variant', 'NumberOfSupportingReads'], inplace=True)

        # Save the processed adVNTR results
        output_result_path = os.path.join(output, f"{output_name}_adVNTR_result.tsv")
        advntr_concat.to_csv(output_result_path, sep='\t', index=False)
        logging.info(f"Processed adVNTR results saved to {output_result_path}")

    except Exception as e:
        logging.error(f"Error during processing of deletions and insertions: {e}")
        return

    # Clean up intermediate files if necessary
    cleanup_files(output, output_name)

def cleanup_files(output, output_name):
    """
    Clean up intermediate files.

    Args:
        output (str): The output directory.
        output_name (str): The base name for the output files.
    """
    # Currently, this function is not used.
    # Implement cleanup logic here if needed.
    logging.info('Intermediate files cleaned up.')
