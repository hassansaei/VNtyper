#!/usr/bin/env python3
# vntyper/modules/advntr/advntr_genotyping.py

import logging
import os
import subprocess as sp
from pathlib import Path

import numpy as np
import pandas as pd

from vntyper.scripts.utils import run_command, load_config

# -------------------------------------------------------------------------
# Configure logging
# -------------------------------------------------------------------------
# You can adjust the logging level here to control the verbosity.
# For production, you might set this to logging.INFO or logging.WARNING.
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
)


def load_advntr_config(config_path=None):
    """
    Loads the adVNTR configuration file.
    """
    if config_path is None:
        # Default path to advntr_config.json
        config_path = os.path.join(os.path.dirname(__file__), "advntr_config.json")
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

    Returns:
        int: Return code indicating success (0) or failure (non-zero).
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

    # ---------------------------------------------------------------------
    # Validate input paths before proceeding
    # ---------------------------------------------------------------------
    if not os.path.isfile(db_file):
        logging.critical(f"VNTR database file not found: {db_file}")
        return 1
    if not os.path.isfile(sorted_bam):
        logging.critical(f"Sorted BAM file not found: {sorted_bam}")
        return 1
    if not os.path.isdir(output):
        logging.warning(f"Output directory does not exist, creating: {output}")
        try:
            os.makedirs(output, exist_ok=True)
        except Exception as e:
            logging.critical(f"Could not create output directory {output}: {e}")
            return 1

    advntr_command = (
        f"{advntr_path} genotype -fs -vid {vid} "
        f"--alignment_file {sorted_bam} -o {output}/{output_name}_adVNTR{output_ext} "
        f"-m {db_file} --working_directory {output} -t {threads} {additional_commands}"
    )

    # Define log file for adVNTR output
    log_file = os.path.join(output, f"{output_name}_advntr.log")

    logging.info("Launching adVNTR genotyping...")
    logging.debug(f"Command: {advntr_command}")

    try:
        # Run the adVNTR command and log output to the specified log file
        if not run_command(advntr_command, log_file, critical=True):
            # run_command() returned False, indicating an error
            logging.error("adVNTR genotyping failed. Check the log for details.")
            # Pipeline will continue, but we return a non-zero exit code
            return 1
    except sp.CalledProcessError as cpe:
        logging.error(f"adVNTR genotyping CalledProcessError: {cpe}")
        return 1
    except Exception as e:
        logging.error(f"adVNTR genotyping encountered an unexpected error: {e}")
        return 1

    logging.info("adVNTR genotyping of MUC1-VNTR completed successfully.")
    return 0


def advntr_processing_del(df):
    """
    Process adVNTR deletions by calculating deletion length, computing the frameshift,
    and filtering variants based on valid frameshift patterns.

    This function:
      - Renames columns for clarity.
      - Calculates the number of deletion ('D') and insertion ('I') characters in the entire
        variant string. Note that multi-allelic strings (e.g. "D55_6&D56_6&D57_6&D58_6&D59_6&D60_6")
        are treated as a single value; the function aggregates the counts.
      - Extracts the insertion length information by looking for the first occurrence of a substring
        starting with "LEN". For example, for "D18_7&D19_7&D20_7&I20_7_C_LEN16", it extracts "LEN16" and then
        isolates "16" as the insertion length.
      - Splits the extracted value on 'LEN', converts the numeric part to an integer,
        and computes the frameshift as the absolute difference between insertion length and deletion count.
      - Defines valid deletion frames using parameters from advntr_settings and filters the DataFrame.

    For example:
      - "D55_6&D56_6&D57_6&D58_6&D59_6&D60_6" results in deletion count = 6, insertion length = 0,
        and frameshift = |0 - 6| = 6.
      - "D18_7&D19_7&D20_7&I20_7_C_LEN16" results in deletion count = 3, insertion count = 1,
        insertion length = 16, and frameshift = |16 - 3| = 13.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only those deletions that pass the frameshift filter.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting deletion processing.")

    # Create a copy of the input DataFrame to avoid modifying the original data.
    df1 = df.copy()
    logger.debug("Copied input DataFrame for deletion processing.")

    # Rename columns for clarity.
    df1.rename(columns={"State": "Variant", "Pvalue\n": "Pvalue"}, inplace=True)
    logger.debug("Renamed columns: 'State' -> 'Variant', 'Pvalue\\n' -> 'Pvalue'.")

    # Calculate deletion and insertion counts.
    df1["Deletion_length"] = df1["Variant"].str.count("D")
    df1["Insertion_length"] = df1["Variant"].str.count("I")
    logger.debug("Calculated 'Deletion_length' and 'Insertion_length'.")

    # Extract insertion length values using regex.
    # Fix: Use [0] to ensure extraction returns a Series rather than a DataFrame.
    df1["Insertion_len"] = df1["Variant"].str.extract("(LEN.*)")[0]
    logger.debug("Extracted 'Insertion_len' values from 'Variant' (as Series).")

    # Fill missing values and split the extracted string on 'LEN'.
    df1["Insertion_len"] = df1["Insertion_len"].fillna("LEN")
    df1[["I", "Insertion_len"]] = df1["Insertion_len"].str.split("LEN", expand=True)
    logger.debug("Split 'Insertion_len' column using 'LEN' as separator.")

    # Replace empty strings with '0' (avoiding downcasting warning) and convert to integer.
    df1["Insertion_len"] = (
        df1["Insertion_len"].astype(str).replace("^$", "0", regex=True)
    )
    df1["Insertion_len"] = (
        pd.to_numeric(df1["Insertion_len"], errors="coerce").fillna(0).astype(int)
    )
    df1["Deletion_length"] = df1["Deletion_length"].astype(int)
    logger.debug("Converted 'Insertion_len' and 'Deletion_length' to integers.")

    # Calculate frameshift as the absolute difference between insertion and deletion lengths.
    df1["frame"] = abs(df1["Insertion_len"] - df1["Deletion_length"]).astype(str)
    logger.debug(f"Computed frameshift values; sample: {df1['frame'].head().tolist()}")

    # Define valid deletion frames based on settings from advntr_settings.
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    del_frame = (np.arange(max_frameshift) * frameshift_multiplier + 2).astype(str)
    logger.debug(f"Valid deletion frames (first 5): {del_frame[:5].tolist()}")

    # Filter DataFrame: keep rows with at least one deletion and a valid frameshift.
    df1 = df1[(df1["Deletion_length"] >= 1) & df1["frame"].isin(del_frame)]
    logger.debug(f"Filtered DataFrame shape after deletion processing: {df1.shape}")

    return df1


def advntr_processing_ins(df):
    """
    Process adVNTR insertions by calculating insertion length, computing the frameshift,
    and filtering variants based on valid frameshift patterns.

    This function:
      - Renames columns for clarity.
      - Calculates the number of deletion ('D') and insertion ('I') characters in the entire
        variant string. Note that multi-allelic strings (e.g. "D18_7&D19_7&D20_7&I20_7_C_LEN16")
        are treated as a single value; the function aggregates the counts.
      - Extracts the insertion length information from the variant string by searching for a substring
        starting with "LEN". In the example "D18_7&D19_7&D20_7&I20_7_C_LEN16", it extracts "LEN16" and isolates "16".
      - Splits the extracted value on 'LEN', converts the numeric part to an integer,
        and computes the frameshift as the absolute difference between insertion length and deletion count.
      - Defines valid insertion frames using parameters from advntr_settings and filters the DataFrame.

    For example:
      - For a variant string without a LEN value, such as "D55_6&D56_6&D57_6&D58_6&D59_6&D60_6",
        the insertion length is considered 0 and the frameshift is the deletion count.
      - For "D18_7&D19_7&D20_7&I20_7_C_LEN16", the deletion count is 3, insertion count is 1,
        the extracted insertion length is 16, and the frameshift is |16 - 3| = 13.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only those insertions that pass the frameshift filter.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting insertion processing.")

    # Create a copy of the input DataFrame to avoid modifying the original data.
    df1 = df.copy()
    logger.debug("Copied input DataFrame for insertion processing.")

    # Rename columns for clarity.
    df1.rename(columns={"State": "Variant", "Pvalue\n": "Pvalue"}, inplace=True)
    logger.debug("Renamed columns: 'State' -> 'Variant', 'Pvalue\\n' -> 'Pvalue'.")

    # Calculate deletion and insertion counts.
    df1["Deletion_length"] = df1["Variant"].str.count("D")
    df1["Insertion_length"] = df1["Variant"].str.count("I")
    logger.debug("Calculated 'Deletion_length' and 'Insertion_length'.")

    # Extract insertion length values using regex.
    # Fix: Use [0] to ensure extraction returns a Series rather than a DataFrame.
    df1["Insertion_len"] = df1["Variant"].str.extract("(LEN.*)")[0]
    logger.debug("Extracted 'Insertion_len' values from 'Variant' (as Series).")

    # Fill missing values and split the extracted string on 'LEN'.
    df1["Insertion_len"] = df1["Insertion_len"].fillna("LEN")
    df1[["I", "Insertion_len"]] = df1["Insertion_len"].str.split("LEN", expand=True)
    logger.debug("Split 'Insertion_len' column using 'LEN' as separator.")

    # Replace empty strings with '0' (avoiding downcasting warning) and convert to integer.
    df1["Insertion_len"] = (
        df1["Insertion_len"].astype(str).replace("^$", "0", regex=True)
    )
    df1["Insertion_len"] = (
        pd.to_numeric(df1["Insertion_len"], errors="coerce").fillna(0).astype(int)
    )
    df1["Deletion_length"] = df1["Deletion_length"].astype(int)
    logger.debug("Converted 'Insertion_len' and 'Deletion_length' to integers.")

    # Calculate frameshift as the absolute difference between insertion and deletion lengths.
    df1["frame"] = abs(df1["Insertion_len"] - df1["Deletion_length"]).astype(str)
    logger.debug(f"Computed frameshift values; sample: {df1['frame'].head().tolist()}")

    # Define valid insertion frames based on settings from advntr_settings.
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    ins_frame = (np.arange(max_frameshift) * frameshift_multiplier + 1).astype(str)
    logger.debug(f"Valid insertion frames (first 5): {ins_frame[:5].tolist()}")

    # Filter DataFrame: keep rows with at least one insertion and a valid frameshift.
    df1 = df1[(df1["Insertion_len"] >= 1) & df1["frame"].isin(ins_frame)]
    logger.debug(f"Filtered DataFrame shape after insertion processing: {df1.shape}")

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

    logging.info("Processing adVNTR result...")

    # Read the output file
    with open(output_path, "r") as file:
        content = file.readlines()

    # Define default column names
    default_columns_advntr_output = [
        "#VID",
        "State",
        "NumberOfSupportingReads",
        "MeanCoverage",
        "Pvalue",
    ]

    # Check if the file is empty or only contains comments
    if not content or all(line.startswith("#") for line in content):
        logging.warning("No pathogenic variant was found with adVNTR!")

        # Write the header and a negative result line
        with open(output_path, "w") as f:
            f.write(f"#Input File: {output_name} \n")
            f.write("#Reference file: None\n")
            f.write("#P-value cutoff: 0.001\n")
            f.write("\t".join(default_columns_advntr_output) + "\n")
            f.write(
                "Negative\t"
                + "\t".join(["None"] * (len(default_columns_advntr_output) - 1))
                + "\n"
            )

        logging.info(
            f"Empty or comment-only output found, written header and negative line to {output_path}"
        )
        return

    # Read the output file again to update the header
    with open(output_path, "r") as file:
        content = file.readlines()

    # Replace #VID with VID in the header line
    content = [
        line.replace("#VID", "VID") if line.startswith("#VID") else line
        for line in content
    ]

    # Write the modified content back to the file
    with open(output_path, "w") as file:
        file.writelines(content)

    try:
        # Load the data using pandas
        logging.info("Loading data into DataFrame...")
        df = pd.read_csv(output_path, sep="\t", comment="#")
        logging.info(f"Data loaded successfully with shape: {df.shape}")
        logging.debug(f"First few rows of the DataFrame:\n{df.head()}")
    except Exception as e:
        logging.error(f"Error loading data into DataFrame: {e}")
        return

    try:
        # Process deletions and insertions
        logging.info("Processing deletions...")
        df_del = advntr_processing_del(df)

        logging.info("Processing insertions...")
        df_ins = advntr_processing_ins(df)

        # Concatenate deletions and insertions
        logging.info("Concatenating deletions and insertions...")
        advntr_concat = pd.concat([df_del, df_ins], axis=0)

        # Keep relevant columns
        default_columns_advntr_final = [
            "VID",
            "Variant",
            "NumberOfSupportingReads",
            "MeanCoverage",
            "Pvalue",
        ]
        advntr_concat = advntr_concat[default_columns_advntr_final]

        # Remove duplicates
        logging.info("Removing duplicates...")
        advntr_concat.drop_duplicates(
            subset=["VID", "Variant", "NumberOfSupportingReads"], inplace=True
        )

        # Save the processed adVNTR results
        output_result_path = os.path.join(output, f"{output_name}_adVNTR_result.tsv")
        advntr_concat.to_csv(output_result_path, sep="\t", index=False)
        logging.info(f"Processed adVNTR results saved to {output_result_path}")
    except Exception as e:
        logging.error(f"Error during processing of deletions and insertions: {e}")
        return

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
    logging.info("Intermediate files cleaned up.")
