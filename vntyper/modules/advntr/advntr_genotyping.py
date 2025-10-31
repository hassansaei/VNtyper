#!/usr/bin/env python3
# vntyper/modules/advntr/advntr_genotyping.py

import logging
import os
import re
import subprocess as sp

import numpy as np
import pandas as pd

from vntyper.scripts.utils import load_config, run_command

# -------------------------------------------------------------------------
# Configure logging
# -------------------------------------------------------------------------
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")


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


def run_advntr(db_file, sorted_bam, output, output_name, config, cwd=None):
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
        if not run_command(advntr_command, log_file, critical=True, cwd=cwd):
            logging.error("adVNTR genotyping failed. Check the log for details.")
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

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only those deletions that pass the frameshift filter.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting deletion processing.")
    df1 = df.copy()
    logger.debug("Copied input DataFrame for deletion processing.")
    df1.rename(columns={"State": "Variant", "Pvalue\n": "Pvalue"}, inplace=True)
    logger.debug("Renamed columns: 'State' -> 'Variant', 'Pvalue\\n' -> 'Pvalue'.")
    df1["Deletion_length"] = df1["Variant"].str.count("D")
    df1["Insertion_length"] = df1["Variant"].str.count("I")
    logger.debug("Calculated 'Deletion_length' and 'Insertion_length'.")
    df1["Insertion_len"] = df1["Variant"].str.extract("(LEN.*)")[0]
    logger.debug("Extracted 'Insertion_len' values from 'Variant' (as Series).")
    df1["Insertion_len"] = df1["Insertion_len"].fillna("LEN")
    df1[["I", "Insertion_len"]] = df1["Insertion_len"].str.split("LEN", expand=True)
    logger.debug("Split 'Insertion_len' column using 'LEN' as separator.")
    df1["Insertion_len"] = df1["Insertion_len"].astype(str).replace("^$", "0", regex=True)
    df1["Insertion_len"] = pd.to_numeric(df1["Insertion_len"], errors="coerce").fillna(0).astype(int)
    df1["Deletion_length"] = df1["Deletion_length"].astype(int)
    logger.debug("Converted 'Insertion_len' and 'Deletion_length' to integers.")
    df1["frame"] = abs(df1["Insertion_len"] - df1["Deletion_length"]).astype(str)
    logger.debug(f"Computed frameshift values; sample: {df1['frame'].head().tolist()}")
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    del_frame = (np.arange(max_frameshift) * frameshift_multiplier + 2).astype(str)
    logger.debug(f"Valid deletion frames (first 5): {del_frame[:5].tolist()}")
    df1 = df1[(df1["Deletion_length"] >= 1) & df1["frame"].isin(del_frame)]
    logger.debug(f"Filtered DataFrame shape after deletion processing: {df1.shape}")
    return df1


def advntr_processing_ins(df):
    """
    Process adVNTR insertions by calculating insertion length, computing the frameshift,
    and filtering variants based on valid frameshift patterns.

    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only those insertions that pass the frameshift filter.
    """
    logger = logging.getLogger(__name__)
    logger.debug("Starting insertion processing.")
    df1 = df.copy()
    logger.debug("Copied input DataFrame for insertion processing.")
    df1.rename(columns={"State": "Variant", "Pvalue\n": "Pvalue"}, inplace=True)
    logger.debug("Renamed columns: 'State' -> 'Variant', 'Pvalue\\n' -> 'Pvalue'.")
    df1["Deletion_length"] = df1["Variant"].str.count("D")
    df1["Insertion_length"] = df1["Variant"].str.count("I")
    logger.debug("Calculated 'Deletion_length' and 'Insertion_length'.")
    df1["Insertion_len"] = df1["Variant"].str.extract("(LEN.*)")[0]
    logger.debug("Extracted 'Insertion_len' values from 'Variant' (as Series).")
    df1["Insertion_len"] = df1["Insertion_len"].fillna("LEN")
    df1[["I", "Insertion_len"]] = df1["Insertion_len"].str.split("LEN", expand=True)
    logger.debug("Split 'Insertion_len' column using 'LEN' as separator.")
    df1["Insertion_len"] = df1["Insertion_len"].astype(str).replace("^$", "0", regex=True)
    df1["Insertion_len"] = pd.to_numeric(df1["Insertion_len"], errors="coerce").fillna(0).astype(int)
    df1["Deletion_length"] = df1["Deletion_length"].astype(int)
    logger.debug("Converted 'Insertion_len' and 'Deletion_length' to integers.")
    df1["frame"] = abs(df1["Insertion_len"] - df1["Deletion_length"]).astype(str)
    logger.debug(f"Computed frameshift values; sample: {df1['frame'].head().tolist()}")
    max_frameshift = advntr_settings.get("max_frameshift", 100)
    frameshift_multiplier = advntr_settings.get("frameshift_multiplier", 3)
    ins_frame = (np.arange(max_frameshift) * frameshift_multiplier + 1).astype(str)
    logger.debug(f"Valid insertion frames (first 5): {ins_frame[:5].tolist()}")
    df1 = df1[(df1["Insertion_len"] >= 1) & df1["frame"].isin(ins_frame)]
    logger.debug(f"Filtered DataFrame shape after insertion processing: {df1.shape}")
    return df1


def load_ru_sequences(ru_fasta_path):
    """
    Load repeat unit (RU) sequences from a FASTA file.

    Args:
        ru_fasta_path (str): Path to the RU FASTA file.

    Returns:
        dict: A dictionary mapping RU identifier (as string) to its sequence.
    """
    ru_dict = {}
    with open(ru_fasta_path) as f:
        current_ru = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_ru and seq_lines:
                    ru_dict[current_ru] = "".join(seq_lines)
                header = line[1:]
                current_ru = header[2:] if header.startswith("RU") else header
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_ru and seq_lines:
            ru_dict[current_ru] = "".join(seq_lines)
    return ru_dict


def annotate_advntr_variants(variant_series, ru_fasta_path):
    """
    Annotate adVNTR variants with RU, POS, REF, and ALT using the RU FASTA file.

    Args:
        variant_series (pd.Series): Series of variant strings (possibly with multiple parts separated by '&').
        ru_fasta_path (str): Path to the RU FASTA file.

    Returns:
        tuple: Four lists corresponding to RU, POS, REF, and ALT annotations.
    """
    ru_dict = load_ru_sequences(ru_fasta_path)
    ru_annotations = []
    pos_annotations = []
    ref_annotations = []
    alt_annotations = []

    ins_pattern = re.compile(r"^I(\d+)_([0-9]+)_([ACGT])_LEN(\d+)$")
    del_pattern = re.compile(r"^D(\d+)_([0-9]+)$")

    for variant in variant_series:
        parts = variant.split("&")
        ru_parts = []
        pos_parts = []
        ref_parts = []
        alt_parts = []
        for part in parts:
            part = part.strip()
            ins_match = ins_pattern.match(part)
            del_match = del_pattern.match(part)
            if ins_match:
                pos_val = int(ins_match.group(1))
                ru_val = ins_match.group(2)
                inserted_base = ins_match.group(3)
                ins_len = int(ins_match.group(4))
                ru_seq = ru_dict.get(ru_val, "")
                ref_base = ru_seq[pos_val - 1] if ru_seq and pos_val - 1 < len(ru_seq) else "."
                alt_val = ref_base + inserted_base * ins_len
                ru_parts.append(ru_val)
                pos_parts.append(str(pos_val))
                ref_parts.append(ref_base)
                alt_parts.append(alt_val)
            elif del_match:
                pos_val = int(del_match.group(1))
                ru_val = del_match.group(2)
                ru_seq = ru_dict.get(ru_val, "")
                if ru_seq and pos_val > 1 and pos_val - 1 < len(ru_seq):
                    prev_base = ru_seq[pos_val - 2]
                    del_base = ru_seq[pos_val - 1]
                    ref_allele = prev_base + del_base
                    alt_allele = prev_base
                else:
                    ref_allele = "."
                    alt_allele = "."
                ru_parts.append(ru_val)
                pos_parts.append(str(pos_val))
                ref_parts.append(ref_allele)
                alt_parts.append(alt_allele)
            else:
                ru_parts.append(".")
                pos_parts.append(".")
                ref_parts.append(".")
                alt_parts.append(".")
        ru_annotations.append(",".join(ru_parts))
        pos_annotations.append(",".join(pos_parts))
        ref_annotations.append(",".join(ref_parts))
        alt_annotations.append(",".join(alt_parts))

    return ru_annotations, pos_annotations, ref_annotations, alt_annotations


def process_advntr_output(output_path, output, output_name, config=None):
    """
    Process the adVNTR output to extract relevant information and generate final results.

    Optionally, if a configuration is provided and it includes a valid
    'reference_data.code_adVNTR_RUs' FASTA file, the function will annotate each
    variant with the affected repeat unit (RU), position (POS), REF and ALT values.

    The final output always contains the columns:
      "VID, Variant, NumberOfSupportingReads, MeanCoverage, Pvalue, RU, POS, REF, ALT, Flag".

    If the VCF data is empty, a negative result is generated immediately with
    'VID' set to "Negative" and all other columns set to "None", and further processing is skipped.

    Args:
        output_path (str): Path to the adVNTR output file.
        output (str): Directory where the final results will be saved.
        output_name (str): Base name for the output files.
        config (dict, optional): Main configuration dictionary.
    """
    if not os.path.exists(output_path):
        logging.error(f"adVNTR output file {output_path} not found!")
        return

    logging.info("Processing adVNTR result...")

    with open(output_path) as file:
        content = file.readlines()

    # Replace header to ensure consistency
    content = [line.replace("#VID", "VID") if line.startswith("#VID") else line for line in content]
    with open(output_path, "w") as file:
        file.writelines(content)

    try:
        logging.info("Loading data into DataFrame...")
        df = pd.read_csv(output_path, sep="\t", comment="#")
        logging.info(f"Data loaded successfully with shape: {df.shape}")
        logging.debug(f"First few rows of the DataFrame:\n{df.head()}")
    except Exception as e:
        logging.error(f"Error loading data into DataFrame: {e}")
        return

    # Immediately check if the loaded DataFrame is empty
    final_columns = [
        "VID",
        "Variant",
        "NumberOfSupportingReads",
        "MeanCoverage",
        "Pvalue",
        "RU",
        "POS",
        "REF",
        "ALT",
        "Flag",
    ]
    if df.empty:
        logging.warning("VCF file is empty. Generating default negative result.")
        advntr_concat = pd.DataFrame(
            [
                {
                    "VID": "Negative",
                    "Variant": "Not applicable",
                    "NumberOfSupportingReads": "Not applicable",
                    "MeanCoverage": "Not applicable",
                    "Pvalue": "Not applicable",
                    "RU": "Not applicable",
                    "POS": "Not applicable",
                    "REF": "Not applicable",
                    "ALT": "Not applicable",
                    "Flag": "Not applicable",
                }
            ]
        )
        output_result_path = os.path.join(output, f"{output_name}_adVNTR_result.tsv")
        advntr_concat = advntr_concat[final_columns]
        advntr_concat.to_csv(output_result_path, sep="\t", index=False)
        logging.info(f"Processed adVNTR results saved to {output_result_path}")
        cleanup_files(output, output_name)
        return

    try:
        logging.info("Processing deletions...")
        df_del = advntr_processing_del(df)

        logging.info("Processing insertions...")
        df_ins = advntr_processing_ins(df)

        logging.info("Concatenating deletions and insertions...")
        advntr_concat = pd.concat([df_del, df_ins], axis=0)

        if advntr_concat.empty:
            logging.warning("No pathogenic variant found after filtering. Generating default negative result.")
            advntr_concat = pd.DataFrame(
                [
                    {
                        "VID": "Negative",
                        "Variant": "Not applicable",
                        "NumberOfSupportingReads": "Not applicable",
                        "MeanCoverage": "Not applicable",
                        "Pvalue": "Not applicable",
                        "RU": "Not applicable",
                        "POS": "Not applicable",
                        "REF": "Not applicable",
                        "ALT": "Not applicable",
                        "Flag": "Not applicable",
                    }
                ]
            )
        else:
            base_columns = [
                "VID",
                "Variant",
                "NumberOfSupportingReads",
                "MeanCoverage",
                "Pvalue",
            ]
            advntr_concat = advntr_concat[base_columns]
            logging.info("Removing duplicates...")
            advntr_concat.drop_duplicates(subset=["VID", "Variant", "NumberOfSupportingReads"], inplace=True)

            # Perform RU-level annotation if possible
            if config:
                ru_fasta_path = config.get("reference_data", {}).get("code_adVNTR_RUs")
                if ru_fasta_path and os.path.exists(ru_fasta_path):
                    logging.info("Annotating variants with RU-level information.")
                    ru_ann, pos_ann, ref_ann, alt_ann = annotate_advntr_variants(
                        advntr_concat["Variant"], ru_fasta_path
                    )
                    advntr_concat["RU"] = ru_ann
                    advntr_concat["POS"] = pos_ann
                    advntr_concat["REF"] = ref_ann
                    advntr_concat["ALT"] = alt_ann

            # Apply flagging rules if available
            flagging_rules = advntr_config.get("flagging_rules", {})
            if flagging_rules:
                logging.info("Applying flagging rules to adVNTR output.")
                from vntyper.scripts.flagging import add_flags

                advntr_concat = add_flags(advntr_concat, flagging_rules)

            # Ensure all final columns are present
            for col in final_columns:
                if col not in advntr_concat.columns:
                    advntr_concat[col] = "Not applicable"

        advntr_concat = advntr_concat[final_columns]
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
    logging.info("Intermediate files cleaned up.")
