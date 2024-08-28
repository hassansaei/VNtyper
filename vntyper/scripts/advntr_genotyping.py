import subprocess as sp
import logging
import pandas as pd
import numpy as np
import os
from pathlib import Path

def run_advntr(db_file_hg19, sorted_bam, output, output_name, config):
    """
    Run adVNTR genotyping using the specified database file and BAM file, fetching settings from the config.

    Args:
        db_file_hg19 (str): Path to the adVNTR VNTR database file.
        sorted_bam (str): Path to the sorted BAM file.
        output (str): Directory where the results will be saved.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary.
    """
    advntr_path = config["tools"]["advntr"]
    advntr_settings = config["advntr_settings"]
    
    # Set the number of threads from the config or default to 8
    threads = advntr_settings.get("threads", 8)
    
    # Retrieve additional command parts from the config, if available
    additional_commands = advntr_settings.get("additional_commands", "-aln")
    
    advntr_command = (
        f"{advntr_path} genotype -fs -vid {advntr_settings['vid']} "
        f"--alignment_file {sorted_bam} -o {output}{output_name}_adVNTR.vcf "
        f"-m {db_file_hg19} --working_directory {output} -t {threads} {additional_commands}"
    )
    
    # Redirect stdout and stderr to log files
    with open(f"{output}{output_name}_advntr_stdout.log", "w") as stdout_log, open(f"{output}{output_name}_advntr_stderr.log", "w") as stderr_log:
        logging.info("Launching adVNTR genotyping!")
        process = sp.Popen(advntr_command, shell=True, stdout=stdout_log, stderr=stderr_log)
        process.wait()
    
    if process.returncode != 0:
        logging.error("adVNTR genotyping failed. Please check the logs.")
        raise RuntimeError("adVNTR genotyping failed.")
    
    logging.info("adVNTR genotyping of MUC1-VNTR done!")


def read_vcf(path):
    """
    Read the header from a VCF file to extract column names.
    
    Args:
        path (str): Path to the VCF file.
    
    Returns:
        list: A list of column names from the VCF file, or an empty list if no header found.
    """
    vcf_names = []
    if not os.path.exists(path):
        logging.error(f"VCF file {path} does not exist.")
        return vcf_names
    
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#VID"):
                vcf_names = [x.strip() for x in line.split('\t')]
                break
    if not vcf_names:
        logging.error(f"No header found in VCF file: {path}")
    return vcf_names

def advntr_processing_del(df, config):
    """
    Process adVNTR deletions by calculating deletion length and filtering by frameshift.
    
    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.
        config (dict): Configuration dictionary containing advntr settings.
    
    Returns:
        pd.DataFrame: Filtered DataFrame containing deletions.
    """
    df1 = df.copy()
    df1.rename(columns={'State': 'Variant', 'Pvalue\n': 'Pvalue'}, inplace=True)
    
    # Calculate deletion and insertion lengths in the variant sequences.
    df1['Deletion_length'] = df1['Variant'].str.count('D')
    df1['Insertion'] = df1['Variant'].str.count('I')
    
    # Extract insertion length and handle missing values.
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1[['I', 'Insertion_len']] = df1['Insertion_len'].str.split('LEN', expand=True)
    df1.Insertion_len = df1.Insertion_len.replace('', 0).astype(int)
    df1.Deletion_length = df1.Deletion_length.astype(int)
    
    # Calculate the frameshift value as the difference between insertion and deletion lengths.
    df1['frame'] = abs(df1.Insertion_len - df1.Deletion_length).astype(str)
    
    # Define frameshift patterns using config settings.
    max_frameshift = config["advntr_settings"]["max_frameshift"]
    frameshift_multiplier = config["advntr_settings"]["frameshift_multiplier"]
    
    del_frame = (np.arange(max_frameshift) * frameshift_multiplier + 2).astype(str)
    
    # Filter for valid deletions with matching frameshift values.
    df1 = df1[(df1['Deletion_length'] >= 1) & df1['frame'].isin(del_frame)]
    
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
    df1 = df.copy()
    df1.rename(columns={'State': 'Variant', 'Pvalue\n': 'Pvalue'}, inplace=True)
    
    # Calculate deletion and insertion lengths in the variant sequences.
    df1['Deletion_length'] = df1['Variant'].str.count('D')
    df1['Insertion'] = df1['Variant'].str.count('I')
    
    # Extract insertion length and handle missing values.
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1[['I', 'Insertion_len']] = df1['Insertion_len'].str.split('LEN', expand=True)
    df1.Insertion_len = df1.Insertion_len.replace('', 0).astype(int)
    df1.Deletion_length = df1.Deletion_length.astype(int)
    
    # Calculate the frameshift value as the difference between insertion and deletion lengths.
    df1['frame'] = abs(df1.Insertion_len - df1.Deletion_length).astype(str)
    
    # Define frameshift patterns using config settings.
    max_frameshift = config["advntr_settings"]["max_frameshift"]
    frameshift_multiplier = config["advntr_settings"]["frameshift_multiplier"]
    
    ins_frame = (np.arange(max_frameshift) * frameshift_multiplier + 1).astype(str)
    
    # Filter for valid insertions with matching frameshift values.
    df1 = df1[(df1['Insertion_len'] >= 1) & df1['frame'].isin(ins_frame)]
    
    return df1

def process_advntr_output(vcf_path, output, output_name, config):
    """
    Process the adVNTR output VCF to extract relevant information and generate final results.
    
    Args:
        vcf_path (str): Path to the adVNTR VCF file.
        output (str): Directory where the final results will be saved.
        output_name (str): Base name for the output files.
        config (dict): Configuration dictionary.
    """
    if not os.path.exists(vcf_path):
        logging.error(f"adVNTR VCF file {vcf_path} not found!")
        return

    logging.info('Processing code-adVNTR result...')
    
    names = read_vcf(vcf_path)
    df = pd.read_csv(vcf_path, comment='#', sep='\s+', header=None, names=names)

    # Define frameshift patterns for both insertions and deletions using config settings.
    max_frameshift = config["advntr_settings"]["max_frameshift"]
    frameshift_multiplier = config["advntr_settings"]["frameshift_multiplier"]
    
    ins_frame = (np.arange(max_frameshift) * frameshift_multiplier + 1).astype(str)
    del_frame = (np.arange(max_frameshift) * frameshift_multiplier + 2).astype(str)

    if df.empty:
        logging.warning('No pathogenic variant was found with code-adVNTR!')
        pre_result_path = output + output_name + '_pre_result.tsv'
        if not os.path.exists(pre_result_path):
            logging.error(f"Pre-result file {pre_result_path} not found!")
            return

        with open(pre_result_path, 'r') as f:
            with open(output + output_name + '_Final_result.tsv', 'w') as f1:
                f1.write(f'## VNtyper_Analysis_for_{output_name} \n')
                f1.write('# Kestrel_Result \n')
                f1.write('\t'.join(columns_kmer) + '\n')
                next(f)
                for line in f:
                    f1.write(line)

        # Clean up intermediate files.
        rm_list_4 = config.get("advntr_cleanup_files", [])
        for db in rm_list_4:
            db_str = output + output_name + db
            rm_command = "rm " + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
        
        logging.info('The final result is saved in *_Final_result.tsv')
        return  # Replaced sys.exit() with return to avoid abrupt script termination
    else:
        # Process deletions and insertions using the respective functions.
        df_del = advntr_processing_del(df, config)
        df_ins = advntr_processing_ins(df, config)

    # Concatenate the processed results for deletions and insertions.
    advntr_concat = pd.concat([df_del, df_ins], axis=0)
    
    # Keep only the relevant columns.
    advntr_concat = advntr_concat[['#VID', 'Variant', 'NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']]

    # Remove duplicate records based on specific columns.
    advntr_concat.drop_duplicates(subset=['#VID', 'Variant', 'NumberOfSupportingReads'], inplace=True)

    # Save the processed adVNTR results to a TSV file.
    output_result_path = os.path.join(output, f"{output_name}_adVNTR_result.tsv")
    advntr_concat.to_csv(output_result_path, sep='\t', index=False)
    logging.info(f"Processed adVNTR results saved to {output_result_path}")

    # Clean up intermediate files, if configured
    rm_list_4 = config.get("advntr_cleanup_files", [])
    for db in rm_list_4:
        db_str = os.path.join(output, f"{output_name}{db}")
        if os.path.exists(db_str):
            rm_command = f"rm {db_str}"
            process = sp.Popen(rm_command, shell=True)
            process.wait()

    logging.info('The final result is saved and intermediate files cleaned up.')
   
