import subprocess as sp
import logging
import pandas as pd
import numpy as np
import os
from pathlib import Path

# Constants for frame shift calculations
MAX_FRAMESHIFT = 100  # The maximum number of frameshift cycles to check
FRAMESHIFT_MULTIPLIER = 3  # The multiplier for generating insertion and deletion frames

def run_advntr(reference, db_file_hg19, sorted_bam, output, output_name):
    """
    Run adVNTR genotyping using the specified reference and BAM file.
    
    Args:
        reference (str): Path to the reference FASTA file.
        db_file_hg19 (str): Path to the adVNTR VNTR database file.
        sorted_bam (str): Path to the sorted BAM file.
        output (str): Directory where the results will be saved.
        output_name (str): Base name for the output files.
    """
    advntr_command = (
        f"advntr genotype -fs -vid 25561 --outfmt vcf "
        f"--alignment_file {sorted_bam} -o {output}{output_name}_adVNTR.vcf "
        f"-m {db_file_hg19} -r {reference} --working_directory {output}"
    )
    
    logging.info("Launching adVNTR genotyping!")
    process = sp.Popen(advntr_command, shell=True)
    process.wait()
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
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#VID"):
                vcf_names = [x.strip() for x in line.split('\t')]
                break
    if not vcf_names:
        logging.error(f"No header found in VCF file: {path}")
    return vcf_names

def advntr_processing_del(df):
    """
    Process adVNTR deletions by calculating deletion length and filtering by frameshift.
    
    Args:
        df (pd.DataFrame): DataFrame containing adVNTR variant data.
    
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
    
    # Filter for valid deletions with matching frameshift values.
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
    df1 = df.copy()  # Renamed from dff to df1 for consistency
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
    
    # Filter for valid insertions with matching frameshift values.
    df1 = df1[(df1['Insertion_len'] >= 1) & df1['frame'].isin(ins_frame)]
    
    return df1

def process_advntr_output(vcf_path, output, output_name):
    """
    Process the adVNTR output VCF to extract relevant information and generate final results.
    
    Args:
        vcf_path (str): Path to the adVNTR VCF file.
        output (str): Directory where the final results will be saved.
        output_name (str): Base name for the output files.
    """
    if not os.path.exists(vcf_path):
        logging.error(f"adVNTR VCF file {vcf_path} not found!")
        return

    logging.info('Processing code-adVNTR result...')
    
    names = read_vcf(vcf_path)
    df = pd.read_csv(vcf_path, comment='#', delim_whitespace=True, header=None, names=names)

    # Define frameshift patterns for both insertions and deletions using constants.
    ins_frame = np.arange(MAX_FRAMESHIFT) * FRAMESHIFT_MULTIPLIER + 1
    ins_frame = ins_frame.astype(str)
    del_frame = np.arange(MAX_FRAMESHIFT) * FRAMESHIFT_MULTIPLIER + 2
    del_frame = del_frame.astype(str)

    if df.empty:
        print('No pathogenic variant was found with code-adVNTR!')
        with open(output + output_name + '_pre_result.tsv', 'r') as f:
            with open(output + output_name + '_Final_result.tsv', 'w') as f1:
                f1.write(f'## VNtyper_Analysis_for_{output_name} \n')
                f1.write('# Kestrel_Result \n')
                f1.write('\t'.join(columns_kmer) + '\n')
                next(f)
                for line in f:
                    f1.write(line)

        # Clean up intermediate files.
        for db in rm_list_4:
            db_str = output + output_name + db
            rm_command = "rm " + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
        
        logging.info('The final result is saved in *_Final_result.tsv')
        return  # Replaced sys.exit() with return to avoid abrupt script termination
    else:
        df_del = advntr_processing_del(df)
        df_ins = advntr_processing_ins(df)

    advntr_concat = pd.concat([df_del, df_ins], axis=0)
    advntr_concat = advntr_concat[['#VID', 'Variant', 'NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']]

    advntr_concat.drop_duplicates(subset=['#VID', 'Variant', 'NumberOfSupportingReads'], inplace=True)
    advntr_concat.to_csv(output + output_name + '_adVNTR_result.tsv', sep='\t', index=False)
    logging.info(f"Processed adVNTR results saved to {output + output_name + '_adVNTR_result.tsv'}")
