#!/usr/bin/env python3
# vntyper/scripts/motif_processing.py

import pandas as pd
from Bio import SeqIO


def load_muc1_reference(reference_file):
    """
    Load the MUC1 VNTR reference motifs from a FASTA file.

    This function parses a FASTA file containing MUC1 VNTR reference motifs and
    returns a pandas DataFrame with two columns: 'Motifs' and 'Motif_sequence'.

    Args:
        reference_file (str): Path to the FASTA file containing reference motifs.

    Returns:
        pd.DataFrame: DataFrame containing motif identifiers and their sequences.
    """
    identifiers = []
    sequences = []

    # Open and parse the FASTA file
    with open(reference_file) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq)

    # Create a DataFrame from the parsed data
    return pd.DataFrame({
        "Motifs": identifiers,
        "Motif_sequence": sequences
    })


def preprocessing_insertion(df, muc1_ref):
    """
    Preprocess insertion variants by merging with the MUC1 reference motifs.

    This function renames specific columns, drops unnecessary ones, merges the
    DataFrame with reference motifs, and adds a 'Variant' column to indicate
    insertion events.

    Args:
        df (pd.DataFrame): DataFrame containing insertion variant data.
        muc1_ref (pd.DataFrame): DataFrame with MUC1 reference motifs.

    Returns:
        pd.DataFrame: Preprocessed DataFrame with insertion variants annotated.
    """
    # Rename the '#CHROM' column to 'Motifs'
    df.rename(columns={'#CHROM': 'Motifs'}, inplace=True)

    # Drop unwanted columns
    columns_to_drop = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    df.drop(columns=columns_to_drop, axis=1, inplace=True)

    # Rename the last column to 'Sample'
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: 'Sample'}, inplace=True)

    # Merge with the MUC1 reference motifs on 'Motifs'
    df = pd.merge(df, muc1_ref, on='Motifs', how='left')

    # Add a 'Variant' column to indicate this is an insertion
    df['Variant'] = 'Insertion'

    return df


def preprocessing_deletion(df, muc1_ref):
    """
    Preprocess deletion variants by merging with the MUC1 reference motifs.

    This function renames specific columns, drops unnecessary ones, merges the
    DataFrame with reference motifs, and adds a 'Variant' column to indicate
    deletion events.

    Args:
        df (pd.DataFrame): DataFrame containing deletion variant data.
        muc1_ref (pd.DataFrame): DataFrame with MUC1 reference motifs.

    Returns:
        pd.DataFrame: Preprocessed DataFrame with deletion variants annotated.
    """
    # Rename the '#CHROM' column to 'Motifs'
    df.rename(columns={'#CHROM': 'Motifs'}, inplace=True)

    # Drop unwanted columns
    columns_to_drop = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    df.drop(columns=columns_to_drop, axis=1, inplace=True)

    # Rename the last column to 'Sample'
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: 'Sample'}, inplace=True)

    # Merge with the MUC1 reference motifs on 'Motifs'
    df = pd.merge(df, muc1_ref, on='Motifs', how='left')

    # Add a 'Variant' column to indicate this is a deletion
    df['Variant'] = 'Deletion'

    return df


def load_additional_motifs(config):
    """
    Load additional motifs for final annotation from specified FASTA files.

    This function reads motifs from 'code-adVNTR_RUs.fa' and
    'MUC1_motifs_Rev_com.fa' files as specified in the configuration,
    merges them into a single DataFrame, and returns it.

    Args:
        config (dict): Configuration dictionary containing paths to reference data.

    Returns:
        pd.DataFrame: DataFrame containing merged additional motifs with IDs and sequences.
    """
    identifiers = []
    sequences = []

    # Retrieve file paths from the configuration
    code_advntr_file = config["reference_data"]["code_adVNTR_RUs"]
    muc1_motifs_rev_com_file = config["reference_data"]["muc1_motifs_rev_com"]

    # Open and parse both FASTA files
    with open(code_advntr_file) as ru_file, open(muc1_motifs_rev_com_file) as motif_file:
        for seq_record in SeqIO.parse(ru_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(str(seq_record.seq.upper()))

        for seq_record in SeqIO.parse(motif_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(str(seq_record.seq.upper()))

    # Create pandas Series from the lists
    s1 = pd.Series(identifiers, name='ID')
    s2 = pd.Series(sequences, name='Sequence')

    # Merge the Series into a single DataFrame
    merged_motifs = pd.DataFrame({
        'Motif': s1,
        'Motif_sequence': s2
    })

    return merged_motifs
