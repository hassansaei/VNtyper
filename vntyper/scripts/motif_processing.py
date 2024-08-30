import pandas as pd
from Bio import SeqIO

def load_muc1_reference(reference_file):
    """
    Loads the MUC1 VNTR reference motifs from a FASTA file.
    """
    identifiers = []
    sequences = []
    with open(reference_file) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq)
    
    return pd.DataFrame({"Motifs": identifiers, "Motif_sequence": sequences})

def preprocessing_insertion(df, muc1_ref):
    """
    Preprocesses insertion variants by merging with the reference motifs.
    """
    # Rename the '#CHROM' column to 'Motifs'
    df.rename(columns={'#CHROM': 'Motifs'}, inplace=True)
    
    # Drop unwanted columns
    df.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    
    # Rename the last column to 'Sample'
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: 'Sample'}, inplace=True)
    
    # Merge with the MUC1 reference motifs
    df = pd.merge(df, muc1_ref, on='Motifs', how='left')
    
    # Add a 'Variant' column to indicate this is an insertion
    df['Variant'] = 'Insertion'
    
    return df

def preprocessing_deletion(df, muc1_ref):
    """
    Preprocesses deletion variants by merging with the reference motifs.
    """
    # Rename the '#CHROM' column to 'Motifs'
    df.rename(columns={'#CHROM': 'Motifs'}, inplace=True)
    
    # Drop unwanted columns
    df.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    
    # Rename the last column to 'Sample'
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: 'Sample'}, inplace=True)
    
    # Merge with the MUC1 reference motifs
    df = pd.merge(df, muc1_ref, on='Motifs', how='left')
    
    # Add a 'Variant' column to indicate this is a deletion
    df['Variant'] = 'Deletion'
    
    return df

# Load the additional motifs for final processing
def load_additional_motifs(config):
    """
    Loads additional motifs from the code-adVNTR_RUs.fa and MUC1_motifs_Rev_com.fa files for final annotation.
    """
    identifiers = []
    sequences = []
    
    # Get the file paths from the config
    code_advntr_file = config["reference_data"]["code_adVNTR_RUs"]
    muc1_motifs_rev_com_file = config["reference_data"]["muc1_motifs_rev_com"]

    with open(code_advntr_file) as RU_file, open(muc1_motifs_rev_com_file) as Motif_file:
        for seq_record in SeqIO.parse(RU_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq.upper())

        for seq_record in SeqIO.parse(Motif_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq.upper())

    s1 = pd.Series(identifiers, name='ID')
    s2 = pd.Series(sequences, name='Sequence')
    merged_motifs = pd.DataFrame(dict(Motif=s1, Motif_sequence=s2))

    return merged_motifs
    