import pandas as pd
from Bio import SeqIO

def process_motifs(reference_VNTR, tools_path):
    # Code for processing motifs from reference files
    pass

def preprocessing_insertion(df, motif_ref):
    df1 = df.copy()
    df1.rename(columns={'#CHROM': 'Motifs'}, inplace=True)
    df1.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    df1 = pd.merge(df1, motif_ref, on='Motifs', how='left')
    df1['Variant'] = 'Insertion'
    return df1

def preprocessing_deletion(df, motif_ref):
    df2 = df.copy()
    df2.rename(columns={'#CHROM': 'Motifs'}, inplace=True)
    df2.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    df2 = pd.merge(df2, motif_ref, on='Motifs', how='left')
    df2['Variant'] = 'Deletion'
    return df2
