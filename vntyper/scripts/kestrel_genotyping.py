import subprocess as sp
import logging
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from vntyper.scripts.file_processing import filter_vcf, filter_indel_vcf, read_vcf

# Construct the Kestrel command based on kmer size and config settings
def construct_kestrel_command(kmer_size, kestrel_path, reference_vntr, output_dir, fastq_1, fastq_2, temp_dir, vcf_out, java_path, java_memory, max_align_states, max_hap_states):
    if not fastq_1 or not fastq_2:
        raise ValueError("FASTQ input files are missing or invalid.")
    
    return (
        f"{java_path} -Xmx{java_memory} -jar {kestrel_path} -k {kmer_size} "
        f"--maxalignstates {max_align_states} --maxhapstates {max_hap_states} "
        f"-r {reference_vntr} -o {vcf_out} "
        f"{fastq_1} {fastq_2} "
        f"--temploc {temp_dir} "
        f"--hapfmt sam -p {temp_dir}/output.sam"
    )

# Kestrel processing logic
def run_kestrel(vcf_path, output_dir, fastq_1, fastq_2, reference_vntr, kestrel_path, temp_dir, kestrel_settings):
    # Extract kestrel settings from the config
    java_path = kestrel_settings.get("java_path", "java")
    java_memory = kestrel_settings.get("java_memory", "15g")
    kmer_sizes = kestrel_settings.get("kmer_sizes", [20, 17, 25, 41])
    max_align_states = kestrel_settings.get("max_align_states", 30)
    max_hap_states = kestrel_settings.get("max_hap_states", 30)
    
    if not fastq_1 or not fastq_2:
        raise ValueError("FASTQ files are not provided for Kestrel.")

    for kmer_size in kmer_sizes:
        kmer_command = construct_kestrel_command(
            kmer_size=kmer_size,
            kestrel_path=kestrel_path,
            reference_vntr=reference_vntr,
            output_dir=output_dir,
            fastq_1=fastq_1,
            fastq_2=fastq_2,
            temp_dir=temp_dir,
            vcf_out=vcf_path,
            java_path=java_path,
            java_memory=java_memory,
            max_align_states=max_align_states,
            max_hap_states=max_hap_states
        )
        
        if vcf_path.is_file():
            logging.info("VCF file already exists, skipping Kestrel run...")
            return
        else:
            logging.info(f"Launching Kestrel with kmer size {kmer_size}...")
            process = sp.Popen(kmer_command, shell=True)
            process.wait()
            logging.info(f"Mapping-free genotyping of MUC1-VNTR with kmer size {kmer_size} done!")

            # Process the VCF if Kestrel run was successful
            if vcf_path.is_file():
                process_kestrel_output(output_dir, vcf_path, reference_vntr)
                break  # If successful, break out of the loop and stop trying other kmer sizes

# Process Kestrel VCF results to generate the final output format
def process_kestrel_output(output_dir, vcf_path, reference_vntr):
    logging.info("Processing Kestrel VCF results...")
    indel_vcf = f"{output_dir}/output_indel.vcf"
    output_ins = f"{output_dir}/output_insertion.vcf"
    output_del = f"{output_dir}/output_deletion.vcf"
    
    # Filter the VCF to extract indels, insertions, and deletions
    filter_vcf(vcf_path, indel_vcf)
    filter_indel_vcf(indel_vcf, output_ins, output_del)
    
    # Read the filtered VCF files into dataframes
    names = read_vcf(vcf_path)
    vcf_insertion = pd.read_csv(output_ins, comment='#', sep='\s+', header=None, names=names)
    vcf_deletion = pd.read_csv(output_del, comment='#', sep='\s+', header=None, names=names)
    
    # Load MUC1 VNTR reference motifs
    MUC1_ref = load_muc1_reference(reference_vntr)

    # Preprocess insertion and deletion dataframes
    insertion_df = preprocessing_insertion(vcf_insertion, MUC1_ref)
    deletion_df = preprocessing_deletion(vcf_deletion, MUC1_ref)
    
    # Combine insertion and deletion dataframes
    combined_df = pd.concat([insertion_df, deletion_df], axis=0)

    # Process and filter results based on frameshifts and confidence scores
    processed_df = process_kmer_results(combined_df)

    # Save the processed dataframe to a TSV file
    final_output_path = os.path.join(output_dir, "kestrel_result.tsv")
    processed_df.to_csv(final_output_path, sep='\t', index=False)
    
    logging.info("Kestrel VCF processing completed.")
    return processed_df

def load_muc1_reference(reference_file):
    identifiers = []
    sequences = []
    with open(reference_file) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq)
    
    return pd.DataFrame({"Motifs": identifiers, "Motif_sequence": sequences})

def preprocessing_insertion(df, muc1_ref):
    df1 = df.copy()
    df1.rename(columns={'#CHROM': 'Motifs'}, inplace=True)
    df1.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    df1 = pd.merge(df1, muc1_ref, on='Motifs', how='left')
    df1['Variant'] = 'Insertion'
    return df1

def preprocessing_deletion(df, muc1_ref):
    df2 = df.copy()
    df2.rename(columns={'#CHROM': 'Motifs'}, inplace=True)
    df2.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    df2 = pd.merge(df2, muc1_ref, on='Motifs', how='left')
    df2['Variant'] = 'Deletion'
    return df2

def process_kmer_results(combined_df):
    # Split and filter based on frameshifts and confidence criteria
    combined_df["ref_len"] = combined_df["REF"].str.len()
    combined_df["alt_len"] = combined_df["ALT"].str.len()
    combined_df["Frame_Score"] = round((combined_df.alt_len - combined_df.ref_len) / 3, 2)
    combined_df['Frame_Score'] = combined_df['Frame_Score'].astype(str).apply(lambda x: x.replace('.0', 'C'))
    combined_df["TrueFalse"] = combined_df['Frame_Score'].str.contains('C', regex=True)
    combined_df["TrueFalse"] = combined_df["TrueFalse"].astype(str)
    combined_df = combined_df[combined_df['TrueFalse'].str.contains("False")]

    # More processing, filtering, and labeling steps can be added here

    return combined_df
