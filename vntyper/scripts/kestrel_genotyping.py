import subprocess as sp
import logging
import pandas as pd
from pathlib import Path

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

# Filter VCF for indels
def filter_vcf(input_path, output_path):
    with open(input_path, "r") as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line[:2] == "##":
                indel_file.write(line)
            else:
                [_, _, _, ref, alt, *_] = line.split("\t")
                if len(ref) == 1 and len(alt) != 1 or len(ref) != 1 and len(alt) == 1:
                    indel_file.write(line)

# Further filter the indel VCF into insertion and deletion files
def filter_indel_vcf(indel_vcf, output_ins, output_del):
    with open(indel_vcf, "r") as vcf_file, open(output_ins, "w") as insertion_file, open(output_del, "w") as deletion_file:
        for line in vcf_file:
            if line[:2] == "##":
                insertion_file.write(line)
                deletion_file.write(line)
            else:
                [_, _, _, ref, alt, *_] = line.split("\t")
                if len(ref) == 1 and len(alt) > 1:
                    insertion_file.write(line)
                else:
                    deletion_file.write(line)

# Read VCF headers
def read_vcf(path):
    vcf_names = None
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    
    if not vcf_names:
        raise ValueError(f"No valid VCF headers found in {path}")
    
    return vcf_names

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
                process_vcf_results(output_dir, vcf_path)
                break  # If successful, break out of the loop and stop trying other kmer sizes

def process_vcf_results(output_dir, vcf_path):
    logging.info("Processing Kestrel VCF results...")
    indel_vcf = f"{output_dir}/output_indel.vcf"
    output_ins = f"{output_dir}/output_insertion.vcf"
    output_del = f"{output_dir}/output_deletion.vcf"
    
    # Filter the VCF to extract indels, insertions, and deletions
    filter_vcf(vcf_path, indel_vcf)
    filter_indel_vcf(indel_vcf, output_ins, output_del)
    
    # Read the filtered VCF files into dataframes
    names = read_vcf(vcf_path)
    vcf_insertion = pd.read_csv(output_ins, comment='#', delim_whitespace=True, header=None, names=names)
    vcf_deletion = pd.read_csv(output_del, comment='#', delim_whitespace=True, header=None, names=names)
    
    # Further processing logic can go here...
    logging.info("VCF processing completed.")
