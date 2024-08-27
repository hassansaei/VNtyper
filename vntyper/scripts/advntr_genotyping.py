import subprocess as sp
import pandas as pd

def run_advntr(reference, db_file_hg19, sorted_bam, output, output_name):
    advntr_command = f"advntr genotype -fs -vid 25561 --outfmt vcf --alignment_file {sorted_bam} -o {output}{output_name}_adVNTR.vcf " \
                     f"-m {db_file_hg19} -r {reference} --working_directory {output}"
    print("Launching adVNTR genotyping!\n")
    process = sp.Popen(advntr_command, shell=True)
    process.wait()
    print("adVNTR genotyping of MUC1-VNTR done!")

def process_advntr_output(vcf_path):
    # Process adVNTR results
    pass
