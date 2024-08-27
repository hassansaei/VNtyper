import pandas as pd

def read_vcf(path):
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    return vcf_names

def filter_vcf(input_path, output_path):
    with open(input_path, "r") as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line.startswith("##"):
                indel_file.write(line)
            else:
                _, _, _, ref, alt, *_ = line.split("\t")
                if len(ref) == 1 and len(alt) != 1 or len(ref) != 1 and len(alt) == 1:
                    indel_file.write(line)

def filter_indel_vcf(indel_vcf, output_ins, output_del):
    with open(indel_vcf, "r") as vcf_file, open(output_ins, "w") as insertion_file, open(output_del, "w") as deletion_file:
        for line in vcf_file:
            if line.startswith("##"):
                insertion_file.write(line)
                deletion_file.write(line)
            else:
                _, _, _, ref, alt, *_ = line.split("\t")
                if len(ref) == 1 and len(alt) > 1:
                    insertion_file.write(line)
                else:
                    deletion_file.write(line)
