#!/usr/bin/env python
# coding: utf-8

#pip3 install regex
#pip3 install biopython
#pip3 install setuptools==58
#pip3 install PyVCF


import subprocess as sp
import pandas as pd
import numpy as np
import sys,os,argparse,tempfile,shutil
from os.path import isfile, join
from pathlib import Path
import timeit
import regex as re

parser = argparse.ArgumentParser(description='Given raw fastq files, this tool genotype MUC1-VNTR using k-mer and profile-HMM based mathods')
parser.add_argument('-ref', '--reference_file', type=str, metavar='Referense', help='FASTA-formatted reference file and indexes', required=True)
parser.add_argument('-r1', '--fastq1', type=str, default=None, help='Fastq file first pair', required=True)
parser.add_argument('-r2', '--fastq2', type=str, default=None, help='Fastq file second pair', required=True)
parser.add_argument('-k', '--Kmer', type=str, default=None, help='Kmer size', required=True)
parser.add_argument('-o', '--output', type=str, default=None, help='Output file name', required=True)
parser.add_argument('-ref_VNTR', '--reference_VNTR', type=str, metavar='Referense', help='MUC1-specific reference file', required=True)
parser.add_argument('-t', '--threads', type=str, default=None, help='Number of threads (CPU)')
parser.add_argument('-p', '--tools_path', type=str, default=None, help='Path to the scripts directory', required=True)
parser.add_argument('-w', '--working_dir', type=str, default=None, help='the path to the output', required=True)


args = parser.parse_args()
tools_path = args.tools_path

if not os.path.exists(args.working_dir + args.output):
    os.mkdir(args.working_dir + args.output)

if not os.path.exists(args.working_dir + args.output + "/kmer_method"):
    os.mkdir(args.working_dir + args.output + "/kmer_method")

if not os.path.exists(args.working_dir + args.output + "/adVNTR"):
    os.mkdir(args.working_dir + args.output + "/adVNTR")
    
if not os.path.exists(args.working_dir + args.output + "/temp"):
    os.mkdir(args.working_dir + args.output + "/temp")

output = args.working_dir + args.output + "/"


welcomeMessage = """=========================================================
Genotyping MUC1-VNTR using mapping free and profile-HMM based methods
v. 1.0
This is free non-commercial software. 
=============================================================================
"""
endMessage = """========================
Thanks for using our pipeline!
========================================
"""

print (welcomeMessage)

start = timeit.default_timer()

# Kmer-based mapping free genotyping of MUC1-VNTR
fastq_1 = args.fastq1
fastq_2 = args.fastq2
reference_VNTR = args.reference_VNTR
vcf_out = output + "kmer_method/" + args.output + ".vcf"
sam_out = output + "kmer_method/" + args.output + ".sam"
vcf_path = Path(vcf_out)

if vcf_path.is_file():
    print ("VCF file already exists...")
else:
    if None not in (fastq_1, fastq_2, reference_VNTR):
        kmer_command = "java -Xmx10g -jar " + tools_path + "kestrel-1.0.1/" + "./kestrel.jar " +  "-k " +  args.Kmer + " -r " + reference_VNTR  + " -o " + vcf_out + " -p " + sam_out + " " + fastq_1 + " " + fastq_2 + " " + " --hapfmt " + "sam" + " && "  + "java -Xmx10g -jar " + tools_path + "picard.jar SortSam VALIDATION_STRINGENCY=SILENT " + "I=" + sam_out  + " " + "O=" + output + "kmer_method/" + args.output + ".bam " + " SORT_ORDER=coordinate"
        print ("Launching Kestrel!\n")
        print ("Kestrel Command: " + kmer_command  + "\n")
        process = sp.Popen(kmer_command , shell=True)
        process.wait()
        print ("Mapping-free genotyping of MUC1-VNTR done!\n")
    else:
        print (args.fastq1 + " is not an expected fastq file... skipped...")


with open(vcf_out, "r") as vcf_file, \
        open(output + "kmer_method/" + args.output + "_indel.vcf", "w") as indel_file:
    for line in vcf_file:
        if line[:2] == "##":
            indel_file.write(line)
        else:
            [_, _, _, ref, alt, *_] = line.split("\t")  
            if len(ref) == 1 and len(alt) != 1 or len(ref) != 1 and len(alt) == 1:
                indel_file.write(line)
                

with open(output + "kmer_method/" + args.output + "_indel.vcf", "r") as vcf_file, \
        open(output + "kmer_method/" + args.output + "_insertion.vcf", "w") as insertion_file, \
            open(output + "kmer_method/" + args.output + "_deletion.vcf", "w") as deletion_file:
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

def read_vcf(path):
    with open(path,'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    f.close()
    return vcf_names

names = read_vcf(vcf_out)
vcf_insertion = pd.read_csv(output + "kmer_method/" + args.output + "_insertion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
vcf_deletion = pd.read_csv(output + "kmer_method/" + args.output + "_deletion.vcf", comment='#', delim_whitespace=True, header=None, names=names)


from Bio import SeqIO
with open(reference_VNTR) as fasta_file: 
    identifiers = []
    seq = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        seq.append(seq_record.seq)
          
s1 = pd.Series(identifiers, name='ID')
s2 = pd.Series(seq, name='Sequence')
MUC1_ref = pd.DataFrame(dict(Motifs=s1, Motif_sequence=s2))

def preprocessing_insertion(df):
    df1=df.copy()
    df1.rename(columns={'#CHROM':'Motifs'}, inplace = True)
    df1.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    df1 = pd.merge(df1, MUC1_ref, on='Motifs', how='left')
    df1['Variant'] = 'Insertion'
    return df1

def preprocessing_deletion(df):
    df2=df.copy()
    df2.rename(columns={'#CHROM':'Motifs'}, inplace = True)
    df2.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1, inplace=True)
    df2 = pd.merge(df2, MUC1_ref, on='Motifs', how='left')
    df2['Variant'] = 'Deletion'
    return df2

insertion = preprocessing_insertion(vcf_insertion)
deletion = preprocessing_deletion(vcf_deletion)
vertical_concat = pd.concat([insertion, deletion], axis=0)

def final_processing(df):
    df3=df.copy()
    df3 = df3.rename(columns=lambda c: 'Depth' if c.endswith('\n') else c)
    df3.drop_duplicates(subset=['POS', 'REF', 'ALT'], inplace=True)
    df3[['Del', 'Estimated depth of all haplotypes supporting the alternate variant', 'Estimated depth of all haplotypes in the variant active region']] = df3['Depth'].str.split(':', expand=True)
    df3 = df3[['Motifs', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence', 'Estimated depth of all haplotypes supporting the alternate variant', 'Estimated depth of all haplotypes in the variant active region']]
    df3["ref_len"] = df3["REF"].str.len()
    df3["alt_len"] = df3["ALT"].str.len()
    df3["Codon_frame"] = (df3.alt_len - df3.ref_len) / 3
    df3['Codon_frame'] = df3['Codon_frame'].astype(str).apply(lambda x: x.replace('.0','C'))
    df3["TrueFalse"]=df3['Codon_frame'].str.contains('C', regex=True)
    df3["TrueFalse"] = df3["TrueFalse"].astype(str)
    df3 = df3[df3['TrueFalse'].str.contains("False")]
    df3.drop(['TrueFalse'], axis=1, inplace=True)

    return df3

Kmer_result = final_processing(vertical_concat)
Kmer_result.to_csv(output + "kmer_method/" + args.output + "_result.tsv", index=False, sep='\t')

#import pdfkit

f = open(output + "kmer_method/" + args.output + "_result.html",'w')
a = Kmer_result.to_html()
f.write(a)
f.close()

#pdfkit.from_file('NPH770.html', 'NPH770.pdf')

stop = timeit.default_timer()
print('Run Time (min): ', (stop - start)/60)  
