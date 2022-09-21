#!/usr/bin/env python
# coding: utf-8

#pip3 install regex
#pip3 install biopython
#pip3 install setuptools==58
#pip3 install PyVCF

from concurrent.futures import thread
import subprocess as sp
from xml.dom.expatbuilder import theDOMImplementation
import pandas as pd
import numpy as np
import sys,os,argparse,tempfile,shutil
from os.path import isfile, join
from pathlib import Path
import timeit
from pysam import SamtoolsError
import regex as re
import logging 

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
parser.add_argument('-m', '--reference_vntr', type=str, default=None, help='adVNTR reference vntr database', required=True)
parser.add_argument("--ignore_advntr", action="store_true", help="Skip adVNTR genotyping of MUC1-VNTR")

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


welcomeMessage = """==============================================
VNtyper: Genotyping MUC1 coding VNTR using mapping-free genotyping (fast) and profile-HMM-based (slow) algorithms
v. 1.0.0
This is free non-commercial software. 
=================================================================
"""
endMessage = """========================
Thanks for using VNtyper!
========================================
"""

print (welcomeMessage)

log_file = args.output + ".log"
logging.basicConfig(filename= log_file , level= logging.DEBUG , encoding= 'utf-8', format="%(asctime)s %(message)s")

logging.info('Fastq file quality control: deduplication, low quality read removal, and size correction')

start = timeit.default_timer()

# Fastq quality control (deduplication, low quality read removal, size correction...)
fastq_1 = args.fastq1
fastq_2 = args.fastq2

# Running fastp
if None not in (fastq_1, fastq_2):
    QC_command = tools_path + "./fastp "+ "--thread " + args.threads + " " + "--in1 " + fastq_1 + " --in2 " + fastq_2 + " --out1 " + args.working_dir + args.output + "/kmer_method/" + args.output + "_R1.fastq.gz" + " " + "--out2 " + args.working_dir + args.output + "/kmer_method/" + args.output + "_R2.fastq.gz" + " --phred64 --compression 6 --disable_adapter_trimming --dedup --dup_calc_accuracy 3 --length_required 120 --html " + output + "temp/" + args.output + ".html"
    print ("Quality control step!\n")
    process = sp.Popen(QC_command , shell=True)
    process.wait()
    print ("QC passed!\n")
    logging.info('Quality control passed..')
else:
    print (args.fastq1 + " is not an expected fastq file... skipped...")
    logging.error('Provided raw file is not the expected fastq file...skipped...')


# Kmer-based mapping free genotyping of MUC1-VNTR
logging.info('Kmer-based mapping free genotyping of MUC1-VNTR...')
reference_VNTR = args.reference_VNTR
vcf_out = output + "kmer_method/" + args.output + ".vcf"
sam_out = output + "kmer_method/" + args.output + ".sam"
vcf_path = Path(vcf_out)

fastq_1 = args.working_dir + args.output + "/kmer_method/" + args.output + "_R1.fastq.gz"
fastq_2 = args.working_dir + args.output + "/kmer_method/" + args.output + "_R2.fastq.gz"

if vcf_path.is_file():
    print ("VCF file already exists...")
else:
    if None not in (fastq_1, fastq_2, reference_VNTR):
        kmer_command = "java -Xmx10g -jar " + tools_path + "kestrel-1.0.1/" + "./kestrel.jar " +  "-k " +  args.Kmer + " -r " + reference_VNTR  + " -o " + vcf_out + " " + fastq_1 + " " + fastq_2 
        print ("Launching Kestrel!\n")
        print ("Kestrel Command: " + kmer_command  + "\n")
        process = sp.Popen(kmer_command , shell=True)
        process.wait()
        print ("Mapping-free genotyping of MUC1-VNTR done!\n")
        logging.info('Mapping-free genotyping of MUC1-VNTR done!')
    else:
        print (args.fastq1 + " is not an expected fastq file... skipped...")

logging.info('VCF file preprocessing...')
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

logging.info('Writing variants to the output file...')
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

stop = timeit.default_timer()
print ('Run Time (min): ', (stop - start)/60)  
#import pdfkit

f = open(output + "kmer_method/" + args.output + "_result.html",'w')
a = Kmer_result.to_html()
f.write(a)
f.close()

logging.info('Removing temporary and intermediate files...')
# Remove temporary intermediate files
rm_list = ["_deletion", "_indel", "_insertion"]
for db in rm_list:
    db_str = args.output + db + ".vcf"
    rm_command =  "rm " + args.working_dir + args.output + "/kmer_method/" + db_str
    process = sp.Popen(rm_command, shell=True)
    process.wait()
else:
    print (db_str + " is not found!")


# VNTR genotyping with adVNTR
if args.ignore_advntr:
    logging.info('MUC1 VNTR genotyping with adVNTR skippied...')
    print (endMessage)
    sys.exit()
else:
    print('launching adVNTR (This will take a while)...')
    logging.info('launching adVNTR (Profile-HMM based VNTR genotyping)...')

start_advntr = timeit.default_timer()

reference = args.reference_file
advntr_vcf_out = output + "adVNTR/" + args.output + ".vcf"
sam_out = output + "adVNTR/" + args.output + ".sam"
advntr_vcf_path = Path(advntr_vcf_out)

if advntr_vcf_path.is_file():
    print ("adVNTR VCF file already exists...")
else:
    if None not in (fastq_1, fastq_2, reference):
        Mapping_command = "bwa mem -t " + args.threads + " " + reference + " " + fastq_1 + " " +  fastq_2 + " -o " + sam_out +  " && "  + "java -Xmx10g -jar " + tools_path + "picard.jar SortSam VALIDATION_STRINGENCY=SILENT " + "I=" + sam_out  + " " + "O=" + output + "adVNTR/" + args.output + "_sorted.bam" + " SORT_ORDER=coordinate" + " && " + "samtools index " + output + "adVNTR/" + args.output + "_sorted.bam"
        print ("Launching BWA!\n")
        print ("BWA Command: " + Mapping_command  + "\n")
        process = sp.Popen(Mapping_command , shell=True)
        process.wait()
        print ("Mapping reads to the Chr1 Done!\n")
    else:
        print (args.fastq1 + " is not an expected fastq file... skipped...")

sorted_bam = output + "adVNTR/" + args.output + "_sorted.bam"
db_file_hg19 = args.reference_vntr

if None not in (sorted_bam , reference):
    adVNTR_command = "singularity exec " + tools_path + "code-adVNTR.sif advntr genotype -fs -vid 25561 --outfmt vcf --alignment_file " + sorted_bam + " -o " + output + "adVNTR/" + args.output + ".vcf" + " " + "-m " + db_file_hg19 + " -r " + args.reference_file + " --working_directory " + output + "adVNTR/"
    print('Launching adVNTR genotyping!\n')
    print('adVNTR command: ' + adVNTR_command)
    process = sp.Popen(adVNTR_command, shell=True)
    process.wait()
    print('adVNTR genotyping of MUC1-VNTR done!')
    logging.info('adVNTR genotyping of MUC1-VNTR done')
else:
    print('Input files are not expected files...skipped')
    logging.error('Input files are not the expected files for adVNTR...skipped')

logging.info('adVNTR result preprocessing...')
path = output + "adVNTR/" + args.output + ".vcf"
def read_vcf(path):
    with open(path,'r') as f:
        for line in f:
            if line.startswith("#VID"):
                vcf_names = [x for x in line.split('\t')]
                break
    f.close()
    return vcf_names

advntr_out = output + "adVNTR/" + args.output + ".vcf"

names = read_vcf(advntr_out)
df = pd.read_csv(advntr_out, comment='#', delim_whitespace=True, header=None, names=names)

def search(regex: str, df, case=False):
    #Search all the text columns of df
    textlikes = df.select_dtypes(include=[object, "object"])
    return df[
        textlikes.apply(
            lambda column: column.str.contains(regex, regex=True, case=case, na=False)
        ).any(axis=1)
    ]

def advntr_processing(df):
    df_ = df.copy()
    df_.rename(columns={'#VID':'#ID', 'State': 'Variant', 'Pvalue\n':'Pvalue'}, inplace = True)
    df_['Deletion_length'] = df_['Variant'].str.count('D').add(0).fillna(0)
    df_['Insertion'] = df_['Variant'].str.count('I').add(0).fillna(0)
    df_[['D', 'I']]= df_['Variant'].str.split('I', expand=True)
    df_.I = df_.I.fillna('LEN')
    df_[['I', 'Insertion_len']]= df_['I'].str.split('LEN', expand=True)
    df_.drop(['I', 'D', 'Insertion'], axis=1, inplace=True)
        
    return df_

if df.empty:
    print ('No variant found..')
else:
    df_ = advntr_processing(df)

IG = search('I22_2_G_LEN1', df_)


if IG.empty:
    print('DataFrame is empty!')
else:
    IG.to_csv(args.output+'I_G.bed', sep='\t', index=False)
    df_.to_csv(args.output+'_advntr.bed', sep='\t', index=False)

# Remove temporary intermediate files
rm_list = ["_deletion", "_indel", "_insertion"]
rm_list_1 = ["_R1", "_R2"]
for db in rm_list:
    db_str = args.output + db + ".vcf"
    rm_command =  "rm " + args.working_dir + args.output + "/kmer_method/" + db_str
    process = sp.Popen(rm_command, shell=True)
    process.wait()
else:
    print (db_str + " is not found!")

for db1 in rm_list_1:
    db1_str = args.output + db1 + ".fastq.gz"
    rm_command =  "rm " + args.working_dir + args.output + "/kmer_method/" + db1_str
    process = sp.Popen(rm_command, shell=True)
    process.wait()
else:
    print (db1_str + " is not found!")

print (endMessage)

stop_advntr = timeit.default_timer()
print('adVNTR Run Time (min): ', (stop_advntr - start_advntr)/60)
