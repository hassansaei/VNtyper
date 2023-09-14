#!/usr/bin/env python3.9 or Above
# coding: utf-8

#Python >= 3.9
#pip3 install regex
#pip3 install pandas
#pip3 install numpy
#pip3 install biopython
#pip3 install setuptools==58
#pip3 install pysam

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
import warnings
warnings.filterwarnings("ignore")

# Parse Arguments 
parser = argparse.ArgumentParser(description='Given raw fastq files, this pipeline genotype MUC1-VNTR using kestrel (Mapping-free genotyping) and Code-adVNTR mathods')
parser.add_argument('-ref', '--reference_file', type=str, metavar='Referense', help='FASTA-formatted reference file and indexes', required=True)
parser.add_argument('-r1', '--fastq1', type=str, default=None, help='Fastq file first pair', required=False)
parser.add_argument('-r2', '--fastq2', type=str, default=None, help='Fastq file second pair', required=False)
parser.add_argument('-o', '--output', type=str, default=None, help='Output file name', required=True)
parser.add_argument('-ref_VNTR', '--reference_VNTR', type=str, metavar='Referense', help='MUC1-specific reference file', required=True)
parser.add_argument('-t', '--threads', type=str, default=None, help='Number of threads (CPU)')
parser.add_argument('-p', '--tools_path', type=str, default=None, help='Path to the VNtyper directory', required=True)
parser.add_argument('-w', '--working_dir', type=str, default=None, help='the path to the output', required=True)
parser.add_argument('-m', '--reference_vntr', type=str, default=None, help='adVNTR reference vntr database', required=False)
parser.add_argument("--ignore_advntr", action="store_true", help="Skip adVNTR genotyping of MUC1-VNTR")
parser.add_argument("--bam", action="store_true", help="BAM file as an input")
parser.add_argument("--fastq", action="store_true", help="Paired-end fastq files as an input")
parser.add_argument('-a', '--alignment', type=str, default=None, help='Alignment File (with an index file .bai)', required=False)


args = parser.parse_args()
tools_path = args.tools_path

if not os.path.exists(args.working_dir + args.output):
    os.mkdir(args.working_dir + args.output) 

if not os.path.exists(args.working_dir + args.output + "/temp"):
    os.mkdir(args.working_dir + args.output + "/temp")

output = args.working_dir + args.output + "/"


welcomeMessage = """
==========================================================================================================
Given alignment (BAM) or raw file (FASTQ), this tool genotypes MUC1 coding-VNTR 
-- For rapid genotyping, BAM files are preferred!
-- User can Skip code-adVNTR genotyping step using --ignore_advntr option (This step will take a while..)
v. 1.0.0
This is free non-commercial software. 
==========================================================================================================
"""
endMessage = """
==============================
Thanks for using VNtyper pipeline!
Contact: hassan.saei@inserm.fr
==============================
"""
start = timeit.default_timer()
print (welcomeMessage)


# General Functions 
def search(regex: str, df, case=False):
    textlikes = df.select_dtypes(include=[object, "object"])
    return df[
        textlikes.apply(
            lambda column: column.str.contains(regex, regex=True, case=case, na=False)
        ).any(axis=1)
    ]
    
def read_vcf(path):
    with open(path,'r') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    f.close()
    return vcf_names

def filter_vcf(input_path, output_path):

    with open(input_path, "r") as vcf_file, open(output_path, "w") as indel_file:
        for line in vcf_file:
            if line[:2] == "##":
                indel_file.write(line)
            else:
                [_, _, _, ref, alt, *_] = line.split("\t")  
                if len(ref) == 1 and len(alt) != 1 or len(ref) != 1 and len(alt) == 1:
                    indel_file.write(line)

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


# Saving log file in temp/ directory
log_file = args.working_dir + args.output + "/temp/" + args.output + ".log"
logging.basicConfig(filename= log_file , level= logging.DEBUG , encoding= 'utf-8', format="%(asctime)s %(message)s")

# Fastq input - Fastq quality control (deduplication, low quality read removal, size correction...)
if args.fastq:
    fastq_1 = args.fastq1
    fastq_2 = args.fastq2
    # Running fastp
    if None not in (fastq_1, fastq_2):
        QC_command = "fastp "+ "--thread " + args.threads + " " + "--in1 " + fastq_1 + " --in2 " + fastq_2 + " --out1 " + output + args.output + "_R1.fastq.gz" + " " + "--out2 " + output + args.output + "_R2.fastq.gz" + " --compression 6 --disable_adapter_trimming --dedup --dup_calc_accuracy 3 --length_required 40 --html " + output + "temp/" + args.output + ".html"
        logging.info('Fastq file quality control: deduplication, low quality read removal, and size correction...')
        print ("Quality control step!\n")
        process = sp.Popen(QC_command , shell=True)
        process.wait()
        print ("QC passed!\n")
        logging.info('Quality control passed..')
    else:
        print (args.fastq1 + " is not an expected fastq file... skipped...")
        logging.error('Provided raw file is not the expected fastq file...skipped...')

# BAM input (BAM clean up and BAM to Fastq conversion...)
# hg19 MUC1 region chr1:155158000-155163000
if args.bam:
    in_bam = args.alignment
    #out_bam = output + args.output + "__mkdup.bam"
    #command_mkdup = "/WORKSPACE/VNtyper/Scripts/./sambamba-0.6.8 markdup -t " + args.threads + " " +  in_bam + " " + out_bam
    command_slice = "/SOFT/./sambamba-0.6.8 slice " + in_bam + " chr1:155158000-155163000 " + " -o " + output + args.output + "_chr1.bam"
    command_flag = "samtools view -u -f 4 -F264 -@ " + args.threads + " " + in_bam +  " > " + output + args.output + "_unmapped1.bam"  + " && " + "samtools view -u -f 8 -F260 -@ " + args.threads + " " + in_bam +  " > " + output + args.output + "_unmapped2.bam" + " && " +  "samtools view -u -f 12 -F256 -@ " + args.threads + " " + in_bam +  " > " + output + args.output + "_unmapped3.bam" 
    command_merge = "/SOFT/./sambamba-0.6.8 merge -t " + args.threads + " " + output + args.output + "_vntyper.bam" + " " + output + args.output + "_chr1.bam" + " " + output + args.output + "_unmapped1.bam" + " " + output + args.output + "_unmapped2.bam" + " " + output + args.output + "_unmapped3.bam"
    command_sort_fastq = "samtools sort -n " + " -@ " + args.threads + " " + output + args.output + "_vntyper.bam" +  " -o " + output + args.output + "_VN.bam" + " && " + "samtools fastq " + " -@ " + args.threads + " " +  output + args.output + "_VN.bam" +  " -1 " + output + args.output + "_R1.fastq.gz" + " -2 " + output + args.output + "_R2.fastq.gz"
    print ('BAM cleanup and converting to fastq...')
    logging.info('BAM file cleanup and converting to fastq...')
    #process = sp.Popen(command_mkdup , shell=True)
    #process.wait()
    process = sp.Popen(command_slice , shell=True)
    process.wait()
    process = sp.Popen(command_flag , shell=True)
    process.wait()
    process = sp.Popen(command_merge , shell=True)
    process.wait()
    process = sp.Popen(command_sort_fastq, shell=True)
    process.wait()
    logging.info('BAM to Fastq conversion finished!')
    print ('BAM to Fastq conversion finished!')
    command_rm = "rm " + output + args.output + "*.bam" + " && " + "rm " + output + args.output + "*.bai"
    process = sp.Popen(command_rm , shell=True)
    process.wait()
    fastq_1 = output + args.output + "_R1.fastq.gz"
    fastq_2 = output + args.output + "_R2.fastq.gz"


# Kestrel algorithm (Kmer frequency-based genotyping of MUC1-VNTR)
logging.info('Kmer-based mapping free genotyping of MUC1-VNTR (Kestrel)...')
reference_VNTR = args.reference_VNTR
vcf_out = output + args.output + ".vcf"
vcf_path = Path(vcf_out)

#fastq_1 = args.working_dir + args.output + "_R1.fastq.gz"
#fastq_2 = args.working_dir + args.output + "_R2.fastq.gz"

kmer_command_17 = "#java -Xmx15g -jar  /usr/local/lib/kestrel-1.0.1/kestrel.jar -k  17  "  + " -r " + reference_VNTR  + " -o " + vcf_out + " " + fastq_1 + " " + fastq_2 + " " + " --temploc " + args.working_dir + args.output + "/temp/" 
kmer_command_20 = "java  -Xmx15g -jar  /usr/local/lib/kestrel-1.0.1/kestrel.jar -k  20  --maxalignstates 30  --maxhapstates 30 "  + " -r " + reference_VNTR  + " -o " + vcf_out + "  " + fastq_1 + " " + fastq_2 + " " + " --temploc " + args.working_dir + args.output + "/temp/"  + " --hapfmt " + " sam " +  " -p " +  args.working_dir + args.output + "/temp/" + args.output + ".sam"
kmer_command_25 = "#java -Xmx15g -jar  /usr/local/lib/kestrel-1.0.1/kestrel.jar -k  25  "  + " -r " + reference_VNTR  + " -o " + vcf_out + "  " + fastq_1 + " " + fastq_2 + " " + " --temploc " + args.working_dir + args.output + "/temp/" 
kmer_command_41 = "#java -Xmx15g -jar  /usr/local/lib/kestrel-1.0.1/kestrel.jar -k  41  "  + " -r " + reference_VNTR  + " -o " + vcf_out + "  " + fastq_1 + " " + fastq_2 + " " + " --temploc " + args.working_dir + args.output + "/temp/" 

# Run Kestrel algorithm
if vcf_path.is_file():
    print ("VCF file already exists...")
    sys.exit()
else:
    if None not in (fastq_1, fastq_2, reference_VNTR):
        print ("Launching Kestrel...with Kmer size 20!\n")
        process = sp.Popen(kmer_command_20 , shell=True)
        process.wait()
        logging.info('Mapping-free genotyping of MUC1-VNTR with kmer size 20 done!')
    else:
        print (args.fastq1 + " is not an expected fastq file... skipped...")

# Kestrel output file processing
logging.info('VCF file preprocessing...')
input_path = vcf_out
output_path = output + args.output + "_indel.vcf"
indel_vcf = output + args.output + "_indel.vcf"
output_ins = output + args.output + "_insertion.vcf"
output_del = output + args.output + "_deletion.vcf"

filter_vcf(input_path, output_path)
filter_indel_vcf(indel_vcf, output_ins, output_del)

names = read_vcf(vcf_out)
vcf_insertion = pd.read_csv(output + args.output + "_insertion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
vcf_deletion = pd.read_csv(output + args.output + "_deletion.vcf", comment='#', delim_whitespace=True, header=None, names=names)

# MUC1 VNTR motif dictionary processing
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

def StepA_processing(df):
    df3=df.copy()
    df3 = df3.rename(columns=lambda c: 'Depth' if c.endswith('\n') else c)
    df3[['Del', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']] = df3['Depth'].str.split(':', expand=True)
    df3 = df3[['Motifs', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion']]
    df3["ref_len"] = df3["REF"].str.len()
    df3["alt_len"] = df3["ALT"].str.len()
    df3["Frame_Score"] = round((df3.alt_len - df3.ref_len) / 3, 2)
    df3['Frame_Score'] = df3['Frame_Score'].astype(str).apply(lambda x: x.replace('.0','C'))
    df3["TrueFalse"]=df3['Frame_Score'].str.contains('C', regex=True)
    df3["TrueFalse"] = df3["TrueFalse"].astype(str)
    df3 = df3[df3['TrueFalse'].str.contains("False")]

    return df3

Kmer_A = StepA_processing(vertical_concat)

def StepB_proccessing(df):
    df4=df.copy()
    df4.drop(['TrueFalse', 'ref_len', 'alt_len'], axis=1, inplace=True)
    df4[['left', 'Right']] = df4.Frame_Score.astype(str).str.split('.', expand=True)
    df4.left.replace('-0', '-1', inplace=True)

    return df4

# Label function 
def conditions(Kestrel_concat):
    if (Kestrel_concat['Depth_Score'] <= 0.00469) or (Kestrel_concat['Estimated_Depth_Variant_ActiveRegion'] <= 200):
        return 'Low_Precision'
    elif (Kestrel_concat['Estimated_Depth_AlternateVariant'] >= 21) and (Kestrel_concat['Depth_Score'] >= 0.00469) and (Kestrel_concat['Depth_Score'] <= 0.00515):
        return 'Low_Precision'
    elif (Kestrel_concat['Estimated_Depth_AlternateVariant'] > 100):
        return 'High_Precision'
    elif (Kestrel_concat['Estimated_Depth_AlternateVariant'] <= 20):
        return 'Low_Precision'
    elif (Kestrel_concat['Estimated_Depth_AlternateVariant'] >= 21) and (Kestrel_concat['Estimated_Depth_AlternateVariant'] < 100) and (Kestrel_concat['Depth_Score'] >= 0.00515):
        return 'High_Precision'
    elif (Kestrel_concat['Estimated_Depth_AlternateVariant'] >= 100) and (Kestrel_concat['Depth_Score'] >= 0.00515):
        return 'High_Precision*'

# Extract good frameshitfs (3n+1 for insertion and 3n+2 for deletion)
def process_kmer(Kmer_A):
    if Kmer_A.empty:
        print ('Proceeding to the next Kmer size...')
        Kestrel_concat = Kmer_A
    else:
        Kmer_B = StepB_proccessing(Kmer_A)
        Ins = Kmer_B[Kmer_B["left"].apply(lambda x: '-' not in x) &  Kmer_B["Right"].apply(lambda y: '33' in y)]
        Del = Kmer_B[Kmer_B["left"].apply(lambda x: '-' in x) &  Kmer_B["Right"].apply(lambda y: '67' in y)]
        Kestrel_concat = pd.concat([Ins, Del], axis=0)
        Kestrel_concat['Estimated_Depth_AlternateVariant'] = Kestrel_concat['Estimated_Depth_AlternateVariant'].astype(int)
        Kestrel_concat['Estimated_Depth_Variant_ActiveRegion'] = Kestrel_concat['Estimated_Depth_Variant_ActiveRegion'].astype(int)
        Kestrel_concat['Depth_Score'] = round(Kestrel_concat.Estimated_Depth_AlternateVariant / Kestrel_concat.Estimated_Depth_Variant_ActiveRegion, 5)
        Kestrel_concat['Depth_Score'] = Kestrel_concat['Depth_Score'].astype(float)
        Kestrel_concat['Confidence'] = Kestrel_concat.apply(conditions, axis=1)
        if Kestrel_concat['ALT'].str.contains(r'\bGG\b').any():
            Kestrel_concat = pd.concat([
                Kestrel_concat[Kestrel_concat['ALT'] != 'GG'],
                Kestrel_concat[Kestrel_concat['ALT'] == 'GG'].loc[Kestrel_concat['Depth_Score'] >= 0.00469]
            ])
        cg_exists = Kestrel_concat['ALT'].str.contains(r'\bCG\b').any()
        tg_exists = Kestrel_concat['ALT'].str.contains(r'\bTG\b').any()
        if cg_exists and tg_exists:
            Kestrel_concat = Kestrel_concat[~Kestrel_concat['ALT'].isin(['CG', 'TG'])]
        Kestrel_concat = Kestrel_concat[Kestrel_concat['Confidence'] != 'Red_Zone']
        Kestrel_concat.drop(['left', 'Right'], axis=1, inplace=True)
        
    return Kestrel_concat

Kestrel_concat = process_kmer(Kmer_A)

# Running Kmer size 17 and variant processing...
if Kestrel_concat.empty:
    #logging.info('Running Kestrel with Kmer size 17...')
    #print ('Running Kestrel with kmer size 17...')
    process = sp.Popen(kmer_command_17 , shell=True)
    process.wait()
    logging.info('Mapping-free genotyping of MUC1-VNTR with kmer size 20 done!')
    filter_vcf(output + args.output + ".vcf", output + args.output + "_indel.vcf")
    filter_indel_vcf(output + args.output + "_indel.vcf", output + args.output + "_insertion.vcf", output + args.output + "_deletion.vcf")
    names = read_vcf(output + args.output + ".vcf")
    vcf_insertion = pd.read_csv(output + args.output + "_insertion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
    vcf_deletion = pd.read_csv(output + args.output + "_deletion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
    insertion = preprocessing_insertion(vcf_insertion)
    deletion = preprocessing_deletion(vcf_deletion)
    vertical_concat = pd.concat([insertion, deletion], axis=0)
    Kmer_A = StepA_processing(vertical_concat)
    Kestrel_concat = process_kmer(Kmer_A)
    which_kmer = '17'
else:
    which_kmer = '20'

# Running Kmer size 25 and variant processing...
if Kestrel_concat.empty:
    #logging.info('Running Kestrel with kmer size 25...')
    #print ('Running Kestrel with kmer size 25...')
    process = sp.Popen(kmer_command_25 , shell=True)
    process.wait()
    logging.info('Mapping-free genotyping of MUC1-VNTR with kmer size 25 done!')
    filter_vcf(output + args.output + ".vcf", output + args.output + "_indel.vcf")
    filter_indel_vcf(output + args.output + "_indel.vcf", output + args.output + "_insertion.vcf", output + args.output + "_deletion.vcf")
    names = read_vcf(output + args.output + ".vcf")
    vcf_insertion = pd.read_csv(output + args.output + "_insertion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
    vcf_deletion = pd.read_csv(output + args.output + "_deletion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
    insertion = preprocessing_insertion(vcf_insertion)
    deletion = preprocessing_deletion(vcf_deletion)
    vertical_concat = pd.concat([insertion, deletion], axis=0)
    Kmer_A = StepA_processing(vertical_concat)
    Kestrel_concat = process_kmer(Kmer_A)
    which_kmer = '25'

# Running Kmer size 41 and variant processing...
if Kestrel_concat.empty:
    #logging.info('Running Kestrel with kmer size 41...')
    #print ('Running Kestrel with kmer size 41...')
    process = sp.Popen(kmer_command_41 , shell=True)
    process.wait()
    logging.info('Mapping-free genotyping of MUC1-VNTR with kmer size 41 done!')
    filter_vcf(output + args.output + ".vcf", output + args.output + "_indel.vcf")
    filter_indel_vcf(output + args.output + "_indel.vcf", output + args.output + "_insertion.vcf", output + args.output + "_deletion.vcf")
    names = read_vcf(output + args.output + ".vcf")
    vcf_insertion = pd.read_csv(output + args.output + "_insertion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
    vcf_deletion = pd.read_csv(output + args.output + "_deletion.vcf", comment='#', delim_whitespace=True, header=None, names=names)
    insertion = preprocessing_insertion(vcf_insertion)
    deletion = preprocessing_deletion(vcf_deletion)
    vertical_concat = pd.concat([insertion, deletion], axis=0)
    Kmer_A = StepA_processing(vertical_concat)
    Kestrel_concat = process_kmer(Kmer_A)
    which_kmer = '41'


# MUC1 VNTR motif and RU processing
from Bio import SeqIO
with open(args.tools_path + 'Files/code-adVNTR_RUs.fa') as RU_file, open(args.tools_path + 'Files/MUC1_motifs_Rev_com.fa') as Motif_file: 
    identifiers = []
    seq = []
    for seq_record1 in SeqIO.parse(RU_file, 'fasta'):  # (generator)
        identifiers.append(seq_record1.id)
        seq.append(seq_record1.seq.upper())
    for seq_record2 in SeqIO.parse(Motif_file, 'fasta'):  # (generator)
        identifiers.append(seq_record2.id)
        seq.append(seq_record2.seq.upper())
        
s1 = pd.Series(identifiers, name='ID')
s2 = pd.Series(seq, name='Sequence')
Merged_motifs = pd.DataFrame(dict(Motif=s1, Motif_sequence=s2))

# Motif correction and annotation
if Kestrel_concat.empty:
    print('No pathogenic variant was found!')
    logging.info('No pathogenic variant was found...')
    which_kmer = ''
else:
    Kestrel_concat[['Motif_left', 'Motif_right']] = Kestrel_concat['Motifs'].str.split('-', expand= True)
    Kestrel_concat['POS'] = Kestrel_concat['POS'].astype(int)
    motif_left = Kestrel_concat[Kestrel_concat['POS'] < 60]
    motif_right = Kestrel_concat[Kestrel_concat['POS'] >= 60]
    motif_left.rename(columns={'Motif_right': 'Motif'}, inplace = True)
    motif_left.drop(['Motif_sequence'], axis=1, inplace=True)
    motif_left = motif_left.merge(Merged_motifs, on='Motif')
    #motif_left.drop_duplicates(subset=['REF', 'ALT'], inplace=True)
    motif_left = motif_left[['Motif', 'Variant', 'POS', 'REF', 'ALT','Motif_sequence', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion', 'Depth_Score', 'Confidence']]
    motif_left = motif_left.sort_values('Depth_Score', ascending=False)
    motif_left = motif_left.drop_duplicates('ALT', keep='first')
    motif_left = motif_left.sort_values('POS', ascending=False)
    motif_left = motif_left.tail(1)
    motif_right.rename(columns={'Motif_left': 'Motif'}, inplace = True)
    motif_right.drop(['Motif_sequence'], axis=1, inplace=True)
    motif_right = motif_right.merge(Merged_motifs, on='Motif')
    motif_right = motif_right[['Motif', 'Variant', 'POS', 'REF', 'ALT','Motif_sequence', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion', 'Depth_Score', 'Confidence']]
    if motif_right['ALT'].str.contains(r'\bGG\b').any():
        motif_right = motif_right.loc[(motif_right['Motif'] != 'Q') & (motif_right['Motif'] != '8') & (motif_right['Motif'] != '9') &
             (motif_right['Motif'] != '7') & (motif_right['Motif'] != '6p') & (motif_right['Motif'] != '6') & (motif_right['Motif'] != 'V') &
              (motif_right['Motif'] != 'J') & (motif_right['Motif'] != 'I') & (motif_right['Motif'] != 'G') & (motif_right['Motif'] != 'E') & (motif_right['Motif'] != 'A')]
        motif_right = motif_right.loc[motif_right['ALT'] == 'GG']
        motif_right = motif_right.sort_values('Depth_Score', ascending=False)
        motif_right = motif_right.drop_duplicates('ALT', keep='first')
        if motif_right['Motif'].str.contains('X').any():
            motif_right = motif_right[motif_right['Motif'] == 'X']
    else:
        motif_right = motif_right.sort_values('Depth_Score', ascending=False)
        motif_right = motif_right.drop_duplicates('ALT', keep='first')
    motif_right.drop_duplicates(subset=['REF', 'ALT'], inplace=True)
    Kestrel_concat = pd.concat([motif_right, motif_left])
    Kestrel_concat = Kestrel_concat.loc[(Kestrel_concat['ALT'] != 'CCGCC') & (Kestrel_concat['ALT'] != 'CGGCG') & (Kestrel_concat['ALT'] != 'CGGCC')]
    Kestrel_concat = Kestrel_concat[(Kestrel_concat['Motif'] != '6') & (Kestrel_concat['Motif'] != '6p') & (Kestrel_concat['Motif'] != '7')]
    Kestrel_concat['POS'] = Kestrel_concat['POS'].astype(int)
    Kestrel_concat.update(Kestrel_concat['POS'].mask(Kestrel_concat['POS'] >= 60, lambda x: x-60))
    Kestrel_concat['Kmer_Size'] = which_kmer


# Save Kestrel processed result
Kestrel_concat.to_csv(output + args.output + '_pre_result.tsv', sep='\t', index=False)
logging.info('Kestrel output processing done...')

stop = timeit.default_timer()
time = (stop - start)/60
print('-------- %s minutes -------' % time)
logging.info('Running time for Kestrel %s minutes.' % time)

# VNTR genotyping with code-adVNTR
columns_kmer = ['Motif', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence', 'Estimated_Depth_AlternateVariant', 'Estimated_Depth_Variant_ActiveRegion', 'Depth_Score', 'Confidence', 'Kmer_Size']

columns_adVNTR = ['#VID', 'Variant','NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']

rm_list_1 = ['_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_indel.vcf']
rm_list_3 = ["_R1.fastq.gz", "_R2.fastq.gz"]

# Writing output to a file
if args.ignore_advntr:
    logging.info('MUC1 VNTR genotyping with adVNTR skippied...')
    if Kestrel_concat.empty:
        with open(output + args.output + '_Final_result.tsv', 'w') as f3:
            f3.write('## VNtyper_Analysis_for_%s \n' % args.output )
            f3.write('# Kestrel_Result version 1.1.1 \n')
            f3.write('\t'.join(columns_kmer) + '\n') 
        for db in rm_list_1:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
        for db in rm_list_3:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
        print (endMessage)
        sys.exit()
    else:
        with open (output + args.output  + '_pre_result.tsv','r') as f2:
            with open(output + args.output + '_Final_result.tsv', 'w') as f3:
                f3.write('## VNtyper_Analysis_for_%s \n' % args.output )
                f3.write('# Kestrel_Result \n')
                f3.write('\t'.join(columns_kmer) + '\n') 
                next(f2)
                for line in f2:
                    f3.write(line)
        for db in rm_list_1:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
        for db in rm_list_3:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
        logging.info('The final result is saved in *_Final_result.tsv')
        print (endMessage)
        sys.exit()
else:
    if not os.path.exists(args.working_dir + args.output + "/adVNTR"):
        os.mkdir(args.working_dir + args.output + "/adVNTR")
    if args.alignment is not None:
        print('launching adVNTR (This will take a while)...')
        logging.info('launching code-adVNTR (Profile-HMM based VNTR genotyping)...')
        sorted_bam = args.alignment
        reference = args.reference_file
        db_file_hg19 = args.reference_vntr
        if None not in (sorted_bam , reference):
            adVNTR_command = "advntr genotype -fs -vid 25561 --outfmt vcf --alignment_file " + sorted_bam + " -o " + output + args.output + "_adVNTR.vcf" + " " + "-m " + db_file_hg19 + " -r " + args.reference_file + " --working_directory " + output
            print('Launching adVNTR genotyping!\n')
            process = sp.Popen(adVNTR_command, shell=True)
            process.wait()
            print('adVNTR genotyping of MUC1-VNTR done!')
            logging.info('adVNTR genotyping of MUC1-VNTR done')
        else:
            print('Input files are not expected files...skipped')
            logging.error('Input files are not the expected files for code-adVNTR...skipped')
    else:
        reference = args.reference_file
        advntr_vcf_out = output  + args.output + "_adVNTR.vcf"
        sam_out = output  + args.output + ".sam"
        advntr_vcf_path = Path(advntr_vcf_out)

        if advntr_vcf_path.is_file():
            print ("adVNTR VCF file already exists...")
        else:
            if None not in (fastq_1, fastq_2, reference):
                Mapping_command = "bwa mem -t " + args.threads + " " + reference + " " + fastq_1 + " " +  fastq_2 + " -o " + sam_out +  " && "  + "java -Xmx10g -jar /usr/local/lib/picard_2.27.4/picard.jar SortSam VALIDATION_STRINGENCY=SILENT " + "I=" + sam_out  + " " + "O=" + output + args.output + "_sorted.bam" + " SORT_ORDER=coordinate" + " && " + "samtools index " + output + args.output + "_sorted.bam"
                print ("Launching BWA!\n")
                process = sp.Popen(Mapping_command , shell=True)
                process.wait()
                print ("Alignment to the Chr1 Done!\n")
            else:
                print (args.fastq1 + " is not an expected fastq file... skipped...")
                logging.info('launching code-adVNTR (Profile-HMM based VNTR genotyping)...')

        sorted_bam = output + args.output + "_sorted.bam"
        reference = args.reference_file
        db_file_hg19 = args.reference_vntr
        if None not in (sorted_bam , reference):
            print('Launching code-adVNTR genotyping!\n')
            adVNTR_command = "advntr genotype -fs -vid 25561 --outfmt vcf --alignment_file " + sorted_bam + " -o " + output + args.output + "_adVNTR.vcf" + " " + "-m " + db_file_hg19 + " -r " + args.reference_file + " --working_directory " + output
            process = sp.Popen(adVNTR_command, shell=True)
            process.wait()
            print('code-adVNTR genotyping of MUC1-VNTR done!')
            logging.info('code-adVNTR genotyping of MUC1-VNTR done')
        else:
           print('Input files are not expected files...skipped') 
           logging.error('Input files are not the expected files for code-adVNTR...skipped')

logging.info('code-adVNTR result preprocessing...')
path = output + args.output + "_adVNTR.vcf"

def read_vcf(path):
    with open(path,'r') as f:
        for line in f:
            if line.startswith("#VID"):
                vcf_names = [x for x in line.split('\t')]
                break
    f.close()
    return vcf_names

advntr_out = output + args.output + "_adVNTR.vcf"

names = read_vcf(advntr_out)
df = pd.read_csv(advntr_out, comment='#', delim_whitespace=True, header=None, names=names)

# Array of ORF frameshifts for both insertions and deletions
ins_frame = np.arange(100) * 3 + 1
ins_frame = ins_frame.astype(str)
del_frame = np.arange(100) * 3 + 2
del_frame = del_frame.astype(str)

# Processing code-adVNTR output
def advntr_processing_del(df):
    
    df1 = df.copy()
    df1.rename(columns={'State': 'Variant', 'Pvalue\n':'Pvalue'}, inplace = True)
    df1['Deletion_length'] = df1['Variant'].str.count('D').add(0).fillna(0)
    df1['Insertion'] = df1['Variant'].str.count('I').add(0).fillna(0)
    df1['Insertion_len'] = df1['Variant'].str.extract('(LEN.*)')
    df1.Insertion_len = df1.Insertion_len.fillna('LEN')
    df1[['I', 'Insertion_len']]= df1['Insertion_len'].str.split('LEN', expand=True)
    df1.Insertion_len = df1.Insertion_len.replace('', 0)
    df1.Deletion_length = df1.Deletion_length.fillna('0')
    df1.Deletion_length = df1.Deletion_length.astype(int)
    df1.Insertion_len = df1.Insertion_len.astype(int)
    df1['frame'] = abs(df1.Insertion_len - df1.Deletion_length)
    df1.frame = df1['frame'].astype(str)
    df1 = df1.loc[(df1['Deletion_length'] >= 1)]
    df1 = df1[df1['frame'].isin(del_frame)]

    return df1

def advntr_processing_ins(df):
    
    dff = df.copy()
    dff.rename(columns={'State': 'Variant', 'Pvalue\n':'Pvalue'}, inplace = True)
    dff['Deletion_length'] = dff['Variant'].str.count('D').add(0).fillna(0)
    dff['Insertion'] = dff['Variant'].str.count('I').add(0).fillna(0)
    dff['Insertion_len'] = dff['Variant'].str.extract('(LEN.*)')
    dff.Insertion_len = dff.Insertion_len.fillna('LEN')
    dff[['I', 'Insertion_len']]= dff['Insertion_len'].str.split('LEN', expand=True)
    dff.Insertion_len = dff.Insertion_len.replace('', 0)
    dff.Deletion_length = dff.Deletion_length.fillna('0')
    dff.Deletion_length = dff.Deletion_length.astype(int)
    dff.Insertion_len = dff.Insertion_len.astype(int)
    dff['frame'] = abs(dff.Insertion_len - dff.Deletion_length)
    dff.frame = dff['frame'].astype(str)
    dff = dff.loc[(dff['Insertion_len'] >= 1)]
    dff = dff[dff['frame'].isin(ins_frame)]
         
    return dff

rm_list_4 = ['_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_indel.vcf', '.sam', '_sorted.bam', '_sorted.bam.bai', '_R1.fastq.gz', '_R2.fastq.gz', '_adVNTR.vcf', '_adVNTR_result.tsv', '.vcf']
rm_list_5 = ['_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_indel.vcf', '_R1.fastq.gz', '_R2.fastq.gz', '_adVNTR.vcf', '_adVNTR_result.tsv', '.vcf']

if df.empty:
    print('No pathogenic varinat was found with code-adVNTR!')
    with open(output + args.output  + '_pre_result.tsv','r') as f:
        with open(output  + args.output + '_Final_result.tsv', 'w') as f1:
            f1.write('## VNtyper_Analysis_for_%s \n' % args.output )
            f1.write('# Kestrel_Result \n')
            f1.write('\t'.join(columns_kmer) + '\n') 
            next(f)
            for line in f:
                f1.write(line)
    for db in rm_list_4:
        db_str = args.output + db
        rm_command =  "rm " + output + db_str
        process = sp.Popen(rm_command, shell=True)
        process.wait()
    logging.info('The final result is saved in *_Final_result.tsv')
    sys.exit()
else:
    df_del = advntr_processing_del(df)
    df_ins = advntr_processing_ins(df)

advntr_concat = pd.concat([df_del, df_ins], axis=0)
advntr_concat = advntr_concat[['#VID', 'Variant', 'NumberOfSupportingReads', 'MeanCoverage', 
                              'Pvalue']]
                              
advntr_concat.drop_duplicates(subset=['#VID', 'Variant', 'NumberOfSupportingReads'], inplace=True)

advntr_concat.to_csv(output + args.output + '_adVNTR_result.tsv', sep='\t', index=False)

# Saving output
with open(output + args.output  + '_pre_result.tsv' ,'r') as f, open(output + args.output +'_adVNTR_result.tsv', 'r') as ff:
    with open(output + args.output + '_Final_result.tsv', 'w') as f1:
        f1.write('## VNtyper_Analysis_for_%s \n' % args.output )
        f1.write('# Kestrel_Result \n')
        f1.write('\t'.join(columns_kmer) + '\n') 
        next(f)
        for line in f:
            f1.write(line)
        f1.write('# Code-adVNTR_Result \n')
        f1.write('\t'.join(columns_adVNTR) + '\n')
        next(ff)
        for  line in ff:
            f1.write(line)
    if args.alignment is not None:
        for db in rm_list_5:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()
    else:
        for db in rm_list_4:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()

print (endMessage)
logging.info('The final result is saved in *_Final_result.tsv')
sys.exit() 
