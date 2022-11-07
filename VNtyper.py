#!/usr/bin/env python
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

# Parse Arguments 
parser = argparse.ArgumentParser(description='Given raw fastq files, this pipeline genotype MUC1-VNTR using kestrel (Mapping-free genotyping) and Code-adVNTR mathods')
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
parser.add_argument('-a', '--alignment', type=str, default=None, help='Alignment File (Indexed BAM file)', required=False)


args = parser.parse_args()
tools_path = args.tools_path

if not os.path.exists(args.working_dir + args.output):
    os.mkdir(args.working_dir + args.output) 

if not os.path.exists(args.working_dir + args.output + "/temp"):
    os.mkdir(args.working_dir + args.output + "/temp")

output = args.working_dir + args.output + "/"


welcomeMessage = """
======================================================================
Given raw fastq files, this tool genotype MUC1 coding-VNTR using kestrel (Mapping-free genotyping) and Code-adVNTR mathods
Kestrel is faster than code-adVNTR.
It is recommended to run both algorithms to decrease False Neagtives.
v. 1.0.0
This is free non-commercial software. 
======================================================================
"""
endMessage = """
==============================
Thanks for using VNtyper pipeline!
==============================
"""
start = timeit.default_timer()
print (welcomeMessage)


# Functions
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

# Saving log file in temp/ directory
log_file = args.working_dir + args.output + "/temp/" + args.output + ".log"
logging.basicConfig(filename= log_file , level= logging.DEBUG , encoding= 'utf-8', format="%(asctime)s %(message)s")

logging.info('Fastq file quality control: deduplication, low quality read removal, and size correction...')

# Fastq quality control (deduplication, low quality read removal, size correction...)
fastq_1 = args.fastq1
fastq_2 = args.fastq2

# Running fastp
if None not in (fastq_1, fastq_2):
    QC_command = tools_path + "./fastp "+ "--thread " + args.threads + " " + "--in1 " + fastq_1 + " --in2 " + fastq_2 + " --out1 " + output + args.output + "_R1.fastq.gz" + " " + "--out2 " + output + args.output + "_R2.fastq.gz" + " --compression 6 --disable_adapter_trimming --dedup --dup_calc_accuracy 3 --length_required 40 --html " + output + "temp/" + args.output + ".html"
    print ("Quality control step!\n")
    process = sp.Popen(QC_command , shell=True)
    process.wait()
    print ("QC passed!\n")
    logging.info('Quality control passed..')
else:
    print (args.fastq1 + " is not an expected fastq file... skipped...")
    logging.error('Provided raw file is not the expected fastq file...skipped...')


# Kestrel algorithm (Kmer frequency-based genotyping of MUC1-VNTR)
logging.info('Kmer-based mapping free genotyping of MUC1-VNTR (Kestrel)...')
reference_VNTR = args.reference_VNTR
vcf_out = output + args.output + ".vcf"
vcf_path = Path(vcf_out)

#fastq_1 = args.working_dir + args.output + "_R1.fastq.gz"
#fastq_2 = args.working_dir + args.output + "_R2.fastq.gz"

# Run Kestrel
if vcf_path.is_file():
    print ("VCF file already exists...")
else:
    if None not in (fastq_1, fastq_2, reference_VNTR):
        kmer_command = "java -Xmx15g -jar " + tools_path + "kestrel-1.0.1/" + "./kestrel.jar " +  "-k " +  args.Kmer + " -r " + reference_VNTR  + " -o " + vcf_out + " " + fastq_1 + " " + fastq_2 + " " + " --temploc " + args.working_dir + args.output + "/temp/" + " " + " --noseqfilter"
        print ("Launching Kestrel!\n")
        process = sp.Popen(kmer_command , shell=True)
        process.wait()
        print ("Mapping-free genotyping of MUC1-VNTR done!\n")
        logging.info('Mapping-free genotyping of MUC1-VNTR done!')
    else:
        print (args.fastq1 + " is not an expected fastq file... skipped...")

# Kestrel output file processing
logging.info('VCF file preprocessing...')
with open(vcf_out, "r") as vcf_file, \
        open(output + args.output + "_indel.vcf", "w") as indel_file:
    for line in vcf_file:
        if line[:2] == "##":
            indel_file.write(line)
        else:
            [_, _, _, ref, alt, *_] = line.split("\t")  
            if len(ref) == 1 and len(alt) != 1 or len(ref) != 1 and len(alt) == 1:
                indel_file.write(line)
                

with open(output + args.output + "_indel.vcf", "r") as vcf_file, \
        open(output + args.output + "_insertion.vcf", "w") as insertion_file, \
            open(output + args.output + "_deletion.vcf", "w") as deletion_file:
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
    df3[['Del', 'Estimated depth of all haplotypes supporting the alternate variant', 'Estimated depth of all haplotypes in the variant active region']] = df3['Depth'].str.split(':', expand=True)
    df3 = df3[['Motifs', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence', 'Estimated depth of all haplotypes supporting the alternate variant', 'Estimated depth of all haplotypes in the variant active region']]
    df3["ref_len"] = df3["REF"].str.len()
    df3["alt_len"] = df3["ALT"].str.len()
    df3["Frame_Score"] = round((df3.alt_len - df3.ref_len) / 3, 2)
    df3['Frame_Score'] = df3['Frame_Score'].astype(str).apply(lambda x: x.replace('.0','C'))
    df3["TrueFalse"]=df3['Frame_Score'].str.contains('C', regex=True)
    df3["TrueFalse"] = df3["TrueFalse"].astype(str)

    return df3

Kmer_A = StepA_processing(vertical_concat)
dupC_in_frame = search('GGGGGGGG', Kmer_A)

if dupC_in_frame.empty:
    print('')
else:
    dupC_in_frame.drop(['TrueFalse', 'ref_len', 'alt_len', 'Frame_Score'], axis=1, inplace=True)


def StepB_proccessing(df):
    df4=df.copy()
    df4 = df4[df4['TrueFalse'].str.contains("False")]
    df4.drop(['TrueFalse', 'ref_len', 'alt_len'], axis=1, inplace=True)
    df4[['left', 'Right']] = df4.Frame_Score.astype(str).str.split('.', expand=True)
    df4.left.replace('-0', '-1', inplace=True)
    
    return df4

Kmer_B = StepB_proccessing(Kmer_A)

# Extract good frameshitfs (3n+1 for insertion and 3n+2 for deletion)
Ins = Kmer_B[Kmer_B["left"].apply(lambda x: '-' not in x) & 
                        Kmer_B["Right"].apply(lambda y: '33' in y)]

Del = Kmer_B[Kmer_B["left"].apply(lambda x: '-' in x) & 
                        Kmer_B["Right"].apply(lambda y: '66' in y)]

Kestrel_concat = pd.concat([Ins, Del], axis=0)
Kestrel_concat.drop(['left', 'Right'], axis=1, inplace=True)

# MUC1 VNTR motif and RU processing
from Bio import SeqIO
with open('Files/code-adVNTR_RUs.fa') as RU_file, open('Files/MUC1_motifs_Rev_com.fa') as Motif_file: 
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

    if dupC_in_frame.empty:
        print ('No Disease-causing varinat found with Kestrel')
    else:
        dupC_in_frame[['Motif_left', 'Motif_right']] = dupC_in_frame['Motifs'].str.split('-', expand= True)
        
        for i in dupC_in_frame.POS:
            if i >= 60:
                dupC_in_frame['POS'] = i-60
                dupC_in_frame.rename(columns={'Motif_right': 'Motif'}, inplace = True)
                dupC_in_frame.drop(['Motif_sequence'], axis=1, inplace=True)
                dupC_in_frame = dupC_in_frame.merge(Merged_motifs, on='Motif')
                dupC_in_frame = dupC_in_frame[['Motif', 'Variant', 'POS', 'REF', 'ALT','Motif_sequence', 'Estimated depth of all haplotypes supporting the alternate variant',
                'Estimated depth of all haplotypes in the variant active region']]
                dupC_in_frame.to_csv(output + args.output + '_dupC_inframe.tsv', sep='\t', index = False)
            else:
                dupC_in_frame['POS'] = i
                dupC_in_frame.rename(columns={'Motif_left': 'Motif'}, inplace = True)
                dupC_in_frame.drop(['Motif_sequence'], axis=1, inplace=True)
                dupC_in_frame = dupC_in_frame.merge(Merged_motifs, on='Motif')
                dupC_in_frame = dupC_in_frame[['Motif', 'Variant', 'POS', 'REF', 'ALT','Motif_sequence', 'Estimated depth of all haplotypes supporting the alternate variant',
                'Estimated depth of all haplotypes in the variant active region']]
                dupC_in_frame.to_csv(output + args.output + '_dupC_inframe.tsv', sep='\t', index = False)
else:
    Kestrel_concat[['Motif_left', 'Motif_right']] = Kestrel_concat['Motifs'].str.split('-', expand= True)

    for i in Kestrel_concat.POS:
        if i >= 60:
            Kestrel_concat['POS'] = i-60
            Kestrel_concat.rename(columns={'Motif_right': 'Motif'}, inplace = True)
            Kestrel_concat.drop(['Motif_sequence'], axis=1, inplace=True)
            Kestrel_concat = Kestrel_concat.merge(Merged_motifs, on='Motif')
            Kestrel_concat = Kestrel_concat[['Motif', 'Variant', 'POS', 'REF', 'ALT','Motif_sequence', 'Estimated depth of all haplotypes supporting the alternate variant',
            'Estimated depth of all haplotypes in the variant active region']]
            if Kestrel_concat['Motif'].str.contains('X').any():
                Kestrel_concat = Kestrel_concat.loc[(Kestrel_concat['Motif'] == 'X') & (Kestrel_concat['ALT'] == 'GG')]
                Kestrel_concat.drop_duplicates(subset=['Motif', 'POS', 'REF', 'ALT'], inplace=True)

        else:
            Kestrel_concat['POS'] = i
            Kestrel_concat.rename(columns={'Motif_left': 'Motif'}, inplace = True)
            Kestrel_concat.drop(['Motif_sequence'], axis=1, inplace=True)
            Kestrel_concat = Kestrel_concat.merge(Merged_motifs, on='Motif')
            Kestrel_concat = Kestrel_concat[['Motif', 'Variant', 'POS', 'REF', 'ALT','Motif_sequence', 'Estimated depth of all haplotypes supporting the alternate variant',
                            'Estimated depth of all haplotypes in the variant active region']]           
            if Kestrel_concat['Motif'].str.contains('X').any():
                Kestrel_concat = Kestrel_concat.loc[(Kestrel_concat['Motif'] == 'X') & (Kestrel_concat['ALT'] == 'GG')]
                Kestrel_concat.drop_duplicates(subset=['Motif', 'POS', 'REF', 'ALT'], inplace=True)
            else:
                Kestrel_concat.drop_duplicates(subset=['Motif', 'POS', 'REF', 'ALT'], inplace=True)

# Save Kestrel processed result
Kestrel_concat.to_csv(output + args.output + '_pre_result.tsv', sep='\t', index=False)


stop = timeit.default_timer()
logging.info('Kestrel output processing done...')

# VNTR genotyping with code-adVNTR
columns_kmer = ['Motif', 'Variant', 'POS', 'REF', 'ALT', 'Motif_sequence (RevCom)', 'Estimated depth of all haplotypes supporting the alternate variant', 
'Estimated depth of all haplotypes in the variant active region']

columns_adVNTR = ['#VID', 'Variant','NumberOfSupportingReads', 'MeanCoverage', 'Pvalue']

rm_list_1 = ['_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_indel.vcf']
rm_list_2 = ['_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_dupC_inframe.tsv', '_indel.vcf']
rm_list_3 = ["_R1.fastq.gz", "_R2.fastq.gz"]

if args.ignore_advntr:
    logging.info('MUC1 VNTR genotyping with adVNTR skippied...')

    if dupC_in_frame.empty:

        with open(output + args.output  + '_pre_result.tsv','r') as f:
            with open(output + args.output + '_Final_result.tsv', 'w') as f1:
                f1.write('## VNtyper_Analysis_for_%s \n' % args.output )
                f1.write('# Kestrel_Result \n')
                f1.write('\t'.join(columns_kmer) + '\n') 
                next(f)
                for line in f:
                    f1.write(line)
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
        with open(output + args.output  + '_dupC_inframe.tsv','r') as f2:
            with open(output + args.output + '_Final_result.tsv', 'w') as f3:
                f3.write('## VNtyper_Analysis_for_%s \n' % args.output )
                f3.write('# Kestrel_Result \n')
                f3.write('\t'.join(columns_kmer) + '\n') 
                next(f2)
                for line in f2:
                    f3.write(line)
        for db in rm_list_2:
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
    if args.alignment is not None:
        print('launching adVNTR (This will take a while)...')
        logging.info('launching code-adVNTR (Profile-HMM based VNTR genotyping)...')
        sorted_bam = args.alignment
        reference = args.reference_file
        db_file_hg19 = args.reference_vntr
        if None not in (sorted_bam , reference):
            adVNTR_command = "singularity exec " + tools_path + "code-adVNTR.sif advntr genotype -fs -vid 25561 --outfmt vcf --alignment_file " + sorted_bam + " -o " + output + args.output + "_adVNTR.vcf" + " " + "-m " + db_file_hg19 + " -r " + args.reference_file + " --working_directory " + output
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
                Mapping_command = "bwa mem -t " + args.threads + " " + reference + " " + fastq_1 + " " +  fastq_2 + " -o " + sam_out +  " && "  + "java -Xmx10g -jar " + tools_path + "picard.jar SortSam VALIDATION_STRINGENCY=SILENT " + "I=" + sam_out  + " " + "O=" + output + args.output + "_sorted.bam" + " SORT_ORDER=coordinate" + " && " + "samtools index " + output + args.output + "_sorted.bam"
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
            adVNTR_command = "singularity exec " + tools_path + "code-adVNTR.sif advntr genotype -fs -vid 25561 --outfmt vcf --alignment_file " + sorted_bam + " -o " + output + args.output + "_adVNTR.vcf" + " " + "-m " + db_file_hg19 + " -r " + args.reference_file + " --working_directory " + output
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
    if (df1.Deletion_length > 0).any() and (df1.Insertion_len > 0).any():
        
        if (df1.Deletion_length > df1.Insertion_len).any():
            
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
    if (dff.Deletion_length > 0).any() and (dff.Insertion_len > 0).any():
            
        if (dff.Insertion_len > dff.Deletion_length).any():
            
            dff = dff[dff['frame'].isin(ins_frame)]
            
    return dff


rm_list_4 = ['_dupC_inframe.tsv', '_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_indel.vcf', '.sam', '_sorted.bam', '_sorted.bam.bai', '_R1.fastq.gz', '_R2.fastq.gz', '_adVNTR.vcf', '_adVNTR_result.tsv', '.vcf']
rm_list_5 = ['_pre_result.tsv', '_insertion.vcf', '_deletion.vcf', '_indel.vcf', '_R1.fastq.gz', '_R2.fastq.gz', '_adVNTR.vcf', '_adVNTR_result.tsv', '.vcf']

if df.empty:

    with open(output + args.output  + '_pre_result.tsv','r') as f:
        with open(output  + args.output + '_Final_result.tsv', 'w') as f1:
            f1.write('## VNtyper_Analysis_for_%s \n' % args.output )
            f1.write('# Kestrel_Result \n')
            f1.write('\t'.join(columns_kmer) + '\n') 
            next(f)
            for line in f:
                f1.write(line)
            f1.write('# Code-adVNTR_Result \n')
            f1.write('\t'.join(columns_adVNTR) + '\n')
            f1.write('No variant was found!')

        for db in rm_list_4:
            db_str = args.output + db
            rm_command =  "rm " + output + db_str
            process = sp.Popen(rm_command, shell=True)
            process.wait()

    sys.exit()  
else:

    df_del = advntr_processing_del(df)
    df_ins = advntr_processing_ins(df)

advntr_concat = pd.concat([df_del, df_ins], axis=0)
advntr_concat = advntr_concat[['#VID', 'Variant', 'NumberOfSupportingReads', 'MeanCoverage', 
                              'Pvalue']]
advntr_concat.drop_duplicates(subset=['#VID', 'Variant', 'NumberOfSupportingReads'], inplace=True)

advntr_concat.to_csv(output + args.output + '_adVNTR_result.tsv', sep='\t', index=False)

if dupC_in_frame.empty:

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
    sys.exit() 

else:
    with open(output + args.output + '_dupC_inframe.tsv','r') as f, open(output + args.output +'_adVNTR_result.tsv', 'r') as ff:
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
    sys.exit()

logging.info('The final result is saved in *_Final_result.tsv')
