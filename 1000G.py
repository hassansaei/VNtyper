#python 3.9

import boto3

bucket = "1000genomes"
folder = "1000G_2504_high_coverage/data"
s3 = boto3.resource("s3")
s3_bucket = s3.Bucket(bucket)
files_in_s3 = [f.key.split(folder + "/")[1] for f in s3_bucket.objects.filter(Prefix=folder).all()]
list = [s.partition('/')[0] for s in files_in_s3]
list_2 = [s.partition('/')[2] for s in files_in_s3]
list_2 = [i for i in list_2 if '.crai' not in i]
list_2 = [s.replace(".final.cram", "") for s in list_2]


import pandas as pd

#list = pd.DataFrame(list)
#list = list.drop_duplicates()
#list.to_csv('1000G_2504_high_coverage.csv', index=False) 

list_final = [ii for n ,ii in enumerate(list) if ii not in list[:n]]
list_final_2 = [ii for n ,ii in enumerate(list_2 ) if ii not in list_2[:n]]
length = len(list_final)
length_2 = len(list_final_2)

from concurrent.futures import thread
import subprocess as sp
import os
import os.path
from pathlib import Path


for i in range(length):
	print(list_final[i])
	if os.path.isfile("/WORKSPACE/1000G/crams/" + list_final_2[i] + "_chr1.bam"):
		print ("File %s exist!" %i)
		#extract_chr1 = "/WORKSPACE/1000G/samtools-1.3.1/./samtools view -@ 12 -T /WORKSPACE/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa -b -o " + " /WORKSPACE/1000G/crams/" + list_final_2[i] + "_chr1.bam" + " " + " /WORKSPACE/1000G/crams/" + list_final_2[i] + ".final.cram" + " chr1:155184000-155194000 " + " && " + " /WORKSPACE/1000G/samtools-1.3.1/./samtools index " + " /WORKSPACE/1000G/crams/" +  list_final_2[i] + "_chr1.bam" 
		#process = sp.Popen(extract_chr1 , shell=True)
		#process.wait()
		#print ("Cram to bam conversion done!\n")
		run_vntyper = "python /WORKSPACE/1000G/crams/VNtyper_1.1.py  -ref /WORKSPACE/hg19/chr1.fa --bam -a " + "/WORKSPACE/1000G/crams/" +  list_final_2[i] + "_chr1.bam " + " -ref_VNTR /SOFT/VNtyper/Files/MUC1-VNTR.fa -t 12 -p /SOFT/VNtyper/ -w " +  "/WORKSPACE/1000G/vntyper/" + " -o " +  list_final_2[i] + " --ignore_advntr" 
		process = sp.Popen(run_vntyper , shell=True)
		process.wait()
		print ("VNtyper analysis done!\n")
	else:
		download_data = "/WORKSPACE/bin/aws s3 cp s3://1000genomes/1000G_2504_high_coverage/data/" + list_final[i] + " " + " /WORKSPACE/1000G/crams/" +  " --recursive "
		process = sp.Popen(download_data , shell=True)
		process.wait()
		print ("AWS cram file donwloaded! number: %s" %i)
		extract_chr1 = "/WORKSPACE/1000G/samtools-1.3.1/./samtools view -@ 12 -T /WORKSPACE/1000G/GRCh38_full_analysis_set_plus_decoy_hla.fa -b -o " + " /WORKSPACE/1000G/crams/" + list_final_2[i] + "_chr1.bam" + " " + " /WORKSPACE/1000G/crams/" + list_final_2[i] + ".final.cram" + " chr1:155184000-155194000 " + " && " + " /WORKSPACE/1000G/samtools-1.3.1/./samtools index " + " /WORKSPACE/1000G/crams/" +  list_final_2[i] + "_chr1.bam" 
		process = sp.Popen(extract_chr1 , shell=True)
		process.wait()
		print ("Cram to bam conversion done!\n")
		run_vntyper = "python /WORKSPACE/1000G/crams/VNtyper_1.1.py  -ref /WORKSPACE/hg19/chr1.fa --bam -a " + "/WORKSPACE/1000G/crams/" +  list_final_2[i] + "_chr1.bam " + " -ref_VNTR /SOFT/VNtyper/Files/MUC1-VNTR.fa -t 12 -p /SOFT/VNtyper/ -w " +  "/WORKSPACE/1000G/vntyper/" + " -o " +  list_final_2[i] + " --ignore_advntr" + " && " + " rm " + " /WORKSPACE/1000G/crams/" + list_final_2[i] + ".final.cram"
		process = sp.Popen(run_vntyper , shell=True)
		process.wait()
		print ("VNtyper analysis done!\n")



