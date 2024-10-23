## VNtyper - A pipeline to genotype MUC1-VNTR 
Genotyping MUC1 coding-VNTR in ADTKD-MUC1 using short-read sequencing (SRS) data. Vntyper pipeline embedded two different variant calling algorithms:
- Mapping free genotyping using kmer frequencies [(Kestrel)](https://github.com/paudano/kestrel)
- Profile-HMM based method [(code-adVNTR)](https://github.com/mehrdadbakhtiari/adVNTR/tree/enhanced_hmm)

## Installation & Requirements
The tool can be downloaded by cloning from the GitHub page (You should have sudo access to the server, otherwise you can use docker/singularity image):

```bashscript
# Make a directory that you want to download VNtyper
mkdir vntyper
git clone https://github.com/hassansaei/VNtyper.git
# Go to the directory in which you downloaded the source code
cd VNtyper
```
The following command will automatically download and install all prerequisites:
```bashscrip
chmod u+x install_prerequisites.sh
bash install_prerequisites.sh or ./install_prerequisites.sh
```
The requirements are as follows:
1. Python >= 3.9 and libraries
    - Pandas ``` pip3 install pandas```
    - numpy ``` pip3 install numpy```
    - regex ``` pip3 install regex```
    - biopython ``` pip3 install biopython ```
    - setuptools==58 ``` pip3 install setuptools```
    - pysam ``` pip3 install pysam ```

3. Install [(BWA)](https://bio-bwa.sourceforge.net/)
4. Download chr1.fa file form [(UCSC genome browser)](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz)
5. Index fasta file with BWA
6. Install Singularity
7. Download [(Kestrel)](https://github.com/paudano/kestrel)
8. Building singularity image for code-adVNTR
9. Download [(VNTR database)](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip) for code-adVNTR
10. The MUC1 VNTR motif dictionary and index files are provided in the File directory

## VNtyper docker image

Docker image is also provided and can be pulled from the docker hub. 
You have to make a directory to store both you inputs and outputs in the host machine.
The instructions for installing docker on Linux can be found [(here)](https://docs.docker.com/desktop/install/linux-install/)

```bashscript
mkdir shared
sudo docker pull saei/vntyper:1.3.0
```

__Run docker with only the kmer method:__

```bashscript
sudo docker run --rm -it -v /PATH to the shared directory/shared:/SOFT/shared saei/vntyper:1.3.0 \
-t 8 --bam  --hg19[or --hg38] -p /SOFT/VNtyper/ -ref  /SOFT/VNtyper/Files/chr1.fa  \
-ref_VNTR /SOFT/VNtyper/Files/MUC1-VNTR.fa \
-a /SOFT/shared/SAAMPLE.bam -t 8 -w /SOFT/shared/ -o SAMPLE_NAME --ignore_advntr
```
__Run docker with both methods:__

```bashscript

sudo docker run --rm -it -v /PATH to the shared directory/shared:/SOFT/shared saei/vntyper:1.3.0 \
-t 8 --bam  --hg19 -p /SOFT/VNtyper/  -ref  /SOFT/VNtyper/Files/chr1.fa  \
-ref_VNTR /SOFT/VNtyper/Files/MUC1-VNTR.fa  -m /SOFT/VNtyper/Files/hg19_genic_VNTRs.db \
-a /SOFT/shared/SAMPLE.bam -t 8 -w /SOFT/shared/ -o SAMPLE_NAME

```

## Execution
Use the following command to see the help for running the tool.

```bashscript
python3 VNtyper.py --help 

usage: VNtyper.py [-h] [--version] -ref Referense [-r1 FASTQ1] [-r2 FASTQ2] -o OUTPUT -ref_VNTR Referense [-t THREADS] -p TOOLS_PATH -w WORKING_DIR [-m REFERENCE_VNTR]
                  [--ignore_advntr] [--bam] [--hg19] [--hg38] [--fastq] [-a ALIGNMENT] [-c CUSTOM_REGION]

Given raw fastq files or BAM files, this pipeline genotype MUC1-VNTR using two independent genotyping methods (Kestrel and code-adVNTR)

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -ref Referense, --reference_file Referense
                        FASTA-formatted reference file and indexes for read alignment
  -r1 FASTQ1, --fastq1 FASTQ1
                        Fastq file first pair
  -r2 FASTQ2, --fastq2 FASTQ2
                        Fastq file second pair
  -o OUTPUT, --output OUTPUT
                        Output file name
  -ref_VNTR Referense, --reference_VNTR Referense
                        MUC1-specific reference file
  -t THREADS, --threads THREADS
                        Number of threads (CPU)
  -p TOOLS_PATH, --tools_path TOOLS_PATH
                        Path to the VNtyper directory
  -w WORKING_DIR, --working_dir WORKING_DIR
                        the path to the output
  -m REFERENCE_VNTR, --reference_vntr REFERENCE_VNTR
                        adVNTR reference vntr database
  --ignore_advntr       Skip adVNTR genotyping of MUC1-VNTR
  --bam                 BAM file as an input
  --hg19                Input bam is from hg19 alignment
  --hg38                Input bam is from hg38 alignment
  --fastq               Paired-end fastq files as an input
  -a ALIGNMENT, --alignment ALIGNMENT
                        Alignment File (with an index file .bai)
  -c CUSTOM_REGION, --custom_region CUSTOM_REGION
                        Custom region for MUC1 gene to be subseted in the provided BAM file

```
[Note] Since the program uses python3.9 logging system, it can not be executed using lower versions of the python.

__Running only kmer-based genotyping:__
```bashscript
python3 VNtyper.py --bam -ref Files/chr1.fa -a SAMPLE.bam -o SAMPLE_NAME -ref_VNTR Files/MUC1-VNTR.fa -t Threads -p VNtyper/ -w WORKING_DIRECTORY --ignore_advntr
```
[Note] This algorithm is far more faster than the second method. 

__Running both genotyping methods:__
```bashscript
python3 VNtyper.py --bam -ref Files/chr1.fa -a SAMPLE.bam -o SAMPLE_NAME -ref_VNTR Files/MUC1-VNTR.fa -t Threads -p VNtyper/  -w WORKING_DIRECTORY -m Files/vntr_data/hg19_genic_VNTRs.db
```
[Note] This algorithm is far more slower than the first method.

## Results from high-coverage 1000G project
We analyzed MUC1 region in 2300 samples from 1000G 30X project. The results from this analysis could be found [(here)](https://e.pcloud.link/publink/show?code=kZxlpjZEVQGRPuTSgRyYq3xOrs9hkWENpRX) 

## Evaluating MUC1 VNTR region coverage using samtools
```bashscript
for f in *.bam; do  samtools depth -b MUC_hg19.bed $f | awk '{sum+=$3} END { print sum/NR}' > $f.coverage; done
```
MUC_hg19.bed is provide. MUC1_hg19.bed could also be replaced by : -r chr1:155160500-155162000

## Sample bam files MUC1 8C positive 
Here we provided five (example_1.bam to example_5.bam) MUC1 8C positive bam files for evaluation.
Link to bam files: [(Bam)](https://e.pcloud.link/publink/show?code=kZGSejZWTuXKX6IQnzyD5yxpUJMNpiONMXk)
Example_1 to 3 from NTI cohort and example_4 and 5 from renome cohort.

## Output
The tool creates a folder for each case in the working directory which is assigned by the user. Inside the folder there is a directory for temporary files and log files, and the final output:
- Temp folder: Fastp QC report (.html) and log file for VNtyper
- The output of VNtyper '*_Final_result.tsv'


The Kestrel output is a VCF file, which is processed by VNtyper and the final result is stored in *_Final_result.tsv. The resulting file contains information for the motifs, varinant types, position of the variant, and its corresponding depth. The output for code-adVNTR is a bed or vcf file with variant information and P-value.



## Citations
1. Saei H, Morinière V, Heidet L, Gribouval O, Lebbah S, Tores F, Mautret-Godefroy M, Knebelmann B, Burtey S, Vuiblet V, Antignac C, Nitschké P, Dorval G. VNtyper enables accurate alignment-free genotyping of MUC1 coding VNTR using short-read sequencing data in autosomal dominant tubulointerstitial kidney disease. iScience. 2023 Jun 17;26(7):107171. https://doi.org/10.1016/j.isci.2023.107171
   
2. Peter A Audano, Shashidhar Ravishankar, Fredrik O Vannberg, Mapping-free variant calling using haplotype reconstruction from k-mer frequencies, Bioinformatics, Volume 34, Issue 10, May 2018, Pages 1659–1665, https://doi.org/10.1093/bioinformatics/btx753
   
3. Park J, Bakhtiari M, Popp B, Wiesener M, Bafna V. Detecting tandem repeat variants in coding regions using code-adVNTR. iScience. 2022 Jul 19;25(8):104785. https://doi.org/10.1016/j.isci.2022.104785 

