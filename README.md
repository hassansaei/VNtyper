## VNtyper - A pipeline to genotype MUC1-VNTR 
Genotyping MUC1 coding-VNTR in medullary-cystic kidney disease type 1  (MCKD1) using short-read seqeuncing (SRS) data. Vntyper pipeline embeded two different varinat calling algorithms:
- Mapping free genotyping using kmer frequencies [(Kestrel)](https://github.com/paudano/kestrel)
- Profile-HMM based method [(code-adVNTR)](https://github.com/mehrdadbakhtiari/adVNTR/tree/enhanced_hmm)

## Installation & Requirements
The tool can be downloaded by cloning from the github page:

```bashscript
# Make a directory that you want to download VNtyper
mkdir vntyper
git clone https://github.com/hassansaei/VNtyper.git
# Go to the directory that you downloaded the source code
cd VNtyper
```
The following command will automatically download and install all prerequisites:
```bashscrip
chmod u+x install_prerequisites.sh
bash install_prerequisites.sh or ./install_prerequisites.sh
```
The requeirments are as follows:
1. Python >= 3.9
2. Install [(BWA)](https://bio-bwa.sourceforge.net/)
3. Download chr1.fa file form [(UCSC genome browser)](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz)
4. Index fasta file with BWA
5. Install Singularity
6. Download [(Kestrel)](https://github.com/paudano/kestrel)
7. Building singularity image for code-adVNTR
8. Download [(VNTR database)](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip) for code-adVNTR
9. The MUC1 VNTR motif dictionary and index files are provided in the File directory


## Execution
Use following command to see the help for running the tool.
```bashscript
python3 VNtyper.py --help 

usage: VNtyper.py [-h] -ref Referense -r1 FASTQ1 -r2 FASTQ2 -k KMER -o OUTPUT
                  -ref_VNTR Referense [-t THREADS] -p TOOLS_PATH -w
                  WORKING_DIR -m REFERENCE_VNTR [--ignore_advntr]

Given raw fastq files, this tool call frameshift variants in MUC1-VNTR using Kestrel and code-adVNTR algorithms

options:
  -h, --help            show this help message and exit
  -ref Referense, --reference_file Referense
                        FASTA-formatted reference file and indexes
  -r1 FASTQ1, --fastq1 FASTQ1
                        Fastq file first pair
  -r2 FASTQ2, --fastq2 FASTQ2
                        Fastq file second pair
  -k KMER, --Kmer KMER  Kmer size
  -o OUTPUT, --output OUTPUT
                        Output file name
  -ref_VNTR Referense, --reference_VNTR Referense
                        MUC1 VNTR reference file
  -t THREADS, --threads THREADS
                        Number of threads (CPU)
  -p TOOLS_PATH, --tools_path TOOLS_PATH
                        Path to the scripts directory
  -w WORKING_DIR, --working_dir WORKING_DIR
                        the path to the output
  -m REFERENCE_VNTR, --reference_vntr REFERENCE_VNTR
                        adVNTR VNTR database
  --ignore_advntr       Skip code-adVNTR genotyping of MUC1-VNTR
```
[Note] Since the program uses python3.9 logging system, it can not be executed using lower versions of the python.

__Running only kmer-based genotyping:__
```bashscript
python3 VNtyper.py -ref Files/chr1.fa -r1 FASTQ1  -r2 FASTQ2 -k 20 -o SAMPLE_NAME -ref_VNTR Files/MUC1-VNTR.fa -t Threads -p Scripts/ -w WORKING_DIRECTORY -m Files/vntr_data/hg19_genic_VNTRs.db --ignore_advntr
```
[Note] This algorithm is far more faster than the second method. 

__Running both genotyping methods:__
```bashscript
python3 VNtyper.py -ref Files/chr1.fa -r1 FASTQ1  -r2 FASTQ2 -k 20 -o SAMPLE_NAME -ref_VNTR Files/MUC1-VNTR.fa -t Threads -p Scripts/ -w WORKING_DIRECTORY -m Files/vntr_data/hg19_genic_VNTRs.db
```
[Note] To decrease the chance of missing any postive cases it is prefered to run both methods.


## Output



## Citation


