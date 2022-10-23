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
1. Install [(BWA)](https://bio-bwa.sourceforge.net/)
2. Download chr1.fa file form [(UCSC genome browser)](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz)
3. Index fasta file with BWA
4. Install Singularity
5. Download [(Kestrel)](https://github.com/paudano/kestrel)
6. Building singularity image for code-adVNTR
7. Download [(VNTR database)](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip) for code-adVNTR
8. The MUC1 VNTR motif dictionary and index files are provided in the File directory


## Execution



## Output



## Citation


