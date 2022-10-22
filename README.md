## VNtyper - A pipeline to genotype MUC1-VNTR 
Genotyping MUC1 coding-VNTR in autosomal dominant tubulointerstitial kidney disease (MUC1-ADTKD) using short-read seqeuncing (SRS) data. Vntyper pipeline embeded two different varinat calling algorithms:
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
The following command will donwload and install all prerequestics:
```bashscrip
chmod u+x install_prerequisites.sh
bash install_prerequisites.sh or ./install_prerequisites.sh
```
The requeirments are as follows:
1. Download and index chr1.fa file form UCSC genome browser using [(BWA)](https://bio-bwa.sourceforge.net/)
2. Download, build and install Singularity
3. Download Kestrel
4. Building singularity image for code-adVNTR
5. Download [(hg19 VNTR database)](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip) for code-adVNTR 


## Execution



## Output



## Citation


