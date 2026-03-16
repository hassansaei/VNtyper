# Scientific Background

## The MUC1 Gene

The *MUC1* gene, located on chromosome 1 (1q22), encodes the mucin-1 glycoprotein. Mucin-1 is a transmembrane protein expressed on the apical surface of epithelial cells throughout the body, including the kidneys. It plays roles in cell signaling, adhesion, and protection of epithelial surfaces. A defining feature of *MUC1* is a large coding Variable Number Tandem Repeat (VNTR) region within exon 2.

## Variable Number Tandem Repeats

Variable Number Tandem Repeats (VNTRs) are stretches of DNA where a short motif -- typically tens to hundreds of base pairs -- is repeated consecutively. The number of repeats varies between individuals, making VNTRs a source of significant genomic diversity. In the case of *MUC1*, the VNTR consists of a ~60 bp repeat unit, and an individual may carry anywhere from roughly 20 to over 125 copies per allele.

## ADTKD-MUC1

Autosomal Dominant Tubulointerstitial Kidney Disease caused by *MUC1* mutations (ADTKD-MUC1) is a hereditary kidney disorder characterized by progressive tubular damage and renal fibrosis, typically leading to end-stage kidney disease in adulthood. The causative mutations are frameshift insertions or deletions within the *MUC1* VNTR that produce a truncated, misfolded protein (termed MUC1-fs). Because the mutation can occur in any one of the many repeat copies, conventional Sanger sequencing and standard next-generation sequencing analysis pipelines have difficulty detecting it.

## The Challenge

Standard short-read aligners such as BWA and Bowtie2 rely on unique mapping positions to place reads in the reference genome. Within highly repetitive VNTR regions, reads map ambiguously to multiple locations, resulting in poor alignment quality and unreliable variant calls. This means that pathogenic frameshift mutations in *MUC1* are routinely missed by genome-wide variant-calling workflows.

## VNtyper's Approach

VNtyper overcomes these limitations by employing a mapping-free, k-mer-based strategy through the Kestrel variant caller. Instead of aligning reads to a reference, Kestrel reconstructs local haplotypes directly from k-mer frequency spectra. This sidesteps the alignment bias inherent in repetitive regions and enables sensitive detection of small insertions and deletions that cause frameshifts. VNtyper further refines the raw calls with empirical scoring, depth-based confidence classification, and optional cross-validation using the alignment-based adVNTR tool.

!!! info "Reference"
    Saei H, Moriniere V, Heidet L, et al. VNtyper enables accurate alignment-free genotyping of MUC1 coding VNTR using short-read sequencing data. *iScience*. 2023;26(7):107171. doi:[10.1016/j.isci.2023.107171](https://doi.org/10.1016/j.isci.2023.107171)
