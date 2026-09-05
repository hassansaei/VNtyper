# Scientific Background

## The MUC1 Gene

The *MUC1* gene, located on chromosome 1q22, encodes mucin-1, a transmembrane glycoprotein expressed on the apical surface of epithelial cells, including renal tubular cells (Kirby et al., *Nat Genet* 2013). A defining feature of *MUC1* is a large coding Variable Number Tandem Repeat (VNTR) in exon 2.

## Variable Number Tandem Repeats

VNTRs are genomic regions where short sequence motifs repeat in tandem, with copy numbers varying across individuals. In *MUC1*, each repeat unit spans 60 bp (encoding 20 amino acids), and alleles carry between 20 and over 125 repeat units (Saei et al., *iScience* 2023; Kirby et al., *Nat Genet* 2013).

## ADTKD-MUC1

Autosomal Dominant Tubulointerstitial Kidney Disease caused by *MUC1* mutations (ADTKD-MUC1) produces progressive tubular damage and renal fibrosis, leading to end-stage kidney disease in adulthood (Bleyer et al., *Kidney Int* 2014). Causative mutations are frameshift insertions or deletions within the VNTR coding sequence that generate a truncated, toxic protein termed MUC1-fs (Dvela-Levitt et al., *Cell* 2019). The most common mutation is a single cytosine duplication (dupC) in the homopolymer tract at position 67 of the repeat unit (Kirby et al., *Nat Genet* 2013).

## The Genotyping Challenge

Standard short-read aligners (BWA, Bowtie2) depend on unique mapping locations. Within the repetitive, GC-rich VNTR, reads map ambiguously across multiple repeat units with zero mapping quality. Conventional variant callers (such as GATK HaplotypeCaller or DeepVariant) routinely miss pathogenic frameshift insertions in this region (Saei et al., *iScience* 2023).

## VNtyper 2 Approach

VNtyper 2 incorporates the Kestrel variant caller (Audano et al., *Bioinformatics* 2018), which reconstructs local haplotypes directly from k-mer frequency spectra without reference alignment. This design avoids mapping bias in repetitive sequence. VNtyper 2 refines raw Kestrel outputs through empirical depth-score thresholds, quality metric checks, and confidence classifications, with optional orthogonal validation via the alignment-based adVNTR caller (Park et al., *iScience* 2022).

!!! info "Key References"
    - Saei H et al. *iScience* 26, 107171 (2023): VNtyper method and validation
    - Popp B, Saei H et al. *medRxiv* (2026): VNtyper 2 and high-speed MUC1 genotyping
    - Kirby A et al. *Nat Genet* 45, 299-303 (2013): ADTKD-MUC1 discovery
    - Audano PA et al. *Bioinformatics* 34, 1659-1665 (2018): Kestrel algorithm
    - Park J et al. *iScience* 25, 104785 (2022): code-adVNTR
    - Dvela-Levitt M et al. *Cell* 179, 1222-1233 (2019): MUC1-fs disease mechanism
