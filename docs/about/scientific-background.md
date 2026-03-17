# Scientific Background

## The MUC1 Gene

The *MUC1* gene, located on chromosome 1q22, encodes mucin-1, a transmembrane glycoprotein expressed on the apical surface of epithelial cells, including kidney tubular cells (Kirby et al., *Nat Genet* 2013). A defining feature of *MUC1* is a large coding Variable Number Tandem Repeat (VNTR) in exon 2.

## Variable Number Tandem Repeats

VNTRs are stretches of DNA where a short sequence motif is repeated in tandem. The number of copies varies between individuals. In *MUC1*, each repeat unit is approximately 60 bp (encoding 20 amino acids), and individuals carry between 20 and over 125 copies per allele (Saei et al., *iScience* 2023; Kirby et al., *Nat Genet* 2013).

## ADTKD-MUC1

Autosomal Dominant Tubulointerstitial Kidney Disease caused by *MUC1* mutations (ADTKD-MUC1) is characterized by progressive tubular damage and renal fibrosis, typically leading to end-stage kidney disease in adulthood (Bleyer et al., *Kidney Int* 2014). The causative mutations are frameshift insertions or deletions within the VNTR coding sequence that produce a truncated, misfolded protein termed MUC1-fs (Dvela-Levitt et al., *Cell* 2019). The most common mutation is a cytosine duplication (dupC) at position 67 of the repeat unit (Kirby et al., *Nat Genet* 2013).

## The Genotyping Challenge

Standard short-read aligners (BWA, Bowtie2) rely on unique mapping positions. Within the highly repetitive VNTR, reads map ambiguously to multiple repeat units, resulting in poor alignment quality and unreliable variant calls. This causes pathogenic frameshift mutations to be routinely missed by conventional variant-calling pipelines (Saei et al., *iScience* 2023).

## VNtyper 2's Approach

VNtyper 2 uses the Kestrel variant caller (Audano et al., *Bioinformatics* 2018), which reconstructs local haplotypes directly from k-mer frequency spectra instead of aligning reads to a reference. This bypasses the alignment bias inherent in repetitive regions. VNtyper 2 further refines raw Kestrel calls with empirical scoring thresholds, depth-based confidence classification, and optional cross-validation via the alignment-based adVNTR tool (Park et al., *iScience* 2022).

!!! info "Key References"
    - Saei H et al. *iScience* 26, 107171 (2023) — VNtyper method and validation
    - Kirby A et al. *Nat Genet* 45, 299–303 (2013) — ADTKD-MUC1 discovery
    - Audano PA et al. *Bioinformatics* 34, 1659–1665 (2018) — Kestrel algorithm
    - Park J et al. *iScience* 25, 104785 (2022) — code-adVNTR
    - Dvela-Levitt M et al. *Cell* 179, 1222–1233 (2019) — MUC1-fs disease mechanism
