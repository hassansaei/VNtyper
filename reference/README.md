# Reference folder for vntyper

## 1. VNTR Reference Generator

This repository contains a Python script `generate_vntr_reference.py` which generates pairwise combinations of VNTR motifs from a given input FASTA file and writes them out to a new FASTA file. The script also uses a JSON configuration file to filter out certain disallowed motif combinations.

## Overview

The script:
1. Reads an input FASTA file that contains multiple contigs (motifs).
2. Generates all pairwise combinations of these contigs, including self-combinations.
3. Concatenates the sequences in a specified order to form the new reference entries.
4. Uses a JSON configuration file to skip (filter out) certain disallowed combinations.

## Files

- **MUC1_motifs_Rev_com.fa**: The input FASTA file containing the original VNTR motifs.  
  Each contig is represented as:
  ```
  >contigName
  ACTGACTG...
  ```

- **filter_config.json**: The JSON configuration file that specifies disallowed combinations.  
  Example format:
  ```json
  {
    "1": ["1", "3", "4", ...],
    "2": ["1", "2", "4", ...],
    ...
  }
  ```
  This means if the first contig is `"1"`, it should not be paired with `"1"`, `"3"`, `"4"`, etc.

- **generate_vntr_reference.py**: The Python script that:
  - Loads the configuration JSON file.
  - Iterates over all contigs from the input FASTA file.
  - Generates all pairwise combinations.
  - Skips any combinations found in the JSON configuration file.
  - Writes the allowed combinations to the output FASTA file.

## Prerequisites

- Python 3.x
- A JSON configuration file with the appropriate format.
- A FASTA file with your source motifs.

## Usage

1. **Prepare Your Input Files**:
   - Ensure `MUC1_motifs_Rev_com.fa` (or your chosen input FASTA) is in the same directory as the script.
   - Place `filter_config.json` in the same directory and confirm it has the correct format.

2. **Run the Script**:
   From the terminal or command prompt:
   ```bash
   python3 generate_vntr_reference.py
   ```
   
   By default, the script is set to:
   - Read from `MUC1_motifs_Rev_com.fa`.
   - Write to `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`.
   - Use `filter_config.json` for filtering.

   If you want to customize paths, open `generate_vntr_reference.py` and modify:
   ```python
   input_file_path = 'MUC1_motifs_Rev_com.fa'
   output_file_path = 'All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa'
   config_file_path = 'filter_config.json'
   ```
   Then run the script again.

3. **Check the Output**:
   After running the script, you will have a new FASTA file (`All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`) containing only the allowed motif combinations.

## Configuration Details

- The `filter_config.json` file maps each "first contig" to a list of "second contigs" that cannot follow it.
- If a combination `first-second` is disallowed, it will be skipped entirely.
- Modify `filter_config.json` as needed to add or remove disallowed combinations. Any combination where the first contig is not listed or the second contig is not in the corresponding array/set will be allowed.

## Example

If you have:
```json
{
  "1": ["1", "3"]
}
```

This means:
- When the first contig is `1`, do not allow `1-1` or `1-3`.
- All other combinations starting with `1` (except `1-1` and `1-3`) are allowed.

## 2. BAM/FASTQ Pseudonymization and Anonymization

### Overview

`pseudonymize.py` is a comprehensive script for properly anonymizing and pseudonymizing BAM and FASTQ files for research use, particularly optimized for creating test datasets from patient data.

**Key Features:**
- **Efficient workflow**: Subsets BAM files to region of interest FIRST, then anonymizes only the small subset
- **Complete anonymization**: Removes ALL identifying information from headers (@RG, @PG, @CO) and read group tags
- **Deterministic pseudonyms**: Generates consistent pseudonyms based on file MD5 hashes
- **Region subsetting**: Extracts MUC1 VNTR region from whole genome/exome BAMs
- **FASTQ generation**: Optionally reverts subsetted BAMs back to paired-end FASTQ files
- **Verification**: Automatically verifies no forbidden strings remain in output

### Workflow

The script implements best practices for efficient BAM anonymization:

1. **Compute MD5 checksums** - Calculate checksums for deterministic pseudonym generation
2. **Subset to MUC1 region** - Extract only chr1:155158000-155163000 (hg19) or equivalent
3. **Anonymize subset** - Anonymize only the small subset (~100x faster than full BAM)
4. **Revert to FASTQ** - Optionally convert anonymized subset to paired-end FASTQ files
5. **Generate manifest** - Create JSON with MD5 checksums for all output files

**Efficiency:** Processing a 2GB BAM takes <30 seconds (subset + anonymize) vs 10+ minutes for full BAM anonymization.

### Usage

#### Basic Usage
```bash
python pseudonymize.py \
  --input-dir /path/to/original/bams \
  --output-dir /path/to/anonymized/output \
  --ref-assembly hg19 \
  --subset-muc1 \
  --revert-fastq
```

#### Key Arguments

**Required:**
- `--input-dir` - Directory containing original BAM files
- `--output-dir` - Output directory for anonymized files

**Subsetting & Conversion:**
- `--subset-muc1` - Create MUC1 region-only subsets (required for output)
- `--revert-fastq` - Convert subset BAMs to paired-end FASTQ files
- `--ref-assembly {hg19,hg38,GRCh37,GRCh38}` - Reference assembly for region coordinates

**Anonymization:**
- `--forbidden-strings FILE` - File with identifiers to remove (one per line, auto-detected if not provided)
- `--anonymize-read-names` - Also anonymize read names (slow, adds ~5-10min per file)

**Performance:**
- `--workers N` - Number of parallel workers (default: all available cores)

**Output Control:**
- `--json-out FILE` - JSON manifest file path (default: pseudonymization_output.json)
- `--json-filter {all,subset}` - Include all files or only subset files in JSON (default: all)
- `--log-level {DEBUG,INFO,WARNING,ERROR}` - Logging verbosity (default: INFO)

### Examples

**1. Full workflow - subset, anonymize, and create FASTQs:**
```bash
python pseudonymize.py \
  --input-dir raw_bams/ \
  --output-dir tests/data_anonymized/ \
  --ref-assembly hg19 \
  --subset-muc1 \
  --revert-fastq \
  --json-filter subset \
  --workers 4
```

**2. With custom forbidden strings:**
```bash
# Create forbidden_strings.txt with sample IDs (one per line)
echo -e "SAMPLE123\nPATIENT456" > forbidden_strings.txt

python pseudonymize.py \
  --input-dir raw_bams/ \
  --output-dir anonymized/ \
  --ref-assembly hg38 \
  --subset-muc1 \
  --forbidden-strings forbidden_strings.txt
```

**3. Multiple reference assemblies (with mapping file):**
```bash
# Create ref_mapping.csv
# filename,reference
# sample1.bam,hg19
# sample2.bam,hg38

python pseudonymize.py \
  --input-dir raw_bams/ \
  --output-dir anonymized/ \
  --ref-mapping-file ref_mapping.csv \
  --subset-muc1 \
  --revert-fastq
```

### Output Files

For each input BAM, the script generates:

**Subset BAM files:**
- `example_XXXX_hg19_subset.bam` - Anonymized subset BAM (MUC1 region only)
- `example_XXXX_hg19_subset.bam.bai` - BAM index

**FASTQ files (if --revert-fastq specified):**
- `example_XXXX_hg19_subset_R1.fastq.gz` - Forward reads
- `example_XXXX_hg19_subset_R2.fastq.gz` - Reverse reads

**Metadata files:**
- `pseudonymization_table.csv` - Mapping of original names to pseudonyms
- `pseudonymization_output.json` - Resource manifest with MD5 checksums

### Anonymization Details

The script performs comprehensive anonymization:

1. **Header cleaning:**
   - Removes all @RG (read group) lines with identifying information
   - Removes all @PG (program) lines containing file paths and commands
   - Removes all @CO (comment) lines that may contain metadata
   - Adds minimal generic @RG header for tool compatibility

2. **Read modifications:**
   - Replaces ALL RG:Z tags in reads with generic value using `samtools addreplacerg --mode overwrite_all`
   - Optionally anonymizes read names (QNAMEs) by hashing flowcell/instrument IDs

3. **Verification:**
   - Checks headers and sample of reads for forbidden strings
   - Raises error if any identifying information remains

### Supported Reference Assemblies

| Assembly | MUC1 Region Coordinates |
|----------|------------------------|
| hg19 | chr1:155158000-155163000 |
| hg38 | chr1:155184000-155194000 |
| GRCh37 | 1:155158000-155163000 |
| GRCh38 | 1:155184000-155194000 |

### Performance Benchmarks

Tested on 7 BAM files (~14 GB total):
- **Processing time**: ~3 minutes total
- **Output size**: ~50 MB (280x reduction)
- **Per-file time**: 10-40 seconds (vs 5-15 minutes for full BAM anonymization)

### Dependencies

- **samtools** (v1.15+) - BAM manipulation
- Python 3.7+
- Standard library: `os`, `re`, `csv`, `json`, `hashlib`, `subprocess`, `argparse`, `logging`, `tempfile`, `concurrent.futures`

### Important Notes

- **Research use only** - Not validated for clinical use
- **Subsetting required** - Use `--subset-muc1` to generate output files
- **Irreversible** - Pseudonymization is one-way; keep original mapping secure
- **Verification** - Always check output files do not contain identifying information
- **Paired reads** - Subsetting uses `samtools view -P` to preserve mate pairs

## 3. Adapted adVNTR references

The referneces in vntr_db_advntr.zip have been created by removing all non-Muc1 entries from the hg9 database and adding the MUC1 reference with the same ID to the hgg38 database where it was missing and changing the start position to the corresponding start position after finding this with UCSC BLAT usng the left sequence.