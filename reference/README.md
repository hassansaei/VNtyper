# Reference folder for vntyper

> **Milestone 5 note:** the files this section describes -
> `MUC1_motifs_Rev_com.fa`, `filter_config.json`, `generate_vntr_reference.py` and
> `vntr_db_advntr.zip` - no longer live in this directory. They moved to
> [`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data)'s `seeds/`
> directory, alongside the release spec that pins their digests. `reference/` (this
> directory) now populates at runtime from that repository's published bundle via
> `vntyper install-references`; only this README and the two pseudonymization scripts
> below stay tracked here. The rest of this section documents the generator's logic and
> history and remains accurate as documentation - only the file locations changed.

## 1. Reference genomes, MUC1 motifs and adVNTR databases

The pairwise MUC1 motif generation (`generate_vntr_reference.py`, `filter_config.json`),
the curation history of the adVNTR databases, and the source-of-truth build tooling for
every published reference bundle now live in
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data) rather than in this
repository. `vntyper install-references` downloads the versioned, checksummed release that
repository publishes by default; see
[Reference Setup](https://hassansaei.github.io/VNtyper/getting-started/reference-setup/)
(or `docs/getting-started/reference-setup.md` in this repository) for the trust model, and
`vntyper/scripts/install_references_config.json` for the exact digests VNtyper pins.

`--from-source` still rebuilds those same references locally from their original upstream
hosts, and it is what the `berntpopp/vntyper-data` build workflow itself runs -- see that
flag in [install-references](https://hassansaei.github.io/VNtyper/cli/install-references/)
or `docs/cli/install-references.md` here.

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

### Multi-Reference Remapping

`pseudonymize.py` supports remapping anonymized FASTQs to multiple reference genomes, which is useful for validating reference consistency and testing pipelines across different assemblies.

#### Configuration File

Create a JSON configuration file (`pseudonymize_config.json`) to define reference paths and alignment settings:

```json
{
  "references": {
    "hg19": {
      "path": "reference/alignment/chr1.hg19.fa",
      "region": "chr1:155158000-155163000"
    },
    "hg38": {
      "path": "reference/alignment/chr1.hg38.fa",
      "region": "chr1:155184000-155194000"
    },
    "GRCh37": {
      "path": "reference/alignment/chr1.GRCh37.fna",
      "region": "1:155158000-155163000"
    },
    "GRCh38": {
      "path": "reference/alignment/chr1.GRCh38.fna",
      "region": "1:155184000-155194000"
    }
  },
  "bwa": {
    "threads": 4,
    "required_index_files": [".amb", ".ann", ".bwt", ".pac", ".sa"]
  },
  "forbidden_strings_file": "forbidden_strings.txt"
}
```

#### Usage with Multi-Reference Remapping

```bash
# Remap to all 4 reference genomes
python pseudonymize.py \
  --input-dir raw_bams/ \
  --output-dir tests/data_anonymized/ \
  --ref-assembly hg19 \
  --subset-muc1 \
  --revert-fastq \
  --remap-to-reference hg19,hg38,GRCh37,GRCh38 \
  --forbidden-strings forbidden_strings.txt \
  --log-file tests/data_anonymized/pseudonymize.log \
  --workers 4
```

#### Output Directory Structure

When using multi-reference remapping, files are organized by type and reference:

```
output_dir/
├── pseudonymize.log              # Comprehensive logging
├── console.log                   # Console output
├── pseudonymization_table.csv    # Original → pseudonym mapping
├── pseudonymization_output.json  # Manifest with MD5 checksums
│
├── example_XXXX_hg19_subset.bam        # Subset BAMs (at root level)
├── example_XXXX_hg19_subset.bam.bai
│
├── fastqs/                              # All FASTQs in subdirectory
│   ├── example_XXXX_hg19_subset_R1.fastq.gz
│   ├── example_XXXX_hg19_subset_R2.fastq.gz
│   └── ...
│
└── remapped/                            # Remapped BAMs organized by mapper/reference
    └── bwa/                             # Mapper subdirectory (extensible to minimap2, bowtie2, etc.)
        ├── hg19/
        │   ├── example_XXXX_hg19_bwa.bam
        │   ├── example_XXXX_hg19_bwa.bam.bai
        │   └── ...
        ├── hg38/
        │   ├── example_XXXX_hg38_bwa.bam
        │   ├── example_XXXX_hg38_bwa.bam.bai
        │   └── ...
        ├── GRCh37/
        │   ├── example_XXXX_GRCh37_bwa.bam
        │   ├── example_XXXX_GRCh37_bwa.bam.bai
        │   └── ...
        └── GRCh38/
            ├── example_XXXX_GRCh38_bwa.bam
            ├── example_XXXX_GRCh38_bwa.bam.bai
            └── ...
```

**Directory organization rationale:**
- **Root level**: Subset BAMs (primary outputs from subsetting workflow)
- **fastqs/**: All FASTQ files in one subdirectory for easy access
- **remapped/\<mapper\>/\<reference\>/**: Organized by mapper AND reference for extensibility
  - Supports future addition of other mappers (minimap2, bowtie2, etc.)
  - Keeps outputs from different references cleanly separated
  - BAM naming: `{sample}_{reference}_{mapper}.bam` (e.g., `example_6449_hg38_bwa.bam`)

#### Multi-Aligner Support

`pseudonymize.py` supports multiple alignment tools for remapping. Use `--aligners` to specify which aligners to use:

**Supported aligners:**
- **BWA** (default) - Standard, well-tested aligner
- **BWA-MEM2** - 1.3-3x faster than BWA, drop-in replacement
- **Minimap2** - Versatile, fast for both short and long reads
- **Bowtie2** - Memory-efficient, very accurate
- **DRAGMAP** - DRAGEN-compatible (currently disabled due to stdout piping issues)

**Usage with specific aligners:**
```bash
# Use only BWA (default)
python pseudonymize.py --input-dir raw/ --output-dir out/ \
  --subset-muc1 --revert-fastq --remap-to-reference hg38

# Use multiple aligners
python pseudonymize.py --input-dir raw/ --output-dir out/ \
  --subset-muc1 --revert-fastq --remap-to-reference hg38 \
  --aligners bwa bwa-mem2 minimap2
```

**Zero-length read filtering:**
- Automatically filters zero-length reads for minimap2 (prevents crashes)
- Uses `seqtk` if available (fast), falls back to Python implementation
- No filtering needed for BWA/BWA-MEM2/Bowtie2 (handle gracefully)

**Aligner configuration:**
Edit `reference/pseudonymize_config.json` to enable/disable aligners or modify settings:
```json
{
  "aligners": {
    "bwa": {
      "enabled": true,
      "executable": "bwa",
      "alignment_command": "bwa mem -t {threads} {ref_path} {r1} {r2}",
      ...
    }
  }
}
```

#### Remapping Features

**Automatic index verification:**
- Checks for required index files for each aligner
- Warns if indices are missing
- Use `install_references.py` to generate indices (see below)

**Efficient parallel alignment:**
- Uses multi-threading for all supported aligners
- Pipes directly to samtools for sorting and compression
- Generates BAI indices automatically

**VNtyper-optimized extraction:**
- Uses VNtyper's efficient offset-based unmapped read extraction
- Preserves both mapped (MUC1 region) and unmapped reads
- Maintains read pairing throughout workflow

**Logging and verification:**
- Comprehensive logging to both file and console
- Reports read counts after each remapping step
- Includes alignment statistics in logs

#### Reference Index Installation with install_references.py

`vntyper install-references -d reference/` (with no flags) fetches the published,
BWA-indexed bundle from `berntpopp/vntyper-data` -- see
[Reference Setup](https://hassansaei.github.io/VNtyper/getting-started/reference-setup/).
That default path only ever produces a BWA index, so remapping with `bwa-mem2`, `minimap2`,
`bowtie2` or DRAGMAP needs `--from-source`, which builds every reference locally instead:

**Basic usage:**
```bash
# Fetch the published bundle (BWA index only)
python vntyper/scripts/install_references.py -d reference/

# Build from source with specific aligners only
python vntyper/scripts/install_references.py -d reference/ --from-source --aligners bwa bwa-mem2 minimap2

# Build from source with more threads (faster)
python vntyper/scripts/install_references.py -d reference/ --from-source --threads 16
```

**Configuration:**
Edit `vntyper/scripts/install_references_config.json` to:
- Enable/disable specific aligners for indexing (`--from-source` only)
- Modify NCBI reference URLs and checksums (`--from-source` only; the default bundle path
  verifies against the `asset_sha256` digests in the same file's `bundle`/genome sections
  instead)
- Adjust index parameters

**Supported index types:**
- **in_place**: BWA, BWA-MEM2 (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa` files)
- **separate_file**: Minimap2 (`.mmi` file)
- **index_base**: Bowtie2 (`.bt2` files with base prefix)
- **index_directory**: DRAGMAP (directory with hash tables)

**Performance:**
- BWA/BWA-MEM2: ~5 minutes per reference
- Minimap2: ~2 minutes per reference
- Bowtie2: ~5 minutes per reference
- DRAGMAP: ~10 minutes per reference (disabled by default)

**Manual indexing (alternative):**
```bash
# BWA
bwa index reference/alignment/chr1.hg19.fa

# BWA-MEM2
bwa-mem2 index reference/alignment/chr1.hg38.fa

# Minimap2
minimap2 -d reference/alignment/chr1.GRCh37.fna.mmi reference/alignment/chr1.GRCh37.fna

# Bowtie2
bowtie2-build reference/alignment/chr1.GRCh38.fna reference/alignment/chr1.GRCh38_bowtie2

# Verify indices
ls -lh reference/alignment/
```

#### Forbidden Strings File

Create a text file listing identifiers to remove during anonymization (one per line):

```
# forbidden_strings.txt - Sample identifiers to anonymize
# Comments starting with # are ignored

NPH1149
NPH1908593
NTI1179
SAMPLE_ID_123
PATIENT_456
```

The script automatically scans BAM headers and read groups to extract forbidden strings if not provided.

#### Performance with Multi-Reference Remapping

**Example workload**: 7 BAM files remapped to 4 references each (28 total alignments)

- **Per-sample subset + anonymize**: 10-40 seconds
- **Per-reference alignment**: 20-30 seconds (BWA MEM with 4 threads)
- **Total time**: ~2-3 hours for 7 samples × 4 references
- **Output size**: ~50MB subset BAMs + ~200MB remapped BAMs

**Optimization tips:**
- Use `--workers N` to parallelize sample processing
- BWA threading (`bwa.threads` in config) speeds up individual alignments
- SSD storage significantly improves I/O performance
- Consider increasing `bwa.threads` to 8-16 for faster alignments

## 3. Adapted adVNTR references

The curation history of the adVNTR databases (`hg19_muc1.db`, `hg38_muc1.db` -- built by
removing all non-MUC1 entries from the hg19 database and adding the MUC1 entry to the hg38
database where it was missing, at the start position found via UCSC BLAT) now lives in
[`berntpopp/vntyper-data`](https://github.com/berntpopp/vntyper-data), which builds and
publishes both databases as part of the `muc1` release bundle. `vntyper install-references`
installs both unconditionally, regardless of which genome assemblies are selected with
`--references`.