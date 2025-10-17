#!/usr/bin/env python3
"""
Properly anonymize and pseudonymize BAM and FASTQ files.

This script:
1. Removes ALL identifying information from BAM headers (@RG, @PG, @CO lines)
2. Replaces RG:Z tags in reads with generic values
3. Adds minimal generic @RG header for tool compatibility
4. Creates deterministic pseudonyms based on file MD5 hashes
5. Optionally subsets to MUC1 region and reverts to FASTQ

Best practices implemented:
- Uses samtools addreplacerg with overwrite_all mode to replace ALL read group tags
- Removes @PG lines containing identifying file paths and commands
- Removes @CO comment lines that may contain metadata
- Adds generic @RG with anonymous sample name
- Verifies anonymization was successful
"""

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

# Required for offset-based unmapped read extraction
try:
    import pysam
except ImportError:
    raise ImportError(
        "pysam is required for offset-based unmapped read extraction.\nInstall it with: pip install pysam"
    )

###############################################################################
# Logging Setup
###############################################################################
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Map reference assemblies to MUC1 subsetting regions
REGION_MAP = {
    "hg19": "chr1:155158000-155163000",
    "hg38": "chr1:155184000-155194000",
    "grch37": "NC_000001.10:155158000-155163000",  # NCBI accession for GRCh37
    "grch38": "NC_000001.11:155184000-155194000",  # NCBI accession for GRCh38
    "hg19_ensembl": "1:155158000-155163000",  # ENSEMBL simple numeric (no chr prefix)
    "hg38_ensembl": "1:155184000-155194000",  # ENSEMBL simple numeric (no chr prefix)
}

###############################################################################
# 1) Helpers for parsing filenames and computing MD5
###############################################################################


def parse_filename(filename):
    r"""Parse filename into core name, read suffix, and extension.

    Identify:
      1) The "core" base name (e.g. 'NPH1908593'),
      2) An optional read-suffix like '_R1' or '_R2' (or '_R\d+'),
      3) The file extension, e.g. '.fastq.gz', '.fq.gz', '.fastq', '.fq', '.bam'.
    """
    possible_exts = [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam"]
    read_suffix = ""
    extension = None

    for ext in possible_exts:
        if filename.endswith(ext):
            extension = ext
            break

    if extension is None:
        return filename, "", ""

    core = filename[: -len(extension)]
    match = re.search(r"(_R\d+)$", core)
    if match:
        read_suffix = match.group(1)
        core = core[: -len(read_suffix)]

    return core, read_suffix, extension


def compute_md5(filepath, chunk_size=1_048_576):
    """Compute the MD5 hex digest of the given file in a memory-efficient way."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()


###############################################################################
# 1b) Unmapped Read Extraction Using BAI Offset
###############################################################################


def read_uint32(f):
    """Read 4 bytes from file 'f' in little-endian format as an unsigned integer."""
    return int.from_bytes(f.read(4), byteorder="little", signed=False)


def read_uint64(f):
    """Read 8 bytes from file 'f' in little-endian format as an unsigned integer."""
    return int.from_bytes(f.read(8), byteorder="little", signed=False)


def get_last_chunk_end(bai_filename):
    """Find the maximum virtual offset among all mapped regions in BAI index.

    This is the key to efficient unmapped read extraction: BAM files store
    mapped reads first, then unmapped reads at the end. The BAI index contains
    virtual offsets for all mapped chunks. By finding the MAX offset, we know
    where mapped reads end and can seek directly to unmapped reads.

    Args:
        bai_filename: Path to BAI index file

    Returns:
        int: Maximum virtual offset among all mapped regions
    """
    max_vo = 0
    with open(bai_filename, "rb") as bai:
        # Read magic (4 bytes) and number of references (4 bytes)
        bai.read(4)  # skip magic
        n_ref = read_uint32(bai)
        for _ in range(n_ref):
            n_bins = read_uint32(bai)
            for _ in range(n_bins):
                read_uint32(bai)  # bin number, not used here
                n_chunks = read_uint32(bai)
                for _ in range(n_chunks):
                    # Each chunk: 8 bytes for chunk_beg, 8 bytes for chunk_end
                    read_uint64(bai)  # chunk_beg, not used
                    chunk_end = read_uint64(bai)
                    if chunk_end > max_vo:
                        max_vo = chunk_end
            # Read number of linear index entries and skip them
            n_intv = read_uint32(bai)
            bai.seek(n_intv * 8, os.SEEK_CUR)
    return max_vo


def extract_unmapped_reads_from_offset(bam_file, bai_file, output_bam):
    """
    Extract unmapped reads from a BAM file using the offset-based approach.

    This is MUCH more efficient than `samtools view -f 4` because:
    1. BAM files store reads sequentially: mapped reads first, unmapped at end
    2. We read the BAI index to find the maximum virtual offset of mapped regions
    3. We seek directly to that offset and read only unmapped reads
    4. This avoids scanning the entire multi-GB file

    This is the exact approach used by VNtyper for efficient unmapped read extraction.

    Args:
        bam_file: Path to input BAM file
        bai_file: Path to corresponding BAI index file
        output_bam: Path for output BAM file (unmapped reads only)

    Raises:
        IOError: If reading/writing fails or BAI file is invalid
    """
    last_vo = get_last_chunk_end(bai_file)
    logging.info(f"Last mapped virtual offset (from BAI): {last_vo}")

    with pysam.AlignmentFile(bam_file, "rb") as inbam:
        # Seek to the computed virtual offset - jump directly to unmapped reads!
        inbam.seek(last_vo)
        with pysam.AlignmentFile(output_bam, "wb", header=inbam.header) as outbam:
            count = 0
            for read in inbam:
                if read.is_unmapped:
                    outbam.write(read)
                    count += 1
            logging.info(f"Extracted {count} unmapped reads to {output_bam}")


###############################################################################
# 1c) Reference Download and BWA Remapping Functions
###############################################################################


# Reference path mapping (assembly name → file path)
REFERENCE_PATHS = {
    "hg19": "reference/alignment/chr1.hg19.fa",
    "hg38": "reference/alignment/chr1.hg38.fa",
    "GRCh37": "reference/alignment/chr1.GRCh37.fna",
    "GRCh38": "reference/alignment/chr1.GRCh38.fna",
    "hg19_ensembl": "reference/alignment/chr1.hg19_ensembl.fa",
    "hg38_ensembl": "reference/alignment/chr1.hg38_ensembl.fa",
}


def get_reference_path(assembly):
    """
    Get the file path for a reference assembly.

    Args:
        assembly: Reference assembly name (hg19, hg38, GRCh37, GRCh38)

    Returns:
        Path: Path to the reference FASTA file

    Raises:
        ValueError: If assembly is not supported
    """
    if assembly not in REFERENCE_PATHS:
        raise ValueError(f"Unsupported reference assembly: {assembly}. Supported: {list(REFERENCE_PATHS.keys())}")
    return REFERENCE_PATHS[assembly]


def check_bwa_index(reference_path):
    """
    Check if BWA index files exist for the reference.

    Args:
        reference_path: Path to reference FASTA file

    Returns:
        bool: True if all BWA index files exist, False otherwise
    """
    required_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    for ext in required_extensions:
        index_file = reference_path + ext
        if not os.path.exists(index_file):
            logging.debug(f"Missing BWA index file: {index_file}")
            return False
    return True


def download_and_index_reference(assembly, script_dir):
    """
    Download and index a reference genome if not already present.

    Uses the install_references.py infrastructure to download and index.

    Args:
        assembly: Reference assembly name (hg19, hg38, GRCh37, GRCh38)
        script_dir: Directory containing pseudonymize.py

    Raises:
        RuntimeError: If download or indexing fails
    """
    logging.info(f"Downloading and indexing {assembly} reference...")
    logging.info("This will take several minutes (~10-15 min for download + indexing)")

    # Determine repository root (parent of script_dir which is reference/)
    repo_root = os.path.dirname(script_dir)

    # Run install_references to download and index
    install_script = os.path.join(repo_root, "vntyper", "scripts", "install_references.py")

    if not os.path.exists(install_script):
        raise RuntimeError(f"install_references.py not found at {install_script}. Cannot auto-download reference.")

    cmd = [
        "python",
        install_script,
        "--output-dir",
        os.path.join(repo_root, "reference"),
        "--skip-indexing",  # We'll index separately to show progress
    ]

    logging.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f"Reference download failed: {result.stderr}")
        raise RuntimeError(f"Failed to download {assembly} reference")

    logging.info(f"Successfully downloaded {assembly} reference")


def index_reference(reference_path):
    """
    Create BWA index for a reference genome.

    Args:
        reference_path: Path to reference FASTA file

    Raises:
        RuntimeError: If indexing fails
    """
    logging.info(f"Creating BWA index for {reference_path}...")
    logging.info("This will take 5-10 minutes for chromosome 1...")

    cmd = ["bwa", "index", reference_path]
    logging.info(f"Running: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f"BWA indexing failed: {result.stderr}")
        raise RuntimeError(f"Failed to index {reference_path}")

    logging.info("BWA indexing completed successfully")


def ensure_reference_ready(assembly, auto_download=True):
    """
    Ensure reference genome is downloaded and indexed.

    Checks if reference exists and has BWA index. If not, optionally downloads
    and indexes it.

    Args:
        assembly: Reference assembly name (hg19, hg38, GRCh37, GRCh38)
        auto_download: If True, auto-download and index if missing

    Returns:
        str: Path to the ready reference file

    Raises:
        FileNotFoundError: If reference missing and auto_download=False
        RuntimeError: If download or indexing fails
    """
    reference_path = get_reference_path(assembly)

    # Check if reference file exists
    if not os.path.exists(reference_path):
        if not auto_download:
            raise FileNotFoundError(
                f"Reference file not found: {reference_path}\n"
                f"Run: python vntyper/scripts/install_references.py "
                f"--output-dir reference"
            )

        # Download reference
        script_dir = os.path.dirname(os.path.abspath(__file__))
        download_and_index_reference(assembly, script_dir)

        # Verify download succeeded
        if not os.path.exists(reference_path):
            raise RuntimeError(f"Reference download completed but file not found: {reference_path}")

    # Check if BWA index exists
    if not check_bwa_index(reference_path):
        logging.info(f"BWA index not found for {assembly}, creating index...")
        index_reference(reference_path)

        # Verify indexing succeeded
        if not check_bwa_index(reference_path):
            raise RuntimeError(f"BWA indexing failed for {reference_path}")

    logging.info(f"Reference ready: {reference_path}")
    return reference_path


def remap_fastq_to_reference(fastq1, fastq2, reference_path, output_bam, threads=8):
    """
    Align paired-end FASTQs to reference genome using BWA MEM.

    Performs the following steps:
    1. BWA MEM alignment
    2. Conversion to BAM format
    3. Coordinate sorting
    4. BAM indexing

    Args:
        fastq1: Path to R1 FASTQ file
        fastq2: Path to R2 FASTQ file
        reference_path: Path to reference FASTA (must be BWA indexed)
        output_bam: Path for output sorted BAM file
        threads: Number of threads for alignment and sorting (default: 8)

    Raises:
        RuntimeError: If alignment, sorting, or indexing fails
    """
    logging.info(f"Remapping {os.path.basename(fastq1)} to {reference_path}...")

    # BWA MEM command piped directly to samtools sort
    # Modern samtools infers BAM format from .bam extension
    full_cmd = f"bwa mem -t {threads} {reference_path} {fastq1} {fastq2} | samtools sort -@ {threads} -o {output_bam} -"

    logging.info(f"Running BWA alignment: {full_cmd}")

    # Stream output to console in real-time (don't buffer)
    result = subprocess.run(full_cmd, shell=True)

    if result.returncode != 0:
        raise RuntimeError(f"Failed to remap FASTQs to {reference_path}")

    # Index the output BAM
    logging.info(f"Indexing {output_bam}...")
    index_cmd = ["samtools", "index", output_bam]
    result = subprocess.run(index_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f"BAM indexing failed: {result.stderr}")
        raise RuntimeError(f"Failed to index {output_bam}")

    logging.info(f"Successfully remapped to {output_bam}")


###############################################################################
# 1d) Multi-Aligner Support Functions
###############################################################################


def load_aligner_config(config_file="reference/pseudonymize_config.json"):
    """
    Load aligner configurations from pseudonymize_config.json.

    Args:
        config_file: Path to configuration file

    Returns:
        dict: Dictionary of aligner configurations
    """
    if not os.path.exists(config_file):
        logging.warning(f"Aligner config file not found: {config_file}")
        return {}

    try:
        with open(config_file) as f:
            config = json.load(f)
        return config.get("aligners", {})
    except Exception as e:
        logging.warning(f"Failed to load aligner config: {e}")
        return {}


def check_aligner_available(executable):
    """
    Check if an aligner executable is available in PATH.

    Args:
        executable: Name or path of executable

    Returns:
        bool: True if available, False otherwise
    """
    try:
        result = subprocess.run(["which", executable], capture_output=True, text=True, check=False)
        return result.returncode == 0
    except Exception:
        return False


def get_available_aligners(aligner_config, user_specified=None):
    """
    Get list of available aligners from configuration.

    Args:
        aligner_config: Dictionary of aligner configurations
        user_specified: List of user-specified aligner names (or None for all enabled)

    Returns:
        dict: Dictionary of available aligner names to their configs
    """
    available = {}

    for aligner_name, aligner_info in aligner_config.items():
        # Skip if not enabled in config
        if not aligner_info.get("enabled", False):
            continue

        # Skip if user specified aligners and this isn't one of them
        if user_specified and aligner_name not in user_specified:
            continue

        # Check if executable is available
        executable = aligner_info.get("executable", aligner_name)
        if check_aligner_available(executable):
            available[aligner_name] = aligner_info
        else:
            if user_specified and aligner_name in user_specified:
                logging.warning(
                    f"Aligner '{aligner_name}' was specified but executable '{executable}' not found in PATH"
                )

    return available


def check_aligner_index(ref_path, aligner_name, aligner_info):
    """
    Check if aligner index files exist for a reference.

    Args:
        ref_path: Path to reference FASTA file
        aligner_name: Name of aligner
        aligner_info: Aligner configuration dictionary

    Returns:
        bool: True if all index files exist, False otherwise
    """
    index_type = aligner_info.get("index_type", "in_place")
    index_files = aligner_info.get("index_files", [])

    if index_type == "in_place":
        # BWA, BWA-MEM2 style: ref.fa.amb, ref.fa.bwt, etc.
        for ext in index_files:
            if not os.path.exists(ref_path + ext):
                return False
        return True

    elif index_type == "separate_file":
        # Minimap2 style: ref.fa.mmi
        pattern = aligner_info.get("index_path_pattern", "{ref_path}.mmi")
        index_path = pattern.format(ref_path=ref_path)
        return os.path.exists(index_path)

    elif index_type == "index_base":
        # Bowtie2 style: ref_bowtie2.1.bt2, ref_bowtie2.2.bt2, etc.
        ref_dir = os.path.dirname(ref_path)
        ref_stem = os.path.splitext(os.path.basename(ref_path))[0]
        pattern = aligner_info.get("index_base_pattern", "{ref_dir}/{ref_stem}_{aligner}")
        index_base = pattern.format(ref_dir=ref_dir, ref_stem=ref_stem, aligner=aligner_name)

        for ext in index_files:
            if not os.path.exists(index_base + ext):
                return False
        return True

    elif index_type == "index_directory":
        # DRAGMAP style: ref_dragmap_index/hash_table.cfg, etc.
        ref_dir = os.path.dirname(ref_path)
        ref_stem = os.path.splitext(os.path.basename(ref_path))[0]
        pattern = aligner_info.get("index_dir_pattern", "{ref_dir}/{ref_stem}_{aligner}_index")
        index_dir = pattern.format(ref_dir=ref_dir, ref_stem=ref_stem, aligner=aligner_name)

        if not os.path.exists(index_dir):
            return False

        for index_file in index_files:
            if not os.path.exists(os.path.join(index_dir, index_file)):
                return False
        return True

    return False


def get_aligner_index_path(ref_path, aligner_name, aligner_info):
    """
    Get the index path/base for an aligner.

    Args:
        ref_path: Path to reference FASTA
        aligner_name: Name of aligner
        aligner_info: Aligner configuration

    Returns:
        str: Index path, index base, or index directory depending on aligner type
    """
    index_type = aligner_info.get("index_type", "in_place")

    if index_type == "in_place":
        # BWA, BWA-MEM2: just use ref_path
        return ref_path

    elif index_type == "separate_file":
        # Minimap2: ref.fa.mmi
        pattern = aligner_info.get("index_path_pattern", "{ref_path}.mmi")
        return pattern.format(ref_path=ref_path)

    elif index_type == "index_base":
        # Bowtie2: ref_bowtie2
        ref_dir = os.path.dirname(ref_path)
        ref_stem = os.path.splitext(os.path.basename(ref_path))[0]
        pattern = aligner_info.get("index_base_pattern", "{ref_dir}/{ref_stem}_{aligner}")
        return pattern.format(ref_dir=ref_dir, ref_stem=ref_stem, aligner=aligner_name)

    elif index_type == "index_directory":
        # DRAGMAP: ref_dragmap_index/
        ref_dir = os.path.dirname(ref_path)
        ref_stem = os.path.splitext(os.path.basename(ref_path))[0]
        pattern = aligner_info.get("index_dir_pattern", "{ref_dir}/{ref_stem}_{aligner}_index")
        return pattern.format(ref_dir=ref_dir, ref_stem=ref_stem, aligner=aligner_name)

    return ref_path


def filter_zero_length_reads_seqtk(fastq1, fastq2, output_fastq1, output_fastq2, min_length=1):
    """
    Filter out zero-length and very short reads using seqtk (fast, industry-standard).

    seqtk is a fast and lightweight tool for FASTQ processing. This function uses
    seqtk to filter out reads shorter than min_length while maintaining pairing.

    Args:
        fastq1: Path to input R1 FASTQ (can be gzipped)
        fastq2: Path to input R2 FASTQ (can be gzipped)
        output_fastq1: Path for filtered R1 FASTQ (gzipped)
        output_fastq2: Path for filtered R2 FASTQ (gzipped)
        min_length: Minimum read length to keep (default: 1, filters only zero-length)

    Raises:
        RuntimeError: If seqtk filtering fails
    """
    # Use seqtk seq -L to filter by minimum length
    # seqtk seq -L <min_len> filters out sequences shorter than min_len
    cmd_r1 = f"seqtk seq -L {min_length} {fastq1} | gzip > {output_fastq1}"
    cmd_r2 = f"seqtk seq -L {min_length} {fastq2} | gzip > {output_fastq2}"

    # Filter R1
    result = subprocess.run(cmd_r1, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"seqtk filtering failed for R1: {result.stderr}")

    # Filter R2
    result = subprocess.run(cmd_r2, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"seqtk filtering failed for R2: {result.stderr}")


def filter_zero_length_reads_python(fastq1, fastq2, output_fastq1, output_fastq2):
    """
    Filter out zero-length reads from paired-end FASTQs using Python (fallback method).

    Some BAM to FASTQ conversions (samtools fastq) can produce zero-length reads,
    which cause minimap2 to crash with assertion failures. This function filters
    out such reads while maintaining pairing.

    Note: This is slower than seqtk but doesn't require external dependencies.

    Args:
        fastq1: Path to input R1 FASTQ (can be gzipped)
        fastq2: Path to input R2 FASTQ (can be gzipped)
        output_fastq1: Path for filtered R1 FASTQ (gzipped)
        output_fastq2: Path for filtered R2 FASTQ (gzipped)

    Returns:
        tuple: (num_kept, num_filtered) - number of read pairs kept and filtered
    """
    import gzip

    # Determine if input is gzipped
    open_func = gzip.open if fastq1.endswith(".gz") else open

    num_kept = 0
    num_filtered = 0

    with (
        open_func(fastq1, "rt") as f1,
        open_func(fastq2, "rt") as f2,
        gzip.open(output_fastq1, "wt") as out1,
        gzip.open(output_fastq2, "wt") as out2,
    ):
        while True:
            # Read 4 lines from each file (one FASTQ record)
            try:
                # R1 record
                r1_header = f1.readline().rstrip()
                r1_seq = f1.readline().rstrip()
                r1_plus = f1.readline().rstrip()
                r1_qual = f1.readline().rstrip()

                # R2 record
                r2_header = f2.readline().rstrip()
                r2_seq = f2.readline().rstrip()
                r2_plus = f2.readline().rstrip()
                r2_qual = f2.readline().rstrip()

                # Check if we've reached EOF
                if not r1_header or not r2_header:
                    break

                # Filter out if either read has zero length
                if len(r1_seq) == 0 or len(r2_seq) == 0:
                    num_filtered += 1
                    continue

                # Write both reads
                out1.write(f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n")
                out2.write(f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n")
                num_kept += 1

            except Exception:
                break

    logging.info(f"  Python filter: kept {num_kept} pairs, filtered {num_filtered} pairs")
    return num_kept, num_filtered


def remap_with_aligner(fastq1, fastq2, reference_path, output_bam, aligner_name, aligner_info, threads=4):
    """
    Remap paired-end FASTQs using specified aligner.

    Args:
        fastq1: Path to R1 FASTQ
        fastq2: Path to R2 FASTQ
        reference_path: Path to reference FASTA
        output_bam: Path for output BAM
        aligner_name: Name of aligner to use
        aligner_info: Aligner configuration dictionary
        threads: Number of threads

    Raises:
        RuntimeError: If alignment fails
    """
    logging.info(f"  Aligning with {aligner_name}...")

    # Check if we need to filter zero-length reads (for minimap2 and other sensitive aligners)
    # This is needed because samtools fastq can produce zero-length reads that crash minimap2
    needs_filtering = aligner_name in ["minimap2"]

    working_fastq1 = fastq1
    working_fastq2 = fastq2
    filtered_fastq1 = None
    filtered_fastq2 = None

    if needs_filtering:
        # Create temporary filtered FASTQs
        temp_dir = tempfile.gettempdir()
        filtered_fastq1 = os.path.join(temp_dir, f"filtered_R1_{os.getpid()}.fastq.gz")
        filtered_fastq2 = os.path.join(temp_dir, f"filtered_R2_{os.getpid()}.fastq.gz")

        logging.debug(f"  Filtering zero-length reads for {aligner_name}...")

        # Try seqtk first (fast), fallback to Python if not available
        seqtk_available = check_aligner_available("seqtk")

        if seqtk_available:
            logging.debug("  Using seqtk for filtering (fast)")
            try:
                filter_zero_length_reads_seqtk(fastq1, fastq2, filtered_fastq1, filtered_fastq2, min_length=1)
                logging.info("  Filtered zero-length reads using seqtk")
            except RuntimeError as e:
                logging.warning(f"  seqtk filtering failed: {e}, falling back to Python method")
                filter_zero_length_reads_python(fastq1, fastq2, filtered_fastq1, filtered_fastq2)
        else:
            logging.debug("  Using Python for filtering (slower, seqtk not available)")
            filter_zero_length_reads_python(fastq1, fastq2, filtered_fastq1, filtered_fastq2)

        working_fastq1 = filtered_fastq1
        working_fastq2 = filtered_fastq2

    try:
        # Get index path for this aligner
        index_path = get_aligner_index_path(reference_path, aligner_name, aligner_info)

        # Prepare command parameters
        index_type = aligner_info.get("index_type", "in_place")
        params = {"ref_path": reference_path, "r1": working_fastq1, "r2": working_fastq2, "threads": threads}

        if index_type == "separate_file":
            params["index_path"] = index_path
        elif index_type == "index_base":
            params["index_base"] = index_path
        elif index_type == "index_directory":
            params["index_dir"] = index_path

        # Format alignment command
        align_cmd_template = aligner_info.get("alignment_command", "")
        try:
            align_cmd = align_cmd_template.format(**params)
        except KeyError as e:
            raise RuntimeError(f"Missing parameter in alignment command for {aligner_name}: {e}")

        # Pipe aligner output directly to samtools sort
        # Modern samtools infers BAM format from .bam extension
        full_cmd = f"{align_cmd} | samtools sort -@ {threads} -o {output_bam} -"

        logging.debug(f"  Command: {full_cmd}")

        # Stream output to console in real-time (don't buffer)
        result = subprocess.run(full_cmd, shell=True)

        if result.returncode != 0:
            raise RuntimeError(f"{aligner_name} alignment failed")

        # Index the BAM
        subprocess.run(["samtools", "index", output_bam], check=True)

        logging.info(f"  ✓ {aligner_name} alignment complete")

    finally:
        # Clean up temporary filtered FASTQs
        if filtered_fastq1 and os.path.exists(filtered_fastq1):
            os.remove(filtered_fastq1)
        if filtered_fastq2 and os.path.exists(filtered_fastq2):
            os.remove(filtered_fastq2)


###############################################################################
# 2) Anonymization Functions
###############################################################################


def anonymize_read_names_in_bam(input_bam, output_bam):
    """Anonymize read names by replacing the flowcell/instrument ID portion.

    Anonymizes read names with a hashed value. Preserves pairing structure.

    Illumina read names have format:
    <instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>

    We hash the first 3 fields to remove identifying information while
    preserving uniqueness and pairing.

    Args:
        input_bam: Path to input BAM
        output_bam: Path to output BAM with anonymized read names
    """
    logging.info(f"Anonymizing read names in {input_bam}")

    # Use samtools view to convert to SAM, process, and convert back
    cmd_view = ["samtools", "view", "-h", input_bam]

    # Process with Python to anonymize QNAMEs
    with (
        subprocess.Popen(cmd_view, stdout=subprocess.PIPE, text=True) as proc_view,
        subprocess.Popen(["samtools", "view", "-b", "-o", output_bam], stdin=subprocess.PIPE, text=True) as proc_write,
    ):
        flowcell_hash_cache = {}

        for line in proc_view.stdout:
            if line.startswith("@"):
                # Header line - pass through unchanged
                proc_write.stdin.write(line)
            else:
                # Alignment line - anonymize QNAME
                fields = line.split("\t")
                if len(fields) < 11:
                    # Malformed line, pass through
                    proc_write.stdin.write(line)
                    continue

                qname = fields[0]

                # Check if it's an Illumina-style read name
                parts = qname.split(":")
                if len(parts) >= 7:
                    # Illumina format: instrument:run:flowcell:lane:tile:x:y
                    # Hash the first 3 fields (instrument:run:flowcell)
                    flowcell_id = ":".join(parts[:3])

                    if flowcell_id not in flowcell_hash_cache:
                        # Create a short hash
                        hash_val = hashlib.md5(flowcell_id.encode()).hexdigest()[:8]
                        flowcell_hash_cache[flowcell_id] = f"ANON{hash_val}"

                    # Replace with anonymized prefix
                    anon_prefix = flowcell_hash_cache[flowcell_id]
                    new_qname = anon_prefix + ":" + ":".join(parts[3:])
                    fields[0] = new_qname

                # Write modified line
                proc_write.stdin.write("\t".join(fields))

        proc_write.stdin.close()
        proc_write.wait()

    # Index the output
    subprocess.run(["samtools", "index", output_bam], check=True)

    logging.info(f"Read name anonymization complete: {output_bam}")


def anonymize_bam_complete(
    input_bam, output_bam, new_sample_id="sample", forbidden_strings=None, anonymize_read_names=False
):
    """Comprehensively anonymize a BAM file.

    Performs the following steps:
    1. Replacing ALL read group tags (RG:Z) with a generic value
    2. Removing @PG (program) headers containing identifying info
    3. Removing @CO (comment) headers
    4. Adding a new generic @RG header
    5. Optionally anonymizing read names (QNAMEs)
    6. Verifying no forbidden strings remain

    Uses a multi-step approach for maximum compatibility:
    - samtools addreplacerg: Replace all RG tags in reads
    - samtools reheader: Remove @PG and @CO lines from header
    - Python script to anonymize read names if needed

    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output anonymized BAM file
        new_sample_id: Generic sample identifier to use
        forbidden_strings: List of strings that should not appear in output
        anonymize_read_names: If True, replace read names with sequential IDs
    """
    if forbidden_strings is None:
        forbidden_strings = []

    logging.info(f"Anonymizing BAM: {input_bam} -> {output_bam}")

    # Check if read names contain forbidden strings
    needs_qname_anonymization = False
    if anonymize_read_names and forbidden_strings:
        logging.info("Checking if read names contain identifying information...")
        sample_reads = subprocess.check_output(["samtools", "view", input_bam], text=True).split("\n")[
            :1000
        ]  # Check first 1000 reads

        for line in sample_reads:
            if not line:
                continue
            qname = line.split("\t")[0]
            for forbidden in forbidden_strings:
                if forbidden in qname:
                    logging.warning(f"Found '{forbidden}' in read name: {qname}")
                    needs_qname_anonymization = True
                    break
            if needs_qname_anonymization:
                break

    # Determine which BAM to use as input for next steps
    working_bam = input_bam

    # Step 0 (optional): Anonymize read names if needed
    if needs_qname_anonymization:
        logging.info("Anonymizing read names (QNAMEs)...")
        temp_bam_qname = output_bam + ".tmp_qname.bam"
        anonymize_read_names_in_bam(input_bam, temp_bam_qname)
        working_bam = temp_bam_qname

    # Step 1: Replace ALL read group tags using addreplacerg
    # This adds a new @RG line and replaces RG:Z tags in all reads
    temp_bam_1 = output_bam + ".tmp1.bam"

    rg_line = f"@RG\\tID:{new_sample_id}\\tSM:{new_sample_id}\\tPL:ILLUMINA\\tLB:{new_sample_id}"

    cmd_addreplacerg = [
        "samtools",
        "addreplacerg",
        "-r",
        rg_line,
        "-m",
        "overwrite_all",  # CRITICAL: overwrite ALL existing RG tags
        "--no-PG",  # Don't add @PG line for this operation
        "-o",
        temp_bam_1,
        working_bam,  # Use working_bam (might be qname-anonymized)
    ]

    logging.debug(f"Running: {' '.join(cmd_addreplacerg)}")
    subprocess.run(cmd_addreplacerg, check=True)

    # Step 2: Remove @PG and @CO lines from header using reheader
    # These contain identifying file paths, commands, and metadata
    temp_bam_2 = output_bam + ".tmp2.bam"

    cmd_reheader = [
        "samtools",
        "reheader",
        "-P",  # Don't add @PG line
        "-c",
        "grep -v '^@PG' | grep -v '^@CO'",
        temp_bam_1,
    ]

    logging.debug(f"Running: {' '.join(cmd_reheader)}")
    with open(temp_bam_2, "wb") as fout:
        subprocess.run(cmd_reheader, stdout=fout, check=True)

    # Step 3: Remove old @RG lines (from original), keep only our new one
    # The addreplacerg keeps old @RG lines in header, we need to remove them
    cmd_filter_rg = ["samtools", "view", "-H", temp_bam_2]

    header_lines = subprocess.check_output(cmd_filter_rg).decode().split("\n")

    # Filter: keep @HD, @SQ, and only the NEW @RG line (has our new_sample_id)
    filtered_header = []
    for line in header_lines:
        if line.startswith("@HD") or line.startswith("@SQ"):
            filtered_header.append(line)
        elif line.startswith("@RG"):
            # Only keep if it has our new sample ID
            if f"SM:{new_sample_id}" in line and f"ID:{new_sample_id}" in line:
                filtered_header.append(line)
        # Skip @PG, @CO, and old @RG lines

    # Write filtered header to temp file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sam", delete=False) as f:
        temp_header = f.name
        f.write("\n".join(filtered_header) + "\n")

    try:
        # Apply the cleaned header
        cmd_final_reheader = ["samtools", "reheader", "-P", temp_header, temp_bam_2]

        with open(output_bam, "wb") as fout:
            subprocess.run(cmd_final_reheader, stdout=fout, check=True)

        # Index the output
        subprocess.run(["samtools", "index", output_bam], check=True)

    finally:
        # Clean up temp files
        if os.path.exists(temp_header):
            os.remove(temp_header)
        if os.path.exists(temp_bam_1):
            os.remove(temp_bam_1)
        if os.path.exists(temp_bam_1 + ".bai"):
            os.remove(temp_bam_1 + ".bai")
        if os.path.exists(temp_bam_2):
            os.remove(temp_bam_2)
        if needs_qname_anonymization and os.path.exists(working_bam):
            os.remove(working_bam)
            if os.path.exists(working_bam + ".bai"):
                os.remove(working_bam + ".bai")

    # Step 4: Verify anonymization
    verify_anonymization(output_bam, forbidden_strings)

    logging.info(f"Successfully anonymized BAM: {output_bam}")


def verify_anonymization(bam_file, forbidden_strings):
    """Verify that a BAM file doesn't contain any forbidden strings.

    Checks both header and a sample of read records.

    Args:
        bam_file: Path to BAM file to verify
        forbidden_strings: List of strings that should not appear

    Raises:
        ValueError: If any forbidden string is found
    """
    if not forbidden_strings:
        logging.info("No forbidden strings specified, skipping verification")
        return

    logging.info(f"Verifying anonymization of {bam_file}")

    # Check header
    header = subprocess.check_output(["samtools", "view", "-H", bam_file]).decode()

    for forbidden in forbidden_strings:
        if forbidden in header:
            raise ValueError(f"ANONYMIZATION FAILED: Found '{forbidden}' in BAM header of {bam_file}")

    # Check sample of reads (first 10000)
    reads_output = subprocess.check_output(["samtools", "view", bam_file], text=True)

    # Sample first 10000 reads
    read_lines = reads_output.split("\n")[:10000]
    reads_text = "\n".join(read_lines)

    for forbidden in forbidden_strings:
        if forbidden in reads_text:
            raise ValueError(f"ANONYMIZATION FAILED: Found '{forbidden}' in read records of {bam_file}")

    logging.info(f"✓ Verification passed: No forbidden strings found in {bam_file}")


###############################################################################
# 3) Subset & Revert Logic
###############################################################################


def subset_only(input_bam, region, subset_bam):
    """Subset the input_bam to the specified region.

    Produces an on-disk subset_bam. Uses -P to keep paired reads together.
    """
    logging.info(f"Subsetting BAM {input_bam} to region {region}. Output: {subset_bam}")
    view_cmd = [
        "samtools",
        "view",
        "-P",  # keep paired reads
        "-b",
        input_bam,
        region,
        "-o",
        subset_bam,
    ]
    logging.debug("Subset command: %s", " ".join(view_cmd))

    subprocess.run(view_cmd, check=True)
    subprocess.run(["samtools", "index", subset_bam], check=True)


def revert_only(subset_bam, out_r1, out_r2):
    """Revert the subset_bam to paired-end FASTQ.

    Uses samtools collate and fastq commands.
    Writes to out_r1/out_r2, auto-compressing if .gz is used.
    """
    logging.info(f"Reverting subset {subset_bam} -> {out_r1}, {out_r2}")
    cmd_str = (
        f"samtools collate -u -O {subset_bam} - | samtools fastq -1 {out_r1} -2 {out_r2} -0 /dev/null -s /dev/null -n"
    )
    logging.debug("Revert command: %s", cmd_str)

    subprocess.run(cmd_str, shell=True, check=True)


def subset_revert_task(final_bam, region, new_core_name, file_ref, output_dir, do_revert, forbidden_strings=None):
    """Create subset BAM and optionally revert to FASTQ.

    Always produces a subset .bam named: <new_core_name>_<file_ref>_subset.bam
    Then optionally reverts that .bam to FASTQ:
      <new_core_name>_<file_ref>_subset_R1.fastq.gz, etc.
    """
    # 1) Create subset .bam
    subset_bam = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset.bam")
    subset_only(final_bam, region, subset_bam)

    # Verify the subset is still anonymized
    if forbidden_strings:
        verify_anonymization(subset_bam, forbidden_strings)

    # 2) If do_revert, produce FASTQs
    if do_revert:
        r1 = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset_R1.fastq.gz")
        r2 = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset_R2.fastq.gz")
        revert_only(subset_bam, r1, r2)


###############################################################################
# 4) Parallel tasks for MD5 and anonymize/copy
###############################################################################


def md5_of_file_task(file_path):
    """Compute MD5 for a single file."""
    return file_path, compute_md5(file_path)


def anonymize_or_copy_task(
    file_path, out_path, extension, new_sample_id, forbidden_strings, anonymize_read_names=False
):
    """Anonymize and index if it's a .bam, else copy it."""
    if extension == ".bam":
        anonymize_bam_complete(file_path, out_path, new_sample_id, forbidden_strings, anonymize_read_names)
    else:
        shutil.copy2(file_path, out_path)
    return file_path, out_path


###############################################################################
# 5) Write JSON with file_resources
###############################################################################


def write_json_resources(output_dir, json_out, filter_mode="all"):
    """Write JSON resource file with file metadata.

    Gathers files in output_dir recursively, but only .bam, .bai, .fastq.gz.
    If filter_mode == "all", we include all such files.
    If filter_mode == "subset", we only include files that also have "_subset" in their name.
    Then compute MD5 for each, and generate:

    {
      "file_resources": [
        {
          "filename": "relative/path/to/file",
          "url": "",
          "md5sum": "..."
        }, ...
      ],
      "unit_tests": {},
      "integration_tests": {}
    }

    If json_out is None, default to 'pseudonymization_output.json' in output_dir.
    """
    if not json_out:
        json_out = os.path.join(output_dir, "pseudonymization_output.json")

    logging.info(f"Generating JSON resource file at {json_out} (filter={filter_mode})")
    valid_exts = (".bam", ".bai", ".fastq.gz")

    # Gather files recursively in output_dir
    file_list = []
    for root, dirs, files in os.walk(output_dir):
        for fname in files:
            full_path = os.path.join(root, fname)
            file_list.append(full_path)

    file_resources = []
    for fpath in sorted(file_list):
        # Check extension
        ext = None
        for ve in valid_exts:
            if fpath.endswith(ve):
                ext = ve
                break
        if ext is None:
            continue

        # If we're in "subset" mode, only keep if "_subset" is in the filename
        base_name = os.path.basename(fpath)
        if filter_mode == "subset" and "_subset" not in base_name:
            continue

        # Compute MD5
        md5sum = compute_md5(fpath)
        rel_path = os.path.relpath(fpath, start=os.getcwd())

        file_resources.append({"filename": rel_path, "url": "", "md5sum": md5sum})

    output_data = {
        "file_resources": file_resources,
        "unit_tests": {},
        "integration_tests": {},
    }

    with open(json_out, "w") as jf:
        json.dump(output_data, jf, indent=2)
    logging.info(f"Wrote JSON resource file: {json_out}")


###############################################################################
# 6) Main logic
###############################################################################


def collect_files(input_dir):
    """Collect all relevant files into a dictionary.

    Returns dict mapping core_name -> list of (path, read_suffix, extension).
    Only processes BAM files - FASTQs will be generated from subsets if requested.
    """
    files_by_core = {}
    for fname in os.listdir(input_dir):
        # ONLY process BAM files - FASTQs are generated via --subset-muc1 --revert-fastq
        if fname.lower().endswith(".bam"):
            full_path = os.path.join(input_dir, fname)
            core, read_suffix, extension = parse_filename(fname)
            files_by_core.setdefault(core, []).append((full_path, read_suffix, extension))
    return files_by_core


def generate_deterministic_name(md5_list):
    """Generate deterministic name from list of MD5 strings.

    Sorts them, concatenates, computes final MD5.
    Returns 'example_' + first 4 hex characters, plus the entire final MD5.
    """
    md5_list_sorted = sorted(md5_list)
    combined = "".join(md5_list_sorted)
    final_md5 = hashlib.md5(combined.encode("utf-8")).hexdigest()
    short_name = "example_" + final_md5[:4]
    return short_name, final_md5


def load_reference_mapping(mapping_file):
    """Load a TSV or CSV with columns [filename, reference].

    Returns a dict mapping basename -> reference.
    """
    logging.info(f"Loading reference mapping file: {mapping_file}")
    ref_map = {}
    delimiter = ","

    with open(mapping_file) as f:
        first_line = f.readline()
        if "\t" in first_line and "," not in first_line:
            delimiter = "\t"

    with open(mapping_file) as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if len(row) < 2:
                continue
            base = os.path.basename(row[0])
            ref_map[base] = row[1]
    return ref_map


def load_forbidden_strings_from_file(filepath):
    """Load forbidden strings from a file.

    Skips comments (lines starting with #) and empty lines.

    Args:
        filepath: Path to file containing forbidden strings (one per line)

    Returns:
        List of forbidden strings
    """
    forbidden = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if line and not line.startswith("#"):
                forbidden.append(line)
    return forbidden


def extract_forbidden_strings(input_dir):
    """Extract potential identifying strings from input files.

    These are strings that should NOT appear in the anonymized output.

    Returns list of forbidden strings based on input filenames.
    """
    forbidden = set()

    for fname in os.listdir(input_dir):
        # Add the base filename (without extension) as forbidden
        core, _, ext = parse_filename(fname)
        if core and core not in ["example"]:  # Don't forbid our own naming
            forbidden.add(core)

            # Also check for sample IDs in the filename
            # Common patterns: NPH1234, SAMPLE_123, etc.
            for part in core.split("_"):
                if len(part) >= 4:  # Only meaningful identifiers
                    forbidden.add(part)

    return list(forbidden)


def pseudonymize_files(
    input_dir,
    output_dir,
    max_workers=None,
    ref_assembly=None,
    ref_mapping_file=None,
    do_subset=False,
    do_revert=False,
    remap_to_reference=None,
    aligners=None,
    json_out=None,
    json_filter="all",
    forbidden_strings_file=None,
    anonymize_read_names=False,
):
    """Pseudonymize and anonymize BAM/FASTQ files.

    Multi-step process:
    - Step 1: Compute MD5 for all files in parallel.
    - Step 2: Anonymize BAMs / copy FASTQs in parallel.
    - Step 3: Write CSV with old_core->new_core_name.
    - Step 4: If do_subset, sequentially subset (and optionally revert).
    - Step 5: If remap_to_reference, remap FASTQs to new reference genome.
    - Step 6: Write a JSON resource file with final outputs (only .bam, .bai, .fastq.gz),
      using filter_mode=all or subset as specified by json_filter.
    """
    # Parse remap_to_reference (comma-separated list)
    remap_references = []
    if remap_to_reference:
        # Split by comma and strip whitespace
        remap_references = [ref.strip() for ref in remap_to_reference.split(",")]

        # Validate all references
        valid_refs = ["hg19", "hg38", "GRCh37", "GRCh38", "hg19_ensembl", "hg38_ensembl"]
        for ref in remap_references:
            if ref not in valid_refs:
                raise ValueError(f"Invalid reference '{ref}'. Valid references: {', '.join(valid_refs)}")

        # Remapping requires FASTQ generation
        if not do_revert:
            raise ValueError(
                "--remap-to-reference requires --revert-fastq to be specified. "
                "FASTQs must be generated before remapping."
            )

        logging.info(f"Will remap FASTQs to {len(remap_references)} reference(s): {', '.join(remap_references)}")

    os.makedirs(output_dir, exist_ok=True)
    files_by_core = collect_files(input_dir)

    # Load forbidden strings (identifiers to remove)
    if forbidden_strings_file and os.path.exists(forbidden_strings_file):
        forbidden_strings = load_forbidden_strings_from_file(forbidden_strings_file)
    else:
        # Auto-detect from input filenames
        forbidden_strings = extract_forbidden_strings(input_dir)

    logging.info(f"Forbidden strings for anonymization: {forbidden_strings}")

    # Load the reference mapping if provided
    reference_mapping = {}
    if ref_mapping_file:
        reference_mapping = load_reference_mapping(ref_mapping_file)

    # ------------------------------------------------
    # 1) MD5 all input files in parallel
    # ------------------------------------------------
    all_files = []
    for core_name, file_list in files_by_core.items():
        for fp, _, _ in file_list:
            all_files.append(fp)

    logging.info(f"Computing MD5 checksums for {len(all_files)} files...")
    file_md5_map = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        md5_futures = [executor.submit(md5_of_file_task, fp) for fp in all_files]
        completed = 0
        for fut in as_completed(md5_futures):
            fp, md5sum = fut.result()
            file_md5_map[fp] = md5sum
            completed += 1
            logging.info(f"MD5 progress: {completed}/{len(all_files)} - {os.path.basename(fp)}: {md5sum[:8]}...")

    # ------------------------------------------------
    # 2) Generate pseudonyms
    # ------------------------------------------------
    mapping_rows = []

    for core_name, file_list in files_by_core.items():
        md5s_for_core = [file_md5_map[fp] for (fp, _, _) in file_list]
        new_core_name, final_md5 = generate_deterministic_name(md5s_for_core)
        mapping_rows.append((core_name, final_md5, new_core_name))

    # ------------------------------------------------
    # 4) Write CSV: pseudonymization_table.csv
    # ------------------------------------------------
    csv_path = os.path.join(output_dir, "pseudonymization_table.csv")
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["old_base_name", "combined_md5", "new_pseudonym", "reference_used"])

        for row in sorted(mapping_rows, key=lambda x: x[0]):
            old_core, final_md5, new_core_name = row

            # Determine final reference
            ref_val = ref_assembly
            for item in files_by_core[old_core]:
                base_fname = os.path.basename(item[0])
                if base_fname in reference_mapping:
                    ref_val = reference_mapping[base_fname]
                    break

            if not ref_val or ref_val.lower() not in REGION_MAP:
                ref_val = "hg19"

            writer.writerow([old_core, final_md5, new_core_name, ref_val])

    logging.info(f"Wrote pseudonymization table: {csv_path}")

    # ------------------------------------------------
    # 3) Subset + Extract Unmapped + Merge + Anonymize + Revert
    # ------------------------------------------------
    if do_subset:
        logging.info("Starting optimized workflow: subset MUC1 → extract unmapped → merge → anonymize → revert")
        for core_name, file_list in files_by_core.items():
            new_core_name = [r[2] for r in mapping_rows if r[0] == core_name][0]

            for fp, read_suffix, extension in file_list:
                if extension == ".bam":
                    # Reference assembly for the file
                    base_fname = os.path.basename(fp)
                    file_ref = reference_mapping.get(base_fname, ref_assembly)
                    if not file_ref or file_ref.lower() not in REGION_MAP:
                        logging.warning(f"Unknown or no reference for {base_fname}, defaulting to hg19 region.")
                        file_ref = "hg19"
                    file_ref = file_ref.lower()

                    region = REGION_MAP[file_ref]

                    # Step 3a: Subset ORIGINAL BAM to MUC1 region
                    temp_muc1_subset = os.path.join(output_dir, f"{new_core_name}_{file_ref}_muc1_temp.bam")
                    logging.info(f"Step 1/5: Subsetting {os.path.basename(fp)} to {region}...")
                    subset_only(fp, region, temp_muc1_subset)

                    # Step 3b: Extract unmapped reads using offset-based approach (VNtyper method)
                    temp_unmapped = os.path.join(output_dir, f"{new_core_name}_{file_ref}_unmapped_temp.bam")
                    logging.info("Step 2/5: Extracting unmapped reads using BAI offset method...")

                    # Ensure BAI file exists
                    bai_file = fp + ".bai"
                    if not os.path.exists(bai_file):
                        logging.info(f"Indexing BAM file: {fp}")
                        subprocess.run(["samtools", "index", fp], check=True)

                    # Extract unmapped reads efficiently
                    extract_unmapped_reads_from_offset(bam_file=fp, bai_file=bai_file, output_bam=temp_unmapped)

                    # Step 3c: Merge MUC1 subset + unmapped reads
                    temp_merged = os.path.join(output_dir, f"{new_core_name}_{file_ref}_merged_temp.bam")
                    logging.info("Step 3/5: Merging MUC1 subset + unmapped reads...")

                    cmd_merge = [
                        "samtools",
                        "merge",
                        "-f",  # force overwrite
                        temp_merged,
                        temp_muc1_subset,
                        temp_unmapped,
                    ]
                    subprocess.run(cmd_merge, check=True)

                    # Index the merged BAM
                    subprocess.run(["samtools", "index", temp_merged], check=True)

                    # Clean up intermediate files
                    os.remove(temp_muc1_subset)
                    if os.path.exists(temp_muc1_subset + ".bai"):
                        os.remove(temp_muc1_subset + ".bai")
                    os.remove(temp_unmapped)

                    # Step 3d: Anonymize the merged BAM
                    final_subset = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset.bam")
                    logging.info("Step 4/5: Anonymizing merged BAM...")
                    anonymize_bam_complete(
                        temp_merged, final_subset, new_core_name, forbidden_strings, anonymize_read_names
                    )

                    # Clean up temp merged file
                    if os.path.exists(temp_merged):
                        os.remove(temp_merged)
                        if os.path.exists(temp_merged + ".bai"):
                            os.remove(temp_merged + ".bai")

                    # Step 3e: Revert to FASTQ if requested
                    if do_revert:
                        # Create fastqs subdirectory
                        fastq_dir = os.path.join(output_dir, "fastqs")
                        os.makedirs(fastq_dir, exist_ok=True)

                        r1 = os.path.join(fastq_dir, f"{new_core_name}_{file_ref}_subset_R1.fastq.gz")
                        r2 = os.path.join(fastq_dir, f"{new_core_name}_{file_ref}_subset_R2.fastq.gz")
                        logging.info("Step 5/5: Reverting merged BAM to FASTQ...")
                        revert_only(final_subset, r1, r2)

                        # Step 3f: Remap FASTQs to multiple references with multiple aligners
                        if remap_references:
                            # Load aligner configuration
                            aligner_config = load_aligner_config()
                            # Default to BWA if no aligners specified
                            selected_aligners = aligners if aligners else ["bwa"]
                            available_aligners = get_available_aligners(
                                aligner_config, user_specified=selected_aligners
                            )

                            if not available_aligners:
                                logging.warning(
                                    "No aligners available for remapping. Skipping remapping step. "
                                    "Install aligners (bwa, bwa-mem2, minimap2, bowtie2) and ensure they are in PATH."
                                )
                            else:
                                logging.info(
                                    f"Step 6/6: Remapping FASTQs to {len(remap_references)} reference(s) "
                                    f"with {len(available_aligners)} aligner(s)..."
                                )
                                logging.info(f"  Aligners: {', '.join(available_aligners.keys())}")

                                for ref_idx, target_ref in enumerate(remap_references, 1):
                                    logging.info(f"  Reference [{ref_idx}/{len(remap_references)}]: {target_ref}")

                                    # Ensure reference is ready (download/index if needed)
                                    try:
                                        reference_path = ensure_reference_ready(target_ref, auto_download=True)
                                    except Exception as e:
                                        logging.error(f"Failed to prepare {target_ref} reference: {e}")
                                        logging.warning(
                                            f"Skipping {target_ref} remapping for {new_core_name}. "
                                            "FASTQs are still available for manual remapping."
                                        )
                                        continue

                                    # Remap with each available aligner
                                    for aligner_idx, (aligner_name, aligner_info) in enumerate(
                                        available_aligners.items(), 1
                                    ):
                                        logging.info(
                                            f"    Aligner [{aligner_idx}/{len(available_aligners)}]: {aligner_name}"
                                        )

                                        # Check if aligner index exists for this reference
                                        if not check_aligner_index(reference_path, aligner_name, aligner_info):
                                            logging.warning(
                                                f"    ✗ {aligner_name} index not found for {target_ref}. "
                                                f"Run: python vntyper/scripts/install_references.py -d reference/ "
                                                f"--aligners {aligner_name}"
                                            )
                                            logging.warning(f"    Skipping {aligner_name} for {target_ref}")
                                            continue

                                        # Create remapped subdirectory organized by aligner/reference
                                        # Structure: remapped/<aligner>/<reference>/
                                        remapped_dir = os.path.join(output_dir, "remapped", aligner_name, target_ref)
                                        os.makedirs(remapped_dir, exist_ok=True)

                                        # Remap FASTQs with this aligner
                                        remapped_bam = os.path.join(
                                            remapped_dir, f"{new_core_name}_{target_ref}_{aligner_name}.bam"
                                        )

                                        try:
                                            remap_with_aligner(
                                                fastq1=r1,
                                                fastq2=r2,
                                                reference_path=reference_path,
                                                output_bam=remapped_bam,
                                                aligner_name=aligner_name,
                                                aligner_info=aligner_info,
                                                threads=max_workers if max_workers else 4,
                                            )

                                            # Get stats
                                            read_count = subprocess.check_output(
                                                ["samtools", "view", "-c", remapped_bam], text=True
                                            ).strip()

                                            logging.info(f"    ✓ {aligner_name} → {target_ref}: {read_count} reads")
                                        except Exception as e:
                                            logging.error(f"    ✗ {aligner_name} alignment to {target_ref} failed: {e}")
                                            logging.warning(
                                                f"    Skipping {aligner_name} remapping to {target_ref} "
                                                f"for {new_core_name}"
                                            )
    else:
        logging.warning(
            "--subset-muc1 not specified. No output files will be generated (full BAM anonymization skipped)."
        )

    # ------------------------------------------------
    # 6) Write JSON with final outputs
    # ------------------------------------------------
    write_json_resources(output_dir, json_out, filter_mode=json_filter)

    logging.info("=" * 80)
    logging.info("ANONYMIZATION COMPLETE")
    logging.info("=" * 80)
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Pseudonymization table: {csv_path}")
    if json_out or os.path.exists(os.path.join(output_dir, "pseudonymization_output.json")):
        logging.info(f"Resource manifest: {json_out or os.path.join(output_dir, 'pseudonymization_output.json')}")


def main():
    """Execute the command-line interface."""
    parser = argparse.ArgumentParser(
        description="Properly anonymize and pseudonymize BAM/FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic anonymization with pseudonyms
  python pseudonymize.py --input-dir raw_data/ --output-dir anonymized/

  # With subsetting to MUC1 region
  python pseudonymize.py --input-dir raw_data/ --output-dir anonymized/ \\
    --ref-assembly hg19 --subset-muc1

  # Full workflow: anonymize, subset, and revert to FASTQ
  python pseudonymize.py --input-dir raw_data/ --output-dir anonymized/ \\
    --ref-assembly hg19 --subset-muc1 --revert-fastq \\
    --json-filter subset --workers 4

  # With custom forbidden strings file
  python pseudonymize.py --input-dir raw_data/ --output-dir anonymized/ \\
    --forbidden-strings forbidden.txt
        """,
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the original BAM/FASTQ files.",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory for anonymized/pseudonymized files.")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: use all available cores).",
    )

    parser.add_argument(
        "--ref-assembly",
        default=None,
        choices=["hg19", "hg38", "GRCh37", "GRCh38", "hg19_ensembl", "hg38_ensembl"],
        help="A single reference assembly to apply to all input files.",
    )
    parser.add_argument(
        "--ref-mapping-file",
        default=None,
        help="Path to TSV/CSV file with columns [filename, reference].",
    )

    parser.add_argument(
        "--subset-muc1",
        action="store_true",
        help="If specified, create a MUC1 region-only subset for each BAM (with -P).",
    )
    parser.add_argument(
        "--revert-fastq",
        action="store_true",
        help="If specified (and --subset-muc1), revert the subset BAM to FASTQ in a second step.",
    )

    parser.add_argument(
        "--remap-to-reference",
        default=None,
        help="Remap anonymized FASTQs to specified reference genome(s). "
        "Accepts comma-separated list (e.g., 'hg19,hg38,GRCh37,GRCh38' or just 'hg38'). "
        "Requires --revert-fastq. Auto-downloads references if missing. "
        "Produces remapped BAM files in output directory.",
    )

    parser.add_argument(
        "--aligners",
        nargs="+",
        default=None,
        metavar="ALIGNER",
        help="Specific aligners to use for remapping (e.g., bwa bwa-mem2 minimap2 bowtie2 dragmap). "
        "If not specified, uses all available aligners from pseudonymize_config.json. "
        "Only works with --remap-to-reference.",
    )

    parser.add_argument(
        "--forbidden-strings",
        default=None,
        help="File containing strings that must not appear in output (one per line). "
        "If not provided, will auto-detect from input filenames.",
    )

    parser.add_argument(
        "--anonymize-read-names",
        action="store_true",
        help="Also anonymize read names (QNAMEs). WARNING: Very slow, adds ~5-10min per file.",
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Set logging level. Default=INFO.",
    )

    parser.add_argument(
        "--log-file",
        default=None,
        help="Write logs to this file in addition to console output.",
    )

    parser.add_argument(
        "--json-out",
        default=None,
        help="Write a JSON file listing all files in the output directory with MD5 sums. "
        "If not provided, defaults to 'pseudonymization_output.json' in output-dir.",
    )

    parser.add_argument(
        "--json-filter",
        choices=["all", "subset"],
        default="all",
        help="Filter mode for the JSON resource file. 'all' => all .bam, .bai, .fastq.gz. "
        "'subset' => only files containing '_subset'. Default=all.",
    )

    args = parser.parse_args()

    # Configure logging
    log_level = args.log_level.upper()
    logging.getLogger().setLevel(log_level)

    # Add file handler if log file specified
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file, mode="w")
        file_handler.setLevel(log_level)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logging.getLogger().addHandler(file_handler)
        logging.info(f"Logging to file: {args.log_file}")

    pseudonymize_files(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        max_workers=args.workers,
        ref_assembly=args.ref_assembly,
        ref_mapping_file=args.ref_mapping_file,
        do_subset=args.subset_muc1,
        do_revert=args.revert_fastq,
        remap_to_reference=args.remap_to_reference,
        aligners=args.aligners,
        json_out=args.json_out,
        json_filter=args.json_filter,
        forbidden_strings_file=args.forbidden_strings,
        anonymize_read_names=args.anonymize_read_names,
    )


if __name__ == "__main__":
    main()
