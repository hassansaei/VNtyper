"""
extract_unmapped_from_offset.py

Extract unmapped reads from a BAM file by seeking to the last mapped chunk using the BAI index.
The BAI index is used to find the maximum virtual offset (chunk_end) among all mapped regions.
"""

import logging
import os

import pysam

logger = logging.getLogger(__name__)


def read_uint32(f, field="value"):
    """
    Reads 4 bytes from file 'f' in little-endian format as an unsigned integer.

    Args:
        f: A binary file-like object positioned at the field to read.
        field (str): Name of the field being read, used only in the error
            message if the read comes up short.

    Returns:
        int: The parsed unsigned 32-bit integer.

    Raises:
        ValueError: If fewer than 4 bytes are available to read. A silent
            short read here would otherwise be zero-filled by
            ``int.from_bytes`` and misread as a legitimate value of 0 -- the
            signature of a truncated or corrupt BAI index.
    """
    data = f.read(4)
    if len(data) != 4:
        msg = f"Truncated BAI index: expected 4 bytes for {field}, found {len(data)}"
        logger.error(msg)
        raise ValueError(msg)
    return int.from_bytes(data, byteorder="little", signed=False)


def read_uint64(f, field="value"):
    """
    Reads 8 bytes from file 'f' in little-endian format as an unsigned integer.

    Args:
        f: A binary file-like object positioned at the field to read.
        field (str): Name of the field being read, used only in the error
            message if the read comes up short.

    Returns:
        int: The parsed unsigned 64-bit integer.

    Raises:
        ValueError: If fewer than 8 bytes are available to read. A silent
            short read here would otherwise be zero-filled by
            ``int.from_bytes`` and misread as a legitimate value of 0 -- the
            signature of a truncated or corrupt BAI index.
    """
    data = f.read(8)
    if len(data) != 8:
        msg = f"Truncated BAI index: expected 8 bytes for {field}, found {len(data)}"
        logger.error(msg)
        raise ValueError(msg)
    return int.from_bytes(data, byteorder="little", signed=False)


def get_last_chunk_end(bai_filename):
    """
    Iterates over the BAI file to find the maximum virtual offset (chunk_end)
    among all mapped regions.

    An index with zero references, or references with zero chunks, is a
    legitimate, fully-parseable BAI and returns 0 without raising. Only a
    short read on a fixed-width field -- the signature of a truncated or
    otherwise structurally corrupt file -- raises.

    Args:
        bai_filename (str): Path to the BAI index file.

    Returns:
        int: The maximum virtual offset (chunk_end) recorded across every
            chunk in the index, or 0 if the index has no references or no
            chunks.

    Raises:
        ValueError: If the BAI file is truncated or otherwise structurally
            invalid, i.e. a fixed-width field runs out of bytes before it can
            be fully read.
    """
    max_vo = 0
    with open(bai_filename, "rb") as bai:
        # Read magic (4 bytes) and number of references (4 bytes)
        bai.read(4)  # skip magic
        n_ref = read_uint32(bai, "n_ref")
        for _ in range(n_ref):
            n_bins = read_uint32(bai, "n_bins")
            for _ in range(n_bins):
                _bin = read_uint32(bai, "bin number")  # bin number, not used here
                n_chunks = read_uint32(bai, "n_chunks")
                for _ in range(n_chunks):
                    # Each chunk: 8 bytes for chunk_beg, 8 bytes for chunk_end
                    _chunk_beg = read_uint64(bai, "chunk_beg")
                    chunk_end = read_uint64(bai, "chunk_end")
                    max_vo = max(max_vo, chunk_end)
            # Read number of linear index entries and skip them
            n_intv = read_uint32(bai, "n_intv")
            bai.seek(n_intv * 8, os.SEEK_CUR)
    return max_vo


def extract_unmapped_reads_from_offset(bam_file, bai_file, output_bam):
    """
    Extract unmapped reads from a BAM file by seeking to the maximum virtual offset
    among mapped chunks (from the BAI file), and then writing only the unmapped reads.

    Args:
        bam_file (str): Path to the input BAM file.
        bai_file (str): Path to the corresponding BAI index file.
        output_bam (str): Path for the output BAM file (unmapped reads).

    Raises:
        OSError: If reading from or writing to disk fails (e.g. a missing file or
            a permission error raised by ``open``/``pysam``).
        ValueError: If the BAI file is truncated or otherwise structurally invalid --
            propagated unchanged from ``get_last_chunk_end``.
    """
    last_vo = get_last_chunk_end(bai_file)
    logger.info(f"Last mapped virtual offset (from BAI): {last_vo}")

    with pysam.AlignmentFile(bam_file, "rb") as inbam:
        # Seek to the computed virtual offset.
        inbam.seek(last_vo)
        with pysam.AlignmentFile(output_bam, "wb", header=inbam.header) as outbam:
            count = 0
            for read in inbam:
                if read.is_unmapped:
                    outbam.write(read)
                    count += 1
            logger.info(f"Extracted {count} unmapped reads to {output_bam}")
