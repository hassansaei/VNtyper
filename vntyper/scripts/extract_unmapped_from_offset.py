"""
extract_unmapped_from_offset.py

Extract unmapped reads from a BAM file by seeking to the last mapped chunk using the BAI index.
The BAI index is used to find the maximum virtual offset (chunk_end) among all mapped regions.
"""

import logging
import os

import pysam

logger = logging.getLogger(__name__)

#: The four bytes every BAI index starts with, per the SAM/BAM specification's
#: description of the BAI index format. Sibling index formats that live beside a
#: BAM under a neighbouring extension start with their own signature --
#: ``CSI\x01`` for a CSI index -- so comparing these four bytes is what
#: distinguishes a BAI from a file that merely happens to parse as one.
#: Confirmed against a `samtools index` output: its first four bytes are
#: ``b"BAI\x01"``.
BAI_MAGIC = b"BAI\x01"

#: Bytes per entry in a reference's linear index: each entry is one uint64
#: virtual offset.
LINEAR_INDEX_ENTRY_SIZE = 8


def read_exact(f, size, field="value"):
    """
    Reads exactly ``size`` bytes from ``f``, or raises.

    Args:
        f: A binary file-like object positioned at the field to read.
        size (int): The number of bytes the field occupies.
        field (str): Name of the field being read, used only in the error
            message if the read comes up short.

    Returns:
        bytes: Exactly ``size`` bytes.

    Raises:
        ValueError: If fewer than ``size`` bytes are available to read.
    """
    data = f.read(size)
    if len(data) != size:
        msg = f"Truncated BAI index: expected {size} bytes for {field}, found {len(data)}"
        logger.error(msg)
        raise ValueError(msg)
    return data


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
    return int.from_bytes(read_exact(f, 4, field), byteorder="little", signed=False)


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
    return int.from_bytes(read_exact(f, 8, field), byteorder="little", signed=False)


def get_last_chunk_end(bai_filename):
    """
    Iterates over the BAI file to find the maximum virtual offset (chunk_end)
    among all mapped regions.

    An index with zero references, or references with zero chunks, is a
    legitimate, fully-parseable BAI and returns 0 without raising. Three
    things make a file structurally invalid and each of them raises:

    * its first four bytes are not ``BAI_MAGIC``;
    * a fixed-width field runs out of bytes before it can be fully read;
    * a reference declares a linear index whose bytes are not all present in
      the file.

    The third needs a check of its own. The linear index is stepped over with
    ``seek`` rather than read, and seeking past the end of a file succeeds, so
    the seek returning normally says nothing about whether those bytes exist.
    For every reference but the last, a missing linear-index tail is caught
    incidentally by the following fixed-width read landing at EOF; for the
    final reference nothing follows, so before this check the walk simply
    ended and ``max_vo`` was returned from a file that stops mid-structure.
    The check compares offsets against the file's real length rather than
    consuming the bytes, so a corrupt ``n_intv`` (up to 2**32-1 entries, ~34
    GB) costs nothing to reject.

    The length check is a lower bound, deliberately: ``samtools index`` writes
    an optional trailing ``n_no_coor`` after the last reference, which this
    parser never reads, so trailing bytes are not an error.

    Args:
        bai_filename (str): Path to the BAI index file.

    Returns:
        int: The maximum virtual offset (chunk_end) recorded across every
            chunk in the index, or 0 if the index has no references or no
            chunks.

    Raises:
        ValueError: If the BAI file is truncated or otherwise structurally
            invalid.
    """
    max_vo = 0
    with open(bai_filename, "rb") as bai:
        # The file's real length, taken once, is what the declared linear-index
        # sizes below are checked against.
        bai.seek(0, os.SEEK_END)
        file_size = bai.tell()
        bai.seek(0, os.SEEK_SET)

        # Read magic (4 bytes) and number of references (4 bytes)
        magic = read_exact(bai, 4, "magic")
        if magic != BAI_MAGIC:
            msg = f"Invalid BAI index: expected magic {BAI_MAGIC!r}, found {magic!r} in {bai_filename}"
            logger.error(msg)
            raise ValueError(msg)

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
            # Read the number of linear index entries, check the file is long
            # enough to hold them, and step over them.
            n_intv = read_uint32(bai, "n_intv")
            declared_bytes = n_intv * LINEAR_INDEX_ENTRY_SIZE
            linear_index_start = bai.tell()
            if linear_index_start + declared_bytes > file_size:
                found = file_size - linear_index_start
                msg = f"Truncated BAI index: expected {declared_bytes} bytes for linear index, found {found}"
                logger.error(msg)
                raise ValueError(msg)
            bai.seek(linear_index_start + declared_bytes, os.SEEK_SET)
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
