#!/usr/bin/env python3
"""
Pseudonymize BAM and FASTQ files by:
- Grouping them by a shared "core" name (e.g., NPH1908593).
- Computing per-file MD5 sums in parallel, then combining to generate a final MD5-based pseudonym.
- Reheader + index BAM, rename/copy FASTQs while preserving _R1, _R2, etc.
- Use parallel workers to speed up I/O and reheader steps.

Generates a CSV with columns: old_base_name,combined_md5,new_pseudonym
where:
 - old_base_name is, e.g., "NPH1908593"
 - combined_md5 is the MD5 of all file-specific MD5s (sorted) for that group
 - new_pseudonym is something like "example_3a5c9f2b" (using first 8 hex chars).
"""

import os
import re
import csv
import shutil
import hashlib
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

###############################################################################
# 1) Helpers for parsing filenames and computing MD5
###############################################################################

def parse_filename(filename):
    """
    Identify:
      1) The "core" base name (e.g. 'NPH1908593'),
      2) An optional read-suffix like '_R1' or '_R2' (or '_R\\d+'),
      3) The file extension, e.g. '.fastq.gz', '.bam', etc.

    Returns (core, read_suffix, extension).

    Examples:
      'NPH1908593.bam' -> ('NPH1908593', '', '.bam')
      'NPH1908593_R1.fastq.gz' -> ('NPH1908593', '_R1', '.fastq.gz')
      'NPH1908593_R2.fastq.gz' -> ('NPH1908593', '_R2', '.fastq.gz')
    """
    possible_exts = [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam"]
    read_suffix = ""
    extension = None

    for ext in possible_exts:
        if filename.endswith(ext):
            extension = ext
            break

    if extension is None:
        # If we can't identify the extension, treat everything as the core
        return filename, "", ""

    core = filename[:-len(extension)]
    match = re.search(r'(_R\d+)$', core)
    if match:
        read_suffix = match.group(1)
        core = core[: -len(read_suffix)]

    return core, read_suffix, extension


def compute_md5(filepath, chunk_size=1_048_576):
    """
    Compute the MD5 hex digest of the given file in a memory-efficient way.
    chunk_size=1,048,576 (1MB) by default.
    """
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()


def run_samtools_reheader(input_bam, output_bam):
    """
    Reheader the BAM to remove @pg and @rg lines, then index it.
    """
    cmd_reheader = [
        'samtools', 'reheader',
        '-P',
        '-c', "grep -v ^@pg | grep -v ^@rg",
        input_bam
    ]
    with open(output_bam, 'w') as fout:
        subprocess.run(cmd_reheader, stdout=fout, check=True)

    subprocess.run(['samtools', 'index', output_bam], check=True)


###############################################################################
# 2) Functions to handle parallel tasks
###############################################################################

def md5_of_file_task(file_path):
    """
    Wrapper to compute MD5 for a single file.
    Returns (file_path, md5sum).
    """
    return file_path, compute_md5(file_path)


def reheader_or_copy_task(file_path, out_path, extension):
    """
    Reheader+index if it's a .bam, else copy it.
    Returns (file_path, out_path).
    """
    if extension == ".bam":
        run_samtools_reheader(file_path, out_path)
    else:
        shutil.copy2(file_path, out_path)
    return (file_path, out_path)


###############################################################################
# 3) Main logic: grouping by core, parallel MD5, parallel reheader/copy
###############################################################################

def collect_files(input_dir):
    """
    Collect all relevant files into a dictionary: core_name -> list of (path, read_suffix, extension)
    """
    valid_exts = (".bam", ".fastq", ".fq", ".fastq.gz", ".fq.gz")
    files_by_core = {}
    for fname in os.listdir(input_dir):
        if fname.lower().endswith(valid_exts):
            full_path = os.path.join(input_dir, fname)
            core, read_suffix, extension = parse_filename(fname)
            files_by_core.setdefault(core, []).append((full_path, read_suffix, extension))
    return files_by_core


def generate_deterministic_name(md5_list):
    """
    Given a list of MD5 strings, sort them, concatenate, compute final MD5,
    and return 'example_' + first 8 hex characters, plus the entire final MD5.

    Returns (short_name, final_md5).
    """
    md5_list_sorted = sorted(md5_list)
    combined = "".join(md5_list_sorted)
    final_md5 = hashlib.md5(combined.encode("utf-8")).hexdigest()
    short_name = "example_" + final_md5[:4]
    return short_name, final_md5


def pseudonymize_files(input_dir, output_dir, max_workers=None):
    """
    - Step 1: Collect files by core.
    - Step 2: In parallel, compute MD5 for each file, gather results in a dict: file_path -> md5.
    - Step 3: For each core, compute final MD5-based pseudonym.
    - Step 4: In parallel, reheader or copy all files to the new pseudonym name.
    - Step 5: Write CSV with (old_core, combined_md5, new_pseudonym).
    """
    os.makedirs(output_dir, exist_ok=True)
    files_by_core = collect_files(input_dir)

    # Create a list of all files to MD5
    all_files = []
    for core_name, file_list in files_by_core.items():
        for (fp, _, _) in file_list:
            all_files.append(fp)

    # --- Step 2: compute MD5 in parallel ---
    file_md5_map = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(md5_of_file_task, fp) for fp in all_files]
        for fut in as_completed(futures):
            fp, md5sum = fut.result()
            file_md5_map[fp] = md5sum

    # We'll build a list of tasks for the reheader/copy in parallel
    parallel_tasks = []
    # We'll store the mapping for CSV
    mapping = []

    # --- Step 3: For each core, determine new pseudonym ---
    for core_name, file_list in files_by_core.items():
        # Gather MD5 sums for this core
        md5s_for_core = [file_md5_map[fp] for (fp, _, _) in file_list]
        new_core_name, final_md5 = generate_deterministic_name(md5s_for_core)

        # Keep track for CSV
        mapping.append((core_name, final_md5, new_core_name))

        # Prepare reheader/copy tasks
        for (fp, read_suffix, extension) in file_list:
            new_filename = new_core_name + read_suffix + extension
            out_path = os.path.join(output_dir, new_filename)

            parallel_tasks.append((fp, out_path, extension))

    # --- Step 4: Reheader or copy all files in parallel ---
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(reheader_or_copy_task, fp, out_path, ext)
            for (fp, out_path, ext) in parallel_tasks
        ]
        for fut in as_completed(futures):
            # We don't need the return value in detail, but let's do .result() to catch exceptions
            _ = fut.result()

    # --- Step 5: Write CSV ---
    csv_path = os.path.join(output_dir, "pseudonymization_table.csv")
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["old_base_name", "combined_md5", "new_pseudonym"])
        for row in sorted(mapping, key=lambda x: x[0]):
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Deterministic pseudonymization for BAM/FASTQ using parallel MD5 and reheader/copy steps."
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing the original BAM/FASTQ files.")
    parser.add_argument("--output-dir", required=True, help="Output directory for pseudonymized files.")
    parser.add_argument("--workers", type=int, default=None, 
                        help="Number of parallel workers (default: use all available cores).")
    args = parser.parse_args()

    pseudonymize_files(args.input_dir, args.output_dir, max_workers=args.workers)

if __name__ == "__main__":
    main()
