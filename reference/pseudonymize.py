#!/usr/bin/env python3
"""
Pseudonymize BAM and FASTQ files by:
- Grouping them by a shared "core" name (e.g. 'NPH1908593'),
- Computing per-file MD5 sums in parallel, then combining to generate a final MD5-based pseudonym,
- Reheader + index BAM, rename/copy FASTQs while preserving _R1, _R2, etc.,
- Use parallel workers to speed up I/O and reheader steps.

Generates a CSV with columns: old_base_name,combined_md5,new_pseudonym
where:
 - old_base_name is, e.g., "NPH1908593"
 - combined_md5 is the MD5 of all file-specific MD5s (sorted) for that group
 - new_pseudonym is something like "example_3a5c9f2b" (using first 8 hex chars).

New features:
 - Accept reference-assembly via --ref-assembly or a TSV/CSV mapping file.
 - Subset BAM to MUC1 region (with -P to keep paired reads). For hg19 => chr1:155158000-155163000; for hg38 => chr1:155184000-155194000.
 - Generate <new_core_name>_<ref>_subset.bam and (optionally) revert it to FASTQ in a second step
   (<new_core_name>_<ref>_subset_R1.fastq.gz, etc.).
 - Leave subsetting/reverting sequential (to avoid pickling issues in parallel).
 - Parallelize only the MD5 and reheader/copy steps.
 - Write a JSON file (defaults to `pseudonymization_output.json` in `--output-dir` if not otherwise specified)
   containing only `.bam`, `.bai`, or `.fastq.gz` files in `--output-dir`.
 - Optionally filter JSON output to only `_subset` files by specifying `--json-filter=subset`.
"""

import os
import re
import csv
import json
import shutil
import hashlib
import argparse
import subprocess
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

###############################################################################
# Logging Setup
###############################################################################
logging.basicConfig(level=logging.INFO)

# Map reference assemblies to MUC1 subsetting regions
REGION_MAP = {
    "hg19": "chr1:155158000-155163000",
    "hg38": "chr1:155184000-155194000"
}

###############################################################################
# 1) Helpers for parsing filenames and computing MD5
###############################################################################

def parse_filename(filename):
    """
    Identify:
      1) The "core" base name (e.g. 'NPH1908593'),
      2) An optional read-suffix like '_R1' or '_R2' (or '_R\\d+'),
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

    core = filename[:-len(extension)]
    match = re.search(r'(_R\d+)$', core)
    if match:
        read_suffix = match.group(1)
        core = core[: -len(read_suffix)]

    return core, read_suffix, extension


def compute_md5(filepath, chunk_size=1_048_576):
    """
    Compute the MD5 hex digest of the given file in a memory-efficient way.
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
        '-c', "grep -v ^@PG | grep -v ^@RG | grep -v ^@CO",
        input_bam
    ]
    logging.debug("Reheader command: %s", " ".join(cmd_reheader))

    with open(output_bam, 'w') as fout:
        subprocess.run(cmd_reheader, stdout=fout, check=True)

    subprocess.run(['samtools', 'index', output_bam], check=True)

###############################################################################
# 2) Subset & Revert Logic
###############################################################################

def subset_only(input_bam, region, subset_bam):
    """
    Subset the input_bam to the specified region, producing an on-disk subset_bam.
    Uses -P to keep paired reads together.
    """
    logging.info(f"Subsetting BAM {input_bam} to region {region}. Output: {subset_bam}")
    view_cmd = [
        'samtools', 'view',
        '-P',  # keep paired reads
        '-b', input_bam,
        region,
        '-o', subset_bam
    ]
    logging.debug("Subset command: %s", " ".join(view_cmd))

    subprocess.run(view_cmd, check=True)
    subprocess.run(['samtools', 'index', subset_bam], check=True)


def revert_only(subset_bam, out_r1, out_r2):
    """
    Revert the subset_bam to paired-end FASTQ using:
      samtools collate -u -O <subset_bam> - | samtools fastq ...
    Writes to out_r1/out_r2, auto-compressing if .gz is used.
    """
    logging.info(f"Reverting subset {subset_bam} -> {out_r1}, {out_r2}")
    cmd_str = (
        f"samtools collate -u -O {subset_bam} - | "
        f"samtools fastq "
        f"-1 {out_r1} "
        f"-2 {out_r2} "
        "-0 /dev/null -s /dev/null -n"
    )
    logging.debug("Revert command: %s", cmd_str)

    subprocess.run(cmd_str, shell=True, check=True)


def subset_revert_task(final_bam, region, new_core_name, file_ref, output_dir, do_revert):
    """
    - Always produce a subset .bam named: <new_core_name>_<file_ref>_subset.bam
    - Then optionally revert that .bam to FASTQ:
      <new_core_name>_<file_ref>_subset_R1.fastq.gz, etc.
    """
    # 1) Create subset .bam
    subset_bam = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset.bam")
    subset_only(final_bam, region, subset_bam)

    # 2) If do_revert, produce FASTQs
    if do_revert:
        r1 = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset_R1.fastq.gz")
        r2 = os.path.join(output_dir, f"{new_core_name}_{file_ref}_subset_R2.fastq.gz")
        revert_only(subset_bam, r1, r2)

###############################################################################
# 3) Parallel tasks for MD5 and reheader/copy
###############################################################################

def md5_of_file_task(file_path):
    """
    Wrapper to compute MD5 for a single file.
    """
    return file_path, compute_md5(file_path)


def reheader_or_copy_task(file_path, out_path, extension):
    """
    Reheader+index if it's a .bam, else copy it.
    """
    if extension == ".bam":
        run_samtools_reheader(file_path, out_path)
    else:
        shutil.copy2(file_path, out_path)
    return file_path, out_path

###############################################################################
# 4) Write JSON with file_resources
###############################################################################

def write_json_resources(output_dir, json_out, filter_mode="all"):
    """
    Gather files in output_dir, but only .bam, .bai, .fastq.gz.
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

    # Gather top-level files in output_dir
    file_list = []
    for fname in os.listdir(output_dir):
        full_path = os.path.join(output_dir, fname)
        if os.path.isfile(full_path):
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

        file_resources.append({
            "filename": rel_path,
            "url": "",
            "md5sum": md5sum
        })

    output_data = {
        "file_resources": file_resources,
        "unit_tests": {},
        "integration_tests": {}
    }

    with open(json_out, "w") as jf:
        json.dump(output_data, jf, indent=2)
    logging.info(f"Wrote JSON resource file: {json_out}")

###############################################################################
# 5) Main logic
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
    and return 'example_' + first 4 hex characters, plus the entire final MD5.
    """
    md5_list_sorted = sorted(md5_list)
    combined = "".join(md5_list_sorted)
    final_md5 = hashlib.md5(combined.encode("utf-8")).hexdigest()
    short_name = "example_" + final_md5[:4]
    return short_name, final_md5


def load_reference_mapping(mapping_file):
    """
    Load a TSV or CSV with columns [filename, reference].
    Returns a dict mapping basename -> reference.
    """
    logging.info(f"Loading reference mapping file: {mapping_file}")
    ref_map = {}
    delimiter = ','

    with open(mapping_file, 'r') as f:
        first_line = f.readline()
        if '\t' in first_line and ',' not in first_line:
            delimiter = '\t'

    with open(mapping_file, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if len(row) < 2:
                continue
            base = os.path.basename(row[0])
            ref_map[base] = row[1]
    return ref_map


def pseudonymize_files(
    input_dir,
    output_dir,
    max_workers=None,
    ref_assembly=None,
    ref_mapping_file=None,
    do_subset=False,
    do_revert=False,
    json_out=None,
    json_filter="all"
):
    """
    - Step 1: Compute MD5 for all files in parallel.
    - Step 2: Reheader/copy files in parallel.
    - Step 3: Write CSV with old_core->new_core_name.
    - Step 4: If do_subset, sequentially subset (and optionally revert).
    - Step 5: Write a JSON resource file with final outputs (only .bam, .bai, .fastq.gz),
      using filter_mode=all or subset as specified by json_filter.
    """
    os.makedirs(output_dir, exist_ok=True)
    files_by_core = collect_files(input_dir)

    # Load the reference mapping if provided
    reference_mapping = {}
    if ref_mapping_file:
        reference_mapping = load_reference_mapping(ref_mapping_file)

    # ------------------------------------------------
    # 1) MD5 all input files in parallel
    # ------------------------------------------------
    all_files = []
    for core_name, file_list in files_by_core.items():
        for (fp, _, _) in file_list:
            all_files.append(fp)

    file_md5_map = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        md5_futures = [executor.submit(md5_of_file_task, fp) for fp in all_files]
        for fut in as_completed(md5_futures):
            fp, md5sum = fut.result()
            file_md5_map[fp] = md5sum

    # ------------------------------------------------
    # 2) Reheader/copy in parallel
    # ------------------------------------------------
    parallel_tasks = []
    mapping_rows = []
    for core_name, file_list in files_by_core.items():
        md5s_for_core = [file_md5_map[fp] for (fp, _, _) in file_list]
        new_core_name, final_md5 = generate_deterministic_name(md5s_for_core)
        mapping_rows.append((core_name, final_md5, new_core_name))

        for (fp, read_suffix, extension) in file_list:
            new_filename = new_core_name + read_suffix + extension
            out_path = os.path.join(output_dir, new_filename)
            parallel_tasks.append((fp, out_path, extension))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        reheader_futures = [
            executor.submit(reheader_or_copy_task, fp, out_path, ext)
            for (fp, out_path, ext) in parallel_tasks
        ]
        for fut in as_completed(reheader_futures):
            _ = fut.result()

    # ------------------------------------------------
    # 3) Write CSV: pseudonymization_table.csv
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

    # ------------------------------------------------
    # 4) Subset/Revert (SEQUENTIAL to avoid pickling issues)
    # ------------------------------------------------
    if do_subset:
        for core_name, file_list in files_by_core.items():
            new_core_name = [r[2] for r in mapping_rows if r[0] == core_name][0]

            for (fp, read_suffix, extension) in file_list:
                if extension == ".bam":
                    new_bam = os.path.join(output_dir, new_core_name + read_suffix + extension)

                    # Reference assembly for the file
                    base_fname = os.path.basename(fp)
                    file_ref = reference_mapping.get(base_fname, ref_assembly)
                    if not file_ref or file_ref.lower() not in REGION_MAP:
                        logging.warning(
                            f"Unknown or no reference for {base_fname}, defaulting to hg19 region."
                        )
                        file_ref = "hg19"
                    file_ref = file_ref.lower()

                    region = REGION_MAP[file_ref]
                    subset_revert_task(new_bam, region, new_core_name, file_ref, output_dir, do_revert)

    # ------------------------------------------------
    # 5) Write JSON with final outputs
    # ------------------------------------------------
    write_json_resources(output_dir, json_out, filter_mode=json_filter)


def main():
    parser = argparse.ArgumentParser(
        description="Deterministic pseudonymization for BAM/FASTQ with optional subsetting/revert to FASTQ."
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing the original BAM/FASTQ files.")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for pseudonymized files.")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers (default: use all available cores).")

    parser.add_argument("--ref-assembly", default=None,
                        help="A single reference assembly to apply to all input files (e.g. hg19/hg38).")
    parser.add_argument("--ref-mapping-file", default=None,
                        help="Path to TSV/CSV file with columns [filename, reference].")

    parser.add_argument("--subset-muc1", action="store_true",
                        help="If specified, create a MUC1 region-only subset for each BAM (with -P).")
    parser.add_argument("--revert-fastq", action="store_true",
                        help="If specified (and --subset-muc1), revert the subset BAM to FASTQ in a second step.")

    parser.add_argument("--log-level", default="INFO",
                        help="Set logging level (DEBUG, INFO, WARNING, ERROR). Default=INFO.")

    # If not specified, we'll default to 'pseudonymization_output.json' in the same output directory.
    parser.add_argument("--json-out", default=None,
                        help="Write a JSON file listing all files in the output directory with MD5 sums. "
                             "If not provided, defaults to 'pseudonymization_output.json' in output-dir.")

    parser.add_argument("--json-filter", choices=["all", "subset"], default="all",
                        help="Filter mode for the JSON resource file. 'all' => all .bam, .bai, .fastq.gz. "
                             "'subset' => only files containing '_subset'. Default=all.")

    args = parser.parse_args()

    logging.getLogger().setLevel(args.log_level.upper())

    pseudonymize_files(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        max_workers=args.workers,
        ref_assembly=args.ref_assembly,
        ref_mapping_file=args.ref_mapping_file,
        do_subset=args.subset_muc1,
        do_revert=args.revert_fastq,
        json_out=args.json_out,
        json_filter=args.json_filter
    )


if __name__ == "__main__":
    main()
