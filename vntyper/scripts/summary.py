"""
vntyper/scripts/summary.py

This module provides functions to record and summarize pipeline steps.
Each step is recorded with start and end times, a command description,
a parsed result (based on file type), and an MD5 checksum for the result file.

Additionally, the summary object includes the vntyper version, input files,
and the pipeline end time.

The module also provides functions to convert the summary into CSV/TSV formats.
These conversion functions now flatten nested data structures (e.g. parsed_result)
so that all available data is expanded into individual columns.
"""

import csv
import hashlib
import json
import logging
import os
from datetime import datetime, timezone

logger = logging.getLogger(__name__)

#: Version of the ``pipeline_summary.json`` layout this module writes.
#:
#: Nothing versioned the summary before, so a consumer could not tell "this run
#: did not record the region" from "this run predates the field" - and the report
#: has to distinguish them, because the only honest rendering of the second is to
#: say the value was not recorded rather than to substitute a configured default
#: (#242). Bump this when a key is added, removed or changes meaning.
SUMMARY_SCHEMA_VERSION = 1


def start_summary(
    version=None,
    input_files=None,
    reference_assembly_requested=None,
    reference_key_used=None,
    reference_path=None,
    reference_source_effective=None,
):
    """
    Initializes a new pipeline summary.

    The four ``reference_*`` fields record how the run's reference was actually
    resolved: a UCSC-family fallback is otherwise invisible to the operator once the
    run has finished, since nothing downstream re-derives it from the reference path
    alone. They are optional and default to `None` so a caller that has not resolved a
    reference (or whose run reads none at all) still gets a normal summary.

    Only a FASTQ run's values are the BWA reference it aligned with (#163). A BAM run
    reads no reference at all, and a CRAM run decodes against whatever the alignment
    plan actually resolved - an explicit ``--reference-fasta``, a different configured
    CRAM candidate, an ambient htslib resolution, or a private M5 cache - which can
    differ entirely from the configured BWA path; naming the BWA selection for either
    would be a false provenance claim (MAJOR 5, milestone-5 PR-2 review). Callers pass
    :func:`vntyper.scripts.pipeline_alignment.resolve_summary_reference_provenance`'s
    result rather than the raw BWA resolution, so ``reference_key_used`` and
    ``reference_path`` are `None` for BAM, and ``reference_source_effective`` for BAM
    and CRAM instead carries the alignment plan's own source label (e.g.
    ``"not-required"``, ``"cli"``, ``"htslib_resolved"``) rather than a coordinate
    system name.

    Args:
        version (str, optional): vntyper version. Defaults to "unknown" if not provided.
        input_files (dict, optional): Dictionary of input files. Defaults to empty dict.
        reference_assembly_requested (str, optional): The ``--reference-assembly`` label
            the run was asked for. Recorded regardless of input type: even a BAM/CRAM
            run that reads no reference file still uses this to select the correct
            coordinate-system region, so it is never a false claim the way a BWA path
            would be.
        reference_key_used (str, optional): For FASTQ, the ``reference_data`` config key
            that supplied the BWA reference actually used. `None` for BAM and CRAM,
            which have no comparable configured key to report.
        reference_path (str, optional): The reference path actually used - the BWA
            reference for FASTQ, the alignment plan's proven reference for CRAM, or
            `None` when the run read no reference file at all (BAM) or the winning CRAM
            candidate was an ambient htslib resolution with no run-local path.
        reference_source_effective (str, optional): For FASTQ, the reference source
            ("ucsc", "ncbi" or "ensembl") the run actually used, which can differ from
            the requested assembly's own source when a UCSC-family fallback was taken.
            For BAM and CRAM, the alignment plan's own source label instead.

    Returns:
        dict: A summary dictionary with its schema version, pipeline start timestamp,
        version, input files, the effective reference selection, a placeholder for the
        region the run resolves later, and an empty steps list.
    """
    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "pipeline_start": datetime.now(timezone.utc).replace(tzinfo=None).isoformat(),
        "version": version if version is not None else "unknown",
        "input_files": input_files if input_files is not None else {},
        "reference_assembly_requested": reference_assembly_requested,
        "reference_key_used": reference_key_used,
        "reference_path": reference_path,
        "reference_source_effective": reference_source_effective,
        # The span the coverage stage actually worked over. It is not known yet -
        # `calculate_alignment_coverage` resolves it mid-run - but the key exists
        # from the start so that schema 1 means "this field is present", and an
        # absent one means the summary predates it rather than that the run
        # declined to say (#242).
        "region_resolved": None,
        "steps": [],
    }


def end_summary(summary):
    """
    Adds the pipeline end timestamp to the summary.

    Args:
        summary (dict): The summary dictionary to update.
    """
    summary["pipeline_end"] = datetime.now(timezone.utc).replace(tzinfo=None).isoformat()


def md5sum(file_path):
    """
    Calculates the MD5 checksum of the given file.

    Args:
        file_path (str): Path to the file.

    Returns:
        str: MD5 hash of the file's contents, or None if an error occurs.
    """
    hash_md5 = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except Exception:
        return None


def parse_tsv(file_path):
    """
    Parses a TSV file into a structured JSON object.
    Extracts comment lines (starting with '#') and data rows.

    Args:
        file_path (str): Path to the TSV file.

    Returns:
        dict: A dictionary with keys 'comments' (list) and 'data' (list of dicts).
    """
    comments = []
    data = []
    header = None

    try:
        with open(file_path, encoding="utf-8") as f:
            for line_number, raw_line in enumerate(f, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    comments.append(line.lstrip("#").strip())
                    continue
                if header is None:
                    header = line.split("\t")
                    continue
                row_values = line.split("\t")
                # Validate the field count per row rather than letting zip() silently
                # truncate to the shorter sequence: a ragged row means a malformed or
                # truncated file, and dropping columns without a word turns corruption
                # into a plausible-looking result. Skipping just the bad row keeps one
                # bad line from discarding the whole file, which is what strict=True
                # would do here (the exception would escape to the handler below).
                if len(row_values) != len(header):
                    logger.warning(
                        f"{file_path}:{line_number}: skipping malformed row - expected "
                        f"{len(header)} fields, found {len(row_values)}"
                    )
                    continue
                data.append(dict(zip(header, row_values, strict=True)))
    except Exception as e:
        comments.append(f"Error parsing TSV file: {e}")

    return {"comments": comments, "data": data}


def parse_csv(file_path):
    """
    Parses a CSV file into a structured JSON object.
    Assumes that any row starting with '#' is a comment.

    Args:
        file_path (str): Path to the CSV file.

    Returns:
        dict: A dictionary with keys 'comments' (list) and 'data' (list of dicts).
    """
    comments = []
    data = []
    header = None

    try:
        with open(file_path, encoding="utf-8") as f:
            reader = csv.reader(f)
            for row in reader:
                if not row:
                    continue
                if row[0].startswith("#"):
                    comment = " ".join(row).lstrip("#").strip()
                    comments.append(comment)
                    continue
                if header is None:
                    header = row
                    continue
                # See the note in parse_tsv: warn and skip the individual malformed row
                # rather than silently truncating it or discarding the whole file.
                if len(row) != len(header):
                    logger.warning(
                        f"{file_path}:{reader.line_num}: skipping malformed row - expected "
                        f"{len(header)} fields, found {len(row)}"
                    )
                    continue
                data.append(dict(zip(header, row, strict=True)))
    except Exception as e:
        comments.append(f"Error parsing CSV file: {e}")

    return {"comments": comments, "data": data}


def parse_json_file(file_path):
    """
    Reads and returns the contents of a JSON file.

    Args:
        file_path (str): Path to the JSON file.

    Returns:
        dict: The parsed JSON data, or an error dict if reading fails.
    """
    try:
        with open(file_path, encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        return {"error": f"Error reading JSON file: {e}"}


def record_step(
    summary,
    step_name,
    result_file,
    file_type,
    command,
    start_time,
    end_time,
    write_summary_path=None,
):
    """
    Records a pipeline step in the summary.

    This function calculates the MD5 checksum of the result file,
    parses the file based on its type (TSV, CSV, or JSON),
    and records the provided command and timing information.
    Optionally, if write_summary_path is provided, the summary is immediately written to that file.

    If `result_file` does not exist, an ERROR is logged and the record gains
    `result_file_missing: True`. The key is absent otherwise, so a normal run's summary is
    unchanged (#212).

    Args:
        summary (dict): The summary dictionary to update.
        step_name (str): Name/description of the pipeline step.
        result_file (str): Path to the result file generated by this step.
        file_type (str): Type of the result file ('tsv', 'csv', or 'json').
        command (str): Command or parameters used in this step.
        start_time (datetime): Start time of the step.
        end_time (datetime): End time of the step.
        write_summary_path (str, optional): File path to write the summary after recording the step.
    """
    record = {
        "step": step_name,
        "start": start_time.isoformat(),
        "end": end_time.isoformat(),
        "command": command,
        "result_file": result_file,
        "file_type": file_type,
        "md5sum": None,
        "parsed_result": None,
    }

    # #212: a step whose result file is absent is a failure, not an empty result.
    # `md5sum` swallows the FileNotFoundError into None and `parse_tsv` turns it into a
    # comment and returns `data: []`, and both the HTML report and cohort mode render an
    # empty `data` as a negative. The key is added only when the file is missing, so a
    # normal run's summary -- and the golden-cohort `pipeline_step_records` artefact --
    # is byte-identical. Recording continues rather than raising: the path is the
    # operator's clue to which stage failed, and the stage itself is where an abort
    # belongs (see `run_kestrel`).
    if not os.path.exists(result_file):
        logger.error(f"Step '{step_name}' recorded a result file that does not exist: {result_file}")
        record["result_file_missing"] = True

    # Calculate MD5 checksum
    record["md5sum"] = md5sum(result_file)

    # Parse the result file based on its type
    try:
        if file_type.lower() == "tsv":
            record["parsed_result"] = parse_tsv(result_file)
        elif file_type.lower() == "csv":
            record["parsed_result"] = parse_csv(result_file)
        elif file_type.lower() == "json":
            record["parsed_result"] = parse_json_file(result_file)
        else:
            record["parsed_result"] = {"error": f"Unsupported file type for result parsing: {file_type}"}
    except Exception as e:
        record["parsed_result"] = {"error": f"Error parsing file: {e}"}

    summary["steps"].append(record)

    if write_summary_path is not None:
        write_summary(summary, write_summary_path)


def write_summary(summary, output_path):
    """
    Writes the summary dictionary to a JSON file.

    Args:
        summary (dict): The summary dictionary.
        output_path (str): Path where the summary JSON will be written.
    """
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=4)


def flatten_dict(d, parent_key="", sep="_"):
    """
    Recursively flattens a nested dictionary.

    For any nested dictionary, keys are concatenated with the given separator.
    For lists, if elements are dictionaries, each is flattened and then joined as a JSON string;
    otherwise, list elements are joined by commas.

    Args:
        d (dict): The dictionary to flatten.
        parent_key (str, optional): The base key string for recursive calls.
        sep (str, optional): Separator between keys.

    Returns:
        dict: A flattened dictionary.
    """
    items = {}
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(flatten_dict(v, new_key, sep=sep))
        elif isinstance(v, list):
            if all(isinstance(i, dict) for i in v):
                # Flatten each dict and join the resulting JSON strings
                flattened_list = [json.dumps(flatten_dict(i, "", sep=sep)) for i in v]
                items[new_key] = "; ".join(flattened_list)
            else:
                items[new_key] = ", ".join(str(i) for i in v)
        else:
            items[new_key] = v
    return items


def convert_summary_to_csv(summary, output_csv_path):
    """
    Converts the summary steps into a CSV file.

    Each row in the CSV corresponds to a pipeline step.
    The step records are flattened so that nested data (e.g. parsed_result)
    are expanded into individual columns.

    Args:
        summary (dict): The summary dictionary.
        output_csv_path (str): Path where the CSV file will be written.
    """
    # Flatten each step record
    flattened_steps = [flatten_dict(step) for step in summary.get("steps", [])]

    # Determine all keys across the flattened records
    all_keys = set()
    for record in flattened_steps:
        all_keys.update(record.keys())
    # Use a sorted list of keys for consistent column order
    all_keys = sorted(all_keys)

    with open(output_csv_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_keys)
        writer.writeheader()
        for record in flattened_steps:
            writer.writerow({key: record.get(key, "") for key in all_keys})


def convert_summary_to_tsv(summary, output_tsv_path):
    """
    Converts the summary steps into a TSV file.

    Each row in the TSV corresponds to a pipeline step.
    The step records are flattened so that nested data (e.g. parsed_result)
    are expanded into individual columns.

    Args:
        summary (dict): The summary dictionary.
        output_tsv_path (str): Path where the TSV file will be written.
    """
    flattened_steps = [flatten_dict(step) for step in summary.get("steps", [])]
    all_keys = set()
    for record in flattened_steps:
        all_keys.update(record.keys())
    all_keys = sorted(all_keys)

    with open(output_tsv_path, "w", newline="", encoding="utf-8") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=all_keys, delimiter="\t")
        writer.writeheader()
        for record in flattened_steps:
            writer.writerow({key: record.get(key, "") for key in all_keys})


# Example usage:
if __name__ == "__main__":
    # This example demonstrates how to create a summary, record a step, and write it out.
    summary = start_summary(version="1.2.3", input_files={"sample": "sample.fastq", "bam": "sample.bam"})

    # Simulate a pipeline step with a sample result file (adjust these values as needed)
    step_name = "Example Step"
    result_file = "example_results.tsv"  # Path to your result file
    file_type = "tsv"  # Could be "tsv", "csv", or "json"
    command = "run_example --option value"
    start_time = datetime.now(timezone.utc).replace(tzinfo=None)
    # Simulate some processing delay
    end_time = datetime.now(timezone.utc).replace(tzinfo=None)

    # Record the step (this will calculate the MD5 and parse the file)
    record_step(
        summary,
        step_name,
        result_file,
        file_type,
        command,
        start_time,
        end_time,
        write_summary_path="pipeline_summary.json",
    )

    # Mark pipeline end
    end_summary(summary)

    # Write the summary to a JSON file
    write_summary(summary, "pipeline_summary.json")

    # Optionally, convert the summary to CSV and TSV formats
    convert_summary_to_csv(summary, "pipeline_summary.csv")
    convert_summary_to_tsv(summary, "pipeline_summary.tsv")
