"""
vntyper/scripts/summary.py

This module provides functions to record and summarize pipeline steps.
Each step is recorded with start and end times, a command description,
a parsed result (based on file type), and an MD5 checksum for the result file.

Additionally, the summary object includes the vntyper version, input files,
and the pipeline end time.
"""

import csv
import json
import hashlib
from datetime import datetime


def start_summary(version=None, input_files=None):
    """
    Initializes a new pipeline summary.

    Args:
        version (str, optional): vntyper version. Defaults to "unknown" if not provided.
        input_files (dict, optional): Dictionary of input files. Defaults to empty dict.

    Returns:
        dict: A summary dictionary with pipeline start timestamp, version, input files, and an empty steps list.
    """
    return {
        "pipeline_start": datetime.utcnow().isoformat(),
        "version": version if version is not None else "unknown",
        "input_files": input_files if input_files is not None else {},
        "steps": [],
    }


def end_summary(summary):
    """
    Adds the pipeline end timestamp to the summary.

    Args:
        summary (dict): The summary dictionary to update.
    """
    summary["pipeline_end"] = datetime.utcnow().isoformat()


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
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    comments.append(line.lstrip("#").strip())
                    continue
                if header is None:
                    header = line.split("\t")
                    continue
                row_values = line.split("\t")
                row_dict = {key: value for key, value in zip(header, row_values)}
                data.append(row_dict)
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
        with open(file_path, "r", encoding="utf-8") as f:
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
                row_dict = {key: value for key, value in zip(header, row)}
                data.append(row_dict)
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
        with open(file_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        return {"error": f"Error reading JSON file: {e}"}


def record_step(
    summary, step_name, result_file, file_type, command, start_time, end_time, write_summary_path=None
):
    """
    Records a pipeline step in the summary.

    This function calculates the MD5 checksum of the result file,
    parses the file based on its type (TSV, CSV, or JSON),
    and records the provided command and timing information.
    Optionally, if write_summary_path is provided, the summary is immediately written to that file.

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
            record["parsed_result"] = {
                "error": f"Unsupported file type for result parsing: {file_type}"
            }
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


def convert_summary_to_csv(summary, output_csv_path):
    """
    Converts the summary steps into a CSV file.

    Each row in the CSV corresponds to a pipeline step.

    Args:
        summary (dict): The summary dictionary.
        output_csv_path (str): Path where the CSV file will be written.
    """
    keys = ["step", "start", "end", "command", "result_file", "file_type", "md5sum"]
    with open(output_csv_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        for step in summary.get("steps", []):
            writer.writerow({key: step.get(key, "") for key in keys})


def convert_summary_to_tsv(summary, output_tsv_path):
    """
    Converts the summary steps into a TSV file.

    Each row in the TSV corresponds to a pipeline step.

    Args:
        summary (dict): The summary dictionary.
        output_tsv_path (str): Path where the TSV file will be written.
    """
    keys = ["step", "start", "end", "command", "result_file", "file_type", "md5sum"]
    with open(output_tsv_path, "w", newline="", encoding="utf-8") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        for step in summary.get("steps", []):
            writer.writerow({key: step.get(key, "") for key in keys})


# Example usage:
if __name__ == "__main__":
    # This example demonstrates how to create a summary, record a step, and write it out.
    summary = start_summary(
        version="1.2.3", input_files={"sample": "sample.fastq", "bam": "sample.bam"}
    )

    # Simulate a pipeline step with a sample result file (adjust these values as needed)
    step_name = "Example Step"
    result_file = "example_results.tsv"  # Path to your result file
    file_type = "tsv"  # Could be "tsv", "csv", or "json"
    command = "run_example --option value"
    start_time = datetime.utcnow()
    # Simulate some processing delay
    end_time = datetime.utcnow()

    # Record the step (this will calculate the MD5 and parse the file)
    record_step(
        summary, step_name, result_file, file_type, command, start_time, end_time, write_summary_path="pipeline_summary.json"
    )

    # Mark pipeline end
    end_summary(summary)

    # Write the summary to a JSON file
    write_summary(summary, "pipeline_summary.json")

    # Optionally, convert the summary to CSV and TSV formats
    convert_summary_to_csv(summary, "pipeline_summary.csv")
    convert_summary_to_tsv(summary, "pipeline_summary.tsv")
