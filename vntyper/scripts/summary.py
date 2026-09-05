"""
vntyper/scripts/summary.py

This module provides functions to record and summarize pipeline steps.
Each step is recorded with start and end times, a command description,
a parsed result (based on file type), and an MD5 checksum for the result file.

Additionally, the summary object includes the vntyper version, input files,
and the pipeline end time.

The module also writes the operator-facing ``pipeline_summary.<csv|tsv>`` (one row per
step, the run's provenance first) and ``pipeline_summary_rows.<csv|tsv>`` (one row per
step, result row and field). The cells are computed by ``summary_flattening.py``; no
cell is JSON text (#119).
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
from collections.abc import Mapping
from datetime import datetime, timezone
from typing import Any

from vntyper.scripts.decision_profile import ResolvedDecisionProfile, load_packaged_decision_profile
from vntyper.scripts.profile_provenance import profile_summary_fields
from vntyper.scripts.summary_flattening import (
    LONG_COLUMNS,
    PARSED_RESULT_PREFIX,
    STEP_COLUMNS,
    long_rows,
    run_columns,
    step_rows,
)

logger = logging.getLogger(__name__)

#: Version of the ``pipeline_summary.json`` layout this module writes.
#:
#: Nothing versioned the summary before, so a consumer could not tell "this run
#: did not record the region" from "this run predates the field" - and the report
#: has to distinguish them, because the only honest rendering of the second is to
#: say the value was not recorded rather than to substitute a configured default
#: (#242). Bump this when a key is added, removed or changes meaning.
#:
#: Schema 2 adds ``sample_name_is_explicit``. Schema 1 recorded ``sample_name`` alone,
#: which cannot distinguish the operator's ``--sample-name`` from a name the CLI
#: derived from an input path - and the report has to, because the second is a
#: ``Path.stem`` it must finish deriving while the first is a name to print verbatim.
#: Schema 3 requires the complete molecular-identity quartet on positive caller rows
#: and records the exact run-local decision-profile snapshot identity.
SUMMARY_SCHEMA_VERSION = 3

#: Packaged caller-selection policy recorded by current summaries.
DEFAULT_DECISION_POLICY = "legacy-selection-v1"


def start_summary(
    version=None,
    input_files=None,
    sample_name=None,
    sample_name_is_explicit=False,
    reference_assembly_requested=None,
    reference_key_used=None,
    reference_path=None,
    reference_source_effective=None,
    advntr_evidence_digest=None,
    decision_profile: ResolvedDecisionProfile | None = None,
    canonical_input_files: dict[str, str] | None = None,
    analysis_settings: dict[str, Any] | None = None,
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
        sample_name (str, optional): What the run calls this sample - the same string
            Kestrel embeds in its outputs and VCF header. Recorded verbatim, including
            the literal ``"sample"`` ``cli_handlers`` falls back to when it can resolve
            nothing: this is the run's own record of what it used, and deciding what to
            do with a value belongs to the consumer, not here. Without it the HTML
            report had no way to reach the operator's ``--sample-name`` and derived a
            different name from the input basename, so one run produced two identities
            (#242).
        sample_name_is_explicit (bool, optional): Whether ``sample_name`` is the
            operator's own ``--sample-name`` (True) or a name ``cli_handlers`` derived
            from an input path (False). The string cannot answer this and guessing from
            its shape was wrong in both directions: a derived ``S1_R1.fastq`` is a
            ``Path.stem`` the report must finish deriving into ``S1``, while an explicit
            ``sample`` is a name to print as given. Recording the provenance beside the
            value is what lets the report do both without changing the string Kestrel
            names its output files with. Defaults to False - a caller that does not say
            did not have an operator's name to record.
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
        advntr_evidence_digest (str, optional): Full canonical digest of the
            run-snapshotted governed adVNTR evidence, or None when adVNTR was not used.
        decision_profile: Profile resolved for this run. Direct compatibility callers
            that have no run context load the packaged profile.

    Returns:
        dict: A summary dictionary with its schema version, decision policy, pipeline
        start timestamp, version, input files, the run's sample name and where it came
        from, the effective reference selection,
        a placeholder for the region the run resolves later, and an empty steps list.
    """
    if advntr_evidence_digest is not None and (
        not isinstance(advntr_evidence_digest, str) or re.fullmatch(r"[0-9a-f]{64}", advntr_evidence_digest) is None
    ):
        message = "adVNTR evidence digest must be 64 lowercase hexadecimal characters or None"
        logger.error(message)
        raise ValueError(message)
    resolved_profile = decision_profile if decision_profile is not None else load_packaged_decision_profile()
    if not isinstance(resolved_profile, ResolvedDecisionProfile):
        message = "summary decision profile must be a resolved decision profile"
        logger.error(message)
        raise ValueError(message)
    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "decision_policy": DEFAULT_DECISION_POLICY,
        "advntr_evidence_digest": advntr_evidence_digest,
        **profile_summary_fields(resolved_profile),
        "pipeline_start": datetime.now(timezone.utc).replace(tzinfo=None).isoformat(),
        "version": version if version is not None else "unknown",
        "input_files": input_files if input_files is not None else {},
        "canonical_input_files": canonical_input_files,
        "analysis_settings": analysis_settings,
        "sample_name": sample_name,
        "sample_name_is_explicit": bool(sample_name_is_explicit),
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
        with open(file_path, encoding="utf-8", newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next_physical_line = 1
            for row_values in reader:
                # csv.reader yields logical records, so a quoted cell may span more
                # than one physical line. Keep the first physical line for diagnostics
                # before advancing past every line this logical record consumed.
                line_number = next_physical_line
                next_physical_line = reader.line_num + 1
                if not row_values or not any(value.strip() for value in row_values):
                    continue
                if row_values[0].startswith("#"):
                    comments.append("\t".join(row_values).lstrip("#").strip())
                    continue
                # These are writer-controlled result files. Standard CSV semantics
                # deliberately make an unclosed quote consume later physical lines;
                # the field-count check below then rejects the malformed logical row.
                if header is None:
                    header = row_values
                    continue
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


def refresh_step(summary, step_name, write_summary_path=None):
    """Re-read a recorded step's result file and update its parsed result and md5.

    A step is normally recorded once, when it finishes. The cross-caller nomenclature
    stage runs *after* the Kestrel step has been recorded and rewrites
    ``kestrel_result.tsv`` in place, so without this the summary -- and the HTML
    report and cohort tables built from it -- keep the pre-reconciliation row: the
    tier before promotion, and an empty adVNTR column. The stored checksum no longer
    matched the file on disk either.

    Args:
        summary (dict): The summary to update in place.
        step_name (str): The recorded step to refresh.
        write_summary_path (str, optional): Rewrite the summary file after updating.

    Returns:
        bool: True when a step was found and refreshed.
    """
    for record in summary.get("steps", []):
        if record.get("step") != step_name:
            continue
        result_file = record.get("result_file")
        if not result_file or not os.path.exists(result_file):
            logger.debug(f"Cannot refresh step '{step_name}': {result_file} is absent.")
            return False
        record["md5sum"] = md5sum(result_file)
        # No handler here: each parser already converts its own failure into an
        # ``error`` entry rather than raising, so wrapping this would add a blind
        # except that can never fire.
        record["parsed_result"] = _parse_by_type(result_file, record.get("file_type", "tsv"))
        if write_summary_path:
            write_summary(summary, write_summary_path)
        return True
    logger.debug(f"Cannot refresh step '{step_name}': it is not in the summary.")
    return False


def _parse_by_type(result_file, file_type):
    """Parse a result file according to its recorded type."""
    if file_type.lower() == "tsv":
        return parse_tsv(result_file)
    if file_type.lower() == "csv":
        return parse_csv(result_file)
    if file_type.lower() == "json":
        return parse_json_file(result_file)
    return {"error": f"Unsupported file type for result parsing: {file_type}"}


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


def _summary_table(summary: Mapping[str, Any]) -> tuple[list[str], list[dict[str, str]]]:
    """Assemble the one-row-per-step table.

    Columns are the run provenance first, then :data:`STEP_COLUMNS`, then every
    ``parsed_result_*`` column any step produced, sorted; a cell a step did not produce
    is blank.

    Args:
        summary: The summary dictionary.

    Returns:
        tuple[list[str], list[dict[str, str]]]: The column names and the rows.
    """
    run = run_columns(summary)
    rows = step_rows(summary)
    parsed_columns = sorted({column for row in rows for column in row if column.startswith(PARSED_RESULT_PREFIX)})
    fieldnames = [*run, *STEP_COLUMNS, *parsed_columns]
    return fieldnames, [{**run, **row} for row in rows]


def _write_delimited(
    output_path: str | os.PathLike[str], fieldnames: list[str], rows: list[dict[str, str]], delimiter: str
) -> None:
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter, restval="")
        writer.writeheader()
        writer.writerows(rows)


def convert_summary_to_csv(summary, output_csv_path):
    """
    Writes the summary as a CSV table with one row per pipeline step.

    Every row starts with the run's ``run_*`` provenance columns, then the step record,
    then the ``parsed_result_*`` cells (see :mod:`vntyper.scripts.summary_flattening`).

    Args:
        summary (dict): The summary dictionary.
        output_csv_path (str): Path where the CSV file will be written.
    """
    fieldnames, rows = _summary_table(summary)
    _write_delimited(output_csv_path, fieldnames, rows, ",")


def convert_summary_to_tsv(summary, output_tsv_path):
    """
    Writes the summary as a TSV table with one row per pipeline step.

    Same table as :func:`convert_summary_to_csv`, tab-delimited.

    Args:
        summary (dict): The summary dictionary.
        output_tsv_path (str): Path where the TSV file will be written.
    """
    fieldnames, rows = _summary_table(summary)
    _write_delimited(output_tsv_path, fieldnames, rows, "\t")


def convert_summary_rows_to_delimited(summary, output_path, delimiter=","):
    """
    Writes every result row of every step as a long table.

    The columns are ``step``, ``row_index``, ``field`` and ``value``; a step whose result
    has several rows (adVNTR, cross-match) is complete here, where the per-step table
    only records its row count.

    Args:
        summary (dict): The summary dictionary.
        output_path (str): Path where the rows file will be written, conventionally
            ``pipeline_summary_rows.csv`` or ``pipeline_summary_rows.tsv``.
        delimiter (str, optional): ``","`` (default) or ``"\\t"``.
    """
    _write_delimited(output_path, list(LONG_COLUMNS), long_rows(summary), delimiter)
