#!/usr/bin/env python3
# vntyper/scripts/generate_report.py

import os
import logging
import subprocess
import json
from datetime import datetime
from pathlib import Path
import re

import pandas as pd
from jinja2 import Environment, FileSystemLoader

from vntyper.scripts.utils import load_config


def load_pipeline_summary(summary_file_path):
    """
    Loads the pipeline summary JSON file generated by the pipeline.

    Args:
        summary_file_path (str or Path): Path to the pipeline summary file.

    Returns:
        dict: The loaded summary dictionary or an empty dict if load fails.
    """
    logging.info("Loading pipeline summary from %s", summary_file_path)
    if not os.path.exists(summary_file_path):
        logging.error("Pipeline summary file not found: %s", summary_file_path)
        return {}
    try:
        with open(summary_file_path, "r") as f:
            summary = json.load(f)
        logging.debug("Pipeline summary loaded successfully.")
        return summary
    except Exception as e:
        logging.error("Failed to load pipeline summary: %s", e)
        return {}


def run_igv_report(
    bed_file, bam_file, fasta_file, output_html, flanking=50, vcf_file=None, config=None
):
    """
    Wrapper around `create_report` IGV command. If config is provided and flanking
    is not explicitly set, we fall back to config's default_values.flanking.
    Skips passing None for track arguments (vcf_file or bam_file).
    """
    logging.debug("run_igv_report called with:")
    logging.debug("  bed_file=%s", bed_file)
    logging.debug("  bam_file=%s", bam_file)
    logging.debug("  fasta_file=%s", fasta_file)
    logging.debug("  output_html=%s", output_html)
    logging.debug("  vcf_file=%s", vcf_file)
    logging.debug("  flanking=%s", flanking)

    if config is not None and flanking == 50:
        flanking = config.get("default_values", {}).get("flanking", 50)
        logging.debug("Flanking region set to %s based on config.", flanking)

    bed_file = str(bed_file) if bed_file else None
    bam_file = str(bam_file) if bam_file else None
    fasta_file = str(fasta_file) if fasta_file else None
    output_html = str(output_html) if output_html else None

    igv_report_cmd = [
        "create_report",
        bed_file,
        "--flanking",
        str(flanking),
        "--fasta",
        fasta_file,
        "--tracks",
    ]
    tracks = []
    if vcf_file:
        tracks.append(str(vcf_file))
    if bam_file:
        tracks.append(str(bam_file))
    if not tracks:
        logging.warning(
            "No valid tracks (VCF or BAM) provided to IGV. The IGV report may be empty."
        )
    igv_report_cmd.extend(tracks)
    igv_report_cmd.extend(["--output", output_html])

    logging.debug(
        "IGV report command: %s", " ".join([str(x) for x in igv_report_cmd if x])
    )
    try:
        logging.info(
            "Running IGV report: %s", " ".join([str(x) for x in igv_report_cmd if x])
        )
        subprocess.run(igv_report_cmd, check=True)
        logging.info("IGV report successfully generated at %s", output_html)
    except subprocess.CalledProcessError as e:
        logging.error("Error generating IGV report: %s", e)
        raise
    except Exception as e:
        logging.error("Unexpected error generating IGV report: %s", e)
        raise


def extract_igv_content(igv_report_html):
    """
    Reads the generated IGV HTML report and extracts the IGV content,
    the tableJson variable, and the sessionDictionary variable from the script.
    Returns empty strings if not found or on error.
    """
    logging.debug("extract_igv_content called with igv_report_html=%s", igv_report_html)
    try:
        with open(igv_report_html, "r") as f:
            content = f.read()

        igv_start = content.find('<div id="container"')
        igv_end = content.find("</body>")

        if igv_start == -1 or igv_end == -1:
            logging.error("Failed to extract IGV content from report.")
            return "", "", ""

        igv_content = content[igv_start:igv_end].strip()

        table_json_start = content.find("const tableJson = ") + len(
            "const tableJson = "
        )
        table_json_end = content.find("\n", table_json_start)
        table_json = content[table_json_start:table_json_end].strip()

        session_dict_start = content.find("const sessionDictionary = ") + len(
            "const sessionDictionary = "
        )
        session_dict_end = content.find("\n", session_dict_start)
        session_dictionary = content[session_dict_start:session_dict_end].strip()

        logging.info(
            "Successfully extracted IGV content, tableJson, and sessionDictionary."
        )
        return igv_content, table_json, session_dictionary
    except FileNotFoundError:
        logging.error("IGV report file not found: %s", igv_report_html)
        return "", "", ""
    except Exception as e:
        logging.error("Unexpected error extracting IGV content: %s", e)
        return "", "", ""


def load_fastp_output(fastp_file):
    """
    Loads fastp JSON output (e.g., output.json) for summary metrics if available.
    Returns an empty dict if file not found or if parsing fails.
    """
    logging.debug("load_fastp_output called with fastp_file=%s", fastp_file)
    if not os.path.exists(fastp_file):
        logging.warning("fastp output file not found: %s", fastp_file)
        return {}
    try:
        with open(fastp_file, "r") as f:
            data = json.load(f)
        logging.debug("fastp output successfully loaded.")
        return data
    except Exception as e:
        logging.error("Failed to load or parse fastp output: %s", e)
        return {}


def load_pipeline_log(log_file):
    """
    Loads the pipeline log content from the specified log_file.
    Returns a placeholder string if not found or on error.
    """
    logging.info("Loading pipeline log from %s", log_file)
    if not log_file:
        logging.warning("No pipeline log file provided; skipping log loading.")
        return "No pipeline log file was provided."
    if not os.path.exists(log_file):
        logging.warning("Pipeline log file not found: %s", log_file)
        return "Pipeline log file not found."
    try:
        with open(log_file, "r") as f:
            content = f.read()
        logging.debug("Pipeline log successfully loaded.")
        return content
    except Exception as e:
        logging.error("Failed to read pipeline log file: %s", e)
        return "Failed to load pipeline log."


def load_report_config():
    """
    Loads the report-specific configuration from 'report_config.json'
    located in the same directory as this script.

    Returns:
        dict: The loaded report configuration dictionary.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, "report_config.json")
    try:
        with open(config_path, "r") as f:
            report_config = json.load(f)
        logging.info("Loaded report config from %s", config_path)
        return report_config
    except Exception as e:
        logging.error("Failed to load report config: %s", e)
        return {}


def compute_algorithm_result(df, logic_config):
    """
    Computes the algorithm result (for Kestrel or adVNTR) based on the provided logic configuration.
    Iterates over each rule in logic_config["rules"]. For each condition, compares the plain text value
    from the DataFrame (using the column name) with the expected value.

    Supported operators are: "==", "!=", "in", and "not in". If expected is a list, membership is checked.
    Returns the rule's "result" if matched; otherwise returns logic_config["default"].

    Args:
        df (pandas.DataFrame): DataFrame containing the results.
        logic_config (dict): Configuration dictionary with rules.

    Returns:
        str: The computed algorithm result.
    """
    if df.empty:
        logging.debug("DataFrame is empty; returning default result.")
        return logic_config.get("default", "none")
    row = df.iloc[0]
    logging.debug("Data row for evaluation: %s", row.to_dict())
    logging.debug("Logic configuration: %s", logic_config)
    for idx, rule in enumerate(logic_config.get("rules", [])):
        logging.debug("Evaluating rule %s: %s", idx, rule)
        conditions = rule.get("conditions", {})
        rule_matches = True
        for col, expected in conditions.items():
            if col not in row:
                logging.debug("Rule %s: Column '%s' not found; rule fails.", idx, col)
                rule_matches = False
                break
            actual = str(row.get(col, "")).strip()
            logging.debug(
                "Rule %s, column '%s': actual='%s', expected='%s'",
                idx,
                col,
                actual,
                expected,
            )
            if isinstance(expected, dict):
                op = expected.get("operator")
                exp_val = expected.get("value")
                if op == "==":
                    if actual != str(exp_val).strip():
                        logging.debug(
                            "Rule %s: Condition '%s == %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                elif op == "!=":
                    if actual == str(exp_val).strip():
                        logging.debug(
                            "Rule %s: Condition '%s != %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                elif op == "in":
                    if not isinstance(exp_val, list):
                        exp_val = [exp_val]
                    if actual not in exp_val:
                        logging.debug(
                            "Rule %s: Condition '%s in %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                elif op == "not in":
                    if not isinstance(exp_val, list):
                        exp_val = [exp_val]
                    if actual in exp_val:
                        logging.debug(
                            "Rule %s: Condition '%s not in %s' not met (actual='%s').",
                            idx,
                            col,
                            exp_val,
                            actual,
                        )
                        rule_matches = False
                        break
                else:
                    logging.debug(
                        "Rule %s: Unsupported operator '%s' for column '%s'.",
                        idx,
                        op,
                        col,
                    )
                    rule_matches = False
                    break
            else:
                if isinstance(expected, list):
                    if actual not in expected:
                        logging.debug(
                            "Rule %s: Condition '%s in %s' not met (actual='%s').",
                            idx,
                            col,
                            expected,
                            actual,
                        )
                        rule_matches = False
                        break
                else:
                    if actual != str(expected):
                        logging.debug(
                            "Rule %s: Condition '%s == %s' not met (actual='%s').",
                            idx,
                            col,
                            expected,
                            actual,
                        )
                        rule_matches = False
                        break
        if rule_matches:
            result = rule.get("result")
            logging.debug("Rule %s PASSED; returning result: %s", idx, result)
            return result
        else:
            logging.debug("Rule %s did not pass.", idx)
    logging.debug("No rule matched; returning default result.")
    return logic_config.get("default", "none")


def build_screening_summary(
    kestrel_df,
    advntr_df,
    advntr_available,
    mean_vntr_coverage,
    mean_vntr_cov_threshold,
    report_config,
):
    """
    Build the detailed screening summary text based on Kestrel and adVNTR data.

    Args:
        kestrel_df (pd.DataFrame): Kestrel results DataFrame.
        advntr_df (pd.DataFrame): adVNTR results DataFrame.
        advntr_available (bool): Whether adVNTR results are available.
        mean_vntr_coverage (float): Mean coverage over the VNTR region.
        mean_vntr_cov_threshold (float): Coverage threshold for the VNTR region.
        report_config (dict): Report configuration containing algorithm logic and summary rules.

    Returns:
        str: A detailed summary text describing the findings or negative result.
    """
    summary_text = ""
    try:
        kestrel_logic = report_config.get("algorithm_logic", {}).get("kestrel", {})
        computed_kestrel = compute_algorithm_result(kestrel_df, kestrel_logic)
        logging.debug("Computed Kestrel result: %s", computed_kestrel)

        advntr_logic = report_config.get("algorithm_logic", {}).get("advntr", {})
        computed_advntr = (
            compute_algorithm_result(advntr_df, advntr_logic)
            if advntr_available
            else "none"
        )
        logging.debug("Computed adVNTR result: %s", computed_advntr)

        quality_metrics_pass = True
        if (
            mean_vntr_coverage is not None
            and mean_vntr_coverage < mean_vntr_cov_threshold
        ):
            quality_metrics_pass = False
        logging.debug("Quality metrics pass: %s", quality_metrics_pass)

        current_conditions = {
            "kestrel_result": computed_kestrel,
            "advntr_result": computed_advntr,
            "quality_metrics_pass": quality_metrics_pass,
        }
        logging.debug("Unified screening conditions: %s", current_conditions)

        unified_rules = report_config.get("screening_summary_rules", [])
        default_message = report_config.get(
            "screening_summary_default",
            "The screening was negative (no valid Kestrel or adVNTR data).",
        )

        def condition_matches(current, rule_value):
            if isinstance(rule_value, list):
                return current in rule_value
            return current == rule_value

        for rule in unified_rules:
            conditions = rule.get("conditions", {})
            if all(
                condition_matches(current_conditions.get(key), conditions.get(key))
                for key in conditions
            ):
                summary_text = rule.get("message", "")
                logging.debug("Unified rule matched: %s", conditions)
                break

        if not summary_text:
            summary_text = default_message
            logging.debug("No unified rule matched; using default screening message.")

    except Exception as ex:
        logging.error("Exception in build_screening_summary: %s", ex)
        summary_text = "No summary available."

    logging.debug("Final screening summary: %s", summary_text)
    return summary_text


def generate_summary_report(
    output_dir,
    template_dir,
    report_file,
    log_file,
    bed_file=None,
    bam_file=None,
    fasta_file=None,
    flanking=50,
    vcf_file=None,
    config=None,
):
    """
    Generates a summary report based on a pipeline summary JSON file.
    Instead of parsing results from subfolders, the report now loads the output summary file (JSON)
    and renders the report based on that content. The main configuration is passed in via 'config'.
    This module additionally loads its own report-specific configuration from report_config.json.
    The header information extracted from the BAM Header Parsing step (including warning,
    alignment_pipeline, assembly_text, and assembly_contig) is included in the report context.

    Args:
        output_dir (str): Output directory for the report.
        template_dir (str): Directory containing the report template.
        report_file (str): Name of the report file.
        log_file (str): Path to the pipeline log file.
        bed_file (str, optional): Path to the BED file for IGV reports.
        bam_file (str, optional): Path to the BAM file for IGV reports.
        fasta_file (str, optional): Path to the reference FASTA file for IGV reports.
        flanking (int, optional): Size of the flanking region for IGV reports.
        vcf_file (str, optional): Path to the sorted/indexed VCF file.
        config (dict): Main configuration dictionary (passed from the pipeline).

    Raises:
        ValueError: If config is not provided.
    """
    logging.debug("---- DEBUG: Entered generate_summary_report ----")
    logging.debug(
        "Called with output_dir=%s, template_dir=%s, report_file=%s",
        output_dir,
        template_dir,
        report_file,
    )
    logging.debug(
        "bed_file=%s, bam_file=%s, fasta_file=%s, flanking=%s, log_file=%s, vcf_file=%s",
        bed_file,
        bam_file,
        fasta_file,
        flanking,
        log_file,
        vcf_file,
    )

    if config is None:
        raise ValueError(
            "Config dictionary must be provided to generate_summary_report"
        )

    # Load the script-specific report configuration.
    report_config = load_report_config()

    # Resolve flanking region from main config if needed.
    if flanking == 50 and config is not None:
        flanking = config.get("default_values", {}).get("flanking", 50)
        logging.debug("Flanking region set to %s based on config.", flanking)

    thresholds = config.get("thresholds", {})
    mean_vntr_cov_threshold = thresholds.get("mean_vntr_coverage", 100)
    dup_rate_cutoff = thresholds.get("duplication_rate", 0.1)
    q20_rate_cutoff = thresholds.get("q20_rate", 0.8)
    q30_rate_cutoff = thresholds.get("q30_rate", 0.7)
    passed_filter_rate_cutoff = thresholds.get("passed_filter_reads_rate", 0.8)

    # Load the pipeline summary JSON.
    summary_file_path = Path(output_dir) / "pipeline_summary.json"
    pipeline_summary = load_pipeline_summary(summary_file_path)

    # Extract input_files and pipeline_version from the summary.
    input_files = pipeline_summary.get("input_files", {})
    pipeline_version = pipeline_summary.get("version", "unknown")

    # Extract header info from BAM Header Parsing step (if available)
    header_info = {}
    for step in pipeline_summary.get("steps", []):
        if step.get("step") == "BAM Header Parsing":
            header_info = step.get("parsed_result", {})
            break

    # Extract individual header elements robustly.
    header_warning = header_info.get("warning", "")
    alignment_pipeline = header_info.get("alignment_pipeline", "")
    assembly_text = header_info.get("assembly_text", "")
    assembly_contig = header_info.get("assembly_contig", "")

    # Extract mean VNTR coverage from the "Coverage Calculation" step.
    mean_vntr_coverage = None
    for step in pipeline_summary.get("steps", []):
        if step.get("step") == "Coverage Calculation":
            coverage_info = step.get("parsed_result", {}).get("data", [])
            if (
                coverage_info
                and isinstance(coverage_info, list)
                and len(coverage_info) > 0
            ):
                try:
                    mean_vntr_coverage = float(coverage_info[0].get("mean", 0))
                    logging.debug(
                        "Mean VNTR coverage extracted: %s", mean_vntr_coverage
                    )
                except Exception as e:
                    logging.error("Error parsing mean VNTR coverage: %s", e)
                    mean_vntr_coverage = None
            break

    # Extract Kestrel and adVNTR data from the summary.
    kestrel_data = []
    advntr_data = []
    advntr_available = False
    for step in pipeline_summary.get("steps", []):
        if step.get("step") == "Kestrel Genotyping":
            kestrel_data = step.get("parsed_result", {}).get("data", [])
        elif step.get("step") == "adVNTR Genotyping":
            advntr_data = step.get("parsed_result", {}).get("data", [])
            advntr_available = True

    if kestrel_data:
        kestrel_df = pd.DataFrame(kestrel_data)
        columns_to_display = {
            "Motifs": "Motif",
            "Variant": "Variant",
            "POS": "Position",
            "REF": "REF",
            "ALT": "ALT",
            "Motif_sequence": "Motif Sequence",
            "Estimated_Depth_AlternateVariant": "Depth (Variant)",
            "Estimated_Depth_Variant_ActiveRegion": "Depth (Region)",
            "Depth_Score": "Depth Score",
            "Confidence": "Confidence",
            "Flag": "Flag",
        }
        existing_cols = [col for col in columns_to_display if col in kestrel_df.columns]
        kestrel_df = kestrel_df[existing_cols]
        kestrel_df = kestrel_df.rename(
            columns={col: columns_to_display[col] for col in existing_cols}
        )
        # Sort the DataFrame by Depth Score in descending order if available.
        if "Depth Score" in kestrel_df.columns:
            try:
                kestrel_df["Depth Score"] = pd.to_numeric(
                    kestrel_df["Depth Score"], errors="coerce"
                )
            except Exception as e:
                logging.warning("Could not convert 'Depth Score' to numeric: %s", e)
            kestrel_df = kestrel_df.sort_values(by="Depth Score", ascending=False)
        # Create a copy of the sorted dataframe (without HTML formatting) for matching.
        kestrel_df_raw = kestrel_df.copy()
        # Now apply color-coding to the Confidence column for display.
        if "Confidence" in kestrel_df.columns:
            kestrel_df["Confidence"] = kestrel_df["Confidence"].apply(
                lambda x: (
                    f'<span style="color:orange;font-weight:bold;">{x}</span>'
                    if x == "Low_Precision"
                    else (
                        f'<span style="color:red;font-weight:bold;">{x}</span>'
                        if x in ["High_Precision", "High_Precision*"]
                        else x
                    )
                )
            )
        logging.debug("Kestrel data extracted from summary and formatted.")
    else:
        kestrel_df = pd.DataFrame()
        kestrel_df_raw = pd.DataFrame()
        logging.warning("No Kestrel data found in pipeline summary.")

    if advntr_data:
        advntr_df = pd.DataFrame(advntr_data)
        advntr_columns = [
            "VID",
            "Variant",
            "NumberOfSupportingReads",
            "MeanCoverage",
            "Pvalue",
            "RU",
            "POS",
            "REF",
            "ALT",
            "Flag",
        ]
        advntr_df = advntr_df[
            [col for col in advntr_columns if col in advntr_df.columns]
        ]
        logging.debug("adVNTR data extracted from summary and formatted.")
    else:
        advntr_df = pd.DataFrame()

    log_content = load_pipeline_log(log_file)

    # IGV report generation (if applicable)
    if bed_file and os.path.exists(bed_file):
        logging.info("Running IGV report for BED file: %s", bed_file)
        igv_report_file = Path(output_dir) / "igv_report.html"
        run_igv_report(
            bed_file,
            bam_file,
            fasta_file,
            igv_report_file,
            flanking=flanking,
            vcf_file=vcf_file,
            config=config,
        )
    else:
        logging.warning(
            "BED file does not exist or not provided. Skipping IGV report generation."
        )
        igv_report_file = None

    if igv_report_file and igv_report_file.exists():
        igv_content, table_json, session_dictionary = extract_igv_content(
            igv_report_file
        )
    else:
        logging.warning("IGV report file not found. Skipping IGV content.")
        igv_content, table_json, session_dictionary = "", "", ""

    fastp_file = Path(output_dir) / "fastq_bam_processing/output.json"
    fastp_data = load_fastp_output(fastp_file)

    if mean_vntr_coverage is not None and mean_vntr_coverage < mean_vntr_cov_threshold:
        coverage_icon = '<span style="color:red;font-weight:bold;">&#9888;</span>'
        coverage_color = "red"
        logging.debug("Mean VNTR coverage is below the threshold.")
    else:
        coverage_icon = '<span style="color:green;font-weight:bold;">&#10004;</span>'
        coverage_color = "green"
        logging.debug("Mean VNTR coverage is above the threshold.")

    duplication_rate = None
    q20_rate = None
    q30_rate = None
    passed_filter_rate = None
    sequencing_str = ""
    fastp_available = False
    if fastp_data:
        fastp_available = True
        summary_fastp = fastp_data.get("summary", {})
        duplication = fastp_data.get("duplication", {})
        filtering_result = fastp_data.get("filtering_result", {})

        duplication_rate = duplication.get("rate", None)
        after_filtering = summary_fastp.get("after_filtering", {})
        before_filtering = summary_fastp.get("before_filtering", {})

        q20_rate = after_filtering.get("q20_rate", None)
        q30_rate = after_filtering.get("q30_rate", None)

        total_reads_before = before_filtering.get("total_reads", 1)
        passed_filter_reads = filtering_result.get("passed_filter_reads", 0)
        if total_reads_before > 0:
            passed_filter_rate = passed_filter_reads / total_reads_before
            logging.debug("Passed filter rate calculated: %.2f", passed_filter_rate)
        else:
            passed_filter_rate = None
            logging.debug(
                "Total reads before filtering is zero; passed filter rate set to None."
            )
        sequencing_str = summary_fastp.get("sequencing", "")
        logging.debug("Sequencing setup: %s", sequencing_str)

    def warn_icon(value, cutoff, higher_better=True):
        if value is None:
            logging.debug("warn_icon called with value=None; returning empty strings.")
            return "", ""
        if higher_better:
            if value < cutoff:
                logging.debug(
                    "Value %s is below the cutoff %s (higher_better=True).",
                    value,
                    cutoff,
                )
                return '<span style="color:red;font-weight:bold;">&#9888;</span>', "red"
            else:
                logging.debug(
                    "Value %s is above or equal to the cutoff %s (higher_better=True).",
                    value,
                    cutoff,
                )
                return (
                    '<span style="color:green;font-weight:bold;">&#10004;</span>',
                    "green",
                )
        else:
            if value > cutoff:
                logging.debug(
                    "Value %s is above the cutoff %s (higher_better=False).",
                    value,
                    cutoff,
                )
                return '<span style="color:red;font-weight:bold;">&#9888;</span>', "red"
            else:
                logging.debug(
                    "Value %s is below or equal to the cutoff %s (higher_better=False).",
                    value,
                    cutoff,
                )
                return (
                    '<span style="color:green;font-weight:bold;">&#10004;</span>',
                    "green",
                )

    dup_icon, dup_color = warn_icon(
        duplication_rate, dup_rate_cutoff, higher_better=False
    )
    q20_icon, q20_color = warn_icon(q20_rate, q20_rate_cutoff, higher_better=True)
    q30_icon, q30_color = warn_icon(q30_rate, q30_rate_cutoff, higher_better=True)
    pf_icon, pf_color = warn_icon(
        passed_filter_rate, passed_filter_rate_cutoff, higher_better=True
    )

    kestrel_html = kestrel_df.to_html(
        table_id="kestrel_table",
        classes="table table-bordered table-striped hover compact order-column table-sm",
        index=False,
        escape=False,
    )
    logging.debug("Kestrel results converted to HTML.")

    if advntr_available:
        if not advntr_df.empty:
            advntr_html = advntr_df.to_html(
                classes="table table-bordered table-striped hover compact table-sm",
                index=False,
            )
            logging.debug("adVNTR results converted to HTML.")
        else:
            advntr_html = "<p>No pathogenic variants identified by adVNTR.</p>"
            logging.debug(
                "adVNTR was performed but no variants identified; adding negative message."
            )
    else:
        advntr_html = "<p>adVNTR genotyping was not performed.</p>"
        logging.debug("adVNTR was not performed; adding message to report.")

    # Extract cross-match summary message from the pipeline summary if available.
    cross_match_message = ""
    for step in pipeline_summary.get("steps", []):
        if step.get("step") == "Cross-Match Variant Comparison":
            data = step.get("parsed_result", {}).get("data", [])
            if any(item.get("Match") == "Yes" for item in data):
                cross_match_message = (
                    "At least one match was found between Kestrel and adVNTR results."
                )
            else:
                cross_match_message = (
                    "No matches were found between Kestrel and adVNTR results."
                )
            break

    env = Environment(loader=FileSystemLoader(template_dir))
    try:
        template = env.get_template("report_template.html")
        logging.debug("Jinja2 template 'report_template.html' loaded successfully.")
    except Exception as e:
        logging.error("Failed to load Jinja2 template: %s", e)
        raise

    # Use the sorted (and raw) kestrel dataframe for matching.
    summary_text = build_screening_summary(
        kestrel_df_raw,
        advntr_df,
        advntr_available,
        mean_vntr_coverage,
        mean_vntr_cov_threshold,
        report_config,
    )
    logging.debug("Summary text generated: %s", summary_text)

    context = {
        "kestrel_highlight": kestrel_html,
        "advntr_highlight": advntr_html,
        "advntr_available": advntr_available,
        "log_content": load_pipeline_log(log_file),
        "igv_content": igv_content,
        "table_json": table_json,
        "session_dictionary": session_dictionary,
        "report_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "input_files": input_files,
        "pipeline_version": pipeline_version,
        "header_warning": header_warning,
        "alignment_pipeline": alignment_pipeline,
        "assembly_text": assembly_text,
        "assembly_contig": assembly_contig,
        "mean_vntr_coverage": (
            mean_vntr_coverage if mean_vntr_coverage is not None else "Not calculated"
        ),
        "mean_vntr_coverage_icon": coverage_icon,
        "mean_vntr_coverage_color": coverage_color,
        "fastp_available": fastp_available,
        "duplication_rate": duplication_rate,
        "duplication_rate_icon": dup_icon,
        "duplication_rate_color": dup_color,
        "q20_rate": q20_rate,
        "q20_icon": q20_icon,
        "q20_color": q20_color,
        "q30_rate": q30_rate,
        "q30_icon": q30_icon,
        "q30_color": q30_color,
        "passed_filter_rate": passed_filter_rate,
        "passed_filter_icon": pf_icon,
        "passed_filter_color": pf_color,
        "sequencing_str": sequencing_str,
        "summary_text": summary_text,
        "cross_match_message": cross_match_message,  # New variable for cross-match summary
    }

    try:
        rendered_html = template.render(context)
        logging.debug("Report template rendered successfully.")
    except Exception as e:
        logging.error("Failed to render the report template: %s", e)
        raise

    report_file_path = Path(output_dir) / report_file
    try:
        with open(report_file_path, "w") as f:
            f.write(rendered_html)
        logging.info("Summary report generated and saved to %s", report_file_path)
    except Exception as e:
        logging.error("Failed to write the summary report: %s", e)
        raise
