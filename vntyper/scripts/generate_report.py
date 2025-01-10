#!/usr/bin/env python3
# vntyper/scripts/generate_report.py

import os
import logging
import subprocess
import json
from datetime import datetime
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader


def load_kestrel_results(kestrel_result_file):
    logging.info(f"Loading Kestrel results from {kestrel_result_file}")
    if not os.path.exists(kestrel_result_file):
        logging.warning(f"Kestrel result file not found: {kestrel_result_file}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(kestrel_result_file, sep='\t', comment='#')
        columns_to_display = {
            'Motif': 'Motif',
            'Variant': 'Variant',
            'POS': 'Position',
            'REF': 'REF',
            'ALT': 'ALT',
            'Motif_sequence': 'Motif Sequence',
            'Estimated_Depth_AlternateVariant': 'Depth (Variant)',
            'Estimated_Depth_Variant_ActiveRegion': 'Depth (Region)',
            'Depth_Score': 'Depth Score',
            'Confidence': 'Confidence'
        }
        df = df[list(columns_to_display.keys())]
        df = df.rename(columns=columns_to_display)
        df['Confidence'] = df['Confidence'].apply(
            lambda x: (
                f'<span style="color:orange;font-weight:bold;">{x}</span>'
                if x == 'Low_Precision'
                else f'<span style="color:red;font-weight:bold;">{x}</span>'
                if x == 'High_Precision'
                else x
            )
        )
        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame()
    except Exception as e:
        logging.error(f"Unexpected error loading Kestrel results: {e}")
        return pd.DataFrame()


def load_advntr_results(advntr_result_file):
    logging.info(f"Loading adVNTR results from {advntr_result_file}")
    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame(), False

    try:
        df = pd.read_csv(advntr_result_file, sep='\t', comment='#')
        return df, True
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame(), False
    except Exception as e:
        logging.error(f"Unexpected error loading adVNTR results: {e}")
        return pd.DataFrame(), False


def run_igv_report(bed_file, bam_file, fasta_file, output_html, flanking=50, vcf_file=None, config=None):
    """
    Wrapper around `create_report` IGV command. If config is provided and flanking
    is not explicitly set, we fallback to config's default_values.flanking.
    We skip passing None for track arguments (vcf_file or bam_file).
    """
    logging.debug("run_igv_report called with:")
    logging.debug(f"  bed_file={bed_file}")
    logging.debug(f"  bam_file={bam_file}")
    logging.debug(f"  fasta_file={fasta_file}")
    logging.debug(f"  output_html={output_html}")
    logging.debug(f"  vcf_file={vcf_file}")
    logging.debug(f"  flanking={flanking}")

    if config is not None and flanking == 50:
        flanking = config.get("default_values", {}).get("flanking", 50)

    # Convert each path or None into string or skip
    bed_file = str(bed_file) if bed_file else None
    bam_file = str(bam_file) if bam_file else None
    fasta_file = str(fasta_file) if fasta_file else None
    output_html = str(output_html) if output_html else None

    # We'll build the IGV command piece by piece, skipping None tracks
    igv_report_cmd = [
        'create_report',
        bed_file,               # The region/bed
        '--flanking', str(flanking),
        '--fasta', fasta_file,
        '--tracks'
    ]

    # Collect tracks in a list, skipping any that are None
    tracks = []
    if vcf_file:
        tracks.append(str(vcf_file))
    if bam_file:
        tracks.append(str(bam_file))

    # If you want to enforce at least one track, you can check here:
    if not tracks:
        logging.warning("No valid tracks (VCF or BAM) provided to IGV. The IGV report may be empty.")

    igv_report_cmd.extend(tracks)

    igv_report_cmd.extend([
        '--output', output_html
    ])

    # Now run the command
    try:
        logging.info(f"Running IGV report: {' '.join([str(x) for x in igv_report_cmd if x])}")
        subprocess.run(igv_report_cmd, check=True)
        logging.info(f"IGV report successfully generated at {output_html}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error generating IGV report: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error generating IGV report: {e}")
        raise


def extract_igv_content(igv_report_html):
    logging.debug(f"extract_igv_content called with igv_report_html={igv_report_html}")
    try:
        with open(igv_report_html, 'r') as f:
            content = f.read()

        igv_start = content.find('<div id="container"')
        igv_end = content.find('</body>')

        if igv_start == -1 or igv_end == -1:
            logging.error("Failed to extract IGV content from report.")
            return "", "", ""

        igv_content = content[igv_start:igv_end].strip()

        table_json_start = content.find('const tableJson = ') + len('const tableJson = ')
        table_json_end = content.find('\n', table_json_start)
        table_json = content[table_json_start:table_json_end].strip()

        session_dict_start = content.find('const sessionDictionary = ') + len('const sessionDictionary = ')
        session_dict_end = content.find('\n', session_dict_start)
        session_dictionary = content[session_dict_start:session_dict_end].strip()

        logging.info("Successfully extracted IGV content, tableJson, and sessionDictionary.")
        return igv_content, table_json, session_dictionary
    except FileNotFoundError:
        logging.error(f"IGV report file not found: {igv_report_html}")
        return "", "", ""
    except Exception as e:
        logging.error(f"Unexpected error extracting IGV content: {e}")
        return "", "", ""


def load_fastp_output(fastp_file):
    logging.debug(f"load_fastp_output called with fastp_file={fastp_file}")
    if not os.path.exists(fastp_file):
        logging.warning(f"fastp output file not found: {fastp_file}")
        return {}

    try:
        with open(fastp_file, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        logging.error(f"Failed to load or parse fastp output: {e}")
        return {}


def generate_summary_report(
    output_dir,
    template_dir,
    report_file,
    bed_file=None,
    bam_file=None,
    fasta_file=None,
    flanking=50,
    input_files=None,
    pipeline_version=None,
    mean_vntr_coverage=None,
    vcf_file=None,
    config=None
):
    """
    Generates a summary report.

    Args:
        output_dir (Path): Output directory for the report.
        template_dir (str): Directory containing the report template.
        report_file (str): Name of the report file.
        bed_file (Path, optional): Path to the BED file for IGV reports.
        bam_file (Path, optional): Path to the BAM file for IGV reports.
        fasta_file (Path, optional): Path to the reference FASTA file for IGV reports.
        flanking (int, optional): Size of the flanking region for IGV reports.
        input_files (dict, optional): Dictionary of input filenames.
        pipeline_version (str, optional): The version of the VNtyper pipeline.
        mean_vntr_coverage (float, optional): Mean coverage over the VNTR region.
        vcf_file (Path, optional): Path to the sorted and indexed VCF file.
        config (dict, optional): Configuration dictionary.

    Raises:
        ValueError: If config is not provided.
    """
    logging.debug("---- DEBUG: Entered generate_summary_report ----")
    logging.debug(f"Called with output_dir={output_dir}, template_dir={template_dir}, report_file={report_file}")
    logging.debug(f"bed_file={bed_file}, bam_file={bam_file}, fasta_file={fasta_file}, flanking={flanking}")
    logging.debug(f"vcf_file={vcf_file}, mean_vntr_coverage={mean_vntr_coverage}")

    if config is None:
        raise ValueError("Config dictionary must be provided to generate_summary_report")

    # Debug checks for bed_file existence
    if bed_file:
        abs_bed_file = os.path.abspath(bed_file)
        logging.debug(f"Absolute bed_file => {abs_bed_file}")
        logging.debug(f"Exists? => {os.path.exists(abs_bed_file)}")

    if flanking == 50 and config is not None:
        flanking = config.get("default_values", {}).get("flanking", 50)

    thresholds = config.get("thresholds", {})
    mean_vntr_cov_threshold = thresholds.get("mean_vntr_coverage", 100)
    dup_rate_cutoff = thresholds.get("duplication_rate", 0.1)
    q20_rate_cutoff = thresholds.get("q20_rate", 0.8)
    q30_rate_cutoff = thresholds.get("q30_rate", 0.7)
    passed_filter_rate_cutoff = thresholds.get("passed_filter_reads_rate", 0.8)

    kestrel_result_file = Path(output_dir) / "kestrel/kestrel_result.tsv"
    advntr_result_file = Path(output_dir) / "advntr/output_adVNTR_result.tsv"
    igv_report_file = Path(output_dir) / "igv_report.html"
    fastp_file = Path(output_dir) / "fastq_bam_processing/output.json"

    # Debug checks
    logging.debug(f"kestrel_result_file => {kestrel_result_file}, exists? {kestrel_result_file.exists()}")
    logging.debug(f"advntr_result_file => {advntr_result_file}, exists? {advntr_result_file.exists()}")
    logging.debug(f"igv_report_file => {igv_report_file}, exists? {igv_report_file.exists()}")
    logging.debug(f"fastp_file => {fastp_file}, exists? {fastp_file.exists()}")

    # Attempt IGV generation only if bed_file is present & real
    if bed_file and os.path.exists(bed_file):
        logging.info(f"Running IGV report for BED file: {bed_file}")
        run_igv_report(
            bed_file,
            bam_file,
            fasta_file,
            igv_report_file,
            flanking=flanking,
            vcf_file=vcf_file,
            config=config
        )
    else:
        logging.warning("BED file does not exist or not provided. Skipping IGV report generation.")
        igv_report_file = None

    kestrel_df = load_kestrel_results(kestrel_result_file)
    advntr_df, advntr_available = load_advntr_results(advntr_result_file)

    # Removed log_content since log_file is no longer handled here
    igv_content, table_json, session_dictionary = ("", "", "")

    if igv_report_file and igv_report_file.exists():
        igv_content, table_json, session_dictionary = extract_igv_content(igv_report_file)
    else:
        logging.warning("IGV report file not found. Skipping IGV content.")

    fastp_data = load_fastp_output(fastp_file)

    if mean_vntr_coverage is not None and mean_vntr_coverage < mean_vntr_cov_threshold:
        coverage_icon = '<span style="color:red;font-weight:bold;">&#9888;</span>'
        coverage_color = 'red'
    else:
        coverage_icon = '<span style="color:green;font-weight:bold;">&#10004;</span>'
        coverage_color = 'green'

    duplication_rate = None
    q20_rate = None
    q30_rate = None
    passed_filter_rate = None
    sequencing_str = ""
    fastp_available = False
    if fastp_data:
        fastp_available = True
        summary = fastp_data.get("summary", {})
        duplication = fastp_data.get("duplication", {})
        filtering_result = fastp_data.get("filtering_result", {})

        duplication_rate = duplication.get("rate", None)
        after_filtering = summary.get("after_filtering", {})
        before_filtering = summary.get("before_filtering", {})
        q20_rate = after_filtering.get("q20_rate", None)
        q30_rate = after_filtering.get("q30_rate", None)

        total_reads_before = before_filtering.get("total_reads", 1)
        passed_filter_reads = filtering_result.get("passed_filter_reads", 0)
        if total_reads_before > 0:
            passed_filter_rate = passed_filter_reads / total_reads_before
        else:
            passed_filter_rate = None
        sequencing_str = summary.get("sequencing", "")

    def warn_icon(value, cutoff, higher_better=True):
        if value is None:
            return "", ""
        if higher_better:
            return (
                ('<span style="color:red;font-weight:bold;">&#9888;</span>', 'red') if value < cutoff
                else ('<span style="color:green;font-weight:bold;">&#10004;</span>', 'green')
            )
        else:
            return (
                ('<span style="color:red;font-weight:bold;">&#9888;</span>', 'red') if value > cutoff
                else ('<span style="color:green;font-weight:bold;">&#10004;</span>', 'green')
            )

    dup_icon, dup_color = warn_icon(duplication_rate, dup_rate_cutoff, higher_better=False)
    q20_icon, q20_color = warn_icon(q20_rate, q20_rate_cutoff, higher_better=True)
    q30_icon, q30_color = warn_icon(q30_rate, q30_rate_cutoff, higher_better=True)
    pf_icon, pf_color = warn_icon(passed_filter_rate, passed_filter_rate_cutoff, higher_better=True)

    kestrel_html = kestrel_df.to_html(
        classes='table table-bordered table-striped hover compact order-column table-sm',
        index=False,
        escape=False
    )
    if advntr_available and not advntr_df.empty:
        advntr_html = advntr_df.to_html(
            classes='table table-bordered table-striped hover compact order-column table-sm',
            index=False
        )
    else:
        advntr_html = None

    env = Environment(loader=FileSystemLoader(template_dir))
    try:
        template = env.get_template('report_template.html')
    except Exception as e:
        logging.error(f"Failed to load Jinja2 template: {e}")
        raise

    summary_text = build_screening_summary(
        kestrel_df, advntr_df, advntr_available, mean_vntr_coverage, mean_vntr_cov_threshold
    )
    context = {
        'kestrel_highlight': kestrel_html,
        'advntr_highlight': advntr_html,
        'advntr_available': advntr_available,
        # 'log_content': log_content,  # Removed
        'igv_content': igv_content,
        'table_json': table_json,
        'session_dictionary': session_dictionary,
        'report_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'input_files': input_files or {},
        'pipeline_version': pipeline_version or "unknown",
        'mean_vntr_coverage': (
            mean_vntr_coverage if mean_vntr_coverage is not None else "Not calculated"
        ),
        'mean_vntr_coverage_icon': coverage_icon,
        'mean_vntr_coverage_color': coverage_color,
        'fastp_available': fastp_available,
        'duplication_rate': duplication_rate,
        'duplication_rate_icon': dup_icon,
        'duplication_rate_color': dup_color,
        'q20_rate': q20_rate,
        'q20_icon': q20_icon,
        'q20_color': q20_color,
        'q30_rate': q30_rate,
        'q30_icon': q30_icon,
        'q30_color': q30_color,
        'passed_filter_rate': passed_filter_rate,
        'passed_filter_icon': pf_icon,
        'passed_filter_color': pf_color,
        'sequencing_str': sequencing_str,
        # 'summary_text': summary_text  # Assuming it's included in the template
    }

    try:
        rendered_html = template.render(context)
    except Exception as e:
        logging.error(f"Failed to render the report template: {e}")
        raise

    report_file_path = Path(output_dir) / report_file
    try:
        with open(report_file_path, 'w') as f:
            f.write(rendered_html)
        logging.info(f"Summary report generated and saved to {report_file_path}")
    except Exception as e:
        logging.error(f"Failed to write the summary report: {e}")
        raise


def build_screening_summary(kestrel_df, advntr_df, advntr_available, mean_vntr_coverage, mean_vntr_cov_threshold):
    """
    Build the short screening summary text based on Kestrel and adVNTR data.

    Args:
        kestrel_df (pd.DataFrame): Kestrel results DataFrame.
        advntr_df (pd.DataFrame): adVNTR results DataFrame.
        advntr_available (bool): Whether adVNTR results are available.
        mean_vntr_coverage (float): Mean coverage over the VNTR region.
        mean_vntr_cov_threshold (float): Coverage threshold for the VNTR region.

    Returns:
        str: A short summary text describing the findings or negative result.
    """
    summary_text = ""
    try:
        import re

        def strip_html_tags(confidence_value):
            return re.sub(r"<[^>]*>", "", confidence_value or "")

        # Check if there's any recognized Confidence
        if not kestrel_df.empty and "Confidence" in kestrel_df.columns:
            kestrel_confidences = kestrel_df["Confidence"].apply(strip_html_tags).dropna().unique().tolist()
            kestrel_valid = any(conf in ("High_Precision", "Low_Precision") for conf in kestrel_confidences)
        else:
            kestrel_valid = False

        advntr_has_data = advntr_available and not advntr_df.empty

        logging.debug(f"advntr_available => {advntr_available}")
        logging.debug(f"advntr_df.empty => {advntr_df.empty}")
        logging.debug(f"advntr_has_data => {advntr_has_data}")
        logging.debug(f"adVNTR columns => {list(advntr_df.columns)}")
        logging.debug(f"adVNTR first row => {advntr_df.head(1).to_dict(orient='records')}")

        # Check coverage status
        coverage_status = (
            "above" if (mean_vntr_coverage is not None and mean_vntr_coverage >= mean_vntr_cov_threshold)
            else "below"
        )

        # Check if adVNTR p-value is available
        if advntr_has_data:
            lower_cols = [c.lower() for c in advntr_df.columns]
            logging.debug(f"Lowercased adVNTR columns => {lower_cols}")
            if "p-value" in lower_cols:
                pval_col = [col for col in advntr_df.columns if col.lower() == "p-value"][0]
                pval_val = advntr_df[pval_col].iloc[0]
                pval_msg = f" adVNTR p-value: {pval_val}."
            elif "pvalue" in lower_cols:
                pval_col = [col for col in advntr_df.columns if col.lower() == "pvalue"][0]
                pval_val = advntr_df[pval_col].iloc[0]
                pval_msg = f" adVNTR p-value: {pval_val}."
            else:
                logging.debug("No p-value or pvalue column detected.")
                pval_msg = ""
        else:
            pval_msg = ""

        # Summarize
        if kestrel_valid or advntr_has_data:
            summary_text = (
                f"{'Kestrel detected a variant.' if kestrel_valid else ''}"
                f"{pval_msg} Coverage is {coverage_status} threshold."
            )
        else:
            summary_text = "The screening was negative (no valid Kestrel or adVNTR data)."

    except Exception as ex:
        logging.error(f"Exception in build_screening_summary: {ex}")
        summary_text = "No summary available."

    return summary_text
