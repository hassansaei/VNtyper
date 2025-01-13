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


def load_kestrel_results(kestrel_result_file):
    """
    Loads Kestrel results (kestrel_result.tsv) and adjusts the 'Confidence' column
    to color-code values in HTML. Now also handles 'High_Precision*' by applying
    the same styling as 'High_Precision'.
    """
    logging.info(f"Loading Kestrel results from {kestrel_result_file}")
    if not os.path.exists(kestrel_result_file):
        logging.warning(f"Kestrel result file not found: {kestrel_result_file}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(kestrel_result_file, sep='\t', comment='#')
        logging.debug(f"Kestrel DataFrame loaded with {len(df)} rows.")

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
        logging.debug("Kestrel DataFrame columns renamed for display.")

        # Handle Low_Precision in orange, High_Precision & High_Precision* in red
        df['Confidence'] = df['Confidence'].apply(
            lambda x: (
                f'<span style="color:orange;font-weight:bold;">{x}</span>'
                if x == 'Low_Precision'
                else f'<span style="color:red;font-weight:bold;">{x}</span>'
                if x in ['High_Precision', 'High_Precision*']
                else x
            )
        )
        logging.debug("Kestrel 'Confidence' column color-coded based on precision levels.")
        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame()
    except Exception as e:
        logging.error(f"Unexpected error loading Kestrel results: {e}")
        return pd.DataFrame()


def load_advntr_results(advntr_result_file):
    """
    Loads adVNTR results (output_adVNTR_result.tsv) if present. Returns
    an empty DataFrame and a False flag if not found or parsing fails.
    """
    logging.info(f"Loading adVNTR results from {advntr_result_file}")
    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame(), False

    try:
        df = pd.read_csv(advntr_result_file, sep='\t', comment='#')
        logging.debug(f"adVNTR DataFrame loaded with {len(df)} rows.")
        return df, True
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame(), False
    except Exception as e:
        logging.error(f"Unexpected error loading adVNTR results: {e}")
        return pd.DataFrame(), False()


def load_pipeline_log(log_file):
    """
    Loads the pipeline log content from the specified log_file.
    Returns a placeholder string if not found or on error.
    """
    logging.info(f"Loading pipeline log from {log_file}")
    # If log_file is None or empty, just return a placeholder.
    if not log_file:
        logging.warning("No pipeline log file provided; skipping log loading.")
        return "No pipeline log file was provided."

    # If log_file is not None, check if it exists:
    if not os.path.exists(log_file):
        logging.warning(f"Pipeline log file not found: {log_file}")
        return "Pipeline log file not found."

    try:
        with open(log_file, 'r') as f:
            content = f.read()
        logging.debug("Pipeline log successfully loaded.")
        return content
    except Exception as e:
        logging.error(f"Failed to read pipeline log file: {e}")
        return "Failed to load pipeline log."


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
        logging.debug(f"Flanking region set to {flanking} based on config.")

    # Convert each path or None into string or skip
    bed_file = str(bed_file) if bed_file else None
    bam_file = str(bam_file) if bam_file else None
    fasta_file = str(fasta_file) if fasta_file else None
    output_html = str(output_html) if output_html else None

    igv_report_cmd = [
        'create_report',
        bed_file,
        '--flanking', str(flanking),
        '--fasta', fasta_file,
        '--tracks'
    ]

    tracks = []
    if vcf_file:
        tracks.append(str(vcf_file))
    if bam_file:
        tracks.append(str(bam_file))
    if not tracks:
        logging.warning("No valid tracks (VCF or BAM) provided to IGV. The IGV report may be empty.")

    igv_report_cmd.extend(tracks)
    igv_report_cmd.extend(['--output', output_html])

    logging.debug(f"IGV report command: {' '.join([str(x) for x in igv_report_cmd if x])}")

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
    """
    Reads the generated IGV HTML report and extracts the IGV content,
    the tableJson variable, and the sessionDictionary variable from the script.
    Returns empty strings if not found or on error.
    """
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
    """
    Loads fastp JSON output (e.g., output.json) for summary metrics if available.
    Returns an empty dict if file not found or if parsing fails.
    """
    logging.debug(f"load_fastp_output called with fastp_file={fastp_file}")
    if not os.path.exists(fastp_file):
        logging.warning(f"fastp output file not found: {fastp_file}")
        return {}

    try:
        with open(fastp_file, 'r') as f:
            data = json.load(f)
        logging.debug("fastp output successfully loaded.")
        return data
    except Exception as e:
        logging.error(f"Failed to load or parse fastp output: {e}")
        return {}


def generate_summary_report(
    output_dir,
    template_dir,
    report_file,
    log_file,
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
        output_dir (str): Output directory for the report.
        template_dir (str): Directory containing the report template.
        report_file (str): Name of the report file.
        log_file (str): Path to the pipeline log file.
        bed_file (str, optional): Path to the BED file for IGV reports.
        bam_file (str, optional): Path to the BAM file for IGV reports.
        fasta_file (str, optional): Path to the reference FASTA file for IGV reports.
        flanking (int, optional): Size of the flanking region for IGV reports.
        input_files (dict, optional): Dictionary of input filenames.
        pipeline_version (str, optional): The version of the VNtyper pipeline.
        mean_vntr_coverage (float, optional): Mean coverage over the VNTR region.
        vcf_file (str, optional): Path to the sorted and indexed VCF file.
        config (dict, optional): Configuration dictionary.

    Raises:
        ValueError: If config is not provided.
    """
    logging.debug("---- DEBUG: Entered generate_summary_report ----")
    logging.debug(f"Called with output_dir={output_dir}, template_dir={template_dir}, report_file={report_file}")
    logging.debug(f"bed_file={bed_file}, bam_file={bam_file}, fasta_file={fasta_file}, flanking={flanking}")
    logging.debug(f"log_file={log_file}, vcf_file={vcf_file}, mean_vntr_coverage={mean_vntr_coverage}")

    if config is None:
        raise ValueError("Config dictionary must be provided to generate_summary_report")

    # Debug checks for bed_file existence
    if bed_file:
        abs_bed_file = os.path.abspath(bed_file)
        logging.debug(f"Absolute bed_file => {abs_bed_file}")
        logging.debug(f"Exists? => {os.path.exists(abs_bed_file)}")

    # Debug checks for log_file existence
    if log_file:
        abs_log_file = os.path.abspath(log_file)
        logging.debug(f"Absolute log_file => {abs_log_file}")
        logging.debug(f"Exists? => {os.path.exists(abs_log_file)}")

    if flanking == 50 and config is not None:
        flanking = config.get("default_values", {}).get("flanking", 50)
        logging.debug(f"Flanking region set to {flanking} based on config.")

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
    log_content = load_pipeline_log(log_file)

    # If we did produce an IGV report, extract the content
    if igv_report_file and igv_report_file.exists():
        igv_content, table_json, session_dictionary = extract_igv_content(igv_report_file)
    else:
        logging.warning("IGV report file not found. Skipping IGV content.")
        igv_content, table_json, session_dictionary = "", "", ""

    fastp_data = load_fastp_output(fastp_file)

    # Coverage status
    if mean_vntr_coverage is not None and mean_vntr_coverage < mean_vntr_cov_threshold:
        coverage_icon = '<span style="color:red;font-weight:bold;">&#9888;</span>'
        coverage_color = 'red'
        logging.debug("Mean VNTR coverage is below the threshold.")
    else:
        coverage_icon = '<span style="color:green;font-weight:bold;">&#10004;</span>'
        coverage_color = 'green'
        logging.debug("Mean VNTR coverage is above the threshold.")

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
            logging.debug(f"Passed filter rate calculated: {passed_filter_rate:.2f}")
        else:
            passed_filter_rate = None
            logging.debug("Total reads before filtering is zero; passed filter rate set to None.")
        sequencing_str = summary.get("sequencing", "")
        logging.debug(f"Sequencing setup: {sequencing_str}")

    def warn_icon(value, cutoff, higher_better=True):
        """
        Returns an HTML icon/string and a color ('red' or 'green'), depending
        on whether 'value' passes 'cutoff' under the 'higher_better' logic.
        """
        if value is None:
            logging.debug("warn_icon called with value=None; returning empty strings.")
            return "", ""
        if higher_better:
            if value < cutoff:
                logging.debug(f"Value {value} is below the cutoff {cutoff} (higher_better=True).")
                return '<span style="color:red;font-weight:bold;">&#9888;</span>', 'red'
            else:
                logging.debug(f"Value {value} is above or equal to the cutoff {cutoff} (higher_better=True).")
                return '<span style="color:green;font-weight:bold;">&#10004;</span>', 'green'
        else:
            if value > cutoff:
                logging.debug(f"Value {value} is above the cutoff {cutoff} (higher_better=False).")
                return '<span style="color:red;font-weight:bold;">&#9888;</span>', 'red'
            else:
                logging.debug(f"Value {value} is below or equal to the cutoff {cutoff} (higher_better=False).")
                return '<span style="color:green;font-weight:bold;">&#10004;</span>', 'green'

    # Evaluate duplication, Q20, Q30, and passed filter rates
    dup_icon, dup_color = warn_icon(duplication_rate, dup_rate_cutoff, higher_better=False)
    q20_icon, q20_color = warn_icon(q20_rate, q20_rate_cutoff, higher_better=True)
    q30_icon, q30_color = warn_icon(q30_rate, q30_rate_cutoff, higher_better=True)
    pf_icon, pf_color = warn_icon(passed_filter_rate, passed_filter_rate_cutoff, higher_better=True)

    # Convert Kestrel & adVNTR data to HTML
    kestrel_html = kestrel_df.to_html(
        classes='table table-bordered table-striped hover compact order-column table-sm',
        index=False,
        escape=False
    )
    logging.debug("Kestrel results converted to HTML.")

    if advntr_available and not advntr_df.empty:
        advntr_html = advntr_df.to_html(
            classes='table table-bordered table-striped hover compact order-column table-sm',
            index=False
        )
        logging.debug("adVNTR results converted to HTML.")
    else:
        advntr_html = None
        logging.debug("No adVNTR results available to convert to HTML.")

    # Prepare Jinja2 template
    env = Environment(loader=FileSystemLoader(template_dir))
    try:
        template = env.get_template('report_template.html')
        logging.debug("Jinja2 template 'report_template.html' loaded successfully.")
    except Exception as e:
        logging.error(f"Failed to load Jinja2 template: {e}")
        raise

    # Summarize the results
    summary_text = build_screening_summary(
        kestrel_df, advntr_df, advntr_available, mean_vntr_coverage, mean_vntr_cov_threshold
    )
    logging.debug(f"Summary text generated: {summary_text}")

    # Build the final context for rendering
    context = {
        'kestrel_highlight': kestrel_html,
        'advntr_highlight': advntr_html,
        'advntr_available': advntr_available,
        'log_content': log_content,  # Re-added pipeline log content
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
        'summary_text': summary_text
    }

    # Render the template
    try:
        rendered_html = template.render(context)
        logging.debug("Report template rendered successfully.")
    except Exception as e:
        logging.error(f"Failed to render the report template: {e}")
        raise

    # Write out the final HTML report
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
    Build the detailed screening summary text based on Kestrel and adVNTR data.

    Args:
        kestrel_df (pd.DataFrame): Kestrel results DataFrame.
        advntr_df (pd.DataFrame): adVNTR results DataFrame.
        advntr_available (bool): Whether adVNTR results are available.
        mean_vntr_coverage (float): Mean coverage over the VNTR region.
        mean_vntr_cov_threshold (float): Coverage threshold for the VNTR region.

    Returns:
        str: A detailed summary text describing the findings or negative result.
    """
    summary_text = ""
    try:
        # Function to strip HTML tags from the confidence values
        def strip_html_tags(confidence_value):
            return re.sub(r"<[^>]*>", "", confidence_value or "")

        # Extract unique confidence values from Kestrel results
        kestrel_confidences = []
        if not kestrel_df.empty and "Confidence" in kestrel_df.columns:
            kestrel_confidences = kestrel_df["Confidence"].apply(strip_html_tags).dropna().unique().tolist()
            logging.debug(f"Kestrel confidences extracted: {kestrel_confidences}")

        # Determine if Kestrel identified any pathogenic variants
        pathogenic_kestrel = any(conf in ("High_Precision", "High_Precision*", "Low_Precision")
                                 for conf in kestrel_confidences)
        logging.debug(f"Pathogenic variants identified by Kestrel: {pathogenic_kestrel}")

        # Determine confidence level
        confidence_level = None
        if "High_Precision" in kestrel_confidences or "High_Precision*" in kestrel_confidences:
            confidence_level = "High_Precision"
            logging.debug("Confidence level determined: High_Precision")
        elif "Low_Precision" in kestrel_confidences:
            confidence_level = "Low_Precision"
            logging.debug("Confidence level determined: Low_Precision")

        # Assess quality metrics
        quality_metrics_pass = True
        if mean_vntr_coverage is not None and mean_vntr_coverage < mean_vntr_cov_threshold:
            quality_metrics_pass = False
            logging.debug("Quality metrics assessment: Failed (coverage below threshold).")
        else:
            logging.debug("Quality metrics assessment: Passed (coverage above threshold).")

        # Scenario 1 & 2: Pathogenic variants identified by Kestrel
        if pathogenic_kestrel and confidence_level == "High_Precision":
            if quality_metrics_pass:
                summary_text += ("Pathogenic frameshift variant identified by Kestrel with high precision, "
                                 "and the VNTR coverage and quality metrics are above the threshold.")
                logging.debug("Scenario 1 applied: High precision with passing quality metrics.")
            else:
                summary_text += ("Pathogenic frameshift variant identified by Kestrel with high precision, "
                                 "but one or more quality metrics are below the threshold.")
                logging.debug("Scenario 2 applied: High precision with failing quality metrics.")

        # Scenario 3: Pathogenic variants identified with low precision
        elif pathogenic_kestrel and confidence_level == "Low_Precision":
            if quality_metrics_pass:
                summary_text += ("Warning: Pathogenic variant identified with low precision confidence. "
                                 "Validation through alternative methods (e.g., SNaPshot for dupC or "
                                 "long-read sequencing for other variants) is recommended.")
                logging.debug("Scenario 3a applied: Low precision with passing quality metrics.")
            else:
                summary_text += ("Warning: Pathogenic variant identified with low precision confidence and low-quality metrics. "
                                 "Validation using alternative methods is strongly recommended.")
                logging.debug("Scenario 3b applied: Low precision with failing quality metrics.")

        # Scenario 4: Check if both Kestrel and adVNTR have positive results
        if advntr_available and not advntr_df.empty and pathogenic_kestrel:
            # Determine if adVNTR also has a positive result
            # Assuming that a positive result means that adVNTR identified at least one variant
            adVNTR_positive = not advntr_df.empty

            if adVNTR_positive:
                summary_text += (" Both Kestrel and adVNTR genotyping methods have identified pathogenic variants and are concordant.")
                logging.debug("Scenario 4 applied: Both methods positive and concordant.")
            else:
                summary_text += (" There is a discrepancy between Kestrel and adVNTR genotyping methods regarding the identification of pathogenic variants.")
                logging.debug("Scenario 4 applied: Discrepancy between methods.")

        # If no pathogenic variants are identified by either method
        if not pathogenic_kestrel and (not advntr_available or advntr_df.empty):
            summary_text = "The screening was negative (no valid Kestrel or adVNTR data)."
            logging.debug("No pathogenic variants identified by either method; negative screening.")

    except Exception as ex:
        logging.error(f"Exception in build_screening_summary: {ex}")
        summary_text = "No summary available."

    logging.debug(f"Final summary_text: {summary_text}")
    return summary_text
