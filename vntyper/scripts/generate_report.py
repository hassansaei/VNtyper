import pandas as pd
import os
import logging
import json
import subprocess
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from pathlib import Path

def load_kestrel_results(kestrel_result_file):
    logging.info(f"Loading Kestrel results from {kestrel_result_file}")
    if not os.path.exists(kestrel_result_file):
        logging.warning(f"Kestrel result file not found: {kestrel_result_file}")
        return pd.DataFrame()  # Return an empty DataFrame if the file is missing

    try:
        df = pd.read_csv(kestrel_result_file, sep='\t', comment='#')
        # Filter and rename columns
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

        # Apply conditional styling to the Confidence column
        df['Confidence'] = df['Confidence'].apply(
            lambda x: f'<span style="color:orange;font-weight:bold;">{x}</span>' if x == 'Low_Precision'
            else f'<span style="color:red;font-weight:bold;">{x}</span>' if x == 'High_Precision'
            else x
        )

        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame()

def load_advntr_results(advntr_result_file):
    logging.info(f"Loading adVNTR results from {advntr_result_file}")
    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame(), False  # Return empty DataFrame and False flag

    try:
        df = pd.read_csv(advntr_result_file, sep='\t', comment='#')
        return df, True  # Return DataFrame and True flag
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame(), False  # Return empty DataFrame and False flag

def load_pipeline_log(log_file):
    logging.info(f"Loading pipeline log from {log_file}")
    if not os.path.exists(log_file):
        logging.warning(f"Pipeline log file not found: {log_file}")
        return "Pipeline log file not found."

    try:
        with open(log_file, 'r') as f:
            return f.read()
    except Exception as e:
        logging.error(f"Failed to read pipeline log file: {e}")
        return "Failed to load pipeline log."

def run_igv_report(bed_file, bam_file, fasta_file, output_html, flanking=50):
    """
    Runs the IGV report generation command using the provided BED, BAM, and FASTA files.

    Args:
        bed_file (str or Path): Path to the BED file.
        bam_file (str or Path): Path to the BAM file.
        fasta_file (str or Path): Path to the reference FASTA file.
        output_html (str or Path): Path to the output HTML file for the IGV report.
        flanking (int): Flanking region for IGV reports.
    """
    # Convert Path objects to strings if needed
    bed_file = str(bed_file)
    bam_file = str(bam_file)
    fasta_file = str(fasta_file)
    output_html = str(output_html)

    igv_report_cmd = [
        'create_report',
        bed_file,
        '--flanking', str(flanking),
        '--fasta', fasta_file,
        '--tracks', bam_file,
        '--output', output_html
    ]
    try:
        logging.info(f"Running IGV report: {' '.join(igv_report_cmd)}")
        subprocess.run(igv_report_cmd, check=True)
        logging.info(f"IGV report successfully generated at {output_html}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error generating IGV report: {e}")
        raise

def extract_igv_content(igv_report_html):
    """
    Extracts the relevant content (variant table and IGV div) as well as the tableJson and sessionDictionary from the IGV report.
    """
    try:
        with open(igv_report_html, 'r') as f:
            content = f.read()

        # Extract the necessary parts: container for IGV visualization and variant table
        igv_start = content.find('<div id="container"')
        igv_end = content.find('</body>')

        if igv_start == -1 or igv_end == -1:
            logging.error("Failed to extract IGV content from report.")
            return "", "", ""

        igv_content = content[igv_start:igv_end].strip()

        # Extract tableJson and sessionDictionary
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

def generate_summary_report(output_dir, template_dir, report_file, log_file, bed_file, bam_file,
                            fasta_file, flanking=50, input_files=None, pipeline_version=None):
    """
    Generates a summary report that includes Kestrel results, adVNTR results, pipeline log,
    and IGV alignment visualizations.

    Args:
        output_dir (str): Output directory for the report.
        template_dir (str): Directory containing the report template.
        report_file (str): Name of the report file.
        log_file (str): Path to the pipeline log file.
        bed_file (str): Path to the BED file for IGV reports.
        bam_file (str): Path to the BAM file for IGV reports.
        fasta_file (str): Path to the reference FASTA file for IGV reports.
        flanking (int): Size of the flanking region for IGV reports.
        input_files (dict): Dictionary of input filenames (e.g., {'fastq1': 'sample_R1.fastq', 'fastq2': 'sample_R2.fastq'}).
        pipeline_version (str): The version of the VNtyper pipeline.
    """
    kestrel_result_file = Path(output_dir) / "kestrel/kestrel_result.tsv"
    advntr_result_file = Path(output_dir) / "advntr/output_adVNTR.tsv"
    igv_report_file = Path(output_dir) / "igv_report.html"  # Generated IGV report file

    # Only run IGV report if the BED file exists
    if bed_file and os.path.exists(bed_file):
        logging.info(f"Running IGV report for BED file: {bed_file}")
        run_igv_report(bed_file, bam_file, fasta_file, igv_report_file, flanking=flanking)
    else:
        logging.warning("BED file does not exist or not provided. Skipping IGV report generation.")
        igv_report_file = None  # No IGV report will be created

    # Load Kestrel results
    kestrel_df = load_kestrel_results(kestrel_result_file)

    # Load adVNTR results and get availability flag
    advntr_df, advntr_available = load_advntr_results(advntr_result_file)

    # Load pipeline log
    log_content = load_pipeline_log(log_file)

    # Extract IGV content, tableJson, and sessionDictionary for embedding, if IGV report exists
    if igv_report_file and os.path.exists(igv_report_file):
        igv_content, table_json, session_dictionary = extract_igv_content(igv_report_file)
    else:
        logging.warning("IGV report file not found. Skipping IGV content.")
        igv_content, table_json, session_dictionary = "", "", ""

    # Convert DataFrames to HTML tables with safe argument to allow HTML content in cells
    kestrel_html = kestrel_df.to_html(
        classes='table table-bordered table-striped hover compact order-column table-sm',
        index=False,
        escape=False
    )

    # Only convert advntr_html if advntr_available is True and DataFrame is not empty
    if advntr_available and not advntr_df.empty:
        advntr_html = advntr_df.to_html(
            classes='table table-bordered table-striped hover compact order-column table-sm',
            index=False
        )
    else:
        advntr_html = None  # Set to None if not available

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('report_template.html')

    # Prepare context for the template
    context = {
        'kestrel_highlight': kestrel_html,
        'advntr_highlight': advntr_html,
        'advntr_available': advntr_available,
        'log_content': log_content,
        'igv_content': igv_content,
        'table_json': table_json,
        'session_dictionary': session_dictionary,
        'report_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'input_files': input_files,               # New
        'pipeline_version': pipeline_version      # New
    }

    # Render the template with the context
    rendered_html = template.render(context)

    report_file_path = Path(output_dir) / report_file
    with open(report_file_path, 'w') as f:
        f.write(rendered_html)

    logging.info(f"Summary report generated and saved to {report_file_path}")
