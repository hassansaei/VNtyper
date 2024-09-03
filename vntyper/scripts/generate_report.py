import pandas as pd
import os
import logging
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
        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame()

def load_advntr_results(advntr_result_file):
    logging.info(f"Loading adVNTR results from {advntr_result_file}")
    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(advntr_result_file, sep='\t', comment='#')
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame()

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

def generate_summary_report(output_dir, template_dir, report_file, log_file):
    kestrel_result_file = Path(output_dir) / "kestrel/kestrel_result.tsv"
    advntr_result_file = Path(output_dir) / "advntr/output_adVNTR.tsv"
    log_file_path = log_file

    kestrel_df = load_kestrel_results(kestrel_result_file)
    advntr_df = load_advntr_results(advntr_result_file)
    log_content = load_pipeline_log(log_file_path)

    kestrel_html = kestrel_df.to_html(classes='table table-bordered table-striped hover compact order-column table-sm', index=False)
    advntr_html = advntr_df.to_html(classes='table table-bordered table-striped hover compact order-column table-sm', index=False)

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('report_template.html')

    rendered_html = template.render(
        kestrel_highlight=kestrel_html,
        advntr_highlight=advntr_html if not advntr_df.empty else "No significant adVNTR variants found.",
        log_content=log_content,
        report_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )

    report_file_path = Path(output_dir) / report_file
    with open(report_file_path, 'w') as f:
        f.write(rendered_html)

    logging.info(f"Summary report generated and saved to {report_file_path}")
