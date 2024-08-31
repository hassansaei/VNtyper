import pandas as pd
import os
import logging
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import base64

# Function to load Kestrel results
def load_kestrel_results(kestrel_result_file):
    logging.info(f"Loading Kestrel results from {kestrel_result_file}")
    if not os.path.exists(kestrel_result_file):
        logging.warning(f"Kestrel result file not found: {kestrel_result_file}")
        return pd.DataFrame()  # Return an empty DataFrame if the file is missing
    
    try:
        return pd.read_csv(kestrel_result_file, sep='\t', comment='#')
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame()  # Return an empty DataFrame if parsing fails

# Function to load adVNTR results
def load_advntr_results(advntr_result_file):
    logging.info(f"Loading adVNTR results from {advntr_result_file}")
    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame()  # Return an empty DataFrame if the file is missing
    
    try:
        return pd.read_csv(advntr_result_file, sep='\t', comment='#')
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame()  # Return an empty DataFrame if parsing fails

# Function to recursively find all files with a specific name in a directory tree
def find_results_files(root_dir, filename):
    result_files = []
    for root, dirs, files in os.walk(root_dir):
        if filename in files:
            result_files.append(Path(root) / filename)
    return result_files

# Function to load results from all matching files across all subdirectories
def load_results_from_dirs(input_dirs, filename, file_loader):
    dfs = []
    for input_dir in input_dirs:
        result_files = find_results_files(input_dir, filename)
        for file in result_files:
            df = file_loader(file)
            if not df.empty:
                dfs.append(df)
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()

# Helper function to encode image as base64
def encode_image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

# Main function to aggregate cohort results
def aggregate_cohort(input_dirs, output_dir, summary_file, config_path=None):
    # Load results using the individual functions
    kestrel_df = load_results_from_dirs(input_dirs, "kestrel_result.tsv", load_kestrel_results)
    advntr_df = load_results_from_dirs(input_dirs, "output_adVNTR.tsv", load_advntr_results)

    # Generate summary report and plot
    generate_cohort_summary_report(output_dir, kestrel_df, advntr_df, summary_file)

# Function to generate the cohort summary report
def generate_cohort_summary_report(output_dir, kestrel_positive, advntr_positive, summary_file):
    # Create plots directory within the output directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Generate summary statistics
    total_kestrel = len(kestrel_positive)
    total_advntr = len(advntr_positive)

    # Plot summary and save the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(x=['Kestrel Positive', 'adVNTR Positive'], y=[total_kestrel, total_advntr], ax=ax)
    plot_path = plots_dir / "cohort_summary_plot.png"
    plt.savefig(plot_path)
    plt.close()

    # Convert plot to base64
    plot_base64 = encode_image_to_base64(plot_path)

    # Load the template
    template_dir = "vntyper/templates"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('cohort_summary_template.html')

    # Render the HTML report
    rendered_html = template.render(
        report_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        kestrel_positive=kestrel_positive.to_html(classes='table table-striped'),
        advntr_positive=advntr_positive.to_html(classes='table table-striped') if not advntr_positive.empty else "No significant adVNTR variants found.",
        plot_base64=plot_base64
    )

    # Save the HTML report
    report_file_path = Path(output_dir) / summary_file
    with open(report_file_path, 'w') as f:
        f.write(rendered_html)

    logging.info(f"Cohort summary report generated and saved to {report_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cohort Summary Report Generator")
    parser.add_argument("-i", "--input_dirs", nargs='+', required=True, help="List of input directories or a text file with directories")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for the summary report")
    parser.add_argument("-s", "--summary_file", default="cohort_summary.html", help="Name of the summary report file")
    args = parser.parse_args()

    input_dirs = []
    for i in args.input_dirs:
        if os.path.isfile(i) and i.endswith('.txt'):
            with open(i, 'r') as f:
                input_dirs.extend([line.strip() for line in f])
        else:
            input_dirs.append(i)

    aggregate_cohort(input_dirs, args.output_dir, args.summary_file)
