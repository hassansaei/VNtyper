import matplotlib
matplotlib.use('Agg')  # Use the Agg backend to avoid Qt-related warnings

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
    sample_id = Path(kestrel_result_file).parents[1].name  # Extract the sample ID
    logging.info(f"Loading Kestrel results from {kestrel_result_file}")
    if not os.path.exists(kestrel_result_file):
        logging.warning(f"Kestrel result file not found: {kestrel_result_file}")
        return pd.DataFrame({'Sample': [sample_id]})  # Return DataFrame with Sample ID only
    
    try:
        df = pd.read_csv(kestrel_result_file, sep='\t', comment='#')
        df['Sample'] = sample_id  # Add sample ID column

        # Define the required columns and their human-readable headers
        columns_to_display = {
            'Sample': 'Sample',
            'Motif': 'Motif',
            'Variant': 'Variant',
            'POS': 'Position',
            'REF': 'REF',
            'ALT': 'ALT',
            'Motif_sequence': 'Motif\nSequence',
            'Estimated_Depth_AlternateVariant': 'Depth\n(Variant)',
            'Estimated_Depth_Variant_ActiveRegion': 'Depth\n(Region)',
            'Depth_Score': 'Depth\nScore',
            'Confidence': 'Confidence'
        }

        # Select only the columns that are present in the DataFrame
        available_columns = {col: columns_to_display[col] for col in columns_to_display if col in df.columns}
        missing_columns = set(columns_to_display) - set(available_columns)

        if missing_columns:
            logging.warning(f"Missing columns in {kestrel_result_file}: {missing_columns}")

        # Rename the columns to human-readable names
        df = df[list(available_columns.keys())]
        df = df.rename(columns=available_columns)

        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame({'Sample': [sample_id]})  # Return DataFrame with Sample ID only

# Function to load adVNTR results
def load_advntr_results(advntr_result_file):
    sample_id = Path(advntr_result_file).parents[1].name  # Extract the sample ID
    logging.info(f"Loading adVNTR results from {advntr_result_file}")
    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']})
    
    try:
        df = pd.read_csv(advntr_result_file, sep='\t', comment='#')
        if df.empty:
            logging.warning(f"adVNTR result file {advntr_result_file} is empty.")
            return pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']})
        
        df['Sample'] = sample_id  # Add sample ID column

        # Return only the sample column and any relevant data (if present)
        return df[['Sample'] + [col for col in df.columns if col != 'Sample']]
    except pd.errors.EmptyDataError as e:
        logging.error(f"adVNTR result file {advntr_result_file} is empty: {e}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']})
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['Failed to parse adVNTR results']})
    except Exception as e:
        logging.error(f"Unexpected error occurred while loading adVNTR results: {e}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['Error loading adVNTR results']})

# Function to recursively find all files with a specific name in a directory tree
def find_results_files(root_dir, filename):
    logging.info(f"Searching for {filename} in directory {root_dir}")
    result_files = []
    for root, dirs, files in os.walk(root_dir):
        logging.debug(f"Checking directory: {root}")
        if filename in files:
            file_path = Path(root) / filename
            logging.info(f"Found {filename} in {file_path}")
            result_files.append(file_path)
    return result_files

# Function to load results from all matching files across all subdirectories
def load_results_from_dirs(input_dirs, filename, file_loader):
    dfs = []
    for input_dir in input_dirs:
        logging.info(f"Processing input directory: {input_dir}")
        result_files = find_results_files(input_dir, filename)
        if not result_files:
            # Handle cases where there are no result files in the directory
            sample_id = Path(input_dir).name
            logging.warning(f"No result files found in {input_dir}")
            dfs.append(pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']}))  # Append DataFrame with Sample ID and message
        else:
            for file in result_files:
                logging.info(f"Attempting to load file: {file}")
                df = file_loader(file)
                dfs.append(df)
        # log the final df
        logging.info(f"Loaded {len(dfs)} results from {input_dir}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

# Helper function to encode image as base64
def encode_image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

# Function to generate and save donut charts
def generate_donut_chart(values, labels, total, title, colors, plot_path):
    fig, ax = plt.subplots(figsize=(6, 6))
    wedgeprops = {'width': 0.3, 'edgecolor': 'black', 'linewidth': 2}

    ax.pie(values, 
           wedgeprops=wedgeprops, 
           startangle=90, 
           colors=colors, 
           labels=labels)

    ax.text(0, 0, f"{total}", ha='center', va='center', fontsize=24)
    ax.set_title(title)

    plt.savefig(plot_path)
    plt.close()

    return encode_image_to_base64(plot_path)

# Main function to aggregate cohort results
def aggregate_cohort(input_dirs, output_dir, summary_file, config_path=None):
    # Load results using the individual functions
    kestrel_df = load_results_from_dirs(input_dirs, "kestrel_result.tsv", load_kestrel_results)
    advntr_df = load_results_from_dirs(input_dirs, "output_adVNTR.tsv", load_advntr_results)

    # Generate summary report and plot
    generate_cohort_summary_report(output_dir, kestrel_df, advntr_df, summary_file)

# Function to generate the cohort summary report
def generate_cohort_summary_report(output_dir, kestrel_df, advntr_df, summary_file):
    # Create plots directory within the output directory
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Generate summary statistics for Kestrel
    kestrel_positive = len(kestrel_df[kestrel_df['Confidence'] != 'Negative'])
    kestrel_negative = len(kestrel_df[kestrel_df['Confidence'] == 'Negative'])
    total_kestrel = kestrel_positive + kestrel_negative

    # Generate summary statistics for adVNTR
    advntr_positive = len(advntr_df[(advntr_df['VID'].notna()) & (advntr_df['VID'] != 'Negative') & (advntr_df['Message'].isna())])
    advntr_negative = len(advntr_df[(advntr_df['VID'] == 'Negative')])
    advntr_no_data = len(advntr_df[advntr_df['Message'].notna()])
    total_advntr = advntr_positive + advntr_negative + advntr_no_data

    # Use color-blind-safe colors
    colors = {
        'positive': '#56B4E9',  # Sky blue
        'negative': '#D55E00',  # Vermillion
        'no_data': '#999999'    # Grey
    }

    # Generate the donut chart for Kestrel Results
    kestrel_plot_path = plots_dir / "kestrel_summary_plot.png"
    kestrel_plot_base64 = generate_donut_chart(
        values=[kestrel_positive, kestrel_negative],
        labels=['Positive', 'Negative'],
        total=total_kestrel,
        title='Kestrel Results',
        colors=[colors['positive'], colors['negative']],
        plot_path=kestrel_plot_path
    )

    # Generate the donut chart for adVNTR Results
    advntr_plot_path = plots_dir / "advntr_summary_plot.png"
    advntr_plot_base64 = generate_donut_chart(
        values=[advntr_positive, advntr_negative, advntr_no_data],
        labels=['Positive', 'Negative', 'No Data'],
        total=total_advntr,
        title='adVNTR Results',
        colors=[colors['positive'], colors['negative'], colors['no_data']],
        plot_path=advntr_plot_path
    )

    # Load the template
    template_dir = "vntyper/templates"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('cohort_summary_template.html')

    # Render the HTML report
    rendered_html = template.render(
        report_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        kestrel_positive=kestrel_df.to_html(classes='table table-bordered table-striped hover compact order-column table-sm', index=False),
        advntr_positive=advntr_df.to_html(classes='table table-bordered table-striped hover compact order-column table-sm', index=False),
        kestrel_plot_base64=kestrel_plot_base64,
        advntr_plot_base64=advntr_plot_base64
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
