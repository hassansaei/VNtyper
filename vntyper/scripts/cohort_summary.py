# vntyper/scripts/cohort_summary.py

import os
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd
import base64
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader
import plotly.graph_objects as go

# Use the Agg backend to avoid Qt-related warnings
matplotlib.use('Agg')


def load_kestrel_results(kestrel_result_file):
    """
    Load and process Kestrel genotyping results.

    Args:
        kestrel_result_file (str or Path): Path to the Kestrel result TSV file.

    Returns:
        pandas.DataFrame: Processed Kestrel results.
    """
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
        available_columns = {
            col: columns_to_display[col]
            for col in columns_to_display if col in df.columns
        }
        missing_columns = set(columns_to_display) - set(available_columns)

        if missing_columns:
            logging.warning(f"Missing columns in {kestrel_result_file}: {missing_columns}")

        # Rename the columns to human-readable names
        df = df[list(available_columns.keys())]
        df = df.rename(columns=available_columns)

        # Apply conditional styling to the Confidence column
        df['Confidence'] = df['Confidence'].apply(
            lambda x: (
                f'<span style="color:orange;font-weight:bold;">{x}</span>'
                if x == 'Low_Precision'
                else f'<span style="color:red;font-weight:bold;">{x}</span>'
                if x == 'High_Precision'
                else f'<span style="color:blue;font-weight:bold;">{x}</span>'
                if x == 'Negative'
                else x
            )
        )

        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame({'Sample': [sample_id]})  # Return DataFrame with Sample ID only
    except Exception as e:
        logging.error(f"Unexpected error loading Kestrel results: {e}")
        return pd.DataFrame({'Sample': [sample_id]})  # Return DataFrame with Sample ID only


def load_advntr_results(advntr_result_file):
    """
    Load and process adVNTR genotyping results.

    Args:
        advntr_result_file (str or Path): Path to the adVNTR result TSV file.

    Returns:
        tuple: (pandas.DataFrame, bool) Processed adVNTR results and availability flag.
    """
    sample_id = Path(advntr_result_file).parents[1].name  # Extract the sample ID
    logging.info(f"Loading adVNTR results from {advntr_result_file}")

    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']}), False

    try:
        df = pd.read_csv(advntr_result_file, sep='\t', comment='#')
        if df.empty:
            logging.warning(f"adVNTR result file {advntr_result_file} is empty.")
            return pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']}), False

        df['Sample'] = sample_id  # Add sample ID column

        # Return only the sample column and any relevant data (if present)
        relevant_columns = ['Sample'] + [col for col in df.columns if col != 'Sample']
        df = df[relevant_columns]

        return df, True  # Return DataFrame and True flag
    except pd.errors.EmptyDataError as e:
        logging.error(f"adVNTR result file {advntr_result_file} is empty: {e}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['No adVNTR results found for this sample']}), False
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['Failed to parse adVNTR results']}), False
    except Exception as e:
        logging.error(f"Unexpected error occurred while loading adVNTR results: {e}")
        return pd.DataFrame({'Sample': [sample_id], 'Message': ['Error loading adVNTR results']}), False


def find_results_files(root_dir, filename):
    """
    Recursively find all files with a specific name in a directory tree.

    Args:
        root_dir (str or Path): Root directory to start the search.
        filename (str): Name of the file to search for.

    Returns:
        list: List of Path objects matching the filename.
    """
    logging.info(f"Searching for {filename} in directory {root_dir}")
    result_files = []
    for root, _, files in os.walk(root_dir):
        logging.debug(f"Checking directory: {root}")
        if filename in files:
            file_path = Path(root) / filename
            logging.info(f"Found {filename} in {file_path}")
            result_files.append(file_path)
    return result_files


def load_results_from_dirs(input_dirs, filename, file_loader):
    """
    Load results from all matching files across multiple directories.

    Args:
        input_dirs (list): List of directories to search.
        filename (str): Name of the result files to search for.
        file_loader (function): Function to load and process each file.

    Returns:
        pandas.DataFrame: Concatenated DataFrame of all loaded results.
    """
    dfs = []
    for input_dir in input_dirs:
        logging.info(f"Processing input directory: {input_dir}")
        result_files = find_results_files(input_dir, filename)
        if not result_files:
            # Handle cases where there are no result files in the directory
            sample_id = Path(input_dir).name
            logging.warning(f"No {filename} files found in {input_dir}")
            dfs.append(pd.DataFrame({
                'Sample': [sample_id],
                'Message': [f'No {filename} results found for this sample']
            }))
        else:
            for file in result_files:
                logging.info(f"Attempting to load file: {file}")
                result = file_loader(file)
                
                # Check if the result is a tuple (as in load_advntr_results)
                if isinstance(result, tuple):
                    df = result[0]  # Extract the DataFrame
                else:
                    df = result
                
                dfs.append(df)
        # Log the number of results loaded from the directory
        logging.info(f"Loaded {len(result_files)} {filename} files from {input_dir}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def encode_image_to_base64(image_path):
    """
    Encode an image file to a base64 string.

    Args:
        image_path (str or Path): Path to the image file.

    Returns:
        str: Base64-encoded string of the image.
    """
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"


def generate_donut_chart(values, labels, total, title, colors, plot_path=None, interactive=False):
    """
    Generate and save an interactive or static donut chart.

    Args:
        values (list): Values for each segment of the donut chart.
        labels (list): Labels for each segment.
        total (int): Total value to display in the center.
        title (str): Title of the chart.
        colors (list): Colors for each segment.
        plot_path (str or Path, optional): Path to save the static plot image.
        interactive (bool): Whether to generate an interactive Plotly chart.

    Returns:
        str: Base64-encoded image string for static charts or HTML string for interactive charts.
    """
    if interactive:
        # Generate an interactive donut chart using Plotly
        fig = go.Figure(go.Pie(
            labels=labels,
            values=values,
            hole=0.6,  # Create the donut shape
            marker=dict(colors=colors, line=dict(color='black', width=2)),  # Add black border
            textinfo='none'  # Disable percentage display inside the chart
        ))

        # Add title and central number inside the donut chart
        fig.update_layout(
            title={
                'text': title,
                'y': 0.95,  # Adjust the vertical positioning of the title
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            annotations=[
                dict(
                    text=f'<b>{total}</b>',
                    x=0.5, y=0.5,
                    font_size=40, showarrow=False
                ),  # Centralized total number
                dict(
                    text=labels[0],
                    x=0.15, y=0.5,
                    font_size=14, showarrow=False
                ),  # Place first label to the left
                dict(
                    text=labels[1],
                    x=0.85, y=0.5,
                    font_size=14, showarrow=False
                )   # Place second label to the right
            ],
            showlegend=False,  # Disable the default legend
            margin=dict(t=50, b=50, l=50, r=50),  # Adjust margins
            height=500, width=500  # Set dimensions
        )

        # Return the interactive plot HTML
        return pio.to_html(fig, full_html=False)

    else:
        # Static chart using matplotlib
        fig, ax = plt.subplots(figsize=(6, 6))
        wedgeprops = {'width': 0.3, 'edgecolor': 'black', 'linewidth': 2}
        ax.pie(values, wedgeprops=wedgeprops, startangle=90, colors=colors, labels=labels)
        ax.text(0, 0, f"{total}", ha='center', va='center', fontsize=24)
        ax.set_title(title)
        plt.savefig(plot_path)
        plt.close()
        return encode_image_to_base64(plot_path)


def generate_cohort_summary_report(output_dir, kestrel_df, advntr_df, summary_file, config):
    """
    Generate the cohort summary report with static and interactive plots.

    Args:
        output_dir (str or Path): Output directory for the report.
        kestrel_df (pandas.DataFrame): DataFrame containing Kestrel results.
        advntr_df (pandas.DataFrame): DataFrame containing adVNTR results.
        summary_file (str): Name of the summary report file.
        config (dict): Configuration dictionary.
    """
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Summary statistics for Kestrel
    kestrel_positive = len(kestrel_df[kestrel_df['Confidence'].str.contains('Low_Precision|High_Precision')])
    kestrel_negative = len(kestrel_df[~kestrel_df['Confidence'].str.contains('Low_Precision|High_Precision')])
    total_kestrel = kestrel_positive + kestrel_negative

    # Handle missing 'Message' column in adVNTR results
    if 'Message' not in advntr_df.columns:
        advntr_df['Message'] = None  # Add a default 'Message' column with None values if it's missing

    # Summary statistics for adVNTR
    advntr_positive = len(advntr_df[(advntr_df['VID'].notna()) &
                                    (advntr_df['VID'] != 'Negative') &
                                    (advntr_df['Message'].isna())])
    advntr_negative = len(advntr_df[advntr_df['VID'] == 'Negative'])
    advntr_no_data = len(advntr_df[advntr_df['Message'].notna()])
    total_advntr = advntr_positive + advntr_negative + advntr_no_data

    colors = {
        'positive': '#56B4E9',  # Sky blue
        'negative': '#D55E00',  # Vermillion
        'no_data': '#999999'     # Grey
    }

    # Static Kestrel plot (matplotlib)
    kestrel_plot_path = plots_dir / "kestrel_summary_plot.png"
    kestrel_plot_base64 = generate_donut_chart(
        values=[kestrel_positive, kestrel_negative],
        labels=['Positive', 'Negative'],
        total=total_kestrel,
        title='Kestrel Results',
        colors=[colors['positive'], colors['negative']],
        plot_path=kestrel_plot_path,
        interactive=False
    )

    # Interactive Kestrel plot (Plotly)
    kestrel_plot_html = generate_donut_chart(
        values=[kestrel_positive, kestrel_negative],
        labels=['Positive', 'Negative'],
        total=total_kestrel,
        title='Kestrel Results',
        colors=[colors['positive'], colors['negative']],
        plot_path=None,
        interactive=True
    )

    # Static adVNTR plot (matplotlib)
    advntr_plot_path = plots_dir / "advntr_summary_plot.png"
    advntr_plot_base64 = generate_donut_chart(
        values=[advntr_positive, advntr_negative, advntr_no_data],
        labels=['Positive', 'Negative', 'No Data'],
        total=total_advntr,
        title='adVNTR Results',
        colors=[colors['positive'], colors['negative'], colors['no_data']],
        plot_path=advntr_plot_path,
        interactive=False
    )

    # Interactive adVNTR plot (Plotly)
    advntr_plot_html = generate_donut_chart(
        values=[advntr_positive, advntr_negative, advntr_no_data],
        labels=['Positive', 'Negative', 'No Data'],
        total=total_advntr,
        title='adVNTR Results',
        colors=[colors['positive'], colors['negative'], colors['no_data']],
        plot_path=None,
        interactive=True
    )

    # Load the template from config
    template_dir = config.get('paths', {}).get('template_dir', 'vntyper/templates')
    env = Environment(loader=FileSystemLoader(template_dir))
    try:
        template = env.get_template('cohort_summary_template.html')
    except Exception as e:
        logging.error(f"Failed to load Jinja2 template: {e}")
        raise

    # Prepare context for the template
    context = {
        'report_date': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'kestrel_positive': kestrel_df.to_html(
            classes='table table-bordered table-striped', index=False, escape=False
        ),
        'advntr_positive': advntr_df.to_html(
            classes='table table-bordered table-striped', index=False, escape=False
        ),
        'kestrel_plot_base64': kestrel_plot_base64,
        'advntr_plot_base64': advntr_plot_base64,
        'kestrel_plot_interactive': kestrel_plot_html,
        'advntr_plot_interactive': advntr_plot_html,
        'interactive': True  # Flag to toggle between static and interactive plots in the template
    }

    # Render the template with the context
    try:
        rendered_html = template.render(context)
    except Exception as e:
        logging.error(f"Failed to render the cohort summary template: {e}")
        raise

    # Write the rendered HTML to the summary report file
    report_file_path = Path(output_dir) / summary_file
    try:
        with open(report_file_path, 'w') as f:
            f.write(rendered_html)
        logging.info(f"Cohort summary report generated and saved to {report_file_path}")
    except Exception as e:
        logging.error(f"Failed to write the cohort summary report: {e}")
        raise


def aggregate_cohort(input_dirs, output_dir, summary_file, config):
    """
    Aggregate outputs from multiple runs into a single summary file.

    Args:
        input_dirs (list): List of directories containing output files to aggregate.
        output_dir (str or Path): Output directory for the aggregated summary.
        summary_file (str): Name of the cohort summary report file.
        config (dict): Configuration dictionary.
    """
    # Load results using the individual functions
    kestrel_df = load_results_from_dirs(
        input_dirs=input_dirs,
        filename="kestrel_result.tsv",
        file_loader=load_kestrel_results
    )
    advntr_df = load_results_from_dirs(
        input_dirs=input_dirs,
        filename="output_adVNTR.tsv",
        file_loader=load_advntr_results
    )

    # Generate summary report and plot
    generate_cohort_summary_report(
        output_dir=output_dir,
        kestrel_df=kestrel_df,
        advntr_df=advntr_df,
        summary_file=summary_file,
        config=config
    )
