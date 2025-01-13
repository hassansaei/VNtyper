#!/usr/bin/env python3

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

matplotlib.use('Agg')


def load_kestrel_results(kestrel_result_file):
    """
    Load and process Kestrel genotyping results.

    The Kestrel output is expected to be a TSV file containing genotyping
    information such as motif, variant, position, confidence, and coverage
    statistics. If the file doesn't exist or cannot be parsed, returns a
    DataFrame with just the 'Sample' column.

    Parameters
    ----------
    kestrel_result_file : str or Path
        Path to the Kestrel result TSV file.

    Returns
    -------
    pandas.DataFrame
        Processed Kestrel results with human-readable column names and
        conditional styling applied to the 'Confidence' column.
    """
    sample_id = Path(kestrel_result_file).parents[1].name
    logging.info(f"Loading Kestrel results from {kestrel_result_file}")

    if not os.path.exists(kestrel_result_file):
        logging.warning(f"Kestrel result file not found: {kestrel_result_file}")
        return pd.DataFrame({'Sample': [sample_id]})

    try:
        df = pd.read_csv(kestrel_result_file, sep='\t', comment='#')
        df['Sample'] = sample_id
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
        available_columns = {
            col: columns_to_display[col]
            for col in columns_to_display if col in df.columns
        }
        missing_columns = set(columns_to_display) - set(available_columns)
        if missing_columns:
            logging.warning(f"Missing columns in {kestrel_result_file}: {missing_columns}")
        if not available_columns:
            logging.warning("No expected columns found in Kestrel results, returning minimal DataFrame.")
            return pd.DataFrame({'Sample': [sample_id]})
        df = df[list(available_columns.keys())]
        df = df.rename(columns=available_columns)
        if 'Confidence' in df.columns:
            df['Confidence'] = df['Confidence'].apply(
                lambda x: (
                    f'<span style="color:orange;font-weight:bold;">{x}</span>'
                    if x == 'Low_Precision'
                    else f'<span style="color:red;font-weight:bold;">{x}</span>'
                    if x in ['High_Precision', 'High_Precision*']  # Updated to include 'High_Precision*'
                    else f'<span style="color:blue;font-weight:bold;">{x}</span>'
                    if x == 'Negative'
                    else x
                )
            )
        return df
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse Kestrel result file: {e}")
        return pd.DataFrame({'Sample': [sample_id]})
    except Exception as e:
        logging.error(f"Unexpected error loading Kestrel results: {e}")
        return pd.DataFrame({'Sample': [sample_id]})


def load_advntr_results(advntr_result_file):
    """
    Load and process adVNTR genotyping results.

    adVNTR results are expected in a TSV format, potentially with columns for
    sample ID, VNTR ID (VID), state, supporting reads, coverage, and p-values.
    This function ensures a full schema is always returned, even if no file
    is found or if the file is empty.

    Parameters
    ----------
    advntr_result_file : str or Path
        Path to the adVNTR result TSV file.

    Returns
    -------
    tuple
        A tuple (pandas.DataFrame, bool) where:
        - DataFrame contains the full set of expected columns. If no data is found,
          returns a DataFrame with NaNs and a "Message" column indicating no output.
        - bool indicating whether valid data was found (True) or not (False).
    """
    full_cols = [
        "Sample",
        "VID",
        "State",
        "NumberOfSupportingReads",
        "MeanCoverage",
        "Pvalue",
        "Message"
    ]

    def empty_advntr_row(sample_id, message="no adVNTR output found."):
        row_data = {col: [float('nan')] for col in full_cols if col not in ["Sample", "Message"]}
        row_data["Sample"] = [sample_id]
        row_data["Message"] = [message]
        return pd.DataFrame(row_data, columns=full_cols)

    sample_id = Path(advntr_result_file).parents[1].name
    logging.info(f"Loading adVNTR results from {advntr_result_file}")

    if not os.path.exists(advntr_result_file):
        logging.warning(f"adVNTR result file not found: {advntr_result_file}")
        return empty_advntr_row(sample_id), False

    try:
        df = pd.read_csv(advntr_result_file, sep='\t', comment='#')
        if df.empty:
            logging.warning(f"adVNTR result file {advntr_result_file} is empty.")
            return empty_advntr_row(sample_id), False
        df['Sample'] = sample_id
        for col in full_cols:
            if col not in df.columns:
                df[col] = float('nan')
        df = df[full_cols]
        return df, True
    except pd.errors.EmptyDataError as e:
        logging.error(f"adVNTR result file {advntr_result_file} is empty: {e}")
        return empty_advntr_row(sample_id), False
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        return empty_advntr_row(sample_id, "Failed to parse adVNTR results"), False
    except Exception as e:
        logging.error(f"Unexpected error occurred while loading adVNTR results: {e}")
        return empty_advntr_row(sample_id, "Error loading adVNTR results"), False


def find_results_files(root_dir, filename):
    """
    Recursively find all files with a specific name in a directory tree.

    Parameters
    ----------
    root_dir : str or Path
        Root directory to start the search.
    filename : str
        Name of the file to search for.

    Returns
    -------
    list
        List of Path objects matching the filename.
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

    If loading adVNTR and no files are found in a given directory, this function
    attempts to return one row per subfolder, treating each as a sample. If there
    are no subfolders, it falls back to one row for the directory itself.

    Parameters
    ----------
    input_dirs : list
        List of directories to search.
    filename : str
        Name of the result files to search for.
    file_loader : function
        Function to load and process each file.

    Returns
    -------
    pandas.DataFrame
        Concatenated DataFrame of all loaded results.
    """
    dfs = []
    loading_advntr = ('advntr' in filename.lower())

    if loading_advntr:
        full_cols = [
            "Sample",
            "VID",
            "State",
            "NumberOfSupportingReads",
            "MeanCoverage",
            "Pvalue",
            "Message"
        ]

        def empty_advntr_row(sample_id, message="no adVNTR output found."):
            row_data = {col: [float('nan')] for col in full_cols if col not in ["Sample", "Message"]}
            row_data["Sample"] = [sample_id]
            row_data["Message"] = [message]
            return pd.DataFrame(row_data, columns=full_cols)

    for input_dir in input_dirs:
        logging.info(f"Processing input directory: {input_dir}")
        result_files = find_results_files(input_dir, filename)

        if not result_files:
            if loading_advntr:
                subfolders = [d for d in Path(input_dir).iterdir() if d.is_dir()]
                if subfolders:
                    for subf in subfolders:
                        sample_id = subf.name
                        dfs.append(empty_advntr_row(sample_id))
                else:
                    sample_id = Path(input_dir).name
                    dfs.append(empty_advntr_row(sample_id))
            else:
                sample_id = Path(input_dir).name
                dfs.append(pd.DataFrame({
                    'Sample': [sample_id],
                    'Message': [f'No {filename} results found for this sample']
                }))
        else:
            for file in result_files:
                logging.info(f"Attempting to load file: {file}")
                result = file_loader(file)
                if isinstance(result, tuple):
                    df = result[0]
                else:
                    df = result
                dfs.append(df)

        logging.info(f"Loaded {len(result_files)} {filename} files from {input_dir}")

    if not dfs:
        logging.warning(f"No data found at all for {filename}. Returning empty DataFrame.")
        return pd.DataFrame()

    return pd.concat(dfs, ignore_index=True)


def encode_image_to_base64(image_path):
    """
    Encode an image file to a base64 string.

    Parameters
    ----------
    image_path : str or Path
        Path to the image file.

    Returns
    -------
    str
        Base64-encoded string of the image.
    """
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"


def generate_donut_chart(values, labels, total, title, colors, plot_path=None, interactive=False):
    """
    Generate and save a donut chart (static or interactive).

    For static plots, matplotlib is used. For interactive plots, Plotly is used.
    This chart visualizes categories as parts of a donut, with the total in the center.

    Parameters
    ----------
    values : list
        Values for each segment of the donut chart.
    labels : list
        Labels for each segment.
    total : int
        Total value displayed in the center of the donut.
    title : str
        Title of the chart.
    colors : list
        Colors for each segment of the donut.
    plot_path : str or Path, optional
        Path to save the static plot image.
    interactive : bool
        Whether to generate an interactive Plotly chart.

    Returns
    -------
    str
        Base64-encoded image string for static charts or HTML string for interactive charts.
    """
    if interactive:
        fig = go.Figure(go.Pie(
            labels=labels,
            values=values,
            hole=0.6,
            marker=dict(colors=colors, line=dict(color='black', width=2)),
            textinfo='none'
        ))
        fig.update_layout(
            title={
                'text': title,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            annotations=[
                dict(text=f'<b>{total}</b>', x=0.5, y=0.5, font_size=40, showarrow=False),
                dict(text=labels[0] if len(labels) > 0 else '', x=0.15, y=0.5, font_size=14, showarrow=False),
                dict(text=labels[1] if len(labels) > 1 else '', x=0.85, y=0.5, font_size=14, showarrow=False)
            ],
            showlegend=False,
            margin=dict(t=50, b=50, l=50, r=50),
            height=500, width=500
        )
        return pio.to_html(fig, full_html=False)
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        wedgeprops = {'width': 0.3, 'edgecolor': 'black', 'linewidth': 2}
        try:
            ax.pie(values, wedgeprops=wedgeprops, startangle=90, colors=colors, labels=labels)
            ax.text(0, 0, f"{total}", ha='center', va='center', fontsize=24)
            ax.set_title(title)
            if plot_path:
                plt.savefig(plot_path)
            else:
                logging.warning("No plot_path provided for static donut chart, chart not saved.")
        except Exception as e:
            logging.error(f"Error generating donut chart: {e}")
        plt.close()
        if plot_path and os.path.exists(plot_path):
            return encode_image_to_base64(plot_path)
        else:
            return ""


def generate_cohort_summary_report(output_dir, kestrel_df, advntr_df, summary_file, config):
    """
    Generate the cohort summary report combining Kestrel and adVNTR results.

    This function creates summary statistics, generates donut charts for each data type,
    and then renders a Jinja2 template to produce a final HTML report.

    Parameters
    ----------
    output_dir : str or Path
        Output directory for the report.
    kestrel_df : pandas.DataFrame
        DataFrame containing Kestrel results.
    advntr_df : pandas.DataFrame
        DataFrame containing adVNTR results.
    summary_file : str
        Name of the summary report file.
    config : dict
        Configuration dictionary containing paths and settings.

    Returns
    -------
    None
        Writes the HTML report to the specified summary file.
    """
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    if 'Confidence' in kestrel_df.columns:
        try:
            kestrel_df_conf = kestrel_df['Confidence'].fillna('')
            # Updated to include 'High_Precision*' in positive counts
            kestrel_positive = len(kestrel_df[kestrel_df_conf.str.contains('Low_Precision|High_Precision\\*?', na=False)])
            kestrel_negative = len(kestrel_df[~kestrel_df_conf.str.contains('Low_Precision|High_Precision\\*?', na=False)])
        except Exception as e:
            logging.error(f"Error processing 'Confidence' values: {e}")
            kestrel_positive = 0
            kestrel_negative = 0
    else:
        logging.warning("No 'Confidence' column found in Kestrel results. Setting positive/negative counts to 0.")
        kestrel_positive = 0
        kestrel_negative = 0

    total_kestrel = kestrel_positive + kestrel_negative

    if 'Message' not in advntr_df.columns:
        logging.warning("No 'Message' column found in adVNTR results, adding it.")
        advntr_df['Message'] = None

    if 'VID' in advntr_df.columns:
        try:
            advntr_positive = len(
                advntr_df[
                    (advntr_df['VID'].notna()) &
                    (advntr_df['VID'] != 'Negative') &
                    (advntr_df['Message'].isna())
                ]
            )
            advntr_negative = len(advntr_df[advntr_df['VID'] == 'Negative'])
            advntr_no_data = len(advntr_df[advntr_df['Message'].notna()])
        except Exception as e:
            logging.error(f"Error processing adVNTR values: {e}")
            advntr_positive = 0
            advntr_negative = 0
            advntr_no_data = len(advntr_df)
    else:
        logging.warning("No 'VID' column found in adVNTR results. Treating all as no data.")
        advntr_positive = 0
        advntr_negative = 0
        advntr_no_data = len(advntr_df)

    total_advntr = advntr_positive + advntr_negative + advntr_no_data

    color_list = config.get("visualization", {}).get("donut_colors", ["#56B4E9", "#D55E00", "#999999"])
    colors = {
        'positive': color_list[0],
        'negative': color_list[1],
        'no_data': color_list[2]
    }

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
    kestrel_plot_html = generate_donut_chart(
        values=[kestrel_positive, kestrel_negative],
        labels=['Positive', 'Negative'],
        total=total_kestrel,
        title='Kestrel Results',
        colors=[colors['positive'], colors['negative']],
        plot_path=None,
        interactive=True
    )

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
    advntr_plot_html = generate_donut_chart(
        values=[advntr_positive, advntr_negative, advntr_no_data],
        labels=['Positive', 'Negative', 'No Data'],
        total=total_advntr,
        title='adVNTR Results',
        colors=[colors['positive'], colors['negative'], colors['no_data']],
        plot_path=None,
        interactive=True
    )

    template_dir = config.get('paths', {}).get('template_dir', 'vntyper/templates')
    env = Environment(loader=FileSystemLoader(template_dir))
    try:
        template = env.get_template('cohort_summary_template.html')
    except Exception as e:
        logging.error(f"Failed to load Jinja2 template: {e}")
        raise

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
        'interactive': True
    }

    try:
        rendered_html = template.render(context)
    except Exception as e:
        logging.error(f"Failed to render the cohort summary template: {e}")
        raise

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

    This function uses load_results_from_dirs to gather Kestrel and adVNTR results
    from multiple directories, merges them into DataFrames, and then generates a
    single cohort summary report.

    Parameters
    ----------
    input_dirs : list
        List of directories containing output files to aggregate. Each directory
        represents a separate sample or run.
    output_dir : str or Path
        Output directory for the aggregated summary report.
    summary_file : str
        Name of the summary report file.
    config : dict
        Configuration dictionary containing paths and other settings.

    Returns
    -------
    None
        Writes the cohort summary report to the specified output directory.
    """
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

    generate_cohort_summary_report(
        output_dir=output_dir,
        kestrel_df=kestrel_df,
        advntr_df=advntr_df,
        summary_file=summary_file,
        config=config
    )
