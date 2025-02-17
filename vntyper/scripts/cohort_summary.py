#!/usr/bin/env python3

import os
import logging
from datetime import datetime
from pathlib import Path
import zipfile
import tempfile
import shutil

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
            logging.warning(
                f"Missing columns in {kestrel_result_file}: {missing_columns}"
            )
        if not available_columns:
            logging.warning(
                "No expected columns found in Kestrel results, returning minimal DataFrame."
            )
            return pd.DataFrame({'Sample': [sample_id]})
        df = df[list(available_columns.keys())]
        df = df.rename(columns=available_columns)
        if 'Confidence' in df.columns:
            df['Confidence'] = df['Confidence'].apply(
                lambda x: (
                    f'<span style="color:orange;font-weight:bold;">{x}</span>'
                    if x == 'Low_Precision'
                    else f'<span style="color:red;font-weight:bold;">{x}</span>'
                    if x in ['High_Precision', 'High_Precision*']
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

    The primary file is expected to be "output_adVNTR_result.tsv". If that file
    does not exist, the sample is marked as Not Performed.
    
    - If the TSV exists and contains data rows, the sample is marked as Positive.
    - If the TSV exists but contains no data rows (only header), the sample is marked as Negative.
    - If no file exists, the sample is marked as Not Performed.

    Parameters
    ----------
    advntr_result_file : str or Path
        Path to the adVNTR result file (primary TSV).

    Returns
    -------
    tuple
        A tuple (pandas.DataFrame, bool) where:
        - The DataFrame contains columns: Sample, VID, Variant, NumberOfSupportingReads, MeanCoverage, Pvalue, AdvntrState.
        - The bool is True if the TSV was found (even if negative), False if no file.
    """
    cols = ["VID", "Variant", "NumberOfSupportingReads", "MeanCoverage", "Pvalue"]
    out_cols = ["Sample"] + cols + ["AdvntrState"]
    sample_id = Path(advntr_result_file).parents[1].name
    logging.info(f"Loading adVNTR results from {advntr_result_file}")

    primary_file = Path(advntr_result_file).with_name("output_adVNTR_result.tsv")
    if not primary_file.exists():
        logging.warning(f"Primary adVNTR TSV not found for sample {sample_id}. Marking as Not Performed.")
        data = {
            "Sample": [sample_id],
            "VID": [""],
            "Variant": ["Not Performed"],
            "NumberOfSupportingReads": [float('nan')],
            "MeanCoverage": [float('nan')],
            "Pvalue": [float('nan')],
            "AdvntrState": ["Not Performed"]
        }
        return pd.DataFrame(data, columns=out_cols), False

    try:
        df = pd.read_csv(primary_file, sep='\t', comment='#')
        # If the file has no data rows (only header), mark as Negative.
        if df.empty or len(df) == 0:
            logging.debug("TSV file exists but has no data rows; marking sample as Negative.")
            data = {
                "Sample": [sample_id],
                "VID": [""],
                "Variant": ["Negative"],
                "NumberOfSupportingReads": [float('nan')],
                "MeanCoverage": [float('nan')],
                "Pvalue": [float('nan')],
                "AdvntrState": ["Negative"]
            }
            return pd.DataFrame(data, columns=out_cols), True
        else:
            df["Sample"] = sample_id
            df["AdvntrState"] = "Positive"
            df = df[out_cols]
            return df, True
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse adVNTR result file: {e}")
        data = {
            "Sample": [sample_id],
            "VID": [""],
            "Variant": ["Negative"],
            "NumberOfSupportingReads": [float('nan')],
            "MeanCoverage": [float('nan')],
            "Pvalue": [float('nan')],
            "AdvntrState": ["Negative"]
        }
        return pd.DataFrame(data, columns=out_cols), True
    except Exception as e:
        logging.error(f"Unexpected error occurred while loading adVNTR results: {e}")
        data = {
            "Sample": [sample_id],
            "VID": [""],
            "Variant": ["Negative"],
            "NumberOfSupportingReads": [float('nan')],
            "MeanCoverage": [float('nan')],
            "Pvalue": [float('nan')],
            "AdvntrState": ["Negative"]
        }
        return pd.DataFrame(data, columns=out_cols), True


def find_results_files(root_dir, filenames):
    """
    Recursively find all specified files in a directory tree.

    This function searches recursively through all subdirectories for the specified filenames.

    Parameters
    ----------
    root_dir : str or Path
        Root directory to start the search.
    filenames : list of str
        List of filenames to search for.

    Returns
    -------
    list
        List of Path objects matching any of the filenames.
    """
    logging.info(f"Searching for {filenames} in directory {root_dir}")
    result_files = []
    for root, _, files in os.walk(root_dir):
        logging.debug(f"Checking directory: {root}")
        for filename in filenames:
            if filename in files:
                file_path = Path(root) / filename
                logging.info(f"Found {filename} in {file_path}")
                result_files.append(file_path)
    return result_files


def load_results_from_dirs(input_dirs, filenames, file_loader):
    """
    Load results from all matching files across multiple directories.

    If loading adVNTR and no files are found in a given directory, this function
    attempts to return one row per subfolder, treating each as a sample. If there
    are no subfolders, it falls back to one row for the directory itself.

    Parameters
    ----------
    input_dirs : list
        List of directories to search.
    filenames : list of str
        List of result filenames to search for.
    file_loader : function
        Function to load and process each file.

    Returns
    -------
    pandas.DataFrame
        Concatenated DataFrame of all loaded results.
    """
    dfs = []  # Ensure dfs is defined
    loading_advntr = any('advntr' in filename.lower() for filename in filenames)
    if loading_advntr:
        cols = ["VID", "Variant", "NumberOfSupportingReads", "MeanCoverage", "Pvalue"]
        def empty_advntr_row(sample_id):
            data = {
                "Sample": [sample_id],
                "VID": [""],
                "Variant": ["Not Performed"],
                "NumberOfSupportingReads": [float('nan')],
                "MeanCoverage": [float('nan')],
                "Pvalue": [float('nan')],
                "AdvntrState": ["Not Performed"]
            }
            return pd.DataFrame(data, columns=["Sample"] + cols + ["AdvntrState"])

    for input_dir in input_dirs:
        logging.info(f"Processing input directory: {input_dir}")
        result_files = find_results_files(input_dir, filenames)
        if not result_files:
            if loading_advntr:
                has_kestrel = (Path(input_dir) / 'kestrel').is_dir()
                if has_kestrel:
                    sample_id = Path(input_dir).name
                    dfs.append(empty_advntr_row(sample_id))
                else:
                    logging.debug(f"Directory {input_dir} does not appear to be a sample directory. Skipping.")
            else:
                sample_id = Path(input_dir).name
                dfs.append(pd.DataFrame({
                    'Sample': [sample_id],
                    'Message': [f'No {filenames} results found for this sample']
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
        logging.info(f"Loaded {len(result_files)} files for {filenames} from {input_dir}")
    if not dfs:
        logging.warning(f"No data found at all for {filenames}. Returning empty DataFrame.")
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
    try:
        with open(image_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"
    except Exception as e:
        logging.error(f"Failed to encode image {image_path}: {e}")
        return ""


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
    if sum(values) == 0:
        logging.warning(f"No data to plot for donut chart '{title}'.")
        return ""
    if interactive:
        fig = go.Figure(go.Pie(
            labels=labels,
            values=values,
            hole=0.6,
            marker=dict(colors=colors, line=dict(color='black', width=2)),
            textinfo='none'
        ))
        fig.update_layout(
            title={'text': title, 'y': 0.95, 'x': 0.5, 'xanchor': 'center', 'yanchor': 'top'},
            annotations=[dict(text=f'<b>{total}</b>', x=0.5, y=0.5, font_size=40, showarrow=False)],
            showlegend=False,
            margin=dict(t=50, b=50, l=50, r=50),
            height=500,
            width=500
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

    This function creates summary statistics, generates donut charts for each
    data type, and then renders a Jinja2 template to produce a final HTML
    report.

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

    # Count advntr states based on the "AdvntrState" column.
    advntr_positive = len(advntr_df[advntr_df["AdvntrState"] == "Positive"]) if "AdvntrState" in advntr_df.columns else 0
    advntr_negative = len(advntr_df[advntr_df["AdvntrState"] == "Negative"]) if "AdvntrState" in advntr_df.columns else 0
    advntr_not_performed = len(advntr_df[advntr_df["AdvntrState"] == "Not Performed"]) if "AdvntrState" in advntr_df.columns else 0

    total_advntr = advntr_positive + advntr_negative + advntr_not_performed

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
        values=[advntr_positive, advntr_negative, advntr_not_performed],
        labels=['Positive', 'Negative', 'Not Performed'],
        total=total_advntr,
        title='adVNTR Results',
        colors=[colors['positive'], colors['negative'], colors['no_data']],
        plot_path=advntr_plot_path,
        interactive=False
    )
    advntr_plot_html = generate_donut_chart(
        values=[advntr_positive, advntr_negative, advntr_not_performed],
        labels=['Positive', 'Negative', 'Not Performed'],
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
        'kestrel_positive': kestrel_df.to_html(classes='table table-bordered table-striped', index=False, escape=False),
        'advntr_positive': advntr_df.to_html(classes='table table-bordered table-striped', index=False, escape=False),
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


def aggregate_cohort(input_paths, output_dir, summary_file, config):
    """
    Aggregate outputs from multiple runs into a single summary file.

    This function processes each input path, which can be either a directory
    or a zip file. Zip files are extracted to temporary directories for
    processing. Only the top-level sample directories within extracted zip
    files are processed to prevent internal subdirectories from being treated
    as separate samples.

    Parameters
    ----------
    input_paths : list
        List of directories or zip files containing output files to aggregate.
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
    temp_dirs = []
    processed_dirs = []

    try:
        for path_str in input_paths:
            path = Path(path_str)
            if not path.exists():
                logging.warning(f"Input path does not exist and will be skipped: {path}")
                continue
            if path.is_dir():
                logging.info(f"Adding directory to processing list: {path}")
                processed_dirs.append(path)
            elif zipfile.is_zipfile(path):
                logging.info(f"Extracting zip file: {path}")
                temp_dir = tempfile.mkdtemp(prefix="cohort_zip_")
                try:
                    with zipfile.ZipFile(path, 'r') as zip_ref:
                        zip_ref.extractall(temp_dir)
                    temp_path = Path(temp_dir)
                    logging.info(f"Extracted zip file to temporary directory: {temp_path}")

                    # **Modification Starts Here**
                    # Identify sample directories within the extracted zip
                    # Assuming that each sample directory contains both 'kestrel/' and 'advntr/' subdirectories
                    sample_dirs = [
                        p for p in temp_path.iterdir()
                        if p.is_dir() and (p / 'kestrel').is_dir() and (p / 'advntr').is_dir()
                    ]

                    if not sample_dirs:
                        logging.warning(f"No sample directories found in extracted zip file: {path}")
                        if (temp_path / 'kestrel').is_dir() and (temp_path / 'advntr').is_dir():
                            processed_dirs.append(temp_path)
                            logging.info(f"Added top-level directory as a single sample from zip file: {path}")
                        else:
                            logging.error(f"Extracted zip file does not contain valid sample directories: {path}")
                    else:
                        processed_dirs.extend(sample_dirs)
                        logging.info(f"Added {len(sample_dirs)} sample directories from zip file: {path}")
                    # **Modification Ends Here**

                    temp_dirs.append(temp_dir)
                except zipfile.BadZipFile as e:
                    logging.error(f"Bad zip file {path}: {e}")
                    shutil.rmtree(temp_dir)
                except Exception as e:
                    logging.error(f"Error extracting zip file {path}: {e}")
                    shutil.rmtree(temp_dir)
            else:
                logging.warning(f"Unsupported file type (not a directory or zip): {path}")

        if not processed_dirs:
            logging.error("No valid input directories or zip files found for cohort aggregation.")
            return

        # Load Kestrel results
        kestrel_df = load_results_from_dirs(
            input_dirs=processed_dirs,
            filenames=["kestrel_result.tsv"],
            file_loader=load_kestrel_results
        )

        # Load adVNTR results with the correct filename (only the primary TSV is used)
        advntr_filenames = ["output_adVNTR_result.tsv"]
        advntr_df = load_results_from_dirs(
            input_dirs=processed_dirs,
            filenames=advntr_filenames,
            file_loader=load_advntr_results
        )

        # Generate the cohort summary report
        generate_cohort_summary_report(
            output_dir=output_dir,
            kestrel_df=kestrel_df,
            advntr_df=advntr_df,
            summary_file=summary_file,
            config=config
        )

    finally:
        # Cleanup temporary directories
        for temp_dir in temp_dirs:
            try:
                shutil.rmtree(temp_dir)
                logging.debug(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logging.error(f"Failed to remove temporary directory {temp_dir}: {e}")
