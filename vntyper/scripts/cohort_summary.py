#!/usr/bin/env python3
"""
vntyper/scripts/cohort_summary.py

This module aggregates outputs from multiple runs into a single cohort summary report.
It exclusively loads the pipeline_summary.json from each sample directory (found at the top level
or in subfolders) to construct the cohort tables, donut plots, and additional statistics
(including runtimes, coverage, versions, assembly, and pipeline information).

Note: This module no longer defines its own CLI parser as these are now defined in the main CLI script.
"""

import base64
import hashlib
import json
import logging
import os
import shutil
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader

matplotlib.use("Agg")


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
            encoded_string = base64.b64encode(image_file.read()).decode("utf-8")
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
        fig = go.Figure(
            go.Pie(
                labels=labels,
                values=values,
                hole=0.6,
                marker=dict(colors=colors, line=dict(color="black", width=2)),
                textinfo="none",
            )
        )
        fig.update_layout(
            title={
                "text": title,
                "y": 0.95,
                "x": 0.5,
                "xanchor": "center",
                "yanchor": "top",
            },
            annotations=[dict(text=f"<b>{total}</b>", x=0.5, y=0.5, font_size=40, showarrow=False)],
            showlegend=False,
            margin=dict(t=50, b=50, l=50, r=50),
            height=500,
            width=500,
        )
        return pio.to_html(fig, full_html=False)
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        wedgeprops = {"width": 0.3, "edgecolor": "black", "linewidth": 2}
        try:
            ax.pie(
                values,
                wedgeprops=wedgeprops,
                startangle=90,
                colors=colors,
                labels=labels,
            )
            ax.text(0, 0, f"{total}", ha="center", va="center", fontsize=24)
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
        with open(config_path) as f:
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
        df (pandas.DataFrame): DataFrame containing the results (generally a single row).
        logic_config (dict): Configuration dictionary with rules.

    Returns:
        str: The computed algorithm result (e.g., 'High_Precision', 'positive', etc.).
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


# --------------------------------------------------------------------------
# New helper functions to correctly aggregate row-level results into
# sample-level categories (Positive, Positive_Flagged, Negative).
# --------------------------------------------------------------------------
def unify_kestrel_result(row_result):
    """
    Convert a row-level Kestrel result (e.g. 'High_Precision', 'Low_Precision_flagged')
    into a broader category: 'Positive', 'Positive_Flagged', or 'Negative'.
    """
    if row_result in ["High_Precision", "Low_Precision"]:
        return "Positive"
    elif row_result in ["High_Precision_flagged", "Low_Precision_flagged"]:
        return "Positive_Flagged"
    else:
        return "Negative"


def unify_advntr_result(row_result):
    """
    Convert a row-level adVNTR result (e.g. 'positive', 'positive flagged')
    into a broader category: 'Positive', 'Positive_Flagged', or 'Negative'.
    """
    if row_result == "positive":
        return "Positive"
    elif row_result == "positive flagged":
        return "Positive_Flagged"
    else:
        return "Negative"


def aggregate_sample_category(results):
    """
    Given a list of final row-level categories for a sample (each in
    {'Positive', 'Positive_Flagged', 'Negative'}), pick the highest category
    following the rule:
      - If there's at least one 'Positive' => 'Positive'
      - Else if there's at least one 'Positive_Flagged' => 'Positive_Flagged'
      - Else => 'Negative'
    """
    if any(r == "Positive" for r in results):
        return "Positive"
    elif any(r == "Positive_Flagged" for r in results):
        return "Positive_Flagged"
    else:
        return "Negative"


def load_pipeline_summary_for_sample(sample_dir):
    """
    Load the pipeline_summary.json from a sample directory and extract Kestrel,
    adVNTR data and additional statistics (runtime, coverage, version, assembly, pipeline).

    For the adVNTR step, the algorithm result will later be computed based on the logic
    defined in report_config.json.

    Parameters
    ----------
    sample_dir : str or Path
        Directory containing the pipeline_summary.json file.

    Returns
    -------
    tuple
        Three elements: (kestrel_data, advntr_data, additional_stats).
        additional_stats is a dict with keys:
          - runtime: pipeline run duration (in seconds)
          - version: pipeline version
          - assembly: assembly text from BAM Header Parsing
          - pipeline: alignment pipeline from BAM Header Parsing
          - coverage: coverage metrics dict (mean, median, stdev, min, max)
    """
    sample_dir = Path(sample_dir)
    summary_path = sample_dir / "pipeline_summary.json"
    if not summary_path.exists():
        logging.warning(f"Pipeline summary file not found in {sample_dir}")
        return [], [], {}
    try:
        with open(summary_path) as f:
            summary = json.load(f)
        kestrel_data = []
        advntr_data = []
        additional_stats = {}

        # Compute runtime from top-level timestamps if available
        pipeline_start = summary.get("pipeline_start")
        pipeline_end = summary.get("pipeline_end")
        if pipeline_start and pipeline_end:
            start_dt = datetime.fromisoformat(pipeline_start)
            end_dt = datetime.fromisoformat(pipeline_end)
            runtime_sec = (end_dt - start_dt).total_seconds()
            additional_stats["runtime"] = f"{runtime_sec:.2f} seconds"
        else:
            additional_stats["runtime"] = "N/A"

        # Pipeline version from top-level field
        additional_stats["version"] = summary.get("version", "N/A")

        # Initialize defaults for assembly, pipeline and coverage
        additional_stats["assembly"] = "N/A"
        additional_stats["pipeline"] = "N/A"
        additional_stats["coverage"] = {}

        for step in summary.get("steps", []):
            if step.get("step") == "Kestrel Genotyping":
                kestrel_data = step.get("parsed_result", {}).get("data", [])
            elif step.get("step") == "adVNTR Genotyping":
                advntr_data = step.get("parsed_result", {}).get("data", [])
            elif step.get("step") == "BAM Header Parsing":
                parsed = step.get("parsed_result", {})
                additional_stats["assembly"] = parsed.get("assembly_text", "N/A")
                additional_stats["pipeline"] = parsed.get("alignment_pipeline", "N/A")
            elif step.get("step") == "Coverage Calculation":
                parsed = step.get("parsed_result", {})
                data_list = parsed.get("data", [])
                if data_list:
                    additional_stats["coverage"] = data_list[0]
        return kestrel_data, advntr_data, additional_stats
    except Exception as e:
        logging.error(f"Error loading pipeline summary from {sample_dir}: {e}")
        return [], [], {}


def generate_cohort_summary_report(output_dir, kestrel_df, advntr_df, summary_file, config, additional_stats_html=""):
    """
    Generate the cohort summary report combining Kestrel and adVNTR results along with
    additional statistics (runtimes, coverage, version, assembly, and pipeline).

    This function creates summary statistics, generates donut charts for each
    data type, and then renders a Jinja2 template to produce a final HTML report.

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
    additional_stats_html : str, optional
        HTML table string containing additional statistics.

    Returns
    -------
    None
        Writes the HTML report to the specified summary file.
    """
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load report-specific config to get algorithm logic
    report_cfg = load_report_config()
    kestrel_logic = report_cfg.get("algorithm_logic", {}).get("kestrel", {})
    advntr_logic = report_cfg.get("algorithm_logic", {}).get("advntr", {})

    # -----------------------------
    # Compute sample-level results
    # -----------------------------
    # For each row in Kestrel, compute row-level result, unify it, then group by sample.
    if not kestrel_df.empty and "Sample" in kestrel_df.columns:
        kestrel_df["__row_result"] = kestrel_df.apply(
            lambda row: compute_algorithm_result(pd.DataFrame([row]), kestrel_logic),
            axis=1,
        )
        kestrel_df["__unified"] = kestrel_df["__row_result"].apply(unify_kestrel_result)
        kestrel_sample_results = kestrel_df.groupby("Sample")["__unified"].apply(list).apply(aggregate_sample_category)
    else:
        kestrel_sample_results = pd.Series(dtype=str)

    # For each row in adVNTR, compute row-level result, unify it, then group by sample.
    if not advntr_df.empty and "Sample" in advntr_df.columns:
        advntr_df["__row_result"] = advntr_df.apply(
            lambda row: compute_algorithm_result(pd.DataFrame([row]), advntr_logic),
            axis=1,
        )
        advntr_df["__unified"] = advntr_df["__row_result"].apply(unify_advntr_result)
        advntr_sample_results = advntr_df.groupby("Sample")["__unified"].apply(list).apply(aggregate_sample_category)
    else:
        advntr_sample_results = pd.Series(dtype=str)

    # -------------------------
    # Count final sample-level
    # -------------------------
    # Kestrel
    kestrel_counts = kestrel_sample_results.value_counts()
    k_pos = kestrel_counts.get("Positive", 0)
    k_pos_flag = kestrel_counts.get("Positive_Flagged", 0)
    k_neg = kestrel_counts.get("Negative", 0)
    total_kestrel = k_pos + k_pos_flag + k_neg

    # adVNTR
    advntr_counts = advntr_sample_results.value_counts()
    a_pos = advntr_counts.get("Positive", 0)
    a_pos_flag = advntr_counts.get("Positive_Flagged", 0)
    a_neg = advntr_counts.get("Negative", 0)
    total_advntr = a_pos + a_pos_flag + a_neg

    # --------------------------------------------------------------------
    # Updated color assignments: Positive=Blue, Flagged=Orange, Negative=Dark Grey
    # --------------------------------------------------------------------
    color_list = ["#FF0000", "#FFA500", "#404040"]  # Exactly 3 colors

    # Generate Kestrel donut chart with 3 categories
    kestrel_plot_path = plots_dir / "kestrel_summary_plot.png"
    kestrel_plot_base64 = generate_donut_chart(
        values=[k_pos, k_pos_flag, k_neg],
        labels=["Positive", "Positive (Flagged)", "Negative"],
        total=total_kestrel,
        title="Kestrel Results",
        colors=color_list,
        plot_path=kestrel_plot_path,
        interactive=False,
    )
    kestrel_plot_html = generate_donut_chart(
        values=[k_pos, k_pos_flag, k_neg],
        labels=["Positive", "Positive (Flagged)", "Negative"],
        total=total_kestrel,
        title="Kestrel Results",
        colors=color_list,
        plot_path=None,
        interactive=True,
    )

    # Generate adVNTR donut chart with 3 categories
    advntr_plot_path = plots_dir / "advntr_summary_plot.png"
    advntr_plot_base64 = generate_donut_chart(
        values=[a_pos, a_pos_flag, a_neg],
        labels=["Positive", "Positive (Flagged)", "Negative"],
        total=total_advntr,
        title="adVNTR Results",
        colors=color_list,
        plot_path=advntr_plot_path,
        interactive=False,
    )
    advntr_plot_html = generate_donut_chart(
        values=[a_pos, a_pos_flag, a_neg],
        labels=["Positive", "Positive (Flagged)", "Negative"],
        total=total_advntr,
        title="adVNTR Results",
        colors=color_list,
        plot_path=None,
        interactive=True,
    )

    # Create a separate copy for HTML formatting so that machine-readable outputs remain plain.
    kestrel_df_html = kestrel_df.copy()
    if "Confidence" in kestrel_df_html.columns:
        kestrel_df_html["Confidence"] = kestrel_df_html["Confidence"].apply(
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

    # Reorder Kestrel DataFrame columns: place Sample first then the remaining columns.
    desired_kestrel_cols = [
        "Sample",
        "Motif",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Motif_sequence",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "Depth_Score",
        "Confidence",
        "Flag",
    ]
    kestrel_columns = [col for col in desired_kestrel_cols if col in kestrel_df_html.columns]
    kestrel_html = kestrel_df_html[kestrel_columns].to_html(
        classes="table table-bordered table-striped hover compact order-column table-sm",
        index=False,
        escape=False,
    )

    # Reorder advntr DataFrame columns: ensure Sample is first.
    desired_advntr_cols = [
        "Sample",
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
    advntr_columns = [col for col in desired_advntr_cols if col in advntr_df.columns]
    advntr_html = advntr_df[advntr_columns].to_html(
        classes="table table-bordered table-striped hover compact order-column table-sm",
        index=False,
        escape=False,
    )

    template_dir = config.get("paths", {}).get("template_dir", "vntyper/templates")
    env = Environment(loader=FileSystemLoader(template_dir))
    try:
        template = env.get_template("cohort_summary_template.html")
    except Exception as e:
        logging.error(f"Failed to load Jinja2 template: {e}")
        raise

    context = {
        "report_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "kestrel_positive": kestrel_html,
        "advntr_positive": advntr_html,
        "kestrel_plot_base64": kestrel_plot_base64,
        "advntr_plot_base64": advntr_plot_base64,
        "kestrel_plot_interactive": kestrel_plot_html,
        "advntr_plot_interactive": advntr_plot_html,
        "additional_stats": additional_stats_html,
        "interactive": True,
    }

    try:
        rendered_html = template.render(context)
    except Exception as e:
        logging.error(f"Failed to render the cohort summary template: {e}")
        raise

    report_file_path = Path(output_dir) / summary_file
    try:
        with open(report_file_path, "w") as f:
            f.write(rendered_html)
        logging.info(f"Cohort summary report generated and saved to {report_file_path}")
    except Exception as e:
        logging.error(f"Failed to write the cohort summary report: {e}")
        raise


def aggregate_cohort(
    input_paths,
    output_dir,
    summary_file,
    config,
    additional_formats="",
    pseudonymize_samples=False,
):
    """
    Aggregate outputs from multiple runs into a single summary file.

    This function processes each input path, which can be either a directory or a zip file.
    Zip files are extracted to temporary directories for processing.
    Instead of parsing individual result files, this version exclusively loads the pipeline_summary.json
    from each sample directory (found either at the top level or recursively in subfolders)
    to construct the cohort tables, donut plots, and additional statistics.

    Additionally, if additional output formats are specified, the aggregated cohort
    data for Kestrel and adVNTR are exported as CSV, TSV, and/or JSON files.

    Additional Parameters
    ---------------------
    pseudonymize_samples : bool or str, optional
        If a value is provided, pseudonymize sample names using the given value as the prefix.
        The pseudonym will be the prefix concatenated with the first 5 characters of the MD5 hash
        of the original sample name.

    Parameters
    ----------
    input_paths : list
        List of directories or zip files containing output files to aggregate.
    output_dir : str or Path
        Output directory for the aggregated summary report.
    summary_file : str
        Name of the cohort summary report file.
    config : dict
        Configuration dictionary containing paths and settings.
    additional_formats : str, optional
        Comma-separated list of additional output formats to generate
        (supported: csv, tsv, json). HTML is always generated.

    Returns
    -------
    None
        Writes the cohort summary report to the specified output directory.
    """
    temp_dirs = []
    processed_dirs = set()  # use a set to avoid duplicate directories
    additional_stats_list = []

    # Identify valid directories/zip files (no changes here)
    for path_str in input_paths:
        path = Path(path_str)
        if not path.exists():
            logging.warning(f"Input path does not exist and will be skipped: {path}")
            continue
        if path.is_dir():
            if (path / "pipeline_summary.json").exists():
                logging.info(f"Found pipeline_summary.json in {path}")
                processed_dirs.add(path)
            else:
                found = False
                for summary_file_path in path.rglob("pipeline_summary.json"):
                    sample_dir = summary_file_path.parent
                    logging.info(f"Found pipeline_summary.json in {sample_dir}")
                    processed_dirs.add(sample_dir)
                    found = True
                if not found:
                    logging.warning(f"No pipeline_summary.json found in directory {path}")
        elif zipfile.is_zipfile(path):
            logging.info(f"Extracting zip file: {path}")
            temp_dir = tempfile.mkdtemp(prefix="cohort_zip_")
            try:
                with zipfile.ZipFile(path, "r") as zip_ref:
                    zip_ref.extractall(temp_dir)
                temp_path = Path(temp_dir)
                if (temp_path / "pipeline_summary.json").exists():
                    logging.info(f"Found pipeline_summary.json in {temp_path}")
                    processed_dirs.add(temp_path)
                else:
                    found = False
                    for summary_file_path in temp_path.rglob("pipeline_summary.json"):
                        sample_dir = summary_file_path.parent
                        logging.info(f"Found pipeline_summary.json in {sample_dir}")
                        processed_dirs.add(sample_dir)
                        found = True
                    if not found:
                        logging.warning(f"No pipeline_summary.json found in extracted zip file: {path}")
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

    # If pseudonymization is requested, build a mapping from original to pseudonym names.
    sample_mapping = {}

    kestrel_list = []
    advntr_list = []
    for sample_dir in processed_dirs:
        original_sample = Path(sample_dir).name
        if pseudonymize_samples:
            # Compute MD5 hash of the original sample name and take first 5 characters.
            hash_suffix = hashlib.md5(original_sample.encode()).hexdigest()[:5]
            pseudonym = f"{pseudonymize_samples}{hash_suffix}"
            sample_mapping[pseudonym] = original_sample
        else:
            pseudonym = original_sample

        logging.info(f"Processing sample directory: {sample_dir} as {pseudonym}")
        k_data, a_data, add_stats = load_pipeline_summary_for_sample(sample_dir)
        if k_data:
            for entry in k_data:
                entry["Sample"] = pseudonym
            kestrel_list.extend(k_data)
        else:
            logging.warning(f"No Kestrel data found in pipeline summary for sample {original_sample}.")
        if a_data:
            for entry in a_data:
                entry["Sample"] = pseudonym
            advntr_list.extend(a_data)
        else:
            logging.warning(f"No adVNTR data found in pipeline summary for sample {original_sample}.")
        if add_stats:
            add_stats["Sample"] = pseudonym
            additional_stats_list.append(add_stats)

    if kestrel_list:
        kestrel_df = pd.DataFrame(kestrel_list)
    else:
        logging.warning("No Kestrel data found in any sample.")
        kestrel_df = pd.DataFrame()
    if advntr_list:
        advntr_df = pd.DataFrame(advntr_list)
    else:
        logging.warning("No adVNTR data found in any sample.")
        advntr_df = pd.DataFrame()

    # Create additional statistics DataFrame and HTML table if any stats were gathered.
    if additional_stats_list:
        additional_stats_df = pd.DataFrame(additional_stats_list)
        # For coverage, flatten the dict (if available)
        if "coverage" in additional_stats_df.columns:
            coverage_df = additional_stats_df["coverage"].apply(pd.Series)
            coverage_df = coverage_df.add_prefix("cov_")
            additional_stats_df = pd.concat([additional_stats_df.drop(columns=["coverage"]), coverage_df], axis=1)
        # Reorder columns to place "Sample" first if it exists
        if "Sample" in additional_stats_df.columns:
            cols = ["Sample"] + [col for col in additional_stats_df.columns if col != "Sample"]
            additional_stats_df = additional_stats_df[cols]
        additional_stats_html = additional_stats_df.to_html(
            classes="table table-bordered table-striped hover compact order-column table-sm",
            index=False,
            escape=False,
        )
    else:
        additional_stats_html = ""

    generate_cohort_summary_report(
        output_dir=output_dir,
        kestrel_df=kestrel_df,
        advntr_df=advntr_df,
        summary_file=summary_file,
        config=config,
        additional_stats_html=additional_stats_html,
    )

    for temp_dir in temp_dirs:
        try:
            shutil.rmtree(temp_dir)
            logging.debug(f"Cleaned up temporary directory: {temp_dir}")
        except Exception as e:
            logging.error(f"Failed to remove temporary directory {temp_dir}: {e}")

    # Generate additional machine-readable cohort summaries if requested
    if additional_formats:
        formats = [fmt.strip().lower() for fmt in additional_formats.split(",") if fmt.strip()]
        if not kestrel_df.empty:
            if "csv" in formats:
                csv_path = Path(output_dir) / "cohort_kestrel.csv"
                kestrel_df.to_csv(csv_path, index=False)
                logging.info(f"Cohort Kestrel CSV written to: {csv_path}")
            if "tsv" in formats:
                tsv_path = Path(output_dir) / "cohort_kestrel.tsv"
                kestrel_df.to_csv(tsv_path, sep="\t", index=False)
                logging.info(f"Cohort Kestrel TSV written to: {tsv_path}")
            if "json" in formats:
                json_path = Path(output_dir) / "cohort_kestrel.json"
                kestrel_df.to_json(json_path, orient="records", indent=4)
                logging.info(f"Cohort Kestrel JSON written to: {json_path}")
        if not advntr_df.empty:
            if "csv" in formats:
                csv_path = Path(output_dir) / "cohort_advntr.csv"
                advntr_df.to_csv(csv_path, index=False)
                logging.info(f"Cohort adVNTR CSV written to: {csv_path}")
            if "tsv" in formats:
                tsv_path = Path(output_dir) / "cohort_advntr.tsv"
                advntr_df.to_csv(tsv_path, sep="\t", index=False)
                logging.info(f"Cohort adVNTR TSV written to: {tsv_path}")
            if "json" in formats:
                json_path = Path(output_dir) / "cohort_advntr.json"
                advntr_df.to_json(json_path, orient="records", indent=4)
                logging.info(f"Cohort adVNTR JSON written to: {json_path}")

    # If pseudonymization was enabled, output the pseudonymization table.
    if pseudonymize_samples and sample_mapping:
        pseudonym_table_path = Path(output_dir) / "pseudonymization_table.tsv"
        try:
            with open(pseudonym_table_path, "w") as pt:
                pt.write("Pseudonym\tOriginal\n")
                for pseudonym, original in sample_mapping.items():
                    pt.write(f"{pseudonym}\t{original}\n")
            logging.info(f"Pseudonymization table written to: {pseudonym_table_path}")
        except Exception as e:
            logging.error(f"Failed to write pseudonymization table: {e}")
