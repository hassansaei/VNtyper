"""
vntyper/scripts/cohort_summary.py

This module aggregates outputs from multiple runs into a single cohort summary report.
It exclusively loads the pipeline_summary.json from each sample directory (found at the top level
or in subfolders) to construct the cohort tables, donut plots, and additional statistics
(including runtimes, coverage, versions, assembly, and pipeline information).

Note: This module no longer defines its own CLI parser as these are now defined in the main CLI script.
"""

import base64
import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader, select_autoescape

# The pipeline-summary step names this cohort consumes moved to cohort_inputs with the
# code that matches them. They are compared by exact string against what pipeline.py
# records and a typo silently drops a report section (AGENTS.md trap 5), so they are
# named from summary_steps there, never spelled out.
from vntyper.scripts.cohort_categories import (
    category_counts,
    sample_categories,
    unify_advntr_result,
    unify_kestrel_result,
)
from vntyper.scripts.cohort_exports import (
    parse_output_formats,
    write_cohort_frame,
    write_pseudonymization_table,
)
from vntyper.scripts.cohort_inputs import (
    DEFAULT_PSEUDONYM_ALGORITHM,
    DEFAULT_PSEUDONYM_LENGTH,
    cleanup_temp_dirs,
    discover_sample_directories,
    duplicate_identity,
    load_pipeline_summary_for_sample,
    pseudonymized_sample_name,
)
from vntyper.scripts.cohort_tables import (
    additional_stats_frame,
    advntr_table_html,
    kestrel_table_html,
    stats_table_html,
)
from vntyper.scripts.output_paths import contained_output_path

logger = logging.getLogger(__name__)

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
        logger.error(f"Failed to encode image {image_path}: {e}")
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
        logger.warning(f"No data to plot for donut chart '{title}'.")
        return ""
    if interactive:
        fig = go.Figure(
            go.Pie(
                labels=labels,
                values=values,
                hole=0.6,
                marker={"colors": colors, "line": {"color": "black", "width": 2}},
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
            annotations=[{"text": f"<b>{total}</b>", "x": 0.5, "y": 0.5, "font_size": 40, "showarrow": False}],
            showlegend=False,
            margin={"t": 50, "b": 50, "l": 50, "r": 50},
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
                logger.warning("No plot_path provided for static donut chart, chart not saved.")
        except Exception as e:
            logger.error(f"Error generating donut chart: {e}")
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
        logger.info("Loaded report config from %s", config_path)
        return report_config
    except Exception as e:
        logger.error("Failed to load report config: %s", e)
        return {}


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
    # Neither frame is modified: the reduction annotates its own copy, which is what
    # keeps its two working columns out of the exports aggregate_cohort writes from
    # these same frames afterwards. See cohort_categories.sample_categories.
    kestrel_sample_results = sample_categories(kestrel_df, kestrel_logic, unify_kestrel_result)
    advntr_sample_results = sample_categories(advntr_df, advntr_logic, unify_advntr_result)

    # -------------------------
    # Count final sample-level
    # -------------------------
    k_pos, k_pos_flag, k_neg, total_kestrel = category_counts(kestrel_sample_results)
    a_pos, a_pos_flag, a_neg, total_advntr = category_counts(advntr_sample_results)

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

    # `Confidence` is the only cell in either table that holds markup this module built
    # (a colour span). Every other column is a sample's own string - the sample name most
    # obviously, but the motif, the flag and the alleles too - so it is escaped. Both
    # tables state that per column rather than per table; see cohort_tables.
    kestrel_html = kestrel_table_html(kestrel_df)
    advntr_html = advntr_table_html(advntr_df)

    template_dir = config.get("paths", {}).get("template_dir", "vntyper/templates")
    # Autoescaping, to parity with the per-sample report (AGENTS.md trap 11): anything
    # marked `|safe` in the template must be a fragment VNtyper built, never a value read
    # from a sample. The four `|safe` uses are the two escaped tables above, the stats
    # table, and Plotly's own figure HTML.
    env = Environment(loader=FileSystemLoader(template_dir), autoescape=select_autoescape(["html", "xml"]))
    try:
        template = env.get_template("cohort_summary_template.html")
    except Exception as e:
        logger.error(f"Failed to load Jinja2 template: {e}")
        raise

    context = {
        "report_date": datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S"),
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
        logger.error(f"Failed to render the cohort summary template: {e}")
        raise

    # `--summary-file` is documented as a name; this is what makes that true.
    report_file_path = contained_output_path(output_dir, summary_file, "--summary-file")
    try:
        with open(report_file_path, "w") as f:
            f.write(rendered_html)
        logger.info(f"Cohort summary report generated and saved to {report_file_path}")
    except Exception as e:
        logger.error(f"Failed to write the cohort summary report: {e}")
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
        The pseudonym is the prefix followed by the leading hex characters of the digest of
        the original sample name; the digest and its width are read from
        ``config["cohort"]["pseudonym"]`` and default to 12 characters of SHA-256.

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

    Raises
    ------
    ValueError
        If two discovered samples would be reported under one name - either because two
        directories share an identity, or because two identities share a pseudonym. A
        cohort that silently merges two patients' genotypes is worse than one that refuses
        to run (#206).
    """
    additional_stats_list = []

    # Identify valid directories/zip files
    processed_dirs, temp_dirs = discover_sample_directories(input_paths)

    if not processed_dirs:
        cleanup_temp_dirs(temp_dirs)
        logger.error("No valid input directories or zip files found for cohort aggregation.")
        return

    # The digest and its width are configuration, not code. `--config-path` replaces the
    # whole config rather than merging it (AGENTS.md trap 2), so every read is a .get()
    # chain against the module defaults and a config predating #206 still runs.
    pseudonym_config = config.get("cohort", {}).get("pseudonym", {})
    pseudonym_algorithm = pseudonym_config.get("algorithm", DEFAULT_PSEUDONYM_ALGORITHM)
    pseudonym_length = pseudonym_config.get("digest_characters", DEFAULT_PSEUDONYM_LENGTH)

    # The reported sample must identify exactly one patient, and that has to hold before
    # any digest is taken: two discovered directories can share a basename, so their
    # identities are already equal and no digest width separates them. sample_categories()
    # groups on the reported sample, so an undetected pair is counted as one (#206).
    collision = duplicate_identity(processed_dirs)
    if collision is not None:
        first, second = collision
        # The sample directory of a zip input is an extraction directory, which tells the
        # operator nothing, so each one is named together with the input it came from.
        msg = (
            f"Duplicate cohort sample identity {first.identity!r}: "
            f"{first.directory} (from {first.origin}) and "
            f"{second.directory} (from {second.origin}) would be reported as one sample. "
            "Give them distinct directory names, or distinct recorded input files for a zipped run."
        )
        cleanup_temp_dirs(temp_dirs)
        logger.error(msg)
        raise ValueError(msg)

    # If pseudonymization is requested, build a mapping from original to pseudonym names.
    sample_mapping = {}

    kestrel_list = []
    advntr_list = []
    try:
        for sample in processed_dirs:
            sample_dir = sample.directory
            original_sample = sample.identity
            if pseudonymize_samples:
                pseudonym = pseudonymized_sample_name(
                    pseudonymize_samples,
                    original_sample,
                    algorithm=pseudonym_algorithm,
                    length=pseudonym_length,
                )
                existing = sample_mapping.get(pseudonym)
                if existing is not None:
                    # #206: sample_mapping is keyed on the pseudonym, so this used to
                    # overwrite silently - two patients' rows became one row, and
                    # sample_categories() counted one result where there were two. The
                    # identities are already known to be distinct, so a repeat here is a
                    # digest collision rather than the same sample seen twice.
                    msg = (
                        f"Pseudonym collision: {original_sample!r} and {existing!r} both map to {pseudonym!r}. "
                        "Widen cohort.pseudonym.digest_characters in the configuration and re-run."
                    )
                    logger.error(msg)
                    raise ValueError(msg)
                sample_mapping[pseudonym] = original_sample
            else:
                pseudonym = original_sample

            logger.info(f"Processing sample directory: {sample_dir} as {pseudonym}")
            k_data, a_data, add_stats = load_pipeline_summary_for_sample(sample_dir)
            if k_data:
                for entry in k_data:
                    entry["Sample"] = pseudonym
                kestrel_list.extend(k_data)
            else:
                logger.warning(f"No Kestrel data found in pipeline summary for sample {original_sample}.")
            if a_data:
                for entry in a_data:
                    entry["Sample"] = pseudonym
                advntr_list.extend(a_data)
            else:
                logger.warning(f"No adVNTR data found in pipeline summary for sample {original_sample}.")
            if add_stats:
                add_stats["Sample"] = pseudonym
                additional_stats_list.append(add_stats)

        if kestrel_list:
            kestrel_df = pd.DataFrame(kestrel_list)
        else:
            logger.warning("No Kestrel data found in any sample.")
            kestrel_df = pd.DataFrame()
        if advntr_list:
            advntr_df = pd.DataFrame(advntr_list)
        else:
            logger.warning("No adVNTR data found in any sample.")
            advntr_df = pd.DataFrame()

        # Create additional statistics DataFrame and HTML table if any stats were gathered.
        # The frame is hoisted so the HTML table and the export below are one value: a CSV
        # that could disagree with the rendered table would be worse than no CSV.
        if additional_stats_list:
            additional_stats_df = additional_stats_frame(additional_stats_list)
            additional_stats_html = stats_table_html(additional_stats_df)
        else:
            additional_stats_df = pd.DataFrame()
            additional_stats_html = ""

        generate_cohort_summary_report(
            output_dir=output_dir,
            kestrel_df=kestrel_df,
            advntr_df=advntr_df,
            summary_file=summary_file,
            config=config,
            additional_stats_html=additional_stats_html,
        )
    finally:
        # In a `finally` because the guards above and the render itself can all raise, and
        # an aborted cohort must not leave its extracted archives behind.
        cleanup_temp_dirs(temp_dirs)

    # Generate additional machine-readable cohort summaries if requested. The render
    # above leaves both frames as they were, so what goes out here is what was read out
    # of the samples' pipeline_summary.json files and nothing else; see cohort_exports.
    if additional_formats:
        formats = parse_output_formats(additional_formats)
        write_cohort_frame(kestrel_df, output_dir, "cohort_kestrel", "Kestrel", formats)
        write_cohort_frame(advntr_df, output_dir, "cohort_advntr", "adVNTR", formats)
        # #172: the statistics frame carries every cov_* column, coverage_qc among them.
        # Until now it reached the HTML table only, so no machine-readable cohort output
        # carried a coverage figure at all.
        write_cohort_frame(additional_stats_df, output_dir, "cohort_stats", "Statistics", formats)

    # If pseudonymization was enabled, output the pseudonymization table.
    if pseudonymize_samples and sample_mapping:
        write_pseudonymization_table(output_dir, sample_mapping)
