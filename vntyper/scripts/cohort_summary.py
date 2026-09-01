"""
vntyper/scripts/cohort_summary.py

This module aggregates outputs from multiple runs into a single cohort summary report.
It exclusively loads the pipeline_summary.json from each sample directory (found at the top level
or in subfolders) to construct the cohort tables, donut plots, and additional statistics
(including runtimes, coverage, versions, assembly, and pipeline information).

Note: This module no longer defines its own CLI parser as these are now defined in the main CLI script.
"""

import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

# The pipeline-summary step names this cohort consumes moved to cohort_inputs with the
# code that matches them. They are compared by exact string against what pipeline.py
# records and a typo silently drops a report section (AGENTS.md trap 5), so they are
# named from summary_steps there, never spelled out.
from vntyper.scripts.cohort_categories import (
    category_counts,
    complete_sample_categories,
    sample_categories,
    samples_without_rows,
    unify_advntr_result,
    unify_kestrel_result,
)
from vntyper.scripts.cohort_exports import (
    parse_output_formats,
    write_cohort_frame,
    write_pseudonymization_table,
)
from vntyper.scripts.cohort_inputs import (
    cleanup_temp_dirs,
    discover_sample_directories,
    duplicate_identity,
    load_pipeline_summary_for_sample,
)
from vntyper.scripts.cohort_profiles import PROFILE_EXPORT_COLUMNS, group_decision_profiles
from vntyper.scripts.cohort_pseudonyms import (
    pseudonym_settings,
    pseudonymized_sample_name,
)
from vntyper.scripts.cohort_tables import (
    additional_stats_frame,
    advntr_table_html,
    kestrel_table_html,
    stats_table_html,
)
from vntyper.scripts.nomenclature import (
    FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT,
    FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,
)
from vntyper.scripts.output_paths import contained_output_path
from vntyper.scripts.report_assets import template_search_paths
from vntyper.scripts.report_formatting import nomenclature_legend

logger = logging.getLogger(__name__)


def generate_donut_chart(values, labels, total, title, colors):
    """Generate one interactive Plotly donut as an HTML fragment.

    The fragment omits Plotly's library. The template embeds that library once for all
    rendered charts, removing one duplicate bundle (about 4.85 MB for a two-donut
    report). The former static PNG path was unreachable in the report and wrote two
    files that no rendered page referenced.

    Args:
        values: Values for each donut segment.
        labels: Labels corresponding to ``values``.
        total: Total displayed in the donut centre.
        title: Chart title.
        colors: Colors corresponding to ``values``.

    Returns:
        The Plotly HTML fragment, or an empty string when every value is zero.
    """
    if sum(values) == 0:
        logger.warning(f"No data to plot for donut chart '{title}'.")
        return ""
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
        title={"text": title, "y": 0.95, "x": 0.5, "xanchor": "center", "yanchor": "top"},
        annotations=[{"text": f"<b>{total}</b>", "x": 0.5, "y": 0.5, "font_size": 40, "showarrow": False}],
        showlegend=False,
        margin={"t": 50, "b": 50, "l": 50, "r": 50},
        height=500,
        width=500,
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


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


def generate_cohort_summary_report(
    output_dir,
    kestrel_df,
    advntr_df,
    summary_file,
    config,
    additional_stats_html="",
    sample_names=None,
    advntr_evidence_provenance=None,
    decision_profile_provenance=None,
):
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
    sample_names : sequence of str, optional
        Every sample the report knows about. Direct callers default to the sorted union
        of sample names present in the two result frames.
    advntr_evidence_provenance : sequence of mapping, optional
        Per-sample run-recorded evidence revision and assertion values.
    decision_profile_provenance : sequence of mapping, optional
        Per-sample verified decision-profile ID, revision, and SHA-256.

    Returns
    -------
    None
        Writes the HTML report to the specified summary file.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load report-specific config to get algorithm logic
    report_cfg = load_report_config()
    kestrel_logic = report_cfg.get("algorithm_logic", {}).get("kestrel", {})
    advntr_logic = report_cfg.get("algorithm_logic", {}).get("advntr", {})

    if sample_names is None:
        named: set[str] = set()
        for frame in (kestrel_df, advntr_df):
            if "Sample" in frame.columns:
                named.update(frame["Sample"].astype(str))
        sample_names = sorted(named)

    # -----------------------------
    # Compute sample-level results
    # -----------------------------
    # Neither frame is modified: the reduction annotates its own copy, which is what
    # keeps its two working columns out of the exports aggregate_cohort writes from
    # these same frames afterwards. See cohort_categories.sample_categories.
    kestrel_sample_results = complete_sample_categories(
        sample_categories(kestrel_df, kestrel_logic, unify_kestrel_result), sample_names
    )
    advntr_sample_results = complete_sample_categories(
        sample_categories(advntr_df, advntr_logic, unify_advntr_result), sample_names
    )

    # -------------------------
    # Count final sample-level
    # -------------------------
    k_pos, k_pos_flag, k_neg, k_unest, total_kestrel = category_counts(kestrel_sample_results)
    a_pos, a_pos_flag, a_neg, a_unest, total_advntr = category_counts(advntr_sample_results)

    # --------------------------------------------------------------------
    # Colors: Positive=Red, Flagged=Orange, Negative=Dark Grey, Unestablished=Light Grey
    # --------------------------------------------------------------------
    color_list = ["#FF0000", "#FFA500", "#404040", "#B0B0B0"]  # Exactly 4 colors

    profile_groups = group_decision_profiles(decision_profile_provenance or ())
    pooled_metrics_suppressed = len(profile_groups) > 1
    profile_group_context: list[dict[str, object]] = []
    if pooled_metrics_suppressed:
        kestrel_plot_html = ""
        advntr_plot_html = ""
        for group in profile_groups:
            group_kestrel = kestrel_sample_results.reindex(group.samples)
            group_advntr = advntr_sample_results.reindex(group.samples)
            gk_pos, gk_flag, gk_neg, gk_unest, gk_total = category_counts(group_kestrel)
            ga_pos, ga_flag, ga_neg, ga_unest, ga_total = category_counts(group_advntr)
            profile_group_context.append(
                {
                    "profile_id": group.profile_id,
                    "revision": group.revision,
                    "sha256": group.sha256,
                    "samples": group.samples,
                    "kestrel_plot": generate_donut_chart(
                        [gk_pos, gk_flag, gk_neg, gk_unest],
                        ["Positive", "Positive (Flagged)", "Negative", "Unestablished"],
                        gk_total,
                        "Kestrel Results",
                        color_list,
                    ),
                    "advntr_plot": generate_donut_chart(
                        [ga_pos, ga_flag, ga_neg, ga_unest],
                        ["Positive", "Positive (Flagged)", "Negative", "Unestablished"],
                        ga_total,
                        "adVNTR Results",
                        color_list,
                    ),
                }
            )
    else:
        # One profile authority permits one cohort-level decision-performance aggregate.
        kestrel_plot_html = generate_donut_chart(
            values=[k_pos, k_pos_flag, k_neg, k_unest],
            labels=["Positive", "Positive (Flagged)", "Negative", "Unestablished"],
            total=total_kestrel,
            title="Kestrel Results",
            colors=color_list,
        )
        advntr_plot_html = generate_donut_chart(
            values=[a_pos, a_pos_flag, a_neg, a_unest],
            labels=["Positive", "Positive (Flagged)", "Negative", "Unestablished"],
            total=total_advntr,
            title="adVNTR Results",
            colors=color_list,
        )
        profile_group_context = [
            {
                "profile_id": group.profile_id,
                "revision": group.revision,
                "sha256": group.sha256,
                "samples": group.samples,
                "kestrel_plot": "",
                "advntr_plot": "",
            }
            for group in profile_groups
        ]

    # `Confidence` is the only cell in either table that holds markup this module built
    # (a colour span). Every other column is a sample's own string - the sample name most
    # obviously, but the motif, the flag and the alleles too - so it is escaped. Both
    # tables state that per column rather than per table; see cohort_tables.
    kestrel_html = kestrel_table_html(kestrel_df)
    advntr_html = advntr_table_html(advntr_df)
    legend = nomenclature_legend(kestrel_df, advntr_df)
    bam_evidence_flags = {FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT, FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT}

    template_dir = config.get("paths", {}).get("template_dir")
    # Autoescaping, to parity with the per-sample report (AGENTS.md trap 11): anything
    # marked `|safe` in the template must be a fragment VNtyper built, never a value read
    # from a sample. The six `|safe` uses are the two escaped tables above, the stats
    # table and Plotly's figure/library fragments. Profile-specific figures use the
    # same VNtyper-built Plotly return values; no sample-derived value is marked safe.
    env = Environment(
        loader=FileSystemLoader(template_search_paths(template_dir, entry_template="cohort_summary_template.html")),
        autoescape=select_autoescape(["html", "xml"]),
    )
    try:
        template = env.get_template("cohort_summary_template.html")
    except Exception as e:
        logger.error(f"Failed to load Jinja2 template: {e}")
        raise

    context = {
        "report_date": datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
        "kestrel_positive": kestrel_html,
        "advntr_positive": advntr_html,
        "kestrel_plot_interactive": kestrel_plot_html,
        "advntr_plot_interactive": advntr_plot_html,
        "plotly_library": get_plotlyjs()
        if (
            kestrel_plot_html
            or advntr_plot_html
            or any(group["kestrel_plot"] or group["advntr_plot"] for group in profile_group_context)
        )
        else "",
        "kestrel_missing": samples_without_rows(kestrel_df, sample_names),
        "advntr_missing": samples_without_rows(advntr_df, sample_names),
        "nomenclature_legend": legend,
        "show_kestrel_bam_semantics": any(entry["term"] in bam_evidence_flags for entry in legend),
        "additional_stats": additional_stats_html,
        "advntr_evidence_provenance": advntr_evidence_provenance or (),
        "decision_profile_groups": profile_group_context,
        "pooled_decision_metrics_suppressed": pooled_metrics_suppressed,
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
        If two discovered samples would be reported under one name and nothing in the
        inputs separates them - either because their inputs share a name as well as their
        local values, so ``qualify_colliding_identities`` cannot tell them apart, or
        because two distinct identities share a pseudonym. A cohort that silently merges
        two patients' genotypes is worse than one that refuses to run (#206). Two samples
        that merely share a local value are *not* an error: they are qualified by the name
        of the input each came through and both are reported.
    """
    additional_stats_list = []

    # Identify valid directories/zip files
    processed_dirs, temp_dirs = discover_sample_directories(input_paths)

    # If pseudonymization is requested, build a mapping from original to pseudonym names.
    sample_mapping = {}

    kestrel_list = []
    advntr_list = []
    cohort_samples: list[str] = []
    # The `try` opens here, immediately after discovery, and not at the sample loop: every
    # zip input has already been extracted into a `tempfile.mkdtemp` directory by this
    # point, so anything raising between the two leaks one directory per archive. That was
    # not hypothetical - reading the digest settings below used to raise `AttributeError`
    # on a config carrying `"cohort": null`, from outside the old `try`.
    try:
        if not processed_dirs:
            logger.error("No valid input directories or zip files found for cohort aggregation.")
            return

        # The digest and its width are configuration, not code, and `--config-path`
        # replaces the whole config rather than merging it (AGENTS.md trap 2) - so the
        # read has to survive a hand-written document, null levels included. It lives in
        # cohort_pseudonyms with the defaults it falls back to.
        pseudonym_algorithm, pseudonym_length = pseudonym_settings(config)

        # The reported sample must identify exactly one patient, and that has to hold
        # before any digest is taken: two discovered samples can share a name, so their
        # identities are already equal and no digest width separates them.
        # sample_categories() groups on the reported sample, so an undetected pair is
        # counted as one (#206).
        #
        # Discovery has already qualified every *shared* name with the namespace of the
        # input it came through, so what survives to here is the residue: two samples whose
        # namespace and whose local value are both equal - two archives with one stem, or
        # two directory inputs with one name. Nothing the caller supplied distinguishes
        # them, so there is no name to qualify with and the run stops.
        collision = duplicate_identity(processed_dirs)
        if collision is not None:
            first, second = collision
            # The sample directory of a zip input is an extraction directory, which tells
            # the operator nothing, so each one is named together with the input it came
            # from - and the inputs are exactly what has to change.
            msg = (
                f"Duplicate cohort sample identity {first.identity!r}: "
                f"{first.directory} (from {first.origin}) and "
                f"{second.directory} (from {second.origin}) would be reported as one sample. "
                "Their inputs share a name, so qualifying by it cannot separate them. "
                "Rename one input, or give the two runs distinct recorded input files."
            )
            logger.error(msg)
            raise ValueError(msg)

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
                    #
                    # This guard aborts where the name guard above now qualifies, and the
                    # asymmetry is the point rather than an inconsistency. A *digest*
                    # collision is between two samples that already have perfectly good
                    # distinct names: nothing about them is ambiguous, the report could
                    # name them today, and the fault is entirely in how many characters of
                    # the digest were configured - so qualifying them would be gratuitous
                    # and widening the digest is the fix. A *name* collision is the
                    # opposite: the names themselves are ambiguous, so there is nothing to
                    # widen and qualification is the only way to proceed at all. Different
                    # situations, different answers.
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
            # Discovery found a pipeline summary for this sample.  A summary that is
            # unreadable deliberately yields empty algorithm rows so the rest of the
            # cohort can render; it must still remain in both report denominators as
            # Unestablished rather than disappearing from the cohort entirely.
            cohort_samples.append(pseudonym)
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
            sample_names=cohort_samples,
            advntr_evidence_provenance=[
                {
                    "sample": stats["Sample"],
                    "revision": stats["advntr_evidence_revision"],
                    "assertion": stats["advntr_evidence_assertion"],
                }
                for stats in additional_stats_list
            ],
            decision_profile_provenance=[
                {"Sample": stats["Sample"], **{column: stats[column] for column in PROFILE_EXPORT_COLUMNS}}
                for stats in additional_stats_list
            ],
        )
    finally:
        # In a `finally` because everything above - the config read, the two identity
        # guards, the per-sample loop and the render - can raise, and an aborted cohort
        # must not leave its extracted archives behind. The `return` on an empty discovery
        # runs it too.
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
