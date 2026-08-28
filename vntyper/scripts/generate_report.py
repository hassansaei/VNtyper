"""
generate_report.py

Module Purpose:
---------------
Renders the per-sample HTML report from ``pipeline_summary.json``.

What is left here is the I/O half: reading the summary and the fastp and log
files, shelling out to ``create_report`` for the IGV panel, assembling the
template context and writing the HTML. Everything that turns a value into a
displayed string lives in :mod:`vntyper.scripts.report_formatting`, and the
screening interpretation lives in :mod:`vntyper.scripts.screening_summary` --
both extracted under AGENTS.md rule 3, which is what made either of them
testable.

Two things about this module are load-bearing and easy to break:

* The five pipeline-summary step names are matched by exact string comparison.
  A typo does not fail; it silently drops a report section (AGENTS.md trap 5),
  so they come from :mod:`vntyper.scripts.summary_steps` and are never spelled
  out here.
* The coverage field names come from
  :data:`vntyper.scripts.coverage_stats.COVERAGE_COLUMNS` (contract C1) by way
  of ``report_formatting``. They are read with ``.get(name, 0)``, so a rename
  makes the report state that an uncovered sample has full coverage.
"""

import json
import logging
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from vntyper.scripts import report_assets
from vntyper.scripts.coverage_qc import COVERAGE_QC_NOT_EVALUATED, evaluate_coverage_qc
from vntyper.scripts.cross_match_presentation import build_cross_match_summary
from vntyper.scripts.fastp_cutoffs import build_fastp_cutoffs, build_fastp_measurement
from vntyper.scripts.igv_report import extract_igv_content, run_igv_report
from vntyper.scripts.output_paths import contained_output_path
from vntyper.scripts.report_formatting import (
    ADVNTR_CELL_FORMATS,
    ADVNTR_DISPLAY_CELL_FORMATS,
    ADVNTR_DISPLAY_COLUMNS,
    ADVNTR_DISPLAY_HEADINGS,
    ADVNTR_ESSENTIAL_COLUMNS,
    EMPTY_SESSION_DICTIONARY,
    EMPTY_TABLE_JSON,
    KESTREL_DISPLAY_CELL_FORMATS,
    KESTREL_DISPLAY_COLUMNS,
    KESTREL_ESSENTIAL_COLUMNS,
    annotate_table_columns,
    confidence_html,
    drop_empty_result_rows,
    escaped_table_html,
    flag_html,
    flagged_row_count,
    folded_record_html,
    format_number_columns,
    js_json_literal,
    nomenclature_legend,
    numeric_headings,
    parse_coverage_stats,
    row_count_statement,
    select_display_columns,
    space_flag_tokens,
    summarise_fastp,
    threshold_icon,
    variant_identity,
)
from vntyper.scripts.report_identity import (
    ASSAY_NAME,
    GENE_SYMBOL,
    NOT_RECORDED,
    REPORT_TITLE_DESCRIPTION,
    REPORT_TITLE_PREFIX,
    RESEARCH_USE_STATEMENT,
    SAMPLE_NAME_EXPLICIT_KEY,
    format_input_files,
    format_region,
    format_run_timestamp,
    input_file_names,
    print_running_header_css,
    recorded_or_not,
    resolve_sample_name,
)
from vntyper.scripts.screening_summary import (
    algorithm_state_text,
    build_screening_summary,
    coverage_not_measured_note,
    coverage_qc_word,
    execution_state,
    load_report_config,
    state_chips,
    supports_subthreshold,
)

# These five names are matched by exact string comparison against what pipeline.py
# records. A typo does not fail - it silently drops a report section (AGENTS.md
# trap 5), so they are named, never spelled out.
from vntyper.scripts.subthreshold import find_note
from vntyper.scripts.summary_steps import (
    STEP_ABSENT,
    STEP_ADVNTR,
    STEP_BAM_HEADER,
    STEP_COVERAGE,
    STEP_KESTREL,
    STEP_READ,
    STEP_UNREADABLE,
    get_step_comments,
    get_step_data,
    get_step_result,
    get_step_state,
)

logger = logging.getLogger(__name__)

#: Shown when a coverage figure could not be computed.
NOT_CALCULATED = "Not calculated"

#: How each results table names itself in its visible/total row-count statement.
KESTREL_ROW_NOUN = "Kestrel"
ADVNTR_ROW_NOUN = "adVNTR"

#: A failed embedded-mode generator page is moved here before its temporary
#: directory is cleaned. The fixed single-segment name keeps the evidence inside the
#: run output and makes repeated failures converge on one predictable artifact.
IGV_FAILURE_DIAGNOSTIC_FILENAME = "igv_report.failed.html"

#: CSS classes on both per-sample results tables. Named once so the two cannot drift.
#:
#: Every Bootstrap and DataTables class is gone with the CDN tags that gave them meaning
#: (#242): ``table-bordered``, ``table-sm``, ``hover``, ``compact`` and ``order-column``
#: styled nothing this repository ships, and leaving them would have left the next reader
#: to work out which of five class names were live. Two remain and both are ours:
#: ``table`` is the hook the shared token layer's cell rules select on, and ``sortable``
#: is what the report's own script looks for when it turns column headings into buttons.
#:
#: **No ``table-striped``.** A stripe competes with the one row treatment in these tables
#: that carries meaning, the flagged-value highlight, and it was the ``#f2f2f2`` under the
#: old inline ``orange`` that measured 1.76:1. That used to be enforced only by not
#: emitting a Bootstrap class - a claim with no rule and no test behind it - so re-adding
#: it passed the whole suite. The row background is now the report's own, and
#: ``tests/unit/test_report_presentation.py`` and ``tests/browser/test_flagged_rows.py``
#: hold it to one value: one in the source, one measured in a real browser. The cohort
#: report keeps its stripes (``cohort_tables.TABLE_CLASSES``), and its own stylesheet is
#: what draws them now.
PER_SAMPLE_TABLE_CLASSES = "table sortable"


def load_pipeline_summary(summary_file_path):
    """
    Loads the pipeline summary JSON file generated by the pipeline.

    Args:
        summary_file_path (str or Path): Path to the pipeline summary file.

    Returns:
        dict: The loaded summary dictionary or an empty dict if load fails.
    """
    logger.info("Loading pipeline summary from %s", summary_file_path)
    if not os.path.exists(summary_file_path):
        logger.error("Pipeline summary file not found: %s", summary_file_path)
        return {}
    try:
        with open(summary_file_path) as f:
            summary = json.load(f)
        logger.debug("Pipeline summary loaded successfully.")
        return summary
    except Exception as e:
        logger.error("Failed to load pipeline summary: %s", e)
        return {}


def load_fastp_output(fastp_file):
    """
    Loads fastp JSON output (e.g., output.json) for summary metrics if available.
    Returns an empty dict if file not found or if parsing fails.
    """
    logger.debug("load_fastp_output called with fastp_file=%s", fastp_file)
    if not os.path.exists(fastp_file):
        logger.warning("fastp output file not found: %s", fastp_file)
        return {}
    try:
        with open(fastp_file) as f:
            data = json.load(f)
        logger.debug("fastp output successfully loaded.")
        return data
    except Exception as e:
        logger.error("Failed to load or parse fastp output: %s", e)
        return {}


def load_pipeline_log(log_file):
    """
    Loads the pipeline log content from the specified log_file.
    Returns a placeholder string if not found or on error.
    """
    logger.info("Loading pipeline log from %s", log_file)
    if not log_file:
        logger.warning("No pipeline log file provided; skipping log loading.")
        return "No pipeline log file was provided."
    if not os.path.exists(log_file):
        logger.warning("Pipeline log file not found: %s", log_file)
        return "Pipeline log file not found."
    try:
        with open(log_file) as f:
            content = f.read()
        logger.debug("Pipeline log successfully loaded.")
        return content
    except Exception as e:
        logger.error("Failed to read pipeline log file: %s", e)
        return "Failed to load pipeline log."


def build_kestrel_frames(kestrel_data):
    """Build the display and matching Kestrel frames from summary rows.

    Two frames come out of one input because they are used for different things:
    the screening summary matches on plain values, while the table shows
    colour-coded HTML. Colour-coding the frame the summary reads would make every
    ``Confidence`` comparison fail against a ``<span>``.

    A run that genotyped and called nothing is *not* a run with a result. Kestrel
    records one there anyway - ``output_empty_result`` writes a row of the literal
    string ``"None"`` with ``Confidence`` set to ``Negative``, so that
    ``kestrel_result.tsv`` has a body - and both frames used to carry it: the table
    rendered ``None None None None None None None None NaN Negative`` and the count
    line called it one Kestrel row (#242). It is dropped here, before either frame
    exists, so the display, the count and the screening match all see the same rows.

    Dropping it cannot move the screening verdict. ``compute_algorithm_result`` reads
    only the first row and every configured Kestrel rule breaks on that row's
    ``Confidence``, so the block's ``default`` - "negative" - is what it returned with
    the placeholder present, and an empty frame returns the same default by the
    shortest path.

    Args:
        kestrel_data (list[dict]): The ``Kestrel Genotyping`` step's rows.

    Returns:
        tuple[pandas.DataFrame, pandas.DataFrame]: The frame to render, and the
        unformatted frame to match on. Both are empty when there are no rows, and when
        the only row is the empty-result placeholder.
    """
    kestrel_data = drop_empty_result_rows(kestrel_data) if kestrel_data else kestrel_data
    if not kestrel_data:
        logger.warning("No Kestrel data found in pipeline summary.")
        return pd.DataFrame(), pd.DataFrame()

    frame = select_display_columns(pd.DataFrame(kestrel_data), KESTREL_DISPLAY_COLUMNS)

    if "Depth Score" in frame.columns:
        try:
            frame["Depth Score"] = pd.to_numeric(frame["Depth Score"], errors="coerce")
        except Exception as e:
            logger.warning("Could not convert 'Depth Score' to numeric: %s", e)
        frame = frame.sort_values(by="Depth Score", ascending=False)

    matching_frame = frame.copy()

    # Every displayed number is produced here, by the formatter its column declares.
    # The browser used to do it, table-agnostically and only when three CDNs
    # resolved, so the same file read differently online and offline (#242). The
    # matching frame is copied first and keeps the unformatted values, because the
    # screening rules compare against them.
    frame = space_flag_tokens(format_number_columns(frame, KESTREL_DISPLAY_CELL_FORMATS))

    if "Confidence" in frame.columns:
        frame["Confidence"] = frame["Confidence"].apply(confidence_html)
    if "Flag" in frame.columns:
        frame["Flag"] = frame["Flag"].apply(flag_html)

    # The remaining cells are escaped by `escaped_table_html`, which renders the frame -
    # once, naming `Confidence` and `Flag` as the columns whose markup we built and which
    # escaped their own values. Escaping here as well would double-escape every other
    # cell, so a motif sequence containing a `<` would reach the reader as `&lt;`.
    logger.debug("Kestrel data extracted from summary and formatted.")
    return frame, matching_frame


def build_advntr_frame(advntr_data):
    """Project the adVNTR rows onto the columns the report shows.

    Args:
        advntr_data (list[dict]): The ``adVNTR Genotyping`` step's rows.

    Returns:
        pandas.DataFrame: The frame to render, empty when there are no rows.
    """
    if not advntr_data:
        return pd.DataFrame()
    frame = pd.DataFrame(advntr_data)
    logger.debug("adVNTR data extracted from summary and formatted.")
    return frame[[col for col in ADVNTR_DISPLAY_COLUMNS if col in frame.columns]]


#: The release in which coverage statistics became region-wide (#171). A summary written
#: by an older version carries a ``mean`` over *covered* positions, which is not comparable
#: with the thresholds applied to it today.
REGION_WIDE_MEAN_SINCE = (2, 0, 8)


def _predates_region_wide_mean(version: str | None) -> bool:
    """Whether a recorded pipeline version is older than the region-wide coverage change.

    Args:
        version (str | None): The ``version`` recorded in ``pipeline_summary.json``.

    Returns:
        bool: True only when the version parses and is strictly older than
        :data:`REGION_WIDE_MEAN_SINCE`. An absent, empty or unparseable version returns
        False, so an unrecognised summary is left exactly as it is rather than having its
        numbers silently rescaled.
    """
    if not version:
        return False
    parts = version.strip().lstrip("v").split(".")[:3]
    try:
        parsed = tuple(int(part) for part in parts)
    except ValueError:
        logger.debug("Unparseable pipeline version %r; leaving coverage figures untouched.", version)
        return False
    return parsed < REGION_WIDE_MEAN_SINCE


def generate_summary_report(
    output_dir,
    template_dir,
    report_file,
    log_file,
    bed_file=None,
    bam_file=None,
    fasta_file=None,
    flanking=50,
    vcf_file=None,
    config=None,
    sample_name=None,
    report_igv=report_assets.DEFAULT_REPORT_IGV,
):
    """
    Generates a summary report based on a pipeline summary JSON file.
    Instead of parsing results from subfolders, the report now loads the output summary file (JSON)
    and renders the report based on that content. The main configuration is passed in via 'config'.
    This module additionally loads its own report-specific configuration from report_config.json.
    The header information extracted from the BAM Header Parsing step (including warning,
    alignment_pipeline, assembly_text, and assembly_contig) is included in the report context.

    Args:
        output_dir (str): Output directory for the report.
        template_dir (str | Path | None): Operator directory containing the report
            entry template, or None to use the installed package templates.
        report_file (str): Name of the report file.
        log_file (str): Path to the pipeline log file.
        bed_file (str, optional): Path to the BED file for IGV reports.
        bam_file (str, optional): Path to the BAM file for IGV reports.
        fasta_file (str, optional): Path to the reference FASTA file for IGV reports.
        flanking (int, optional): Size of the flanking region for IGV reports.
        vcf_file (str, optional): Path to the sorted/indexed VCF file.
        config (dict): Main configuration dictionary (passed from the pipeline).
        sample_name (str, optional): What the report calls its sample, in the
            title, the heading and the header block. Reachable through
            ``vntyper report --sample-name``; the pipeline passes nothing and the
            name is derived from the summary's own ``input_files`` instead - see
            :func:`~vntyper.scripts.report_identity.resolve_sample_name`.
        report_igv (str, optional): How the report carries its alignment browser -
            one of :data:`~vntyper.scripts.report_assets.REPORT_IGV_MODES`.
            ``embedded`` (the default) writes the vendored, gzipped igv.js into the
            document, so the archived file is a complete alignment browser needing
            neither a second file nor a network. ``sidecar`` leaves it out and
            points the reader at the self-contained ``igv_report.html`` beside it.
            ``off`` produces no alignment browser at all. Reachable as
            ``--report-igv`` on both ``vntyper pipeline`` and ``vntyper report``.

    Raises:
        ValueError: If config is not provided, an operator template directory lacks
            the report entry template, ``report_igv`` is not a recognised mode, or
            the vendored asset fails its digest check.
    """
    logger.debug("---- DEBUG: Entered generate_summary_report ----")
    logger.debug(
        "Called with output_dir=%s, template_dir=%s, report_file=%s",
        output_dir,
        template_dir,
        report_file,
    )
    logger.debug(
        "bed_file=%s, bam_file=%s, fasta_file=%s, flanking=%s, log_file=%s, vcf_file=%s",
        bed_file,
        bam_file,
        fasta_file,
        flanking,
        log_file,
        vcf_file,
    )

    if config is None:
        raise ValueError("Config dictionary must be provided to generate_summary_report")

    # Load the script-specific report configuration.
    report_config = load_report_config()

    # Resolve flanking region from main config if needed.
    if flanking == 50:
        flanking = config.get("default_values", {}).get("flanking", 50)
        logger.debug("Flanking region set to %s based on config.", flanking)

    thresholds = config.get("thresholds", {})
    if not isinstance(thresholds, dict):
        message = "Config thresholds must be a dictionary."
        logger.error(message)
        raise ValueError(message)
    mean_vntr_cov_threshold = thresholds.get("mean_vntr_coverage", 100)
    percent_vntr_uncovered_threshold = thresholds.get("percent_vntr_uncovered", 50.0)
    fastp_cutoffs = build_fastp_cutoffs(thresholds)

    # Load the pipeline summary JSON.
    summary_file_path = Path(output_dir) / "pipeline_summary.json"
    pipeline_summary = load_pipeline_summary(summary_file_path)

    # Extract input_files and pipeline_version from the summary.
    input_files = pipeline_summary.get("input_files", {})
    pipeline_version = pipeline_summary.get("version", "unknown")

    # Who this report is about (#242). Every report was titled "Summary Report",
    # so two of them were indistinguishable in two browser tabs and a printed one
    # carried no identity at all. Three levels: an explicit `--sample-name` here,
    # then the name the *run* recorded - the same string Kestrel embedded, so the
    # report and the VCF cannot disagree - then a derivation from `input_files`,
    # which is what every summary written before the field existed still uses.
    # `sample_name_is_explicit` says which of the two shapes the recorded name has:
    # an operator's name to print verbatim, or a `Path.stem` to finish deriving.
    resolved_sample_name = resolve_sample_name(
        sample_name,
        *input_file_names(input_files),
        recorded=pipeline_summary.get("sample_name"),
        recorded_is_explicit=pipeline_summary.get(SAMPLE_NAME_EXPLICIT_KEY),
    )
    logger.info("Report sample name resolved to %r.", resolved_sample_name)

    # How the run's reference was actually resolved (#163) - the BWA reference for
    # FASTQ, or whatever the alignment plan proved for BAM/CRAM, never a BWA path a
    # BAM/CRAM run never opened (MAJOR 5, milestone-5 PR-2 review; see
    # `pipeline_alignment.resolve_summary_reference_provenance`). A summary written
    # before this was recorded simply omits these, and the template renders nothing for
    # them - see AGENTS.md trap 5 on step names for the same "absent key, no section"
    # pattern.
    reference_assembly_requested = pipeline_summary.get("reference_assembly_requested")
    reference_key_used = pipeline_summary.get("reference_key_used")
    reference_path = pipeline_summary.get("reference_path")
    reference_source_effective = pipeline_summary.get("reference_source_effective")

    # Header info from the BAM Header Parsing step, whose parsed_result is a flat
    # object rather than {"data": [...]}.
    header_info = get_step_result(pipeline_summary, STEP_BAM_HEADER)
    # The template puts the word "Warning:" in front of this, because that word is what
    # carries the meaning for a reader who does not see the colour. The recorded string
    # begins with `WARNING: ` of its own (`fastq_bam_processing.py`), so the one
    # genuinely urgent sentence in the report rendered as "Warning: WARNING: ...".
    # Stripped here rather than in the step that records it: the summary is a stored
    # artefact that other tooling reads, and this is a presentation decision about one
    # rendering of it.
    header_warning = re.sub(r"^\s*WARNING:\s*", "", header_info.get("warning", "") or "")
    alignment_pipeline = header_info.get("alignment_pipeline", "")
    assembly_text = header_info.get("assembly_text", "")
    assembly_contig = header_info.get("assembly_contig", "")

    # VNTR coverage statistics, keyed by the frozen coverage schema (contract C1).
    coverage_rows = get_step_data(pipeline_summary, STEP_COVERAGE)
    coverage = parse_coverage_stats(coverage_rows)
    mean_vntr_coverage = coverage["mean"]
    percent_vntr_uncovered = coverage["percent_uncovered"]

    # A summary written before 2.0.8 carries no `coverage_qc` column, and its `mean` is the
    # mean over *covered* positions rather than over the region (#171) - so judging it with
    # the new rule would be too lenient: an old mean of 150 at 40% uncovered passes a 100x
    # threshold, while the region-wide mean it stands for is 90 and fails.
    #
    # The mean is recoverable exactly rather than approximately. `percent_uncovered` was
    # already correct before 2.0.8, and with S = sum(depths), c = covered, T = region:
    # (S/c) * (1 - (T-c)/T) = (S/c)(c/T) = S/T. The golden-cohort gate confirmed this
    # identity holds on 61 of 61 cases. So correct the number and judge the corrected one.
    # Two signals are required, and the default is to change nothing. A missing
    # `coverage_qc` column alone is not proof of an old summary - a hand-built or
    # third-party one may simply omit it, and silently rescaling its mean would be worse
    # than the problem being fixed. So also require the recorded version to parse to
    # something older than 2.0.8; an absent or unreadable version leaves the number alone.
    legacy_summary = (
        bool(coverage_rows) and "coverage_qc" not in coverage_rows[0] and _predates_region_wide_mean(pipeline_version)
    )
    mean_for_qc = mean_vntr_coverage
    if legacy_summary and mean_vntr_coverage is not None and percent_vntr_uncovered is not None:
        mean_for_qc = round(mean_vntr_coverage * (1 - percent_vntr_uncovered / 100), 2)
        logger.info(
            "Coverage summary predates 2.0.8; correcting its covered-bases mean %.2f to the "
            "region-wide %.2f before applying the QC rule (#171).",
            mean_vntr_coverage,
            mean_for_qc,
        )

    # Recomputed rather than read out of the stored `coverage_qc`, so `vntyper report` still
    # works against a summary written before #172. Given the same thresholds both sides
    # evaluate the *published* two-decimal figures - serialised with `:.2f` - so they agree.
    # Supplying different thresholds to `vntyper report` deliberately yields a different
    # verdict from the one stored at run time; that is the caller asking for one.
    coverage_qc = evaluate_coverage_qc(
        mean_for_qc,
        percent_vntr_uncovered,
        mean_vntr_cov_threshold,
        percent_vntr_uncovered_threshold,
    )

    # Three states, not two, for **both** algorithms. "Recorded" is not "produced a
    # readable result": when `kestrel_result.tsv` is absent, `record_step` flags the step
    # `result_file_missing` and still records it with an empty `data` list (#212) - the
    # same shape a run that genotyped and called nothing produces. Asking only whether the
    # step is present therefore renders a failed stage as a negative, which is a claim the
    # run never established.
    #
    # `get_step_state` distinguished the three from the start, but its result used to reach
    # only the Kestrel *section*: the screening computation and the whole adVNTR side still
    # asked whether the step was present, so an unreadable adVNTR stage produced "No
    # pathogenic variants identified by adVNTR" and `adVNTR: negative` while the Kestrel
    # section two headings up said "this is not a negative".
    kestrel_state = get_step_state(pipeline_summary, STEP_KESTREL)
    kestrel_df, kestrel_df_raw = build_kestrel_frames(get_step_data(pipeline_summary, STEP_KESTREL))

    # #266. The note is a `##` banner line on `kestrel_result.tsv`, so it arrives through
    # `parsed_result["comments"]` and never through `data` -- it cannot become a row here.
    #
    # `supports_subthreshold` gates *both* uses of it: the state promotion below, and the
    # rendering. A `report_config.json` written before #266 carries neither the
    # `non_finding_results` declaration nor the eight sentences that explain the state, and
    # its `negative` message says "No variant detected by either genotyping method" -- so
    # rendering the line under that configuration would make the report contradict itself.
    # One predicate for both keeps them from diverging.
    kestrel_subthreshold_note = find_note(get_step_comments(pipeline_summary, STEP_KESTREL))
    if kestrel_subthreshold_note and not supports_subthreshold(report_config):
        logger.warning(
            "A below-reporting-floor note was recorded, but this report_config.json does not "
            "declare the subthreshold state; the note is withheld so the report cannot "
            "contradict its own screening message. See issue #266."
        )
        kestrel_subthreshold_note = None

    # `advntr_available` is now "produced a result to match on", not "has a step". That is
    # what makes `advntr_result` the `NOT_PERFORMED` token rather than the block's
    # "negative" default for a stage that ran and lost its result; the execution axis
    # beside it is what keeps the two apart in everything the reader sees.
    advntr_state = get_step_state(pipeline_summary, STEP_ADVNTR)
    advntr_available = advntr_state == STEP_READ
    advntr_df = build_advntr_frame(get_step_data(pipeline_summary, STEP_ADVNTR))

    pipeline_log_content = load_pipeline_log(log_file)

    # IGV report generation (if applicable). `off` skips it outright: the mode says
    # the run wants no alignment browser, and running `create_report` anyway would
    # spend the time and write the file for nothing.
    if report_igv not in report_assets.REPORT_IGV_MODES:
        msg = f"Unknown --report-igv mode {report_igv!r}; expected one of {', '.join(report_assets.REPORT_IGV_MODES)}."
        logger.error(msg)
        raise ValueError(msg)

    temporary_igv_dir = None
    igv_operation_error: Exception | None = None
    if report_igv == report_assets.REPORT_IGV_OFF:
        logger.info("--report-igv off: no alignment browser is produced for this run.")
        igv_report_file = None
        igv_content, table_json, session_dictionary = "", "", ""
    elif bed_file and os.path.exists(bed_file):
        logger.info("Running IGV report for BED file: %s", bed_file)
        if report_igv == report_assets.REPORT_IGV_EMBEDDED:
            temporary_igv_dir = tempfile.TemporaryDirectory(prefix=".vntyper-igv-", dir=output_dir)
            igv_report_file = Path(temporary_igv_dir.name) / "igv_report.html"
        else:
            igv_report_file = Path(output_dir) / "igv_report.html"
        try:
            run_igv_report(
                bed_file,
                bam_file,
                fasta_file,
                igv_report_file,
                flanking=flanking,
                vcf_file=vcf_file,
                config=config,
                report_igv=report_igv,
            )
            if not igv_report_file.exists():
                msg = "The IGV generator completed without writing its expected report."
                logger.error(msg)
                raise ValueError(msg)

            igv_content, table_json, session_dictionary = extract_igv_content(igv_report_file)
            try:
                if not igv_content or not table_json.strip() or not session_dictionary.strip():
                    raise ValueError("one or more required fragments are empty")
                table_value = json.loads(table_json.removesuffix(";").strip())
                session_value = json.loads(session_dictionary.removesuffix(";").strip())
                if not isinstance(table_value, dict) or not isinstance(session_value, dict):
                    raise ValueError("tableJson and sessionDictionary must be JSON objects")
            except ValueError as error:
                msg = "The generated IGV report could not be validated; the summary report was not written."
                logger.error(msg)
                raise ValueError(msg) from error
            table_json = js_json_literal(table_json, EMPTY_TABLE_JSON)
            session_dictionary = js_json_literal(session_dictionary, EMPTY_SESSION_DICTIONARY)
        except Exception as error:
            igv_operation_error = error
            if temporary_igv_dir is not None and igv_report_file.exists():
                try:
                    diagnostic_file = contained_output_path(
                        output_dir,
                        IGV_FAILURE_DIAGNOSTIC_FILENAME,
                        "IGV failure diagnostic",
                    )
                    igv_report_file.replace(diagnostic_file)
                    logger.error("Preserved failed IGV generator output at %s", diagnostic_file)
                except (OSError, ValueError) as preservation_error:
                    logger.error(
                        "Could not preserve failed IGV generator output at %s: %s",
                        IGV_FAILURE_DIAGNOSTIC_FILENAME,
                        preservation_error,
                    )
            raise
        finally:
            if temporary_igv_dir is not None:
                try:
                    temporary_igv_dir.cleanup()
                except Exception as cleanup_error:
                    if igv_operation_error is None:
                        raise
                    logger.error(
                        "Could not clean temporary IGV directory after %s: %s",
                        type(igv_operation_error).__name__,
                        cleanup_error,
                    )
    else:
        logger.warning("BED file does not exist or not provided. Skipping IGV report generation.")
        igv_report_file = None
        igv_content, table_json, session_dictionary = "", "", ""

    # Whether this run has an alignment session at all. The template branches on it to
    # author the right "there is no alignment view here" sentence *in the markup*, which
    # is what makes that sentence true for a reader with scripting off - and it is also
    # what decides whether the 497 KB payload is written at all. A report with no
    # alignments carries none of it.
    session_literal = js_json_literal(session_dictionary, EMPTY_SESSION_DICTIONARY)
    igv_session_available = session_literal.strip() not in ("", EMPTY_SESSION_DICTIONARY)
    igv_payload = report_assets.igv_payload(report_igv) if igv_session_available else None
    if igv_payload is None:
        logger.info(
            "No igv.js payload written into the report (mode=%s, alignment session=%s).",
            report_igv,
            igv_session_available,
        )

    fastp = summarise_fastp(load_fastp_output(Path(output_dir) / "fastq_bam_processing/output.json"))

    coverage_icon, coverage_color = threshold_icon(mean_vntr_coverage, mean_vntr_cov_threshold, higher_better=True)
    uncovered_icon, uncovered_color = threshold_icon(
        percent_vntr_uncovered, percent_vntr_uncovered_threshold, higher_better=False
    )
    duplication_rate = build_fastp_measurement(fastp.duplication_rate, "duplication_rate")
    q20_rate = build_fastp_measurement(fastp.q20_rate, "q20_rate")
    q30_rate = build_fastp_measurement(fastp.q30_rate, "q30_rate")
    passed_filter_rate = build_fastp_measurement(fastp.passed_filter_rate, "passed_filter_rate")
    dup_icon, dup_color = threshold_icon(
        duplication_rate.value, fastp_cutoffs.duplication_rate.value, higher_better=False
    )
    q20_icon, q20_color = threshold_icon(q20_rate.value, fastp_cutoffs.q20_rate.value)
    q30_icon, q30_color = threshold_icon(q30_rate.value, fastp_cutoffs.q30_rate.value)
    pf_icon, pf_color = threshold_icon(passed_filter_rate.value, fastp_cutoffs.passed_filter_rate.value)
    duplication_rate_display = duplication_rate.display
    q20_rate_display = q20_rate.display
    q30_rate_display = q30_rate.display
    passed_filter_rate_display = passed_filter_rate.display

    # "" for an empty frame, which is what the template's authored empty states hang
    # on. This used to call `to_html` directly, bypassing the helper written for
    # exactly this: an empty frame produced a headerless, bodyless `<table>` - a stray
    # empty box under the heading - and every cell was escaped a step earlier instead.
    kestrel_html = escaped_table_html(
        kestrel_df,
        classes=PER_SAMPLE_TABLE_CLASSES,
        html_columns=("Confidence", "Flag"),
        table_id="kestrel_table",
    )
    # Per-column alignment, face and column help, added to the rendered markup rather
    # than carried in the frame: they are presentation over a column, and pandas has no
    # way to attach either to one. Every number in both results tables was centred and
    # the `.num`, `.mono-cell` and `.col-grow` rules the stylesheet has carried since
    # #242 were emitted by nothing at all.
    kestrel_html = annotate_table_columns(
        kestrel_html,
        list(kestrel_df.columns),
        numeric=numeric_headings(KESTREL_DISPLAY_CELL_FORMATS),
        essential=KESTREL_ESSENTIAL_COLUMNS,
        caption="Kestrel variant calls, with the MUC1 name reconciled from both callers",
    )
    kestrel_row_summary = (
        row_count_statement(len(kestrel_df), flagged_row_count(kestrel_df_raw), noun=KESTREL_ROW_NOUN)
        if kestrel_html
        else ""
    )
    logger.debug("Kestrel results converted to HTML.")

    advntr_row_summary = ""
    advntr_folded_record = ""
    if not advntr_available or advntr_df.empty:
        # Every adVNTR state that is not a table is worded by the template, the way the
        # Kestrel side already does it. Two of the four used to be built here instead:
        # the not-performed paragraph was unreachable - the template's
        # `{% if advntr_available and advntr_highlight %}` cannot render it, and the
        # branch beneath that guard printed its own, differently worded line - and the
        # empty-result paragraph *was* reached, and was rendered inside the bordered
        # `.table-container` as if it were a table. One state, one sentence, one place.
        advntr_html = ""
        logger.debug("adVNTR produced no rows to tabulate; the template states which state it is in.")
    else:
        # Rendered from a copy: `advntr_df` itself is what the screening rules match
        # on (`VID`, `Flag`), so it has to keep the unformatted values.
        advntr_display = space_flag_tokens(format_number_columns(advntr_df, ADVNTR_CELL_FORMATS))
        # Headings, on the copy that is rendered. Formatting runs first because
        # `ADVNTR_CELL_FORMATS` is keyed by the source names, which is also why this
        # cannot simply move into `build_advntr_frame`.
        advntr_display = advntr_display.rename(columns=ADVNTR_DISPLAY_HEADINGS)
        if "Flag" in advntr_display.columns:
            advntr_display["Flag"] = advntr_display["Flag"].apply(flag_html)
        advntr_html = annotate_table_columns(
            escaped_table_html(
                advntr_display,
                classes=PER_SAMPLE_TABLE_CLASSES,
                html_columns=("Flag",),
            ),
            list(advntr_display.columns),
            numeric=numeric_headings(ADVNTR_DISPLAY_CELL_FORMATS),
            essential=ADVNTR_ESSENTIAL_COLUMNS,
            caption="adVNTR variant calls, with the MUC1 name reconciled from both callers",
        )
        advntr_row_summary = row_count_statement(len(advntr_df), flagged_row_count(advntr_df), noun=ADVNTR_ROW_NOUN)
        advntr_folded_record = folded_record_html(advntr_display, ADVNTR_ESSENTIAL_COLUMNS, noun=ADVNTR_ROW_NOUN)
        logger.debug("adVNTR results converted to HTML.")

    cross_match_message, cross_match_is_positive, cross_match_is_assessable = build_cross_match_summary(
        pipeline_summary, report_config
    )

    # Autoescaping is on: everything reaching the report from a sample - file
    # names, BAM header fields, motif sequences, log lines - is attacker-influenced
    # and the report is a file people forward. The fragments we build ourselves are
    # marked `|safe` at their interpolation points in the template, and the two
    # results tables and the Confidence spans are escaped before they get there.
    env = Environment(
        loader=FileSystemLoader(
            report_assets.template_search_paths(template_dir, entry_template="report_template.html")
        ),
        autoescape=select_autoescape(["html"]),
    )
    try:
        template = env.get_template("report_template.html")
        logger.debug("Jinja2 template 'report_template.html' loaded successfully.")
    except Exception as e:
        logger.error("Failed to load Jinja2 template: %s", e)
        raise

    # Match on the sorted, unformatted frame: the displayed one has HTML in its
    # Confidence cells and would match nothing.
    screening = build_screening_summary(
        kestrel_df_raw,
        advntr_df,
        advntr_available,
        coverage_qc,
        report_config,
        kestrel_execution=execution_state(kestrel_state),
        advntr_execution=execution_state(advntr_state),
        kestrel_subthreshold=kestrel_subthreshold_note is not None,
    )
    logger.debug("Summary text generated: %s", screening.text)

    def shown(value):
        """Render a coverage figure, or say it was not calculated."""
        return value if value is not None else NOT_CALCULATED

    def _rounded(value):
        """Two decimal places, server-side.

        The cohort template runs `applyRounding` over every table it initialises and
        strips trailing zeroes, so a figure left to the browser renders differently
        online and offline (#242). Rounding here makes the printed value a property of
        the run rather than of the reader's network.
        """
        return None if value is None else round(value, 2)

    # The printed running header, built here so that the template interpolates one
    # opaque fragment into its `<style>` rather than five values. Every value is put
    # through the CSS escaper by `print_running_header_css`, and there is no other way
    # into that stylesheet: `tests/unit/test_report_presentation.py` fails on any other
    # expression inside a `<style>`, and on any `<style>` expression that is not this
    # one marked `|safe`.
    #
    # These three locals are computed once and used twice - here and in the context
    # below - because the margin box and the document's own header block must not be
    # able to disagree about what the report is of.
    assembly_declared_text = recorded_or_not(reference_assembly_requested)
    run_time_text = format_run_timestamp(pipeline_summary.get("pipeline_start"))
    running_header_css = print_running_header_css(
        sample_name=resolved_sample_name,
        assay_name=ASSAY_NAME,
        assembly=assembly_declared_text,
        pipeline_version=str(pipeline_version),
        run_time=run_time_text,
    )

    context = {
        "kestrel_highlight": kestrel_html,
        # The two facts the empty states branch on, both derived from the computed
        # state rather than from the shape of the data - which is what conflated them.
        # `vntyper report` renders a supplied summary (#207) that need not carry the
        # step at all, and a run whose Kestrel stage failed carries one that produced
        # nothing; rendering either as "found nothing" asserts something the run never
        # established. (adVNTR still asks only whether its step is recorded. Its
        # section already words the absent case as a disjunction rather than as a
        # negative, so the same conflation there is not the same defect - but it is the
        # same shape, and it is worth closing separately.)
        "kestrel_step_recorded": kestrel_state != STEP_ABSENT,
        "kestrel_result_unreadable": kestrel_state == STEP_UNREADABLE,
        # The same three facts on the adVNTR side, which had only two branches for them.
        # Its absent-case wording used to be a disjunction ("was not performed or no
        # adVNTR results are available") that covered the failure by being vague about
        # it; each state now says which one it is.
        "advntr_step_recorded": advntr_state != STEP_ABSENT,
        "advntr_result_unreadable": advntr_state == STEP_UNREADABLE,
        # The visible/total statement beside each table, counted here from the frame.
        # DataTables' own "Showing 1 to 3 of 3 entries" footer is switched off: it is
        # a second count that contradicts this one whenever a CDN fails (#242). It is
        # empty when there is no table: the sentence exists to say that nothing was
        # withheld from the rows below it, so counting a non-result defeats its purpose.
        "kestrel_row_summary": kestrel_row_summary,
        "kestrel_subthreshold_note": kestrel_subthreshold_note,
        # The folded columns, printed under each table because a nineteen-column table
        # does not fit A4 and silently lost ten of them - the then-120 bp motif pair
        # among them - off the right edge of every sheet.
        "kestrel_folded_record": folded_record_html(kestrel_df, KESTREL_ESSENTIAL_COLUMNS, noun=KESTREL_ROW_NOUN),
        "advntr_folded_record": advntr_folded_record,
        "advntr_highlight": advntr_html,
        "advntr_row_summary": advntr_row_summary,
        "advntr_available": advntr_available,
        "log_content": pipeline_log_content,
        # Deprecated custom-template API, retained through VNtyper 2.x. The shipped
        # template now builds its alignment panel from the structured JSON/session
        # fields below, but configured legacy templates may still splice this fragment.
        # See report_context_contract.DEPRECATED_KEYS.
        "igv_content": igv_content,
        # Interpolated straight into a <script> block, so they must parse even
        # when there is no IGV report at all.
        "table_json": js_json_literal(table_json, EMPTY_TABLE_JSON),
        "session_dictionary": session_literal,
        # The alignment browser. `igv_payload` is the base64 of the vendored, gzipped
        # igv.js and is None unless this run has a session *and* the mode is
        # `embedded`; the template writes the `IGV_GZ_B64` constant only when it is
        # there, so a sample with no alignments carries no library. It is interpolated
        # into a double-quoted JavaScript string under autoescaping rather than marked
        # `|safe`: base64's alphabet contains none of the five characters Jinja2
        # escapes, so escaping is a no-op over it, and
        # `test_the_embedded_payload_survives_the_round_trip_into_the_document`
        # decodes it back out of the rendered file and checks it against the pinned
        # digest rather than assuming that.
        "igv_payload": igv_payload,
        "igv_mode": report_igv,
        "igv_version": report_assets.IGV_VERSION,
        "igv_session_available": igv_session_available,
        # One line for the Provenance block: which library, which digest, and where it
        # is. Built in the pure module because choosing the wording is presentation
        # logic over a computed state (AGENTS.md trap 11).
        "igv_provenance": report_assets.igv_provenance(report_igv),
        # Two timestamps, labelled, because they are different facts. `report_date`
        # is `datetime.now()` at render, so re-running `vntyper report` over an
        # archived run restamped the only date on the page and the artefact
        # silently claimed to be newer than the analysis (#242). Both carry a zone:
        # the run time is UTC and this one is the rendering machine's local time,
        # so printing them unqualified beside each other invites the reader to
        # subtract them.
        "report_date": datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
        "run_date": run_time_text,
        # The mapping is rendered by iterating it, not by branching on its shape:
        # the template tested for `fastq1 and fastq2` or `bam`, so the other two
        # shapes `resolve_pipeline_input` emits rendered an empty line.
        "input_files_text": format_input_files(input_files),
        "sample_name": resolved_sample_name,
        "assay_name": ASSAY_NAME,
        "report_title_prefix": REPORT_TITLE_PREFIX,
        # The two halves of the same title. `<title>` and the printed running header
        # take the whole string above, because neither can carry markup; the `<h1>`
        # takes these, so the gene symbol can be set in italics the way a gene symbol
        # is. Both are interpolated under autoescaping - nothing here is marked safe.
        "report_gene_symbol": GENE_SYMBOL,
        "report_title_description": REPORT_TITLE_DESCRIPTION,
        # Printed on every archived page. The printed record is the artefact that
        # outlives the HTML, and it is the one that gets filed and forwarded.
        "research_use_statement": RESEARCH_USE_STATEMENT,
        # The only value in either template that is interpolated into a stylesheet,
        # and the only one that may be: it is a fragment VNtyper built, with every
        # sample-derived value already escaped for CSS *and* for the raw text element
        # a `<style>` is.
        "print_running_header_css": running_header_css,
        "pipeline_version": pipeline_version,
        # The provenance block. Two of these four rows are *not* new summary keys:
        # `assembly_declared` is the `reference_assembly_requested` `start_summary`
        # has always written, and `assembly_detected` is the BAM-header step's
        # `assembly_text`. Recording either again under a second name would be the
        # divergent-source problem `cli_report.py`'s docstring warns about. Every
        # row says `not recorded by this run` when its source is absent - never
        # `config["default_values"]["reference_assembly"]`, which would mislabel any
        # `--reference-assembly` override and cannot reconstruct `--custom-regions`.
        "summary_schema_version": recorded_or_not(pipeline_summary.get("schema_version")),
        # The words a provenance value renders as when the run recorded nothing for it.
        # Passed in rather than restated in the template, so the masthead can mark those
        # values as absent instead of drawing them as if they were facts - and so the
        # phrase itself exists in exactly one place (`report_identity.NOT_RECORDED`).
        "not_recorded": NOT_RECORDED,
        "assembly_declared": assembly_declared_text,
        "assembly_detected": recorded_or_not(assembly_text),
        "region_resolved": format_region(pipeline_summary.get("region_resolved")),
        # `reference_assembly_requested` and `assembly_text` are deliberately absent:
        # they reach the template as `assembly_declared` and `assembly_detected`
        # above. Passing the raw values too would put the same fact in the context
        # twice under two names, which is how the template came to have two places
        # that could disagree about the assembly.
        "reference_key_used": reference_key_used,
        "reference_path": reference_path,
        "reference_source_effective": reference_source_effective,
        "header_warning": header_warning,
        "alignment_pipeline": alignment_pipeline,
        "assembly_contig": assembly_contig,
        "mean_vntr_coverage": shown(coverage["mean"]),
        "median_vntr_coverage": shown(coverage["median"]),
        "stdev_vntr_coverage": shown(coverage["stdev"]),
        "min_vntr_coverage": shown(coverage["min"]),
        "max_vntr_coverage": shown(coverage["max"]),
        "vntr_region_length": shown(coverage["region_length"]),
        "vntr_uncovered_bases": shown(coverage["uncovered_bases"]),
        "percent_vntr_uncovered": shown(coverage["percent_uncovered"]),
        # #222. Two figures with different standing, so they are labelled differently in the
        # template rather than sitting beside the window mean as if interchangeable:
        # the flank depth is comparable between assemblies without qualification, while the
        # array sum is comparable only under the counting policy printed with it.
        "vntr_flank_mean_depth": shown(_rounded(coverage["vntr_flank_mean_depth"])),
        "vntr_flank_bases": shown(coverage["vntr_flank_bases"]),
        "vntr_array_length": shown(coverage["vntr_array_length"]),
        "vntr_array_depth_sum": shown(coverage["vntr_array_depth_sum"]),
        "vntr_array_depth_sum_per_unit_length": shown(_rounded(coverage["vntr_array_depth_sum_per_unit_length"])),
        "depth_sum_reference_length": shown(coverage["depth_sum_reference_length"]),
        "depth_counting_policy": shown(coverage["depth_counting_policy"]),
        "coverage_qc": coverage_qc.status,
        "coverage_qc_text": coverage_qc_word(coverage_qc.status),
        # Let templates test measuredness without restating this presentation phrase.
        "not_calculated": NOT_CALCULATED,
        # Whether there was anything to judge. `quality_metrics_pass` stays True for a run
        # with no coverage step at all - the screening axis is unchanged by #172's
        # honesty fix - so the chip would otherwise be painted "passing" beside a status
        # that says the gate was never evaluated.
        "coverage_qc_measured": coverage_qc.status != COVERAGE_QC_NOT_EVALUATED,
        "percent_vntr_uncovered_icon": uncovered_icon,
        # The six raw colour names are deprecated custom-template API retained through
        # VNtyper 2.x. Shipped templates use semantic tokens and the accessible icon
        # fragments instead; see report_context_contract.DEPRECATED_KEYS.
        "percent_vntr_uncovered_color": uncovered_color,
        "mean_vntr_coverage_icon": coverage_icon,
        "mean_vntr_coverage_color": coverage_color,
        "fastp_available": fastp.available,
        # Preserve raw values for custom templates; shipped HTML renders these
        # Decimal-derived presentation fields so the reader sees icon operands.
        "duplication_rate": fastp.duplication_rate,
        "duplication_rate_display": duplication_rate_display,
        "duplication_rate_cutoff": fastp_cutoffs.duplication_rate.label,
        "duplication_rate_icon": dup_icon,
        "duplication_rate_color": dup_color,
        "q20_rate": fastp.q20_rate,
        "q20_rate_display": q20_rate_display,
        "q20_rate_cutoff": fastp_cutoffs.q20_rate.label,
        "q20_icon": q20_icon,
        "q20_color": q20_color,
        "q30_rate": fastp.q30_rate,
        "q30_rate_display": q30_rate_display,
        "q30_rate_cutoff": fastp_cutoffs.q30_rate.label,
        "q30_icon": q30_icon,
        "q30_color": q30_color,
        "passed_filter_rate": fastp.passed_filter_rate,
        "passed_filter_rate_display": passed_filter_rate_display,
        "passed_filter_rate_cutoff": fastp_cutoffs.passed_filter_rate.label,
        "passed_filter_icon": pf_icon,
        "passed_filter_color": pf_color,
        "sequencing_str": fastp.sequencing,
        # The configured message as the ordered parts it was authored in, each rendered as
        # an element of its own and autoescaped with everything else. It used to be one
        # `{{ summary_text|safe }}`, marked safe for no reason but the `<br>` separators
        # inside it - so the whole sentence was exempt from escaping to get a line break.
        # `report_config.json` now carries the split beside the verbatim message, and
        # `screening_summary.render_segments` is what pins that the two still agree.
        "screening_segments": screening.rendered_segments,
        # Additive and config-driven. Empty for a measured gate and for a config written
        # before the key existed; configured screening messages are never reworded.
        "coverage_not_measured_note": coverage_not_measured_note(report_config, coverage_qc),
        # The state as words, computed in the pure module: a chip is the most compressed
        # thing the report says about a stage, so it is the one most easily misread as a
        # verdict, and none of these words may be one the run did not reach.
        "state_chips": state_chips(
            screening,
            report_config,
            cross_match_available=cross_match_is_assessable and bool(cross_match_message),
            cross_match_is_positive=cross_match_is_positive,
        ),
        # The full computed screening state, not just `is_positive` - so a report
        # whose state matched no rule (`emphasis == "indeterminate"`) is distinguishable
        # from a genuine negative rather than rendering identically (#242 I2).
        "screening_state": {
            # Deprecated custom-template API retained through VNtyper 2.x. The shipped
            # provenance line consumes the honest execution-aware `*_state_text`
            # values below instead; see report_context_contract.DEPRECATED_KEYS.
            "kestrel_result": screening.kestrel_result,
            "advntr_result": screening.advntr_result,
            "quality_metrics_pass": screening.quality_metrics_pass,
            "matched_rule": screening.matched_rule,
            "emphasis": screening.emphasis,
            # The binding report-context contract exposes execution separately from
            # algorithm result: performed-and-empty is a negative result, while a
            # stage that did not run or failed established no such result.
            "execution": {
                "kestrel": screening.kestrel_execution,
                "advntr": screening.advntr_execution,
            },
            # What the provenance line prints for each algorithm. Built here rather than
            # branched on in the template because the choice is presentation logic over
            # computed state, and that lives in the pure modules (AGENTS.md trap 11) -
            # and because printing `kestrel_result` unconditionally is the defect: a
            # stage that produced nothing computes the block's "negative" default.
            "kestrel_state_text": algorithm_state_text(screening.kestrel_execution, screening.kestrel_result),
            "advntr_state_text": algorithm_state_text(screening.advntr_execution, screening.advntr_result),
        },
        "cross_match_message": cross_match_message,
        "cross_match_is_positive": cross_match_is_positive,
        # What the run named, at the top of the report rather than in column eleven of a
        # nineteen-column table. Built from the *unformatted* frames: the displayed
        # Kestrel frame carries markup in its `Confidence` and `Flag` cells, and the
        # nomenclature columns are read as text. It is None when no row carries a name,
        # which is what the masthead branches on - a run that named nothing must not
        # render an empty allele panel.
        "variant_identity": variant_identity(kestrel_df_raw, advntr_df),
        # Every coded nomenclature value these rows actually use, in words. A tier is a
        # bare letter and a flag is a kebab-case token; both were printed with nothing
        # anywhere in the artefact saying what they meant. Only the terms present are
        # listed, and they are printed rather than put behind a hover, because this
        # report is read on paper as often as on a screen.
        "nomenclature_legend": nomenclature_legend(kestrel_df_raw, advntr_df),
    }

    try:
        rendered_html = template.render(context)
        logger.debug("Report template rendered successfully.")
    except Exception as e:
        logger.error("Failed to render the report template: %s", e)
        raise

    # `--report-file` is documented as a name; this is what makes that true.
    report_file_path = contained_output_path(output_dir, report_file, "--report-file")
    try:
        with open(report_file_path, "w") as f:
            f.write(rendered_html)
        logger.info("Summary report generated and saved to %s", report_file_path)
    except Exception as e:
        logger.error("Failed to write the summary report: %s", e)
        raise
