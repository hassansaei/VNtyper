#!/usr/bin/env python3
"""
Plot VNtyper summary CSV as scatter plots (no lines), including analysis time.
Generates separate plots for vntyper kestrel and advntr modules.

Usage:
  python plot_vntyper_summary.py --input-csv vntyper_summary.csv --output-png plot.png
"""

import argparse
import ast
import logging
import os
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot scatter plots from a vntyper_summary.csv (no lines), including analysis time. "
        "Separate figures are produced for vntyper kestrel and advntr results."
    )
    parser.add_argument(
        "--input-csv",
        required=True,
        help="Path to the vntyper_summary.csv file.",
    )
    parser.add_argument(
        "--output-png",
        default="vntyper_summary_plot.png",
        help="Base path to output PNG file; module name will be appended (e.g. _vntyper.png, _advntr.png).",
    )
    return parser.parse_args()


def confidence_to_color(conf_str: str) -> str:
    """
    Map confidence string to color:
      - High_Precision  => red
      - Low_Precision   => orange
      - Negative        => black
      - else            => gray
    """
    if not isinstance(conf_str, str):
        return "gray"
    conf_str_lower = conf_str.lower()
    if "high_precision" in conf_str_lower:
        return "red"
    elif "low_precision" in conf_str_lower:
        return "orange"
    elif "negative" in conf_str_lower:
        return "black"
    else:
        return "gray"


def parse_advntr_result(adv_str):
    """
    Parse the advntr_result string.

    If the string (after stripping) equals "Negative" (case-insensitive), return zeros.
    Otherwise, attempt to parse the string as a list of dictionaries and extract:
      - NumberOfSupportingReads
      - MeanCoverage
      - Pvalue
    If parsing fails, return zeros.
    """
    if not isinstance(adv_str, str):
        return 0, 0, 0
    s = adv_str.strip()
    if s.lower() == "negative":
        return 0, 0, 0
    try:
        parsed = ast.literal_eval(s)
        if isinstance(parsed, list) and len(parsed) > 0 and isinstance(parsed[0], dict):
            nsr = parsed[0].get("NumberOfSupportingReads", 0)
            mean_cov = parsed[0].get("MeanCoverage", 0)
            pval = parsed[0].get("Pvalue", 0)
            return nsr, mean_cov, pval
    except Exception:
        return 0, 0, 0
    return 0, 0, 0


def get_advntr_call(row):
    """
    Determine advntr call based on advntr_result.
    Returns "Positive" if a valid result is found (i.e. Pvalue is not None or zero),
    otherwise "Negative".
    """
    adv_str = row.get("advntr_result", "")
    if not isinstance(adv_str, str):
        return "Negative"
    s = adv_str.strip()
    if s.lower() == "negative":
        return "Negative"
    try:
        parsed = ast.literal_eval(s)
        if isinstance(parsed, list) and len(parsed) > 0 and isinstance(parsed[0], dict):
            pval = parsed[0].get("Pvalue", None)
            if pval is None or pval == 0:
                return "Negative"
            return "Positive"
    except Exception:
        return "Negative"
    return "Negative"


def advntr_color_func(row):
    """
    For advntr module, use green for Positive and black for Negative.
    """
    call = row.get("advntr_call", "Negative")
    return "green" if call.lower() == "positive" else "black"


def build_long_form(df, metrics, module_label, value_column="value", color_func=None):
    """
    Build a long-form DataFrame from a given metrics list.

    Each row in the output contains:
      - file_analyzed, method, confidence, module, metric_name, x_value, y_value, color

    If color_func is provided, it is used to compute the color.
    """
    long_rows = []
    for _, row in df.iterrows():
        if color_func is not None:
            color_str = color_func(row)
        else:
            color_str = confidence_to_color(row.get("confidence", ""))
        try:
            y_val = float(row.get(value_column, float("nan")))
        except (ValueError, TypeError):
            y_val = float("nan")
        for metric in metrics:
            raw_val = row.get(metric, "")
            if pd.isna(raw_val):
                continue
            if isinstance(raw_val, str) and "," in raw_val:
                for val_str in raw_val.split(","):
                    val_str = val_str.strip()
                    try:
                        x_val = float(val_str)
                    except ValueError:
                        continue
                    long_rows.append(
                        {
                            "file_analyzed": row.get("file_analyzed", ""),
                            "method": row.get("method", ""),
                            "confidence": row.get("confidence", ""),
                            "module": module_label,
                            "metric_name": metric,
                            "x_value": x_val,
                            "y_value": y_val,
                            "color": color_str,
                        }
                    )
            else:
                try:
                    x_val = float(raw_val)
                except (ValueError, TypeError):
                    continue
                long_rows.append(
                    {
                        "file_analyzed": row.get("file_analyzed", ""),
                        "method": row.get("method", ""),
                        "confidence": row.get("confidence", ""),
                        "module": module_label,
                        "metric_name": metric,
                        "x_value": x_val,
                        "y_value": y_val,
                        "color": color_str,
                    }
                )
    return pd.DataFrame(long_rows)


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # 1) Load the CSV
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as ex:
        logging.error(f"Could not read CSV {args.input_csv}: {ex}")
        return
    if df.empty:
        logging.warning(f"{args.input_csv} is empty, no plot to generate.")
        return

    # 2) For rows with negative confidence, fill NaNs for selected metrics with 0.0.
    common_metrics = ["analysis_time_minutes"]
    vntyper_metrics = [
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "Depth_Score",
    ]
    for m in vntyper_metrics + common_metrics:
        df.loc[
            df["confidence"].fillna("").str.contains("Negative", case=False)
            & df[m].isna(),
            m,
        ] = 0.0

    # 3) Ensure "value" (the downsample fraction or coverage) is numeric.
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    # 4) Parse advntr_result into new columns.
    (
        df["advntr_NumberOfSupportingReads"],
        df["advntr_MeanCoverage"],
        df["advntr_Pvalue"],
    ) = zip(
        *df.apply(lambda row: parse_advntr_result(row.get("advntr_result", "")), axis=1)
    )
    # 5) Compute a new column "advntr_call" for advntr results.
    df["advntr_call"] = df.apply(get_advntr_call, axis=1)

    # 6) Build long-form DataFrames for each module.
    # vntyper module uses its original metric columns.
    vntyper_long = build_long_form(
        df,
        vntyper_metrics + common_metrics,
        module_label="vntyper kestrel",
        value_column="value",
    )
    # advntr module uses the parsed advntr metrics plus analysis_time_minutes.
    advntr_metrics = [
        "advntr_NumberOfSupportingReads",
        "advntr_MeanCoverage",
        "advntr_Pvalue",
    ]
    advntr_long = build_long_form(
        df,
        advntr_metrics + common_metrics,
        module_label="advntr",
        value_column="value",
        color_func=advntr_color_func,
    )

    if vntyper_long.empty and advntr_long.empty:
        logging.warning("No numeric data found to plot. Exiting.")
        return

    def plot_module(long_df, module_label):
        metrics_to_plot = list(long_df["metric_name"].unique())
        num_metrics = len(metrics_to_plot)
        cols = 2
        rows_subplot = (num_metrics + cols - 1) // cols
        fig, axes = plt.subplots(
            nrows=rows_subplot, ncols=cols, figsize=(18, 10), sharey=False
        )
        axes = axes.flatten()

        for idx, metric in enumerate(metrics_to_plot):
            ax = axes[idx]
            sub = long_df[long_df["metric_name"] == metric]
            if sub.empty:
                ax.set_title(f"{metric} (No data)")
                ax.set_xlim(left=0)
                ax.set_ylim(bottom=0)
                continue

            ax.scatter(
                sub["x_value"],
                sub["y_value"],
                c=sub["color"],
                s=100,
                alpha=0.7,
                edgecolors="w",
            )
            if metric == "analysis_time_minutes":
                ax.set_xlabel("Analysis Time (minutes)")
                ax.set_ylabel("Fraction or Coverage (from 'value')")
                ax.set_title("Coverage/Fraction vs Analysis Time")
            else:
                ax.set_xlabel(metric)
                ax.set_ylabel("Fraction or Coverage (from 'value')")
                ax.set_title(metric)
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)

        for j in range(idx + 1, len(axes)):
            fig.delaxes(axes[j])

        # Adapt legend based on module.
        from matplotlib.lines import Line2D

        if module_label == "advntr":
            legend_elements = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="Positive",
                    markerfacecolor="green",
                    markersize=10,
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="Negative",
                    markerfacecolor="black",
                    markersize=10,
                ),
            ]
        else:
            legend_elements = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="High_Precision",
                    markerfacecolor="red",
                    markersize=10,
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="Low_Precision",
                    markerfacecolor="orange",
                    markersize=10,
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="Negative",
                    markerfacecolor="black",
                    markersize=10,
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label="Other",
                    markerfacecolor="gray",
                    markersize=10,
                ),
            ]
        fig.legend(
            handles=legend_elements,
            loc="lower center",
            ncol=len(legend_elements),
            title="Confidence" if module_label != "advntr" else "Advntr Call",
            bbox_to_anchor=(0.5, 0.05),
        )

        plt.suptitle(
            f"VNtyper Scatter Plots ({module_label})\n(X = depth metric or Analysis Time, Y = fraction/coverage)",
            fontsize=20,
        )
        plt.tight_layout(rect=[0, 0.1, 1, 0.95])
        base, ext = os.path.splitext(args.output_png)
        output_file = f"{base}_{module_label.replace(' ', '_')}{ext}"
        plt.savefig(output_file, dpi=300)
        logging.info(f"Plot saved to {output_file}")
        plt.close(fig)

    if not vntyper_long.empty:
        plot_module(vntyper_long, "vntyper kestrel")
    else:
        logging.info("No vntyper kestrel data available to plot.")

    if not advntr_long.empty:
        plot_module(advntr_long, "advntr")
    else:
        logging.info("No advntr data available to plot.")


if __name__ == "__main__":
    main()
