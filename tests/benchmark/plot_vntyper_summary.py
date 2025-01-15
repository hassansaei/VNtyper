#!/usr/bin/env python3

"""
Plot VNtyper summary CSV as scatter plots (no lines), including analysis time.

Usage:
  python plot_vntyper_summary.py --input-csv vntyper_summary.csv --output-png plot.png
"""

import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot scatter plots from a vntyper_summary.csv (no lines), including analysis time."
    )
    parser.add_argument(
        "--input-csv",
        required=True,
        help="Path to the vntyper_summary.csv file.",
    )
    parser.add_argument(
        "--output-png",
        default="vntyper_summary_plot.png",
        help="Path to output PNG file.",
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


def main():
    args = parse_args()

    # Configure logging
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

    # 2) If confidence=="Negative" and the metric is NaN, replace with 0.0
    metrics = [
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "Depth_Score",
        "analysis_time_minutes",  # New metric added
    ]
    negative_mask = df["confidence"].fillna("").str.contains("Negative", case=False)
    for m in metrics:
        if m == "analysis_time_minutes":
            # Assuming analysis_time_minutes should not be NaN based on script, but handle just in case
            df.loc[negative_mask & df[m].isna(), m] = 0.0
        else:
            df.loc[negative_mask & df[m].isna(), m] = 0.0

    # 3) Convert "value" to numeric (Y-axis)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    # 4) Build a long-form DF with columns: [file_analyzed, method, confidence, metric_name, x_value, y_value, color]
    long_rows = []
    for _, row in df.iterrows():
        color_str = confidence_to_color(row.get("confidence", ""))
        y_val = row.get("value", float("nan"))

        for metric_name in metrics:
            if metric_name == "analysis_time_minutes":
                # Plot analysis_time_minutes separately
                analysis_time = row.get(metric_name, float("nan"))
                if pd.isna(analysis_time):
                    continue
                long_rows.append({
                    "file_analyzed": row.get("file_analyzed", ""),
                    "method": row.get("method", ""),
                    "confidence": row.get("confidence", ""),
                    "color": color_str,
                    "metric_name": metric_name,
                    "x_value": 1,  # Dummy x_value since analysis time is a single value per run
                    "y_value": analysis_time,
                })
                continue

            raw_str = str(row.get(metric_name, "")).strip()
            if not raw_str or raw_str.lower() == "nan":
                continue
            # e.g., "6,1593" => split by comma
            for val_str in raw_str.split(","):
                val_str = val_str.strip()
                try:
                    x_val = float(val_str)
                except ValueError:
                    continue
                long_rows.append({
                    "file_analyzed": row.get("file_analyzed", ""),
                    "method": row.get("method", ""),
                    "confidence": row.get("confidence", ""),
                    "color": color_str,
                    "metric_name": metric_name,
                    "x_value": x_val,
                    "y_value": y_val,
                })

    long_df = pd.DataFrame(long_rows)
    if long_df.empty:
        logging.warning("No numeric data found to plot. Exiting.")
        return

    # 5) Create 2x2 subplots, one per metric
    num_metrics = len(metrics)
    cols = 2
    rows = (num_metrics + 1) // cols
    fig, axes = plt.subplots(
        nrows=rows, ncols=cols, figsize=(18, 10), sharey=False
    )
    axes = axes.flatten()  # Flatten in case of multiple rows

    for idx, metric_name in enumerate(metrics):
        ax = axes[idx]
        sub = long_df[long_df["metric_name"] == metric_name]
        if sub.empty:
            ax.set_title(f"{metric_name} (No data)")
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
            continue

        if metric_name == "analysis_time_minutes":
            # Since x_value is dummy (1), plot all points along x=1
            ax.scatter(
                sub["x_value"],
                sub["y_value"],
                c=sub["color"],
                s=100,
                alpha=0.7,
                edgecolors='w'
            )
            ax.set_xlabel("Run")
            ax.set_ylabel("Analysis Time (minutes)")
            ax.set_title("Analysis Time")
            ax.set_xlim(0.5, 1.5)
        else:
            # Regular scatter plot for other metrics
            ax.scatter(
                sub["x_value"],
                sub["y_value"],
                c=sub["color"],
                s=100,
                alpha=0.7,
                edgecolors='w'
            )
            ax.set_xlabel(metric_name)
            ax.set_ylabel("Fraction or Coverage (from 'value')")
            ax.set_title(metric_name)
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)

    # If there are unused subplots, hide them
    for idx in range(len(metrics), len(axes)):
        fig.delaxes(axes[idx])

    # 6) Create a custom legend below the plots
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='High_Precision',
               markerfacecolor='red', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Low_Precision',
               markerfacecolor='orange', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Negative',
               markerfacecolor='black', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Other',
               markerfacecolor='gray', markersize=10),
    ]
    fig.legend(
        handles=legend_elements,
        loc='lower center',
        ncol=4,
        title='Confidence',
        bbox_to_anchor=(0.5, 0.05)
    )

    plt.suptitle("VNtyper Scatter Plots (X = depth metrics, Y = fraction/coverage/time, NO lines)", fontsize=20)
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Adjust to make space for the legend
    plt.savefig(args.output_png, dpi=300)
    logging.info(f"Plot saved to {args.output_png}")


if __name__ == "__main__":
    main()
