#!/usr/bin/env python3

"""
Plot VNtyper summary CSV as scatter plots (no lines).

Usage:
  python plot_vntyper_summary.py --input-csv vntyper_summary.csv --output-png plot.png
"""

import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="Plot scatter plots from a vntyper_summary.csv (no lines).")
    parser.add_argument("--input-csv", required=True, help="Path to the vntyper_summary.csv file.")
    parser.add_argument("--output-png", default="vntyper_summary_plot.png", help="Path to output PNG file.")
    return parser.parse_args()


def confidence_to_color(conf_str: str) -> str:
    """
    Map confidence string -> color:
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
        "Depth_Score"
    ]
    negative_mask = df["confidence"].fillna("").str.contains("Negative", case=False)
    for m in metrics:
        df.loc[negative_mask & df[m].isna(), m] = 0.0

    # 3) Convert "value" to numeric (Y-axis)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    # 4) Build a long-form DF with columns: [file_analyzed, metric_name, x_value, y_value, color]
    long_rows = []
    for _, row in df.iterrows():
        color_str = confidence_to_color(row.get("confidence", ""))
        y_val = row.get("value", float("nan"))

        for metric_name in metrics:
            raw_str = str(row.get(metric_name, "")).strip()
            if not raw_str or raw_str.lower() == "nan":
                continue
            # e.g. "6,1593" => split by comma
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

    # 5) Create 3 subplots, one per metric
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5), sharey=True)
    metric_list = metrics  # same order as above

    for ax, metric_name in zip(axes, metric_list):
        sub = long_df[long_df["metric_name"] == metric_name]
        if sub.empty:
            ax.set_title(f"{metric_name} (No data)")
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
            continue

        # 6) Plot each file_analyzed as individual scatter points
        #    (NO lines connecting points)
        for file_name, group in sub.groupby("file_analyzed"):
            # Each row -> one point
            for _, r in group.iterrows():
                ax.scatter(r["x_value"], r["y_value"], c=r["color"], s=40)

        ax.set_xlabel(metric_name)
        ax.set_ylabel("Fraction or Coverage (from 'value')")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    plt.suptitle("VNtyper Scatter Plots (X = depth metrics, Y = fraction/coverage, NO lines)")
    plt.tight_layout()
    plt.savefig(args.output_png, dpi=150)
    logging.info(f"Plot saved to {args.output_png}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )
    main()
