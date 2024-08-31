import pandas as pd
import os
import logging
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

def load_results_from_dirs(input_dirs, file_name):
    """Load the results from the specified file across multiple directories."""
    dfs = []
    for dir in input_dirs:
        file_path = Path(dir) / file_name
        if file_path.exists():
            df = pd.read_csv(file_path, sep='\t')
            dfs.append(df)
        else:
            logging.warning(f"File not found: {file_path}")
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()

def aggregate_cohort(input_dirs, output_dir, config_path, summary_file):
    # Load the results from multiple directories
    kestrel_df = load_results_from_dirs(input_dirs, "kestrel/kestrel_result.tsv")
    advntr_df = load_results_from_dirs(input_dirs, "advntr/output_adVNTR.tsv")

    # Filter for positive cases
    kestrel_positive = kestrel_df[kestrel_df['Confidence'].str.contains('High_Precision', na=False)]
    advntr_positive = advntr_df[advntr_df['Pvalue'] < 0.05]

    # Generate summary statistics and plots
    generate_cohort_summary_report(output_dir, kestrel_positive, advntr_positive, summary_file)

def generate_cohort_summary_report(output_dir, kestrel_positive, advntr_positive, summary_file):
    # Load the template
    template_dir = "vntyper/templates"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('cohort_summary_template.html')

    # Generate summary statistics
    total_kestrel = len(kestrel_positive)
    total_advntr = len(advntr_positive)

    # Plot summary
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(x=['Kestrel Positive', 'adVNTR Positive'], y=[total_kestrel, total_advntr], ax=ax)
    plot_path = Path(output_dir) / "cohort_summary_plot.png"
    plt.savefig(plot_path)
    plt.close()

    # Render the HTML report
    rendered_html = template.render(
        report_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        kestrel_positive=kestrel_positive.to_html(classes='table table-striped'),
        advntr_positive=advntr_positive.to_html(classes='table table-striped'),
        plot_path=str(plot_path)
    )

    # Save the HTML report
    report_file_path = Path(output_dir) / summary_file
    with open(report_file_path, 'w') as f:
        f.write(rendered_html)

    logging.info(f"Cohort summary report generated and saved to {report_file_path}")
