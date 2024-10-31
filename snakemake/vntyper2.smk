import os

# ----------------------------------------------------------------------------------- #
# Load BAM file paths
bam_files = [line.strip() for line in open('bams.txt')]

# Create a dictionary to store the jobs
jobs = {}
for bam in bam_files:
    bam_basename = os.path.basename(bam).replace(".bam", "")
    bam_file = bam
    output_dir = f"results/{bam_basename}/"
    log_file = f"logs/{bam_basename}.log"
    jobs[bam_basename] = {"bam": bam, "bam_basename": bam_basename, "output_dir": output_dir, "log": log_file}

# ----------------------------------------------------------------------------------- #
# Helper functions
def get_mem_from_threads(wildcards, threads):
    """Calculate memory allocation based on the number of threads."""
    return threads * 2200  # 2.2 GB per thread

# ----------------------------------------------------------------------------------- #
# Define the rules
rule all:
    input:
        expand("logs/{bam_basename}.log", bam_basename=[jobs[key]['bam_basename'] for key in jobs.keys()])

rule run_vntyper_pipeline:
    input:
        bam=lambda wildcards: jobs[wildcards.bam_basename]['bam'],
    output:
        log="logs/{bam_basename}.log",
    params:
        output_dir="results/{bam_basename}/",
        config_path="vntyper/config.json"
    threads: 8
    resources:
        mem_mb=get_mem_from_threads,
        time="72:00:00",
        tmpdir=os.environ.get('TMPDIR', '/tmp')
    conda:
        "vntyper"
    shell:
        """
        vntyper pipeline --bam {input.bam} --thread {threads} --reference-assembly hg38 \
                         --fast-mode --keep-intermediates -o {params.output_dir} \
                         --ignore-advntr --config-path {params.config_path} &> {output.log}
        """
