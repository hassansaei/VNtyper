#!/bin/bash
#
#SBATCH --job-name=sm_vntyper_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=2200M
#SBATCH --output=slurm_logs/%x-%j.log

# Set up temporary directory within your home directory
export TMPDIR=$HOME/scratch/tmp
export TMPDIR=$(mktemp -d)

# Ensure the temporary directory is cleaned up upon script exit
trap "rm -rf $TMPDIR" EXIT

# Create the log directory if it doesn't exist
mkdir -p slurm_logs

# Set default SBATCH output logging
export SBATCH_DEFAULTS=" --output=slurm_logs/%x-%j.log"

# Start date and time logging
date

# Run the Snakemake workflow with the specified Snakefile
srun snakemake -s vntyper2.smk --use-conda -j150

# End date and time logging
date
