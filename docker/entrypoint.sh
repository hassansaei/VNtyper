#!/bin/bash
set -e

# Initialize conda
source /opt/conda/etc/profile.d/conda.sh

# Determine which command to run based on the first argument
if [ "$1" = "vntyper" ]; then
    # Activate the vntyper environment and execute the command
    conda activate vntyper
    exec vntyper "${@:2}"
elif [ "$1" = "advntr" ]; then
    # Activate the envadvntr environment and execute the command
    conda activate envadvntr
    exec advntr "${@:2}"
else
    # If the command is not recognized, execute it directly
    exec "$@"
fi
