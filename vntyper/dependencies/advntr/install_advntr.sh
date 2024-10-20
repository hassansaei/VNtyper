#!/bin/bash

# install_advntr.sh
# A helper script to automate the installation of adVNTR with optional conda environment activation

set -e  # Exit immediately if a command exits with a non-zero status

# Default configuration
INSTALL_DIR="$PWD/adVNTR"
REFERENCE_VNTR_URL="https://github.com/mehrdadbakhtiari/adVNTR/releases/download/0.1/vntr_data_recommended_loci.zip"
OVERWRITE=false
CONDA_ENV=""

# Function to display help message
function display_help() {
    echo "Usage: bash install_advntr.sh [options]"
    echo ""
    echo "Options:"
    echo "  -e, --env              Name of the conda environment to activate (optional)."
    echo "  -d, --install-dir      Directory where adVNTR will be installed (default: $PWD/adVNTR)."
    echo "  -o, --overwrite        Overwrite the installation directory if it exists."
    echo "  -h, --help             Display this help message."
    echo ""
    echo "If a conda environment is specified, the script will attempt to activate it."
    echo "If not specified, the script assumes the conda environment is already activated."
    exit 0
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -e|--env)
            CONDA_ENV="$2"
            shift
            ;;
        -d|--install-dir)
            INSTALL_DIR="$2"
            shift
            ;;
        -o|--overwrite)
            OVERWRITE=true
            ;;
        -h|--help)
            display_help
            ;;
        *)
            echo "Unknown parameter passed: $1"
            display_help
            ;;
    esac
    shift
done

echo "Installation settings:"
echo "  Install directory: $INSTALL_DIR"
echo "  Overwrite if exists: $OVERWRITE"
if [ -n "$CONDA_ENV" ]; then
    echo "  Conda environment to activate: $CONDA_ENV"
else
    echo "  No conda environment specified. Assuming environment is already activated."
fi

# Attempt to activate the conda environment if specified
if [ -n "$CONDA_ENV" ]; then
    # Check if conda is installed
    if ! command -v conda &> /dev/null; then
        echo "Error: Conda is not installed. Please install Conda before running this script."
        exit 1
    fi

    # Check if the conda environment exists
    if conda env list | grep -qw "$CONDA_ENV"; then
        echo "Activating conda environment: $CONDA_ENV"
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    else
        echo "Error: Conda environment '$CONDA_ENV' does not exist. Please create it before running this script."
        exit 1
    fi
fi

# Check if adVNTR directory exists
if [ -d "$INSTALL_DIR" ]; then
    if [ "$OVERWRITE" = true ]; then
        echo "Overwriting existing installation directory: $INSTALL_DIR"
        rm -rf "$INSTALL_DIR"
    else
        echo "Installation directory already exists: $INSTALL_DIR"
        echo "Use the -o or --overwrite option to overwrite it."
        exit 1
    fi
fi

# Clone the adVNTR repository from the enhanced_hmm branch
echo "Cloning adVNTR repository into $INSTALL_DIR..."
git clone https://github.com/mehrdadbakhtiari/adVNTR.git --branch enhanced_hmm "$INSTALL_DIR"

cd "$INSTALL_DIR"

# Install adVNTR
echo "Installing adVNTR..."
python setup.py install

echo "adVNTR installation completed successfully in $INSTALL_DIR."

# Provide usage instructions
echo ""
if [ -n "$CONDA_ENV" ]; then
    echo "To use adVNTR, activate the conda environment with:"
    echo "conda activate $CONDA_ENV"
else
    echo "Note: Ensure that your conda environment is activated before using adVNTR."
fi
