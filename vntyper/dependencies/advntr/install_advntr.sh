#!/bin/bash
# install_advntr.sh
# A helper script to automate the installation of adVNTR with optional conda environment activation

set -e  # Exit immediately if any command exits with a non-zero status

# Default configuration file (can be overridden with -c/--config)
CONFIG_FILE="install_advntr.cfg"

# If a configuration file exists in the current directory, source it.
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
fi

# Set defaults if variables are not defined in the config.
INSTALL_DIR=${INSTALL_DIR:-"$PWD/adVNTR"}
OVERWRITE=false
CONDA_ENV=${CONDA_ENV:-""}
GIT_REPO=${GIT_REPO:-"https://github.com/berntpopp/adVNTR.git"}
GIT_BRANCH=${GIT_BRANCH:-"enhanced_hmm"}

# Function to display help message
function display_help() {
    echo "Usage: bash install_advntr.sh [options]"
    echo ""
    echo "Options:"
    echo "  -e, --env              Name of the conda environment to activate (optional)."
    echo "  -d, --install-dir      Directory where adVNTR will be installed (default: \$INSTALL_DIR)."
    echo "  -o, --overwrite        Overwrite the installation directory if it exists."
    echo "  -c, --config           Path to configuration file (default: install_advntr.cfg)."
    echo "  -h, --help             Display this help message."
    echo ""
    echo "Config file variables:"
    echo "  CONDA_ENV   : Name of the conda environment (e.g., envadvntr)."
    echo "  INSTALL_DIR : Directory for installation."
    echo "  GIT_REPO    : Git repository URL for adVNTR."
    echo "  GIT_BRANCH  : Git branch to clone."
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
        -c|--config)
            CONFIG_FILE="$2"
            shift
            if [ -f "$CONFIG_FILE" ]; then
                source "$CONFIG_FILE"
            else
                echo "Configuration file $CONFIG_FILE not found."
                exit 1
            fi
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
echo "  Git repository: $GIT_REPO"
echo "  Git branch: $GIT_BRANCH"
if [ -n "$CONDA_ENV" ]; then
    echo "  Conda environment to activate: $CONDA_ENV"
else
    echo "  No conda environment specified. Assuming environment is already activated."
fi

# Attempt to activate the conda environment if specified
if [ -n "$CONDA_ENV" ]; then
    if ! command -v conda &> /dev/null; then
        echo "Error: Conda is not installed. Please install Conda before running this script."
        exit 1
    fi

    if conda env list | grep -qw "$CONDA_ENV"; then
        echo "Activating conda environment: $CONDA_ENV"
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    else
        echo "Error: Conda environment '$CONDA_ENV' does not exist. Please create it before running this script."
        exit 1
    fi
fi

# Check if the installation directory exists
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

# Clone the adVNTR repository from the specified branch.
echo "Cloning adVNTR repository from $GIT_REPO (branch: $GIT_BRANCH) into $INSTALL_DIR..."
git clone "$GIT_REPO" --branch "$GIT_BRANCH" "$INSTALL_DIR"

cd "$INSTALL_DIR"

# Install adVNTR
echo "Installing adVNTR..."
python setup.py install

echo "adVNTR installation completed successfully in $INSTALL_DIR."
echo ""
if [ -n "$CONDA_ENV" ]; then
    echo "To use adVNTR, activate the conda environment with:"
    echo "conda activate $CONDA_ENV"
else
    echo "Ensure that your conda environment is activated before using adVNTR."
fi
