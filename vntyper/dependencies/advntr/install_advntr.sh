#!/bin/bash
# install_advntr.sh
# A helper script to automate the installation of adVNTR with optional conda environment activation

set -e  # Exit immediately if any command exits with a non-zero status

# Default configuration file (can be overridden with -c/--config).
#
# Resolved **beside this script**, not in the working directory. The image copies the whole
# directory to /tmp/advntr/ and runs `bash /tmp/advntr/install_advntr.sh` from elsewhere, so
# a bare relative name meant the shipped config was never sourced -- the build silently used
# the script's own fallbacks, and the GIT_COMMIT pin would have had no effect on the image
# it exists to pin (#254). A config that is only read when you happen to be standing in the
# right directory is worse than no config.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/install_advntr.cfg"

# Resolve -c/--config *before* sourcing anything, so an explicit config replaces the shipped
# one instead of layering on top of it.
#
# This became load-bearing when the default above stopped being a bare relative name. That
# name usually matched nothing, so `-c mine.cfg` was in practice the only file sourced. The
# script-adjacent default always exists, so sourcing it first and letting `-c` override
# afterwards would leak every value the caller's file does not happen to set -- GIT_COMMIT
# above all. A `-c` that sets only GIT_BRANCH would then clone the caller's branch and check
# out this repository's pinned commit on top of it, silently building a tree they did not
# choose, or abort reporting that *their* pin is missing when they set no pin at all. The
# shipped config is a complete set of values, not a base layer.
config_override=""
config_flag_seen=false
for ((i = 1; i <= $#; i++)); do
    case "${!i}" in
        -c|--config)
            config_flag_seen=true
            next=$((i + 1))
            config_override="${!next}"
            ;;
    esac
done
if [ "$config_flag_seen" = true ]; then
    if [ -z "$config_override" ] || [ ! -f "$config_override" ]; then
        echo "Configuration file $config_override not found."
        exit 1
    fi
    CONFIG_FILE="$config_override"
fi

# Source exactly one configuration file: the explicit one, or the shipped default.
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
fi

# Set defaults if variables are not defined in the config.
INSTALL_DIR=${INSTALL_DIR:-"$PWD/adVNTR"}
OVERWRITE=false
CONDA_ENV=${CONDA_ENV:-""}
GIT_REPO=${GIT_REPO:-"https://github.com/berntpopp/adVNTR.git"}
GIT_BRANCH=${GIT_BRANCH:-"main"}
GIT_COMMIT=${GIT_COMMIT:-""}

# Function to display help message
function display_help() {
    echo "Usage: bash install_advntr.sh [options]"
    echo ""
    echo "Options:"
    echo "  -e, --env              Name of the conda environment to activate (optional)."
    echo "  -d, --install-dir      Directory where adVNTR will be installed (default: \$INSTALL_DIR)."
    echo "  -o, --overwrite        Overwrite the installation directory if it exists."
    echo "  -c, --config           Path to configuration file. Replaces the config shipped"
    echo "                         beside this script rather than layering over it, so it"
    echo "                         must set every variable it wants to differ from the"
    echo "                         built-in fallbacks below."
    echo "  -h, --help             Display this help message."
    echo ""
    echo "Config file variables:"
    echo "  CONDA_ENV   : Name of the conda environment (e.g., envadvntr)."
    echo "  INSTALL_DIR : Directory for installation."
    echo "  GIT_REPO    : Git repository URL for adVNTR."
    echo "  GIT_BRANCH  : Git branch to clone."
    echo "  GIT_COMMIT  : Exact revision to check out after cloning. This, not GIT_BRANCH,"
    echo "                determines the tree that is built. Empty means the branch tip."
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
            # Already resolved and sourced above, before the defaults were applied. Consume
            # its value here so the rest of the parse stays aligned.
            shift
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
# Printed because it, not the branch, decides the tree that gets built. A pin nobody can see
# in the build log is a pin nobody checks.
echo "  Git commit pin: ${GIT_COMMIT:-<none: branch tip>}"
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

# Check out the pinned revision, if one is configured.
#
# The branch above is only a clone hint; this is what determines the tree. VNtyper records
# upstream line numbers as evidence for why adVNTR runs at `-t 1`, so "whatever the branch
# points at today" is not good enough (#254).
#
# A configured commit that is not present must abort. Silently continuing would leave the
# branch tip checked out while the build reported success -- exactly the mutable-pin
# failure this exists to prevent.
if [ -n "$GIT_COMMIT" ]; then
    echo "Checking out pinned revision $GIT_COMMIT..."
    if ! git checkout --quiet "$GIT_COMMIT"; then
        echo "ERROR: pinned revision $GIT_COMMIT is not present in $GIT_REPO (branch $GIT_BRANCH)." >&2
        echo "       Refusing to build against the branch tip instead: the pin exists because" >&2
        echo "       VNtyper cites exact line numbers from this revision." >&2
        exit 1
    fi
    echo "adVNTR revision: $(git rev-parse HEAD)"
fi

# Install adVNTR
# Note: CFLAGS workaround for GCC 14+ which treats -Wincompatible-pointer-types
# as a hard error. pomegranate's Cython-generated C code passes long int* where
# scipy's BLAS wrapper expects int*. These flags are safely ignored on older GCC.
echo "Installing adVNTR..."
CFLAGS="${CFLAGS:-} -Wno-error=incompatible-pointer-types -Wno-error=int-conversion" python setup.py install

echo "adVNTR installation completed successfully in $INSTALL_DIR."
echo ""
if [ -n "$CONDA_ENV" ]; then
    echo "To use adVNTR, activate the conda environment with:"
    echo "conda activate $CONDA_ENV"
else
    echo "Ensure that your conda environment is activated before using adVNTR."
fi
