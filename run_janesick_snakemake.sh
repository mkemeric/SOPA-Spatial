#!/bin/bash
# Janesick Breast Cancer Pipeline — Snakemake parallel execution.
# Requires conda/mamba. Uses 16 cores by default (configurable via profile).
#
# Usage:
#   ./run_janesick_snakemake.sh [DATA_PATH] [OUTPUT_PATH]
#
# Defaults:
#   DATA_PATH   = ../data/outs/
#   OUTPUT_PATH = results/janesick.zarr

set -e

DATA_PATH="${1:-../data/outs/}"
SDATA_PATH="${2:-results/janesick.zarr}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKFLOW_DIR="$SCRIPT_DIR/sopa/workflow"
CONFIG_FILE="$SCRIPT_DIR/configs/janesick_sopa.yaml"
PROFILE_DIR="$SCRIPT_DIR/profiles/parallel"

echo "========================================"
echo "Janesick Pipeline (Snakemake)"
echo "========================================"
echo "  Input:    $DATA_PATH"
echo "  Output:   $SDATA_PATH"
echo "  Config:   $CONFIG_FILE"
echo "  Profile:  $PROFILE_DIR"
echo ""

# Check conda/mamba is available
if command -v mamba >/dev/null 2>&1; then
    echo "Using mamba for conda environments"
elif command -v conda >/dev/null 2>&1; then
    echo "Using conda for environments"
else
    echo "❌ ERROR: Neither conda nor mamba found."
    echo "   Install Miniforge: https://github.com/conda-forge/miniforge#install"
    echo ""
    echo "   Quick install:"
    echo "     curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-\$(uname)-\$(uname -m).sh"
    echo "     bash Miniforge3-\$(uname)-\$(uname -m).sh -b -p \$HOME/miniforge3"
    echo "     eval \"\$(\$HOME/miniforge3/bin/conda shell.bash hook)\""
    echo ""
    echo "   Or use the CLI pipeline instead: ./run_janesick_pipeline.sh"
    exit 1
fi

# Activate venv/conda if needed for snakemake itself
if [ -d "$SCRIPT_DIR/.venv" ] && [ -z "$VIRTUAL_ENV" ]; then
    source "$SCRIPT_DIR/.venv/bin/activate"
fi

# Verify snakemake is available
if ! command -v snakemake >/dev/null 2>&1; then
    echo "❌ ERROR: snakemake not found. Install with:"
    echo "   pip install 'snakemake>=8.0'"
    exit 1
fi

echo ""
echo "Running Snakemake with $(grep 'cores:' "$PROFILE_DIR/config.yaml" | awk '{print $2}') cores..."
echo ""

snakemake \
    --snakefile "$WORKFLOW_DIR/Snakefile" \
    --configfile "$CONFIG_FILE" \
    --config data_path="$DATA_PATH" sdata_path="$SDATA_PATH" \
    --profile "$PROFILE_DIR" \
    "$@"

# ── Run SPATCH custom modules + visualizations ───────────────
echo ""
echo "── Running SPATCH custom modules ──"
spatch run "$SDATA_PATH" \
    --config "$SCRIPT_DIR/configs/janesick_breast_cancer.yaml" \
    --save
echo "   ✓ SPATCH modules complete"

echo ""
echo "========================================"
echo "Pipeline Complete!"
echo "========================================"
echo ""
echo "  SpatialData: $SDATA_PATH"
echo "  Figures:     results/janesick_breast_cancer/figures/"
echo ""
