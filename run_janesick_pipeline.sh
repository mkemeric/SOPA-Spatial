#!/bin/bash
# Janesick Breast Cancer Pipeline — runs sopa CLI steps in sequence.
# Usage: ./run_janesick_pipeline.sh [DATA_PATH] [OUTPUT_PATH]
#
# Defaults:
#   DATA_PATH   = /mnt/shared/janesick/input
#   OUTPUT_PATH = results/janesick.zarr

set -e

DATA_PATH="${1:-/mnt/shared/janesick/input}"
SDATA_PATH="${2:-results/janesick.zarr}"
EXPLORER_PATH="${SDATA_PATH%.zarr}.explorer"

echo "========================================"
echo "Janesick Breast Cancer Pipeline"
echo "========================================"
echo "  Input:    $DATA_PATH"
echo "  Output:   $SDATA_PATH"
echo "  Explorer: $EXPLORER_PATH"
echo ""

# Activate venv if present and not already active
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ -d "$SCRIPT_DIR/.venv" ] && [ -z "$VIRTUAL_ENV" ]; then
    source "$SCRIPT_DIR/.venv/bin/activate"
    echo "Activated venv: $VIRTUAL_ENV"
fi

mkdir -p "$(dirname "$SDATA_PATH")"

# ── Step 1: Convert Xenium data to SpatialData ─────────────────
echo ""
echo "── Step 1: Converting Xenium data to SpatialData ──"
if [ -d "$SDATA_PATH" ]; then
    echo "   SpatialData already exists at $SDATA_PATH, skipping conversion."
else
    sopa convert "$DATA_PATH" \
        --sdata-path "$SDATA_PATH" \
        --technology xenium
    echo "   ✓ Conversion complete"
fi

# ── Step 2: Patchify for parallel segmentation ──────────────────
echo ""
echo "── Step 2: Patchifying image ──"
sopa patchify image "$SDATA_PATH" \
    --patch-width-pixel 2048 \
    --patch-overlap-pixel 100
echo "   ✓ Patchification complete"

# ── Step 3: Cellpose segmentation ───────────────────────────────
echo ""
echo "── Step 3: Running Cellpose segmentation ──"
sopa segmentation cellpose "$SDATA_PATH" \
    --diameter 35 \
    --channels DAPI \
    --model-type cyto2
echo "   ✓ Segmentation complete"

# ── Step 4: Resolve boundary conflicts ──────────────────────────
echo ""
echo "── Step 4: Resolving boundary conflicts ──"
sopa resolve cellpose "$SDATA_PATH"
echo "   ✓ Conflict resolution complete"

# ── Step 5: Aggregate transcripts per cell ──────────────────────
echo ""
echo "── Step 5: Aggregating transcripts ──"
sopa aggregate "$SDATA_PATH" \
    --min-transcripts 10
echo "   ✓ Aggregation complete"

# ── Step 6: Xenium Explorer export ──────────────────────────────
echo ""
echo "── Step 6: Exporting to Xenium Explorer ──"
sopa explorer write "$SDATA_PATH" \
    --output-path "$EXPLORER_PATH"
echo "   ✓ Explorer export complete"

# ── Step 7: Generate report ─────────────────────────────────────
echo ""
echo "── Step 7: Generating report ──"
sopa report "$SDATA_PATH" "$EXPLORER_PATH/analysis_summary.html"
echo "   ✓ Report generated"

echo ""
echo "========================================"
echo "Pipeline Complete!"
echo "========================================"
echo ""
echo "  SpatialData: $SDATA_PATH"
echo "  Explorer:    $EXPLORER_PATH"
echo "  Report:      $EXPLORER_PATH/analysis_summary.html"
echo ""
echo "Next: run SPATCH custom modules (diffusion_analysis, cell_shape_metrics):"
echo "  python3 -c \""
echo "    import spatialdata as sd"
echo "    from spatch_modules.runner import run_custom_pipeline"
echo "    sdata = sd.read_zarr('$SDATA_PATH')"
echo "    results = run_custom_pipeline(sdata, 'configs/janesick_breast_cancer.yaml')"
echo "  \""
echo ""
