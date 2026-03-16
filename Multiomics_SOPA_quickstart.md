# Multiomics SOPA + SPATCH Quickstart

This guide walks through two complete spatial transcriptomics pipelines
using **Sopa** (core processing) and **SPATCH** (custom analysis modules):

1. **Janesick Breast Cancer** — Xenium (10x Genomics) → Sopa CLI → SPATCH
2. **CODEX COAD** — Colon adenocarcinoma paired ST + proteomics → SPATCH

Both pipelines produce parquet-based analysis outputs and publication-quality
figures. A final section covers viewing and interacting with results in
Jupyter.

---

## 1. Environment Setup

### 1.1 Clone the Repository

```bash
git clone https://github.com/mkemeric/SOPA-Spatial.git
cd /mnt/shared/SOPA-Spatial
```

### 1.2 Run the Setup Script

```bash
source setup_environment.sh
```

This will:
- Create a `spatch` conda environment (Python 3.12)
- Install a **frozen, validated dependency stack** from
  `requirements-spatial-stack.txt` (spatialdata, sopa, zarr, dask, etc.)
- Install analysis packages (scanpy, cellpose, snakemake, etc.)
- Install `spatch_modules` in editable mode
- Register a **SPATCH** Jupyter kernel
- Clone the sopa repo (shallow) for Snakemake workflow files
- Export `PYTHONWARNINGS` and `SOPA_PARALLELIZATION_BACKEND=dask`

If your home directory is small (<15 GB free), the script auto-redirects
conda/pip storage to a larger volume. Override with:

```bash
export SPATCH_STORAGE=/mnt/user
source setup_environment.sh
```

> **Note:** Use `source` (not `bash`) so the environment stays active in
> your current shell. On subsequent runs, the script detects the existing
> environment and activates instantly.

### 1.3 Activate the Environment

Every time you open a new terminal session:

```bash
source setup_environment.sh
```

This is equivalent to `conda activate spatch` but also sets the required
environment variables (`PYTHONWARNINGS`, `SOPA_PARALLELIZATION_BACKEND`).

### 1.4 Verify Installation

```bash
python3 setup_local_imports.py
```
---

## 2. Pipeline 1: Janesick Breast Cancer (Xenium)

This pipeline processes the Janesick et al. (Nature Communications, 2023)
FFPE Human Breast Cancer dataset through the full Sopa pipeline, then
runs SPATCH custom analysis modules.

**Data location:** `/mnt/shared/janesick/input/outs`
**Output zarr:** `results/janesick.zarr`
**SPATCH output:** `results/janesick_breast_cancer/`

### 2.1 Verify Source Data

```bash
ls /mnt/shared/janesick/input/outs/experiment.xenium
```

### 2.2 Sopa Pipeline (Convert → Patchify → Segment → Aggregate)

Enable Dask parallelism for the session:

```bash
export SOPA_PARALLELIZATION_BACKEND=dask
```

**Step 1 — Convert** Xenium output to SpatialData format:

```bash
sopa convert /mnt/shared/janesick/input/outs \
  --sdata-path results/janesick.zarr \
  --technology xenium
```

Takes ~2 minutes. Creates `results/janesick.zarr`.

**Step 2 — Patchify** the image into overlapping tiles:

```bash
sopa patchify image results/janesick.zarr \
  --patch-width-pixel 2048 \
  --patch-overlap-pixel 100
```

Creates ~266 patches. Takes ~1 minute.

**Step 3 — Segment** cells using Cellpose (GPU required for practical
runtimes):

```bash
sopa segmentation cellpose results/janesick.zarr \
  --diameter 35 \
  --channels DAPI \
  --model-type cyto2 \
  --gpu
```

> **Performance:**
> - GPU + Dask: ~15 minutes (266 patches)
> - GPU only: ~42 minutes
> - CPU only: ~40+ hours — not recommended

**Step 4 — Resolve** boundary conflicts between adjacent patches:
    Look for segmentation info block:  
```[INFO] (sopa.segmentation.stainings) Added 280492 cell boundaries in sdata['cellpose_boundaries']
[INFO] (sopa.cli.utils) Segmentation is already resolved. Don't run `sopa resolve`.
```
If the above info lines are posted no need to run the below commnad
    
```bash
sopa resolve cellpose results/janesick.zarr
```
If this drops a stack trace.  Likely the info block above was over looked.  If it is there just proceed to the next step.  

**Step 5 — Aggregate** transcript counts per cell:

```bash
sopa aggregate results/janesick.zarr --min-transcripts 10
```

Takes ~30 seconds. Result: ~218,990 cells × 313 genes.

### 2.3 SPATCH Custom Modules

The Janesick config runs four modules in order:

1. **dapi_tissue_mask** — tissue boundary from DAPI channel, tags cells `in_tissue`
2. **cell_shape_metrics** — morphological metrics (area, circularity, eccentricity, solidity, aspect ratio)
3. **diffusion_analysis** — quantifies transcript diffusion in/out of tissue
4. **pipeline_visualizations** — scanpy preprocessing + 6 publication figures

Run all modules:

```bash
spatch run results/janesick.zarr \
  --config configs/janesick_breast_cancer.yaml
```

Or equivalently via sopa:

```bash
sopa spatch run results/janesick.zarr \
  --config configs/janesick_breast_cancer.yaml
```

Takes ~2–5 minutes. Outputs:
- **Parquet files:** `results/janesick_breast_cancer/spatch/`
- **Figures:** `results/janesick_breast_cancer/figures/`

To run a single module:

```bash
spatch run results/janesick.zarr \
  -c configs/janesick_breast_cancer.yaml \
  -m dapi_tissue_mask
```

### 2.4 Convenience Script

Runs the full Sopa pipeline + SPATCH modules in one command:

```bash
bash run_janesick_pipeline.sh
```

This executes all steps from Section 2.2 and 2.3 sequentially, including
Dask parallelism and GPU segmentation.

### 2.5 Expected Outputs

**Figures** (in `results/janesick_breast_cancer/figures/`):
- `spatial_scatter.png` — cells on tissue, colored by Leiden cluster
- `umap_leiden.png` — UMAP embedding colored by cluster
- `cluster_composition.png` — bar chart of cell counts per cluster
- `marker_heatmap.png` — mean expression of 8 marker genes per cluster
- `neighborhood_enrichment.png` — spatial co-localization matrix
- `cell_shape_distributions.png` — violin plots of morphology metrics

**Metrics** (typical values):
- ~218,966 cells tagged in-tissue
- ~218,983 cells with shape metrics (mean circularity ~0.86)
- 313 genes analyzed for diffusion

---

## 3. Pipeline 2: CODEX COAD (Colon Adenocarcinoma)

This pipeline uses the SPATCH benchmarking dataset: paired spatial
transcriptomics (405K cells × 5001 genes) and CODEX proteomics
(266K cells × 16 proteins) from colon adenocarcinoma tissue.

Unlike Janesick, this data does **not** go through sopa's convert/segment
steps — it arrives as pre-processed h5ad files that are assembled into
SpatialData format by a preparation script.

**Raw data:** `/mnt/shared/data/codex/`
**Output zarr:** `results/codex.zarr`
**SPATCH output:** `results/codex_coad/`

### 3.1 Verify Source Data

```bash
ls /mnt/shared/data/codex/transcriptome/adata.h5ad
ls /mnt/shared/data/codex/proteome/adata_codex.h5ad
ls /mnt/shared/data/codex/segmentation_mask/cell_boundaries.csv
```

The dataset contains:
- `transcriptome/adata.h5ad` — 405K cells × 5001 genes
- `proteome/adata_codex.h5ad` — 266K cells × 16 proteins
- `segmentation_mask/cell_boundaries.csv` — cell boundary polygons

### 3.2 Load Data into SpatialData

The `codex_loader` module reads the h5ad files and boundary CSV, builds
cell polygons, maps `high_quality` → `in_tissue`, and assembles
everything into a SpatialData zarr:

```
spatch run ~/mike/results/codex.zarr -c ~/SOPA-Spatial/configs/codex_config.yaml -m codex_loader
```
spatch CLI above is reccommened: It is an implementation of custom dataloading  

```bash
spatch run results/codex.zarr \
  -c configs/codex_config.yaml \
  -m codex_loader
```

This auto-detects the h5ad multiomics layout under `/mnt/shared/data/codex`
(configured in `codex_config.yaml`). Takes ~2–5 minutes.

> **Note:** If the zarr already exists, delete it first:
> ```bash
> rm -rf results/codex.zarr
> ```
>
> **Legacy:** The standalone `scripts/prepare_codex_sdata.py` still works
> but is superseded by the `codex_loader` module.

### 3.3 SPATCH Custom Modules

The COAD config runs six modules in order:

1. **cell_shape_metrics** — morphological analysis of 405K cells
2. **diffusion_analysis** — uses `high_quality` flag as tissue proxy
3. **gene_protein_correlation** — cross-modal Spearman correlation between
   ST gene expression and CODEX protein levels at multiple spatial
   resolutions (100–500 µm)
4. **annotation_consensus** — CellTypist automated cell type annotation
   using the `Immune_All_Low.pkl` pre-trained model
5. **spatial_cluster** — CellCharter-based spatial clustering with scVI
   latent embedding (5–10 clusters)
6. **pipeline_visualizations** — preprocessing + 5 figure types

Run all modules:

```bash
spatch run results/codex.zarr \
  --config configs/codex_coad.yaml
```

Takes ~10–20 minutes (annotation and clustering are the slowest steps).

To run specific modules:

```bash
# Run just the cross-modal correlation analysis
spatch run results/codex.zarr \
  -c configs/codex_coad.yaml \
  -m gene_protein_correlation

# Run multiple specific modules
spatch run results/codex.zarr \
  -c configs/codex_coad.yaml \
  -m cell_shape_metrics -m diffusion_analysis
```

### 3.4 Expected Outputs

**Figures** (in `results/codex_coad/figures/`):
- `spatial_scatter.png` — cells on tissue, colored by Leiden cluster
- `umap_leiden.png` — UMAP embedding
- `cluster_composition.png` — cell counts per cluster
- `marker_heatmap.png` — COAD marker expression (EPCAM, CDX2, CD3E, etc.)
- `cell_shape_distributions.png` — morphology by cluster

**Parquet files** (in `results/codex_coad/spatch/`):
- Cell shape metrics for 405K cells
- Diffusion ratios per gene
- Gene-protein correlation (13 pairs, mean Spearman r ~0.48)
- CellTypist annotations (98 cell types)
- CellCharter spatial clusters (5–10 clusters)

---

## 4. Viewing Results in Jupyter

### 4.1 Start Jupyter

From a terminal with the `spatch` conda environment active:

```bash
source setup_environment.sh   # activates env + sets exports
jupyter lab
```

Or, if you are inside a Jupyter container, select the **SPATCH** kernel
from the kernel picker (not the default "Python 3" kernel). If the
kernel is missing:

```bash
source setup_environment.sh
python3 -m ipykernel install --user --name spatch --display-name "SPATCH"
```

### 4.2 Notebook Setup Cell

Every notebook should start with a setup cell that resolves the
project root, adds it to `sys.path`, and changes the working directory
so all relative paths (e.g. `results/janesick.zarr`) resolve correctly:

```python
import sys, os
from pathlib import Path

PROJECT_ROOT = Path('..').resolve()
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

from setup_local_imports import setup_notebook_environment
setup_notebook_environment()
```

> **Why `os.chdir`?** Jupyter sets the working directory to the
> notebook's folder (`notebooks/`), but data lives under
> `results/` at the project root. The `os.chdir` ensures all
> relative paths work without modification.

You should see ✓ marks for all core packages. If anything is missing,
re-run `source setup_environment.sh` in a terminal.

### 4.3 Pre-built Notebooks

Two ready-to-run notebooks are in `notebooks/`:

- **`01_spatch_workflow_example.ipynb`** — End-to-end SPATCH workflow:
  list modules, load data, QC, filtering, PCA/UMAP/Leiden, spatial
  analysis, cell shape metrics, and save results.
- **`02_janesick_breast_cancer_analysis.ipynb`** — Detailed Janesick
  dataset analysis: marker gene scoring, differential expression,
  neighborhood enrichment, SPATCH modules, and parquet export.

Both notebooks include a smart data-path fallback: they try
`results/janesick.zarr` first, then fall back to the shared mount at
`/mnt/shared/hannah-aichelman/janesick.zarr` if the local zarr is
incomplete.

To run: open either notebook in JupyterLab, select the **SPATCH**
kernel, and **Run All Cells**.

### 4.4 Browse Figures

In JupyterLab, navigate to `results/janesick_breast_cancer/figures/` or
`results/codex_coad/figures/` in the file browser and double-click any
PNG to view it.

### 4.5 Display Figures in a Notebook

```python
from IPython.display import Image, display
from pathlib import Path

# Janesick figures
fig_dir = Path("results/janesick_breast_cancer/figures")
for png in sorted(fig_dir.glob("*.png")):
    print(png.name)
    display(Image(filename=str(png), width=800))
```

Replace the path with `results/codex_coad/figures` for COAD figures.

### 4.6 Load and Explore the SpatialData Object

```python
import spatialdata as sd

# Load Janesick dataset
sdata = sd.read_zarr("results/janesick.zarr")
print(sdata)

# Inspect the cell × gene table
adata = sdata.tables["table"]
print(f"Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")
print(adata.obs.columns.tolist())
```

```python
# Load COAD dataset
sdata_coad = sd.read_zarr("results/codex.zarr")
print(sdata_coad)

# Transcriptome table
adata_st = sdata_coad.tables["table"]
print(f"ST: {adata_st.n_obs:,} cells × {adata_st.n_vars:,} genes")

# CODEX protein table
adata_codex = sdata_coad.tables["codex_table"]
print(f"CODEX: {adata_codex.n_obs:,} cells × {adata_codex.n_vars:,} proteins")
```

### 4.7 Load SPATCH Analysis Results from Parquet

Each SPATCH module saves its tabular output as parquet. Load them
directly into pandas:

```python
import pandas as pd

# Cell shape metrics
shapes = pd.read_parquet("results/janesick_breast_cancer/spatch/cell_shape_metrics.parquet")
print(shapes.describe())

# Diffusion analysis
diffusion = pd.read_parquet("results/janesick_breast_cancer/spatch/diffusion_analysis.parquet")
print(diffusion.head())
```

For COAD:

```python
# Gene-protein correlation
corr = pd.read_parquet("results/codex_coad/spatch/gene_protein_correlation.parquet")
print(corr.sort_values("spearman_r", ascending=False))

# CellTypist annotations
annot = pd.read_parquet("results/codex_coad/spatch/annotation_consensus.parquet")
print(annot["cell_type_consensus"].value_counts().head(20))
```

### 4.8 Interactive Visualization with Scanpy

```python
import scanpy as sc
import matplotlib.pyplot as plt

sdata = sd.read_zarr("results/janesick.zarr")
adata = sdata.tables["table"]

# If preprocessing was already done by pipeline_visualizations,
# the object will have PCA, UMAP, and Leiden clusters ready
sc.pl.umap(adata, color="leiden", show=True)
sc.pl.spatial(adata, color="leiden", spot_size=10, show=True)
```

### 4.9 Run Individual SPATCH Modules Interactively

```python
from spatch_modules.runner import run_single_module
import spatialdata as sd

sdata = sd.read_zarr("results/janesick.zarr")

# Run cell shape metrics
sdata = run_single_module(
    sdata, "cell_shape_metrics",
    output_dir="results/janesick_breast_cancer/",
    boundaries_key="cellpose_boundaries",
    table_key="table",
)

# Inspect results
print(sdata.tables["table"].obs[["area", "circularity", "eccentricity"]].describe())
```

### 4.10 Compare Janesick vs COAD Side-by-Side

```python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

axes[0].imshow(mpimg.imread("results/janesick_breast_cancer/figures/umap_leiden.png"))
axes[0].set_title("Janesick Breast Cancer")
axes[0].axis("off")

axes[1].imshow(mpimg.imread("results/codex_coad/figures/umap_leiden.png"))
axes[1].set_title("CODEX COAD")
axes[1].axis("off")

plt.tight_layout()
plt.show()
```

---

## 5. Troubleshooting

**"No module named 'spatialdata'" or import errors:**
Activate the environment: `conda activate spatch`. If using Jupyter,
switch to the **SPATCH** kernel. Re-run `bash setup_environment.sh` if
needed.

**"Zarr directory already exists":**
Delete the existing output: `rm -rf results/janesick.zarr` (or
`results/codex.zarr`) and re-run.

**Segmentation is extremely slow (Janesick):**
Ensure `--gpu` is passed to cellpose and `SOPA_PARALLELIZATION_BACKEND=dask`
is exported. Without GPU, expect 40+ hours.

**CellTypist model download fails (COAD):**
CellTypist downloads the `Immune_All_Low.pkl` model on first run. If
behind a proxy, pre-download:
```bash
python -c "import celltypist; celltypist.models.download_models(model='Immune_All_Low.pkl')"
```

**spatial_cluster fails with CUDA errors (COAD):**
The COAD config sets `use_gpu: false` by default for spatial clustering.
If you have a GPU and want to use it, set `use_gpu: true` in
`configs/codex_coad.yaml` and ensure PyTorch CUDA is installed.

**SPATCH kernel not in Jupyter:**
Re-register:
```bash
conda activate spatch
python3 -m ipykernel install --user --name spatch --display-name "SPATCH"
```

**Disk space issues:**
The environment needs ~15 GB. Use `SPATCH_STORAGE` to redirect to a
larger volume before running setup.

---

## 6. Quick Reference

```bash
# Activate environment
conda activate spatch

# ── Janesick (full pipeline) ──
bash run_janesick_pipeline.sh

# ── Janesick (SPATCH modules only, if sopa already ran) ──
spatch run results/janesick.zarr -c configs/janesick_breast_cancer.yaml

# ── COAD (load data + SPATCH modules) ──
spatch run results/codex.zarr -c configs/codex_config.yaml -m codex_loader
spatch run results/codex.zarr -c configs/codex_coad.yaml

# ── List available modules ──
spatch list

# ── Describe a module ──
spatch describe gene_protein_correlation

# ── Run a single module ──
spatch run results/codex.zarr -c configs/codex_coad.yaml -m annotation_consensus
```
