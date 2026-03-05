# SPATCH / Sopa — Quick Start Guide

This guide walks you through setting up and running the SPATCH spatial
transcriptomics pipeline on the Janesick breast-cancer dataset. It assumes
you are working **inside a Jupyter notebook container** and will access
things from either the built-in **terminal** or a **notebook**.

---

## 1. Prerequisites

You need:
- A running Jupyter environment (JupyterHub, JupyterLab, or Notebook)
- Terminal access (every Jupyter environment provides one — look for
  **File → New → Terminal** in JupyterLab, or the **Terminal** button on
  the classic Notebook home page)
- Xenium output data (the `outs/` directory from a 10x Xenium run)
- **GPU access** for cell segmentation (CPU works but is ~50× slower)
- **~15 GB free disk space** for the conda environment and package cache.
  The setup script auto-detects small home directories and redirects
  storage to a larger volume (see [Addendum A](#addendum-a-storage-redirection)).

---

## 2. Install (one-time setup)

Open a **terminal** inside your Jupyter environment and run:

```bash
git clone https://github.com/mkemeric/SOPA-Spatial.git
cd SOPA-Spatial
bash setup_environment.sh
```

**What this does:**
- Detects your Python environment (conda, venv, or system Python)
- If conda/mamba is available, creates a `spatch` conda environment
- Installs all dependencies: `sopa`, `spatialdata`, `scanpy`, `cellpose`,
  `snakemake`, and more
- Installs SPATCH custom analysis modules (`spatch_modules`)
- Registers a **SPATCH** Jupyter kernel so notebooks can find everything
- Skips packages that are already installed (safe to re-run)

**After setup completes**, if the script created a conda environment you
need to activate it each time you open a new terminal:

```bash
conda activate spatch
```

### Verify

From the terminal:

Run the full diagnostic:

```bash
python3 setup_local_imports.py
```

---

## 3. Prepare Your Data

The pipeline expects a Xenium output directory. For the Janesick breast cancer
dataset, the data is located at:

```bash
/mnt/shared/janesick/input/outs
```

No symlink or copying is needed - you can reference this path directly in all
commands below.

Verify the data is accessible:

```bash
ls /mnt/shared/janesick/input/outs/experiment.xenium
```

You should see the file listed without errors.

---

## 4. Configure SPATCH Analysis

The SPATCH config file controls which custom analysis modules run after
the sopa pipeline and how they are parameterized. The default config for
the Janesick dataset is:

```
configs/janesick_breast_cancer.yaml
```

This config enables four modules that align with the sopa pipeline output:

- **dapi_tissue_mask** — uses the DAPI channel (channel 0) from the Xenium
  images to generate a tissue boundary and tag each cell as in-tissue or
  out-of-tissue.  Parameters: `kernel_width`, `dilate_iterations`,
  `erode_iterations`.
- **cell_shape_metrics** — computes morphology metrics (area, circularity,
  eccentricity, solidity, aspect ratio) from the cell boundary polygons
  produced by cellpose segmentation.  Automatically discovers the boundary
  shapes key (e.g. `cellpose_boundaries`).
- **diffusion_analysis** — compares transcript counts inside vs outside the
  tissue mask to quantify per-gene signal diffusion.  Requires
  `dapi_tissue_mask` to run first.  Parameters: `buffer_distance_um`,
  `batch_size`.
- **pipeline_visualizations** — runs scanpy preprocessing if needed
  (normalize, PCA, UMAP, Leiden) and generates six figure types.
  Parameters: `marker_genes`, `leiden_resolution`, `output_dir`,
  `figure_dpi`.

**Customizing the config:**

To disable a module, set `enabled: false` on its step.  To change
marker genes for the heatmap, edit the `marker_genes` list under the
`pipeline_visualizations` step.  Output paths default to
`results/janesick_breast_cancer/figures/` — change `output_dir` if needed.

For a different dataset, copy the config and adjust `data_path`,
`boundaries_key`, channel indices, and marker genes to match your data.

---

## 5. Run the Pipeline

There are three ways to run the pipeline. **Pick whichever fits your
workflow** — all three produce the same results.

### Option A: Sopa CLI (step-by-step from the terminal)

**Best for:** Learning the pipeline, debugging, or running individual
steps manually.

The sopa CLI processes spatial transcriptomics data in four steps:

**Open a terminal** and make sure you are in the project directory with
the correct environment:

```bash
cd SOPA-Spatial
conda activate spatch   # if using conda

# Enable Dask parallel backend for the session (speeds up segmentation)
export SOPA_PARALLELIZATION_BACKEND=dask
```

**Step 1 — Convert**
store that all downstream tools understand):

```bash
sopa convert /mnt/shared/janesick/input/outs \
  --sdata-path results/janesick.zarr \
  --technology xenium
```

This reads the Xenium input directory and writes a `results/janesick.zarr`
folder. Expect this to take 1–3 minutes.

**Step 2 — Patchify** the tissue image into overlapping tiles so cellpose
can segment them in parallel:

```bash
sopa patchify image results/janesick.zarr \
  --patch-width-pixel 2048 \
  --patch-overlap-pixel 100
```

This creates ~250–300 patches (depending on tissue size). Takes ~1 minute.

**Step 3 — Segment** cells using the cellpose deep-learning model:

```bash
sopa segmentation cellpose results/janesick.zarr \
  --diameter 35 \
  --channels DAPI \
  --gpu
```

> **⚠️ GPU is strongly recommended.** With `--gpu` + Dask: ~15 minutes
> total for ~266 patches. Without Dask: ~42 minutes. Without GPU:
> ~600 sec/patch (~40+ hours). If you do not have a GPU, remove
> `--gpu` but expect a very long run.

**Step 4 — Aggregate** transcript counts per cell (filters out cells with
fewer than 10 transcripts):

```bash
sopa aggregate results/janesick.zarr --min-transcripts 10
```

This takes ~30 seconds. When it finishes, `results/janesick.zarr` contains
the full processed dataset ready for analysis.

**Step 5 — Run SPATCH custom modules** (tissue mask, morphology, diffusion
analysis, and visualizations):

```bash
sopa spatch run results/janesick.zarr \
  --config configs/janesick_breast_cancer.yaml
```

This runs four modules in order: `dapi_tissue_mask` → `cell_shape_metrics`
→ `diffusion_analysis` → `pipeline_visualizations`. Each module's tabular
output is saved as a parquet file under `results/janesick_breast_cancer/spatch/`.
Figures are saved to `results/janesick_breast_cancer/figures/`. Takes ~2–5 minutes.

> **Tip:** The convenience script `run_janesick_pipeline.sh` runs all five
> steps above in a single command:
> ```bash
> bash run_janesick_pipeline.sh
> ```

---

### Option B: Snakemake (fully automated pipeline)

**Best for:** Running the complete pipeline hands-off in a single command.
Snakemake automatically chains all four steps (convert → patchify →
segment → aggregate) and handles parallelism.

**Open a terminal:**

```bash
cd SOPA-Spatial
conda activate spatch   # if using conda
```

**Run the pipeline:**

```bash
# Enable Dask parallel backend for faster segmentation
export SOPA_PARALLELIZATION_BACKEND=dask

snakemake \
  --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=/mnt/shared/janesick/input/outs sdata_path=results/janesick.zarr \
  --cores 4
```

**What each flag means:**
- `--snakefile` — path to sopa's built-in workflow definition
- `--configfile` — pipeline parameters (cell diameter, channels, etc.)
- `--config data_path=... sdata_path=...` — input data and output location
- `--cores 4` — number of parallel jobs (increase if you have more CPUs)

The config file (`configs/janesick_sopa.yaml`) already includes `gpu: true`
for cellpose. Snakemake will print each step as it runs. With Dask + GPU,
the full pipeline takes roughly 20–30 minutes.

> **Note:** Do not add `--use-conda` — all dependencies are already
> installed in your environment. Adding it will cause errors.

**To do a dry-run first** (shows what would happen without running anything):

```bash
snakemake \
  --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=/mnt/shared/janesick/input/outs sdata_path=results/janesick.zarr \
  --cores 4 --dry-run
```

**After Snakemake finishes**, run the SPATCH custom modules (tissue mask,
morphology, diffusion analysis, and visualizations):

```bash
sopa spatch run results/janesick.zarr \
  --config configs/janesick_breast_cancer.yaml
```

Figures are saved to `results/janesick_breast_cancer/figures/`.

> **Tip:** The convenience script `run_janesick_snakemake.sh` wraps both
> Snakemake and SPATCH modules into a single command:
> ```bash
> bash run_janesick_snakemake.sh
> ```

---

### Option C: Jupyter Notebook (interactive analysis)

**Best for:** Exploring the data interactively, visualizing results,
and running SPATCH custom analysis modules.

**Step 1 — Open the notebook:**

In JupyterLab, navigate to `SOPA-Spatial/notebooks/` and open
`02_janesick_breast_cancer_analysis.ipynb`.

**Step 2 — Select the correct kernel:**

Click the kernel name in the top-right corner of the notebook and switch
to **SPATCH**. This kernel was registered by `setup_environment.sh` and
has all the required packages.

> If you don't see the SPATCH kernel, run this in a terminal:
> ```bash
> conda activate spatch
> python3 -m ipykernel install --user --name spatch --display-name "SPATCH"
> ```
> Then refresh the notebook page.

**Step 3 — Add the environment setup cell:**

Insert a new cell at the very top of the notebook and run:

```python
from setup_local_imports import setup_notebook_environment
setup_notebook_environment()
```

This verifies all imports and configures matplotlib/scanpy for the
notebook. You should see a list of packages with ✓ marks.

**Step 4 — Run cells sequentially.**

The notebook walks through:
1. Loading Xenium data with `spatialdata_io`
2. Quality control (transcript counts, gene counts, mitochondrial %)
3. Filtering and normalization
4. Clustering (PCA, UMAP, Leiden)
5. Spatial analysis with squidpy
6. SPATCH custom modules (cell shape metrics, diffusion analysis)

Each cell has comments explaining what it does. Run them in order from
top to bottom.

---

## 6. Run SPATCH Custom Modules

After the sopa pipeline completes (via any of the three options above),
you can run the SPATCH custom analysis modules. These add tissue masking,
cell morphology metrics, diffusion analysis, and publication-quality
figures on top of the standard sopa output.

The Janesick config runs four modules in order:
**dapi_tissue_mask** → **cell_shape_metrics** → **diffusion_analysis** →
**pipeline_visualizations**.
See [Addendum B](#addendum-b-spatch-module-details) for module descriptions
and expected output.

### SPATCH CLI (recommended)

The `spatch` and `sopa spatch` CLI commands are installed by
`setup_environment.sh` and are available in every terminal session:

```bash
# Run all enabled SPATCH modules from the config
spatch run results/janesick.zarr --config configs/janesick_breast_cancer.yaml

# Same command via sopa (works identically)
sopa spatch run results/janesick.zarr --config configs/janesick_breast_cancer.yaml

# Run a single module
spatch run results/janesick.zarr -c configs/janesick_breast_cancer.yaml -m dapi_tissue_mask

# Specify a custom output directory for parquet files
spatch run results/janesick.zarr -c configs/janesick_breast_cancer.yaml -o ./results/janesick_breast_cancer/

# List available modules
spatch list

# Show module details
spatch describe pipeline_visualizations
```

**Key flags:**
- `--module / -m <name>` — run only one module (config still read from YAML)
- `--output-dir / -o` — directory for parquet outputs (defaults from config)
- `--config / -c` — path to YAML config file

### Container / agent execution

Each module can run in its own container with externalized data storage.
The `--module` flag enables step-by-step orchestration.  Each module's
tabular output is persisted as a parquet file under `{output_dir}/spatch/`,
so no zarr write-back is needed between steps:

```bash
# Full SPATCH pipeline in one container
docker exec worker sopa spatch run /data/results.zarr \
  --config /data/configs/janesick_breast_cancer.yaml

# Or chain individual modules across containers
docker exec worker sopa spatch run /data/results.zarr \
  -c /data/config.yaml -m dapi_tissue_mask

docker exec worker sopa spatch run /data/results.zarr \
  -c /data/config.yaml -m cell_shape_metrics

docker exec worker sopa spatch run /data/results.zarr \
  -c /data/config.yaml -m diffusion_analysis

docker exec worker sopa spatch run /data/results.zarr \
  -c /data/config.yaml -m pipeline_visualizations
```

The zarr directory (read-only) plus the `spatch/` parquet directory on
the mounted volume provide shared state between containers.  Each module
loads prior outputs from parquet at startup.

### Convenience scripts

The shell scripts wrap the full sopa + SPATCH pipeline into a single command:

```bash
# Sopa CLI steps + SPATCH modules:
bash run_janesick_pipeline.sh

# Snakemake + SPATCH modules:
bash run_janesick_snakemake.sh
```

### Python API (notebooks / interactive)

If you already have a processed zarr from sopa, run SPATCH from Python:

```python
import spatialdata as sd
from spatch_modules.runner import run_custom_pipeline

sdata = sd.read_zarr("results/janesick.zarr")
results = run_custom_pipeline(sdata, "configs/janesick_breast_cancer.yaml")
```

Or run a single module interactively:

```python
from spatch_modules import run_single_module

sdata = run_single_module(sdata, "cell_shape_metrics", output_dir="results/")
```

### Generated figures

The `pipeline_visualizations` module saves six PNG files to
`results/janesick_breast_cancer/figures/`:

- `spatial_scatter.png` — cells on tissue coordinates, colored by Leiden cluster
- `umap_leiden.png` — UMAP embedding colored by Leiden cluster
- `cluster_composition.png` — bar chart of cell counts per cluster
- `marker_heatmap.png` — mean expression of 8 marker genes per cluster
- `neighborhood_enrichment.png` — spatial co-localization matrix (squidpy)
- `cell_shape_distributions.png` — violin plots of morphology metrics by cluster

### Viewing figures in Jupyter

**File browser (quickest):** In JupyterLab, navigate to
`results/janesick_breast_cancer/figures/` in the left sidebar and
double-click any PNG to open it in a viewer tab.

**In a notebook cell:**

```python
from IPython.display import Image, display
from pathlib import Path

fig_dir = Path("results/janesick_breast_cancer/figures")
for png in sorted(fig_dir.glob("*.png")):
    print(png.name)
    display(Image(filename=str(png), width=800))
```

**Side-by-side comparison** (e.g. spatial vs UMAP):

```python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

fig, axes = plt.subplots(1, 2, figsize=(18, 8))
for ax, name in zip(axes, ["spatial_scatter.png", "umap_leiden.png"]):
    ax.imshow(mpimg.imread(f"results/janesick_breast_cancer/figures/{name}"))
    ax.set_title(name.replace(".png", "").replace("_", " ").title())
    ax.axis("off")
plt.tight_layout()
plt.show()
```

---

## 7. Troubleshooting

**"No module named 'spatialdata'" or other import errors:**
Make sure you selected the **SPATCH** kernel (notebooks) or ran
`conda activate spatch` (terminal). If that doesn't help, re-run
`bash setup_environment.sh`.

**"Zarr directory already exists":**
The output zarr from a previous run is still there. Delete it first:
```bash
rm -rf results/janesick.zarr
```

**Segmentation is extremely slow:**
You are running without a GPU. Add `--gpu` to the cellpose command
(Option A) or verify `gpu: true` is in `configs/janesick_sopa.yaml`
(Option B). If no GPU is available, expect 40+ hours for a full slide.

**Kernel not found / SPATCH kernel missing:**
Re-register it:
```bash
conda activate spatch
python3 -m ipykernel install --user --name spatch --display-name "SPATCH"
```

See [Addendum C](#addendum-c-additional-troubleshooting) for less common
issues (disk quota, externally-managed environments, OpenCV, boundary keys).

---

## 8. Reference

- **Full setup docs**: `docs/ENVIRONMENT_SETUP.md`
- **Janesick pipeline details**: `docs/JANESICK_PIPELINE_SETUP.md`
- **Sopa config format**: `configs/janesick_sopa.yaml`
- **SPATCH module config**: `configs/janesick_breast_cancer.yaml`
- **Sopa documentation**: https://qupath.github.io/sopa/
- **SpatialData documentation**: https://spatialdata.scverse.org/

---

## Addendum A: Storage Redirection

Many HPC / container environments have a small home directory quota
(e.g. 10 GB). `setup_environment.sh` automatically detects this and
redirects conda environments, pip cache, and Jupyter data to a larger
volume (`/mnt/user`, `/scratch`, or `$TMPDIR`).

To force a specific location, set `SPATCH_STORAGE` before running setup:

```bash
export SPATCH_STORAGE=/mnt/user   # or any path with ≥15 GB free
bash setup_environment.sh
```

The script writes a `.condarc` so future `conda activate` sessions
automatically use the redirected paths. It also sets `PIP_CACHE_DIR`
and `JUPYTER_DATA_DIR` for the session.

---

## Addendum B: SPATCH Module Details

The custom pipeline (`configs/janesick_breast_cancer.yaml`) runs four
modules in order:

1. **dapi_tissue_mask** — generates a tissue boundary polygon from the
   DAPI image and tags each cell as `in_tissue = 0|1`. This is a
   prerequisite for diffusion analysis. (`opencv-python-headless` is
   installed automatically by `setup_environment.sh`.)
2. **cell_shape_metrics** — computes morphological metrics (area,
   circularity, eccentricity, solidity, aspect ratio) from cell boundary
   polygons. Auto-discovers the boundary shapes key from sopa output
   (e.g. `cellpose_boundaries`).
3. **diffusion_analysis** — compares in-tissue vs out-of-tissue transcript
   counts to quantify signal diffusion per gene.
4. **pipeline_visualizations** — runs scanpy preprocessing (normalize,
   PCA, UMAP, Leiden) if not already present, then generates six figure
   types to `results/janesick_breast_cancer/figures/`:
   `spatial_scatter`, `umap_leiden`, `cluster_composition`,
   `marker_heatmap`, `neighborhood_enrichment`, `cell_shape_distributions`.
   Uses a headless matplotlib backend — no display required.

**Expected output** (Janesick breast cancer dataset):
- ~218,966 cells tagged in-tissue by `dapi_tissue_mask`
- ~218,983 cells with shape metrics (mean circularity ~0.86)
- 313 genes analyzed for diffusion
- 6 PNG figures in `results/janesick_breast_cancer/figures/`

**Dask parallelism note:** The `SOPA_PARALLELIZATION_BACKEND=dask`
setting parallelizes **segmentation** across patches. The other steps
(convert, patchify, aggregate) already use Dask internally for lazy I/O
and memory-efficient processing — no extra configuration needed.

---

## Addendum C: Additional Troubleshooting

**"externally-managed-environment" error:**
The setup script handles this automatically. If you see it when running
pip manually, activate the environment first: `conda activate spatch`
or `source .venv/bin/activate`.

**"Disk quota exceeded" or "No space left on device":**
Your home directory is too small for the full environment (~15 GB).
See [Addendum A](#addendum-a-storage-redirection) for the
`SPATCH_STORAGE` override.

**"opencv-python-headless is required" (dapi_tissue_mask):**
This is installed automatically by `setup_environment.sh`. If you see
this error, re-run the setup script or install manually:
```bash
conda activate spatch
pip install opencv-python-headless
```

**cell_shape_metrics skipped — no cell_boundaries:**
Different segmentation methods store boundaries under different keys
(e.g. `cellpose_boundaries`). The module auto-discovers any shapes key
containing "boundaries". If it still fails, check your zarr's available
shapes keys and set `boundaries_key` in
`configs/janesick_breast_cancer.yaml` accordingly.
