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
- **~15 GB free disk space** for the conda environment and package cache

> **Storage note:** Many HPC / container environments have a small home
> directory quota (e.g. 10 GB). The setup script automatically detects
> this and redirects conda environments, pip cache, and Jupyter data to
> a larger volume (`/mnt/user`, `/scratch`, or `$TMPDIR`). To force a
> specific location, set `SPATCH_STORAGE` before running the installer:
>
> ```bash
> export SPATCH_STORAGE=/mnt/user   # or any path with ≥15 GB free
> bash setup_environment.sh
> ```

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

```bash
python3 -c 'import sopa; import spatch_modules; print("All imports OK!")'
```

Or run the full diagnostic:

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

## 4. Run the Pipeline

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

> **Note on parallelism:** The `SOPA_PARALLELIZATION_BACKEND=dask` setting
> parallelizes **segmentation** across patches. The other steps (convert,
> patchify, aggregate) already use Dask internally for lazy I/O and
> memory-efficient processing — no extra configuration needed.

**Step 1 — Convert** raw Xenium data into the SpatialData format (a zarr
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

## 5. Run SPATCH Custom Modules

After the sopa pipeline completes (via any of the three options above),
you can run the SPATCH custom analysis modules. These work in both
terminal Python and notebooks.

```python
import spatialdata as sd
from spatch_modules.runner import run_custom_pipeline

sdata = sd.read_zarr("results/janesick.zarr")
results = run_custom_pipeline(sdata, "configs/janesick_breast_cancer.yaml")
```

The pipeline config (`configs/janesick_breast_cancer.yaml`) runs three
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

**Expected output** (Janesick breast cancer dataset):
- ~218,966 cells tagged in-tissue by `dapi_tissue_mask`
- ~218,983 cells with shape metrics (mean circularity ~0.86)
- 313 genes analyzed for diffusion

You can also run individual modules interactively:

```python
from spatch_modules import run_single_module

sdata = run_single_module(sdata, "cell_shape_metrics", output_dir="results/")
```

---

## 6. Troubleshooting

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

**"externally-managed-environment" error:**
The setup script handles this automatically. If you see it when running
pip manually, activate the environment first: `conda activate spatch`
or `source .venv/bin/activate`.

**"Disk quota exceeded" or "No space left on device":**
Your home directory is too small for the full environment (~15 GB).
Set `SPATCH_STORAGE` to a path with more space and re-run:
```bash
export SPATCH_STORAGE=/mnt/user   # adjust to your system
bash setup_environment.sh
```
The script writes a `.condarc` so future `conda activate` sessions
automatically use the redirected paths.

**Kernel not found / SPATCH kernel missing:**
Re-register it:
```bash
conda activate spatch
python3 -m ipykernel install --user --name spatch --display-name "SPATCH"
```

**"opencv-python-headless is required" (dapi_tissue_mask):**
The DAPI tissue mask module needs OpenCV. Install it:
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

---

## 7. Reference

- **Full setup docs**: `docs/ENVIRONMENT_SETUP.md`
- **Janesick pipeline details**: `docs/JANESICK_PIPELINE_SETUP.md`
- **Sopa config format**: `configs/janesick_sopa.yaml`
- **SPATCH module config**: `configs/janesick_breast_cancer.yaml`
- **Sopa documentation**: https://qupath.github.io/sopa/
- **SpatialData documentation**: https://spatialdata.scverse.org/
