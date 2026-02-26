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

The pipeline expects a Xenium `outs/` directory. Create a symlink or copy
it so the project can find it:

```bash
cd SOPA-Spatial
mkdir -p data
ln -sf /path/to/your/xenium/outs data/outs
```

Replace `/path/to/your/xenium/outs` with the actual path on your system.
For example: `ln -sf /home/user/shared/janesick/input/outs data/outs`

Verify it worked:

```bash
ls data/outs/experiment.xenium
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
```

**Step 1 — Convert** raw Xenium data into the SpatialData format (a zarr
store that all downstream tools understand):

```bash
sopa convert data/outs/ \
  --sdata-path results/janesick.zarr \
  --technology xenium
```

This reads the Xenium `outs/` directory and writes a `results/janesick.zarr`
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

> **⚠️ GPU is strongly recommended.** With `--gpu`: ~10–12 sec/patch
> (~42 minutes total). Without `--gpu`: ~600 sec/patch (~40+ hours).
> If you do not have a GPU, remove `--gpu` but expect a very long run.

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
snakemake \
  --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=data/outs/ sdata_path=results/janesick.zarr \
  --cores 4
```

**What each flag means:**
- `--snakefile` — path to sopa's built-in workflow definition
- `--configfile` — pipeline parameters (cell diameter, channels, etc.)
- `--config data_path=... sdata_path=...` — input data and output location
- `--cores 4` — number of parallel jobs (increase if you have more CPUs)

The config file (`configs/janesick_sopa.yaml`) already includes `gpu: true`
for cellpose. Snakemake will print each step as it runs. The full pipeline
takes roughly 45–60 minutes with a GPU.

> **Note:** Do not add `--use-conda` — all dependencies are already
> installed in your environment. Adding it will cause errors.

**To do a dry-run first** (shows what would happen without running anything):

```bash
snakemake \
  --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=data/outs/ sdata_path=results/janesick.zarr \
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

This runs:
- **diffusion_analysis** — quantifies transcript diffusion from in-tissue to out-of-tissue regions
- **cell_shape_metrics** — computes morphological metrics (area, circularity, eccentricity, etc.)

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

**Kernel not found / SPATCH kernel missing:**
Re-register it:
```bash
conda activate spatch
python3 -m ipykernel install --user --name spatch --display-name "SPATCH"
```

---

## 7. Reference

- **Full setup docs**: `docs/ENVIRONMENT_SETUP.md`
- **Janesick pipeline details**: `docs/JANESICK_PIPELINE_SETUP.md`
- **Sopa config format**: `configs/janesick_sopa.yaml`
- **SPATCH module config**: `configs/janesick_breast_cancer.yaml`
- **Sopa documentation**: https://qupath.github.io/sopa/
- **SpatialData documentation**: https://spatialdata.scverse.org/
