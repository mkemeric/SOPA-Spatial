# Environment Quick Start

## Install

```bash
git clone https://github.com/mkemeric/SOPA-Spatial.git
cd SOPA-Spatial
./setup_environment.sh
```

This installs `sopa`, `spatialdata`, `spatialdata-io`, analysis dependencies, and `spatch_modules` in editable mode. On conda/mamba systems a `spatch` environment is created; on macOS with Homebrew Python, a `.venv/` is created automatically.

## Verify

```bash
python3 -c 'import sopa; import spatch_modules; print("All imports OK!")'
```

Or run the full diagnostic:

```bash
python3 setup_local_imports.py
```

## Run the Janesick Pipeline

### Option A: Sopa CLI (standard pipeline)

```bash
# 1. Convert Xenium data to SpatialData
sopa convert ../data/outs/ --sdata-path results/janesick.zarr --technology xenium

# 2. Create image patches (required before segmentation)
sopa patchify image results/janesick.zarr --patch-width-pixel 2048 --patch-overlap-pixel 100

# 3. Segment with cellpose (--gpu strongly recommended, ~50x faster)
sopa segmentation cellpose results/janesick.zarr --diameter 35 --channels DAPI --gpu

# 4. Aggregate transcript counts per cell
sopa aggregate results/janesick.zarr --min-transcripts 10
```

> **GPU note:** Without `--gpu`, cellpose runs ~600s/patch (~40+ hours for a full slide).
> With `--gpu`, it drops to ~10-12s/patch (~42 minutes).

### Option B: Sopa Snakemake (full automated pipeline)

```bash
snakemake --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=../data/outs/ sdata_path=results/janesick.zarr \
  --cores 4
```

> **Note:** Do not use `--use-conda` — all dependencies are already
> installed in the current environment via `setup_environment.sh`.

### Option C: Notebook (interactive exploration + SPATCH modules)

Open `notebooks/02_janesick_breast_cancer_analysis.ipynb` and run cells sequentially.
Select the **SPATCH** kernel (registered by `setup_environment.sh`).

Add this cell at the top of the notebook:

```python
from setup_local_imports import setup_notebook_environment
setup_notebook_environment()
```

## Run SPATCH Custom Modules

After the sopa pipeline completes, run the custom analysis modules:

```python
import spatialdata as sd
from spatch_modules.runner import run_custom_pipeline

sdata = sd.read_zarr("results/janesick.zarr")
results = run_custom_pipeline(sdata, "configs/janesick_breast_cancer.yaml")
```

This runs the two Janesick-specific SPATCH modules:
- **diffusion_analysis**: Quantifies transcript diffusion from in-tissue to out-of-tissue regions
- **cell_shape_metrics**: Computes morphological metrics (area, circularity, eccentricity, etc.)

## Troubleshooting

**Import errors?** Re-run `./setup_environment.sh`.

**`externally-managed-environment` error?** The script handles this automatically by creating a `.venv/`. If you see this error running pip manually, activate the venv first: `source .venv/bin/activate`.

## More Information

- **Full setup docs**: `docs/ENVIRONMENT_SETUP.md`
- **Janesick pipeline details**: `docs/JANESICK_PIPELINE_SETUP.md`
- **Sopa config format**: `configs/janesick_sopa.yaml`
- **SPATCH module config**: `configs/janesick_breast_cancer.yaml`
