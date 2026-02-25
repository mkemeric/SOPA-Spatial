# Environment Quick Start

## Install

```bash
cd /Users/mike/projects/University_of_Rochester/multiomics
./setup_environment.sh
```

This installs stock `sopa`, `spatialdata`, and `spatialdata-io` from PyPI, plus `spatch_modules` in editable mode. On macOS with Homebrew Python, a `.venv/` is created automatically.

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
# Convert Xenium data to SpatialData
sopa convert ../data/outs/ --sdata-path results/janesick.zarr --technology xenium

# Segment, aggregate, preprocess
sopa segmentation cellpose results/janesick.zarr --diameter 35 --channels DAPI
sopa aggregate results/janesick.zarr --min-transcripts 10
```

### Option B: Sopa Snakemake (full automated pipeline)

```bash
snakemake --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=../data/outs/ sdata_path=results/janesick.zarr \
  --cores 4
```

### Option C: Notebook (interactive exploration + SPATCH modules)

Open `notebooks/02_janesick_breast_cancer_analysis.ipynb` and run cells sequentially.

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
