# SPATCH Pipeline on SpatialData + Sopa

A production-ready spatial transcriptomics analysis pipeline that rebuilds the original SPATCH research scripts on the **SpatialData** (Nature Methods, 2025) and **Sopa** (Nature Communications, 2024) platforms.

## Overview

This project migrates the ~4,000 lines of custom SPATCH analysis code to a maintainable pipeline built on established open-source infrastructure, reducing custom code to ~400 lines across four plug-in modules while preserving all original analytical capabilities.

### Key Benefits

- **88% code reduction**: From ~4,000 lines to ~400 lines of custom code
- **Unified data model**: All data stored in a single Zarr-backed SpatialData object
- **Multi-platform support**: Native loaders for Xenium, Visium HD, CosMx, Stereo-seq, and MERSCOPE
- **Extensible modules**: Easy-to-write custom analysis modules with a uniform interface
- **Production orchestration**: Ready for Snakemake (HPC) or Nextflow (cloud) deployment

## Installation

### Using conda (recommended)

```bash
# Create environment
conda env create -f environment.yml
conda activate spatch

# Install spatch-modules package
pip install -e .
```

### Using pip

```bash
pip install spatch-modules[all]
```

## Quick Start

### 1. Load spatial data

```python
import spatialdata_io as sdio
import spatialdata as sd

# Built-in loaders for most platforms
sdata = sdio.xenium("/path/to/xenium_output/")
# or: sdata = sdio.visium_hd("/path/to/visium_hd_output/")
# or: sdata = sdio.cosmx("/path/to/cosmx_output/")
```

### 2. Run Sopa segmentation and aggregation

```python
import sopa

sopa.make_image_patches(sdata)
sopa.segmentation.cellpose(sdata, channels="DAPI", diameter=30)
sopa.resolve_conflicts(sdata)
sopa.aggregate(sdata)
```

### 3. Run SPATCH custom modules

```python
from spatch_modules import run_single_module

# Diffusion analysis
sdata = run_single_module(sdata, "diffusion_analysis",
    table_key="table",
    in_tissue_col="in_tissue"
)

# Cell shape metrics  
sdata = run_single_module(sdata, "cell_shape_metrics",
    boundaries_key="cell_boundaries"
)
```

### 4. Run full pipeline from config

```python
from spatch_modules import run_custom_pipeline

results = run_custom_pipeline(sdata, "configs/xenium_example.yaml")
```

## Custom Modules

The `spatch_modules` package provides four built-in analysis modules:

| Module | Description |
|--------|-------------|
| `codex_loader` | Load Akoya CODEX/PhenoCycler data into SpatialData format |
| `diffusion_analysis` | Quantify transcript diffusion from in-tissue to out-of-tissue regions |
| `gene_protein_correlation` | Multi-resolution correlation between ST gene expression and CODEX protein levels |
| `cell_shape_metrics` | Compute morphological metrics (area, circularity, eccentricity) from cell boundaries |

### List available modules

```bash
spatch-modules list
```

### Module details

```bash
spatch-modules describe diffusion_analysis
```

## Creating Custom Modules

Create a new module by subclassing `SpatchModule`:

```python
# user_modules/my_analysis.py
from spatch_modules import SpatchModule, ModuleResult, register

@register
class MyAnalysis(SpatchModule):
    name = "my_analysis"
    version = "1.0.0"
    description = "My custom spatial analysis"
    category = "analysis"
    requires = ["tables/table"]
    produces = ["tables/my_results"]

    def run(self, sdata, **kwargs):
        # Your analysis logic here
        return ModuleResult(sdata=sdata, metrics={"result": 42})
```

Place the file in the `user_modules/` directory and it will be automatically discovered.

## Configuration

Pipeline configuration uses YAML files. See `configs/` for examples:

- `spatch_full.yaml` - Complete pipeline with all options documented
- `xenium_example.yaml` - Minimal Xenium workflow
- `multimodal_xenium_codex.yaml` - Multi-modal ST + CODEX analysis

### Example configuration

```yaml
technology: xenium

segmentation:
  method: cellpose
  cellpose:
    diameter: 30
    channels: [DAPI]

custom_modules:
  steps:
    - module: diffusion_analysis
      enabled: true
      config:
        in_tissue_col: in_tissue
    
    - module: cell_shape_metrics
      enabled: true
      config:
        boundaries_key: cell_boundaries
```

## Project Structure

```
multiomics/
├── spatch_modules/          # Custom module library
│   ├── base.py              # SpatchModule ABC + ModuleResult
│   ├── registry.py          # Module discovery and instantiation
│   ├── runner.py            # YAML-driven pipeline execution
│   ├── cli.py               # Command-line interface
│   └── builtin/             # Shipped modules
│       ├── codex_loader.py
│       ├── diffusion_analysis.py
│       ├── gene_protein_correlation.py
│       └── cell_shape_metrics.py
├── configs/                 # Pipeline configurations
├── samplesheets/            # Input manifests
├── user_modules/            # User-contributed modules
├── notebooks/               # Jupyter notebooks
├── tests/                   # Module tests
├── docs/                    # Documentation
├── sopa/                    # Sopa library (cloned)
├── spatialdata/             # SpatialData library (cloned)
├── pyproject.toml           # Package metadata
├── environment.yml          # Conda environment
└── README.md
```

## Multi-Modal Analysis (ST + CODEX)

For paired spatial transcriptomics and protein imaging data:

1. **Load both datasets**:
   ```python
   sdata_st = sdio.xenium("/path/to/xenium/")
   
   from spatch_modules import get_module
   codex_loader = get_module("codex_loader", data_path="/path/to/codex/")
   sdata_codex = codex_loader.run(None).sdata
   ```

2. **Register datasets** using napari-spatialdata landmarks

3. **Run correlation analysis**:
   ```python
   sdata = run_single_module(sdata, "gene_protein_correlation",
       gene_table_key="table",
       protein_table_key="codex_table",
       resolution_um=[100, 200, 300, 400, 500]
   )
   ```

## Workflow Orchestration

### Snakemake (HPC)

```bash
snakemake \
  --config data_path=/data/xenium_001 \
  --configfile=configs/xenium_example.yaml \
  --profile slurm \
  --cores 64
```

### Nextflow (Cloud)

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheets/experiment.csv \
  --outdir results/
```

## License

MIT License

## Citation

If you use this pipeline, please cite:

- **SpatialData**: Marconato et al., Nature Methods (2025)
- **Sopa**: Blampey et al., Nature Communications (2024)
- **SPATCH**: [Original SPATCH publication]

## Contributing

Contributions are welcome! Please see `CONTRIBUTING.md` for guidelines.
