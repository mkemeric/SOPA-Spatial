# SPATCH Pipeline on SpatialData + Sopa

A production-ready spatial transcriptomics analysis pipeline that rebuilds the original SPATCH research scripts on the **SpatialData** (Nature Methods, 2025) and **Sopa** (Nature Communications, 2024) platforms.

## Overview

This project uses stock **Sopa** and **SpatialData** from PyPI for the core spatial transcriptomics pipeline (data loading, segmentation, aggregation, preprocessing), plus a lightweight `spatch_modules` package (~400 lines) for custom analysis extensions.

### Key Benefits

- **Config-driven**: YAML configs drive the entire pipeline — no custom code needed for standard workflows
- **Stock packages**: sopa, spatialdata, scanpy, squidpy all installed from PyPI
- **Sopa pipeline**: Snakemake workflow + CLI for end-to-end processing
- **Extensible**: Custom SPATCH modules plug in for domain-specific analysis
- **Multi-platform**: Native support for Xenium, Visium HD, CosMx, MERSCOPE, and more

## Installation

```bash
cd /path/to/multiomics
./setup_environment.sh
```

This installs stock sopa + spatialdata from PyPI and `spatch_modules` in editable mode. See `ENVIRONMENT_QUICKSTART.md` for details.

## Quick Start

### 1. Run the sopa pipeline (CLI)

```bash
# Convert Xenium data to SpatialData
sopa convert /path/to/xenium_output/ --sdata-path results/output.zarr --technology xenium

# Segment with Cellpose
sopa segmentation cellpose results/output.zarr --diameter 30 --channels DAPI

# Aggregate transcripts per cell
sopa aggregate results/output.zarr --min-transcripts 10
```

### 2. Or use Snakemake (fully automated)

```bash
snakemake --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=/path/to/xenium_output/ \
  --cores 4
```

### 3. Run SPATCH custom modules

```python
import spatialdata as sd
from spatch_modules.runner import run_custom_pipeline

sdata = sd.read_zarr("results/output.zarr")
results = run_custom_pipeline(sdata, "configs/janesick_breast_cancer.yaml")
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

- `janesick_sopa.yaml` - Sopa-native config for Janesick breast cancer data
- `janesick_breast_cancer.yaml` - Full config including SPATCH custom modules
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
├── spatch_modules/          # Custom SPATCH analysis modules
│   ├── base.py              # SpatchModule ABC + ModuleResult
│   ├── registry.py          # Module discovery and instantiation
│   ├── runner.py            # YAML-driven pipeline execution
│   ├── cli.py               # Command-line interface
│   └── builtin/             # Built-in modules
│       ├── codex_loader.py
│       ├── diffusion_analysis.py
│       ├── gene_protein_correlation.py
│       └── cell_shape_metrics.py
├── configs/                 # Pipeline YAML configurations
│   ├── janesick_sopa.yaml           # Sopa-native config for Janesick data
│   ├── janesick_breast_cancer.yaml  # Full config including SPATCH modules
│   ├── spatch_full.yaml             # Complete annotated reference config
│   └── xenium_example.yaml          # Minimal Xenium example
├── notebooks/               # Jupyter notebooks
├── user_modules/            # User-contributed modules (auto-discovered)
├── tests/                   # Module tests
├── docs/                    # Documentation
├── pyproject.toml           # spatch_modules package metadata
├── environment.yml          # Conda environment spec
├── setup_environment.sh     # Automated setup script
└── setup_local_imports.py   # Notebook environment helper
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

Sopa ships a Snakemake workflow that handles the full pipeline. See `configs/janesick_sopa.yaml` for the config format.

```bash
snakemake --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=/data/xenium_001 \
  --cores 4
```

For HPC, add a Snakemake profile (e.g. `--profile slurm`).

## License

MIT License

## Citation

If you use this pipeline, please cite:

- **SpatialData**: Marconato et al., Nature Methods (2025)
- **Sopa**: Blampey et al., Nature Communications (2024)
- **SPATCH**: [Original SPATCH publication]

## Contributing

Contributions are welcome! Please see `CONTRIBUTING.md` for guidelines.
