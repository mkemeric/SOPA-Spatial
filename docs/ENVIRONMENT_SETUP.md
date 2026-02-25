# Environment Setup Guide

This guide explains how to install the SPATCH multiomics environment. The project uses stock `sopa` and `spatialdata` from PyPI for the core pipeline, plus a lightweight `spatch_modules` package for custom analysis extensions.

## Quick Start

```bash
cd /Users/mike/projects/University_of_Rochester/multiomics
./setup_environment.sh
```

The script installs:
1. `sopa`, `spatialdata`, `spatialdata-io`, `spatialdata-plot` from PyPI
2. `scanpy`, `squidpy`, `cellpose`, and other analysis dependencies
3. `spatch_modules` (local, editable mode) for custom SPATCH extensions

On macOS with Homebrew Python (PEP 668), a `.venv/` is created automatically. Python 3.12 is recommended (3.13+ may not be supported by upstream packages yet).

## Manual Setup

```bash
python3 -m pip install sopa spatialdata spatialdata-io spatialdata-plot
python3 -m pip install scanpy squidpy cellpose matplotlib seaborn jupyter
python3 -m pip install -e . --no-deps   # spatch_modules
```

## Using in Notebooks

Add this cell at the top of your notebook:

```python
from setup_local_imports import setup_notebook_environment
setup_notebook_environment()
```

This verifies all packages are importable and configures matplotlib/scanpy defaults.

## Verification

```bash
python3 setup_local_imports.py
```

Expected output:
```
======================================================================
Verifying package imports
======================================================================
  spatialdata          ✓  v0.x.x
  spatialdata_io       ✓  v0.x.x
  sopa                 ✓  v2.x.x
  spatch_modules       ✓  v1.0.0
  scanpy               ✓  v1.x.x
  squidpy              ✓  v1.x.x
======================================================================

All packages OK.
```

## Troubleshooting

### `ModuleNotFoundError`

Re-run `./setup_environment.sh`.

### `externally-managed-environment` (PEP 668)

The setup script handles this automatically. If running pip manually, activate the venv first:
```bash
source .venv/bin/activate
```

### Jupyter kernel doesn't see packages

Register the kernel:
```bash
python3 -m ipykernel install --user --name spatch --display-name "Python (SPATCH)"
```

Then select "Python (SPATCH)" kernel in Jupyter.

### Dependency conflicts

Start fresh:
```bash
python3 -m pip install --upgrade sopa spatialdata spatialdata-io
python3 -m pip check
```

## Project Structure

```
multiomics/
├── spatch_modules/          # Custom SPATCH analysis modules (editable install)
│   ├── base.py              # SpatchModule ABC + ModuleResult
│   ├── registry.py          # Module discovery and instantiation
│   ├── runner.py            # YAML-driven pipeline execution
│   ├── cli.py               # Command-line interface
│   └── builtin/             # Built-in modules
├── configs/                 # Pipeline YAML configurations
├── notebooks/               # Jupyter notebooks
├── user_modules/            # User-contributed modules (auto-discovered)
├── pyproject.toml           # spatch_modules package metadata
├── environment.yml          # Conda environment spec
├── setup_environment.sh     # Automated setup script
└── setup_local_imports.py   # Notebook environment helper
```

## Next Steps

1. See `ENVIRONMENT_QUICKSTART.md` for running the Janesick pipeline
2. See `docs/JANESICK_PIPELINE_SETUP.md` for detailed pipeline documentation
3. See `docs/JANESICK_RESULTS_EXPLORATION.md` for exploring results
