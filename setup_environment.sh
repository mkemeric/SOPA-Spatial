#!/bin/bash
# Setup script for SPATCH multiomics environment
# Installs stock sopa + spatialdata from PyPI, plus spatch_modules in editable mode.

set -e  # Exit on error

echo "========================================"
echo "SPATCH Environment Setup"
echo "========================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Working directory: $SCRIPT_DIR"
echo ""

# ── Determine Python and environment strategy ──────────────────
# Priority: conda/mamba > existing venv > system python + auto-venv

PYTHON=python3
PIP="$PYTHON -m pip"
USING_CONDA=false

# Option 1: If conda/mamba is available, create or reuse a conda env.
# This is the preferred path for remote/HPC environments and enables
# Snakemake's --use-conda for parallel pipeline execution.
if command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; then
    CONDA_CMD=$(command -v mamba 2>/dev/null || echo conda)
    ENV_NAME="spatch"

    if conda env list 2>/dev/null | grep -q "^$ENV_NAME "; then
        echo "Activating existing conda environment: $ENV_NAME"
    else
        echo "Creating conda environment: $ENV_NAME (Python 3.12)"
        $CONDA_CMD create -y -n "$ENV_NAME" python=3.12
    fi

    # Activate (works in both bash and zsh)
    eval "$(conda shell.$(basename $SHELL) hook)"
    conda activate "$ENV_NAME"
    PYTHON=python3
    PIP="$PYTHON -m pip"
    USING_CONDA=true
    echo "Conda env: $CONDA_DEFAULT_ENV"

# Option 2: Already inside a venv or conda env — use it as-is.
elif [[ -n "${VIRTUAL_ENV}" ]] || [[ -n "${CONDA_DEFAULT_ENV}" ]]; then
    echo "Using active environment: ${VIRTUAL_ENV:-$CONDA_DEFAULT_ENV}"

# Option 3: System python. Check version and PEP 668.
else
    # Prefer Python 3.12 (many scientific packages are <3.13 today)
    DETECTED_VER=$($PYTHON -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")' 2>/dev/null || echo "0.0")
    if [[ "$DETECTED_VER" == 3.* ]] && [[ ${DETECTED_VER#3.} -ge 13 ]]; then
      if command -v python3.12 >/dev/null 2>&1; then
        echo "Detected Python $DETECTED_VER; switching to python3.12 for package compatibility"
        PYTHON=python3.12
        PIP="$PYTHON -m pip"
      else
        echo "⚠️  Detected Python $DETECTED_VER, which is too new for some dependencies."
        echo "   Install Python 3.12: brew install python@3.12  (or install conda/mamba)"
        exit 1
      fi
    fi

    # If externally-managed (Homebrew macOS), auto-create a venv.
    EXTERNALLY_MANAGED=$($PYTHON -c "import sysconfig, os; marker=os.path.join(sysconfig.get_path('stdlib'), 'EXTERNALLY-MANAGED'); print('yes' if os.path.exists(marker) else 'no')" 2>/dev/null || echo "no")

    if [[ "$EXTERNALLY_MANAGED" == "yes" ]]; then
        VENV_DIR="$SCRIPT_DIR/.venv"
        if [ ! -d "$VENV_DIR" ]; then
            echo "Detected externally-managed Python (PEP 668) — creating venv at .venv/"
            $PYTHON -m venv "$VENV_DIR"
        else
            echo "Using existing venv at .venv/"
        fi
        source "$VENV_DIR/bin/activate"
        PYTHON=python3
        PIP="$PYTHON -m pip"
        echo "Activated venv: $VIRTUAL_ENV"
    fi
fi

echo "Checking Python installation..."
$PYTHON --version

echo ""
echo "========================================"
echo "Step 1: Installing core spatial packages"
echo "========================================"
echo "Installing sopa, spatialdata, and spatialdata-io from PyPI..."
$PIP install \
    "sopa>=2.1" \
    "spatialdata>=0.6" \
    "spatialdata-io>=0.2" \
    "spatialdata-plot>=0.2.14,<0.3"
echo "✓ Core spatial packages installed"

echo ""
echo "========================================"
echo "Step 2: Installing analysis dependencies"
echo "========================================"
$PIP install \
    "scanpy>=1.10" \
    "squidpy>=1.5" \
    "cellpose>=3.0" \
    "snakemake>=8.0" \
    matplotlib \
    seaborn \
    scikit-learn \
    scipy \
    pyyaml \
    jupyter \
    ipykernel
echo "✓ Analysis packages installed"

echo ""
echo "========================================"
echo "Step 3: Installing SPATCH modules"
echo "========================================"
if [ -f "pyproject.toml" ]; then
    echo "Installing spatch_modules in editable mode..."
    $PIP install -e . --no-deps
    echo "✓ spatch_modules installed"
else
    echo "⚠️  WARNING: pyproject.toml not found"
fi

echo ""
echo "========================================"
echo "Step 4: Verifying installation"
echo "========================================"
$PYTHON -c "
import sys
import importlib

packages = [
    ('spatialdata', 'sd'),
    ('spatialdata_io', 'sdio'),
    ('sopa', 'sopa'),
    ('spatch_modules', None),
    ('scanpy', 'sc'),
    ('squidpy', 'sq')
]

print('\nPackage verification:')
print('-' * 70)
all_ok = True
for pkg_name, alias in packages:
    try:
        mod = importlib.import_module(pkg_name)
        version = getattr(mod, '__version__', '?')
        print(f'{pkg_name:20s} ✓  v{version}')
    except ImportError as e:
        print(f'{pkg_name:20s} ✗  {str(e)}')
        all_ok = False

print('-' * 70)
if all_ok:
    print('All packages imported successfully!')
else:
    print('Some packages failed to import.')
    sys.exit(1)
"

echo ""
echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo ""
if [[ "$USING_CONDA" == "true" ]]; then
    echo "Conda environment: $CONDA_DEFAULT_ENV"
    echo "  Activate with:  conda activate spatch"
elif [[ -n "$VIRTUAL_ENV" ]]; then
    echo "Virtual environment: $VIRTUAL_ENV"
    echo "  Activate with:  source .venv/bin/activate"
fi
echo ""
echo "Run the Janesick pipeline:"
echo "  CLI (local):      ./run_janesick_pipeline.sh"
if [[ "$USING_CONDA" == "true" ]]; then
    echo "  Snakemake (16x):  ./run_janesick_snakemake.sh"
fi
echo ""
echo "Start Jupyter:"
echo "  jupyter notebook"
echo ""
