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

# ── Storage detection ───────────────────────────────────────────
# Many HPC / container environments give users a small home quota
# (e.g. 10 GB) that fills up fast with conda envs + pip cache.
# This block checks available space and redirects heavy-weight
# directories to a larger volume when necessary.
#
# Override:  export SPATCH_STORAGE=/path/to/big/disk
#            before running this script.

MIN_HOME_GB=15   # minimum free space required on $HOME

home_free_kb=$(df -Pk "$HOME" 2>/dev/null | awk 'NR==2{print $4}')
home_free_gb=$(( ${home_free_kb:-0} / 1048576 ))

if [[ -n "$SPATCH_STORAGE" ]]; then
    STORAGE_ROOT="$SPATCH_STORAGE"
    echo "Using SPATCH_STORAGE override: $STORAGE_ROOT"
elif [[ $home_free_gb -lt $MIN_HOME_GB ]]; then
    # Auto-detect a larger mount
    for candidate in /mnt/user /scratch "$TMPDIR"; do
        if [[ -d "$candidate" && -w "$candidate" ]]; then
            cand_free_kb=$(df -Pk "$candidate" 2>/dev/null | awk 'NR==2{print $4}')
            cand_free_gb=$(( ${cand_free_kb:-0} / 1048576 ))
            if [[ $cand_free_gb -ge $MIN_HOME_GB ]]; then
                STORAGE_ROOT="$candidate"
                break
            fi
        fi
    done

    if [[ -n "$STORAGE_ROOT" ]]; then
        echo "⚠️  Home directory has only ${home_free_gb} GB free (need ${MIN_HOME_GB} GB)."
        echo "   Redirecting conda/pip storage → $STORAGE_ROOT"
    else
        echo "⚠️  Home directory has only ${home_free_gb} GB free and no suitable"
        echo "   alternate volume was found (/mnt/user, /scratch, \$TMPDIR)."
        echo "   Set SPATCH_STORAGE=/path/to/big/disk and re-run, or free space in \$HOME."
        exit 1
    fi
else
    STORAGE_ROOT=""   # home is big enough — use defaults
fi

if [[ -n "$STORAGE_ROOT" ]]; then
    # Conda environments & package cache
    export CONDA_ENVS_PATH="$STORAGE_ROOT/conda_envs"
    export CONDA_PKGS_DIRS="$STORAGE_ROOT/conda_pkgs"
    mkdir -p "$CONDA_ENVS_PATH" "$CONDA_PKGS_DIRS"

    # Persist for future shells (idempotent)
    CONDARC="$HOME/.condarc"
    if ! grep -q "$STORAGE_ROOT/conda_envs" "$CONDARC" 2>/dev/null; then
        cat > "$CONDARC" <<YAML
envs_dirs:
  - $STORAGE_ROOT/conda_envs
pkgs_dirs:
  - $STORAGE_ROOT/conda_pkgs
YAML
        echo "   Wrote $CONDARC"
    fi

    # Pip cache
    export PIP_CACHE_DIR="$STORAGE_ROOT/pip_cache"
    mkdir -p "$PIP_CACHE_DIR"

    # Jupyter data (kernelspecs, etc.)
    export JUPYTER_DATA_DIR="$STORAGE_ROOT/jupyter_data"
    mkdir -p "$JUPYTER_DATA_DIR"

    echo "   CONDA_ENVS_PATH = $CONDA_ENVS_PATH"
    echo "   CONDA_PKGS_DIRS = $CONDA_PKGS_DIRS"
    echo "   PIP_CACHE_DIR   = $PIP_CACHE_DIR"
    echo "   JUPYTER_DATA_DIR= $JUPYTER_DATA_DIR"
    echo ""
fi

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

# ── Helper: filter out already-satisfied packages ───────────────
# Uses importlib.metadata to check installed versions against specs.
# Returns only the package specs that are missing or outdated.
filter_missing() {
    local specs=("$@")
    local missing=()
    for spec in "${specs[@]}"; do
        if ! $PYTHON -c "
import sys
from importlib.metadata import version, PackageNotFoundError
from packaging.requirements import Requirement
from packaging.version import Version
req = Requirement('''$spec''')
try:
    v = Version(version(req.name))
except PackageNotFoundError:
    sys.exit(1)
if req.specifier and not req.specifier.contains(v):
    sys.exit(1)
" 2>/dev/null; then
            missing+=("$spec")
        fi
    done
    # Return the list via stdout, one per line
    printf '%s\n' "${missing[@]}"
}

pip_install_missing() {
    local label="$1"; shift
    local specs=("$@")
    local to_install=()
    while IFS= read -r line; do
        [[ -n "$line" ]] && to_install+=("$line")
    done < <(filter_missing "${specs[@]}")

    if [ ${#to_install[@]} -eq 0 ]; then
        echo "✓ $label — all packages already satisfied, skipping."
    else
        echo "Installing ${#to_install[@]} package(s) for $label..."
        printf '  %s\n' "${to_install[@]}"
        $PIP install "${to_install[@]}"
        echo "✓ $label installed"
    fi
}

echo ""
echo "========================================"
echo "Step 1: Installing core spatial packages"
echo "========================================"
# Pin numpy<2 — many compiled deps (pyarrow, scikit-image, etc.) are
# built against the numpy 1.x ABI and crash with numpy 2.x.
pip_install_missing "NumPy compatibility" "numpy>=1.24,<2"

pip_install_missing "Core spatial packages" \
    "spatialdata>=0.6" \
    "spatialdata-io>=0.2" \
    "spatialdata-plot>=0.2.14,<0.3"

echo ""
echo "========================================"
echo "Step 2: Installing analysis dependencies"
echo "========================================"
pip_install_missing "Analysis packages" \
    "scanpy>=1.10" \
    "squidpy>=1.5" \
    "cellpose>=3.0" \
    "snakemake>=8.0" \
    "boto3>=1.35.40" \
    opencv-python-headless \
    matplotlib \
    seaborn \
    scikit-learn \
    scikit-misc \
    scipy \
    leidenalg \
    pyyaml \
    jupyter \
    ipykernel

echo ""
echo "========================================"
echo "Step 3a: Installing sopa from PyPI"
echo "========================================"
pip_install_missing "sopa" "sopa>=2.1"

# Clone the upstream sopa repo (shallow) for Snakemake workflow files.
# The workflow/ directory is not distributed in the PyPI wheel.
if [ ! -d "sopa/workflow/Snakefile" ]; then
    echo "Cloning sopa repo (shallow) for Snakemake workflow files..."
    rm -rf sopa
    git clone --depth 1 https://github.com/gustaveroussy/sopa.git sopa
    echo "✓ sopa workflow files available at sopa/workflow/"
else
    echo "✓ sopa workflow files already present"
fi

echo ""
echo "========================================"
echo "Step 3b: Installing SPATCH modules"
echo "========================================"
if [ -f "pyproject.toml" ]; then
    # Use pip show (not import) — import can succeed via CWD on sys.path
    # but Jupyter kernels run from a different directory.
    if $PIP show spatch_modules >/dev/null 2>&1; then
        echo "✓ spatch_modules already pip-installed, skipping."
    else
        echo "Installing spatch_modules in editable mode..."
        $PIP install -e . --no-deps
        echo "✓ spatch_modules installed"
    fi
else
    echo "⚠️  WARNING: pyproject.toml not found"
fi

# ── Step 3c: Fix entry-point shebangs ──────────────────────
# When conda envs are created in a non-default path (via storage
# redirection), pip may bake shebang lines that don't resolve.
# This fixes any broken shebangs in the env's bin/ directory.
if [[ -n "$STORAGE_ROOT" ]]; then
    ENV_PREFIX=$($PYTHON -c "import sys; print(sys.prefix)")
    CORRECT_PYTHON="$ENV_PREFIX/bin/python3"
    if [[ -x "$CORRECT_PYTHON" ]]; then
        fixed=0
        for script in "$ENV_PREFIX"/bin/*; do
            [[ -f "$script" && -x "$script" ]] || continue
            head_line=$(head -1 "$script" 2>/dev/null)
            case "$head_line" in
                "#!"*python*)
                    shebang_path=${head_line#\#\!}
                    shebang_path=${shebang_path%% *}   # strip any args
                    if [[ ! -x "$shebang_path" ]]; then
                        sed -i "1s|.*|#!${CORRECT_PYTHON}|" "$script"
                        fixed=$((fixed + 1))
                    fi
                    ;;
            esac
        done
        if [[ $fixed -gt 0 ]]; then
            echo "✓ Fixed $fixed entry-point shebang(s) in $ENV_PREFIX/bin/"
        fi
    fi
fi

# ── Step 3d: Register Jupyter kernel ────────────────────────
# Ensures the current environment is available as a Jupyter kernel
# so notebooks can find all installed packages.
if command -v jupyter >/dev/null 2>&1; then
    KERNEL_NAME="spatch"
    echo ""
    echo "========================================"
    echo "Step 3c: Registering Jupyter kernel"
    echo "========================================"
    $PYTHON -m ipykernel install --user --name "$KERNEL_NAME" --display-name "SPATCH" 2>/dev/null \
        && echo "✓ Jupyter kernel '$KERNEL_NAME' registered" \
        || echo "⚠️  Could not register Jupyter kernel (non-fatal)"
fi

echo ""
echo "========================================"
echo "Step 4: Verifying installation"
echo "========================================"
$PYTHON -c "
import sys
import importlib
import importlib.metadata
import subprocess

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
        version = importlib.metadata.version(pkg_name)
        print(f'{pkg_name:20s} \u2713  v{version}')
    except importlib.metadata.PackageNotFoundError:
        # Bare 'import' can give false positives for local directories
        # (namespace packages). Only trust pip metadata.
        print(f'{pkg_name:20s} \u2717  not installed (pip show {pkg_name} failed)')
        all_ok = False

print('-' * 70)

# Extra check: verify sopa CLI entry point was registered
import shutil
if shutil.which('sopa'):
    print('sopa CLI:            \u2713  ' + shutil.which('sopa'))
else:
    print('sopa CLI:            \u2717  not found in PATH')
    all_ok = False

print('-' * 70)
if all_ok:
    print('All packages installed and verified!')
else:
    print('Some packages failed verification.')
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
