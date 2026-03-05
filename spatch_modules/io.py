"""
Parquet I/O for SPATCH module outputs.

Instead of writing back to the zarr (which fails with Dask locks),
each module's tabular output is saved as a separate parquet file.
At visualization time (or any downstream step), these files are
loaded and joined onto the base table from the zarr.

File layout inside ``{output_dir}/spatch/``::

    dapi_tissue_mask.parquet       # obs columns: in_tissue
    cell_shape_metrics.parquet     # obs columns: area_um2, circularity, …
    table_diffusion_metrics.parquet  # standalone gene-level table
    manifest.json                  # bookkeeping: what ran, when

Convention:
  - Files without the ``table_`` prefix contain **obs-level** columns,
    indexed by the same cell IDs as ``sdata.tables["table"]``.
  - Files with the ``table_`` prefix are **standalone tables**
    (e.g. per-gene metrics) and are not joined onto obs.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SPATCH_SUBDIR = "spatch"
MANIFEST_FILE = "manifest.json"


# ── directory helpers ────────────────────────────────────────────────

def spatch_output_dir(output_dir: str | Path) -> Path:
    """Return the ``spatch/`` subdirectory, creating it if needed."""
    d = Path(output_dir) / SPATCH_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# ── save ─────────────────────────────────────────────────────────────

def save_obs_columns(
    obs: pd.DataFrame,
    columns: list[str],
    module_name: str,
    output_dir: str | Path,
) -> Path | None:
    """Save specific obs columns to a parquet file keyed by *module_name*.

    Args:
        obs: The full ``adata.obs`` DataFrame.
        columns: Column names that this module produced.
        module_name: Used as the parquet filename stem.
        output_dir: Top-level output directory (the ``spatch/`` subdir
            is appended automatically).

    Returns:
        Path to the written file, or ``None`` if *columns* was empty.
    """
    if not columns:
        return None
    d = spatch_output_dir(output_dir)
    path = d / f"{module_name}.parquet"
    obs[columns].to_parquet(path, engine="pyarrow")
    return path


def save_table(
    table_obs: pd.DataFrame,
    table_name: str,
    output_dir: str | Path,
) -> Path:
    """Save a standalone table (e.g. diffusion_metrics) as parquet.

    The file is prefixed with ``table_`` so the loader can distinguish
    it from obs-level parquets.
    """
    d = spatch_output_dir(output_dir)
    path = d / f"table_{table_name}.parquet"
    table_obs.to_parquet(path, engine="pyarrow")
    return path


# ── load ─────────────────────────────────────────────────────────────

def load_obs_parquets(output_dir: str | Path) -> pd.DataFrame | None:
    """Load all obs-level parquet files and left-join them.

    Returns:
        A DataFrame of extra columns (indexed by cell ID) ready to be
        joined onto ``sdata.tables["table"].obs``, or ``None`` if no
        parquet files are present.
    """
    d = Path(output_dir) / SPATCH_SUBDIR
    if not d.exists():
        return None

    frames: list[pd.DataFrame] = []
    for path in sorted(d.glob("*.parquet")):
        # Skip standalone-table files
        if path.stem.startswith("table_"):
            continue
        try:
            frames.append(pd.read_parquet(path, engine="pyarrow"))
        except Exception:
            continue

    if not frames:
        return None

    result = frames[0]
    for df in frames[1:]:
        new_cols = [c for c in df.columns if c not in result.columns]
        if new_cols:
            result = result.join(df[new_cols], how="left")
    return result


def load_table_parquet(
    table_name: str,
    output_dir: str | Path,
) -> pd.DataFrame | None:
    """Load a standalone table parquet (e.g. ``diffusion_metrics``)."""
    path = Path(output_dir) / SPATCH_SUBDIR / f"table_{table_name}.parquet"
    if not path.exists():
        return None
    return pd.read_parquet(path, engine="pyarrow")


def list_saved_outputs(output_dir: str | Path) -> dict[str, list[str]]:
    """Return a summary of saved SPATCH outputs.

    Returns:
        ``{"obs": ["dapi_tissue_mask", ...], "tables": ["diffusion_metrics", ...]}``
    """
    d = Path(output_dir) / SPATCH_SUBDIR
    obs_names: list[str] = []
    table_names: list[str] = []

    if d.exists():
        for path in sorted(d.glob("*.parquet")):
            if path.stem.startswith("table_"):
                table_names.append(path.stem[len("table_"):])
            else:
                obs_names.append(path.stem)

    return {"obs": obs_names, "tables": table_names}


# ── manifest ─────────────────────────────────────────────────────────

def update_manifest(
    module_name: str,
    output_dir: str | Path,
    saved_files: list[str],
) -> None:
    """Append an entry to ``manifest.json`` for bookkeeping."""
    d = spatch_output_dir(output_dir)
    manifest_path = d / MANIFEST_FILE

    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {"modules": {}}

    manifest["modules"][module_name] = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "files": saved_files,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))
