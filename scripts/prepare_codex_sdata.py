#!/usr/bin/env python3
"""Prepare a SpatialData zarr from the SPATCH CODEX dataset.

Reads:
  - transcriptome/adata.h5ad   → tables["table"]
  - proteome/adata_codex.h5ad  → tables["codex_table"]
  - segmentation_mask/cell_boundaries.csv → shapes["cell_boundaries"]

Writes:
  - {output_path}  (SpatialData zarr directory)

Usage:
    python scripts/prepare_codex_sdata.py /mnt/shared/data/codex results/codex.zarr
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import spatialdata as sd
from spatialdata.models import TableModel, ShapesModel
from spatialdata.transformations import Identity


def load_cell_boundaries(csv_path: str | Path, max_cells: int | None = None) -> "GeoDataFrame":
    """Convert vertex CSV to a GeoDataFrame of Shapely polygons.

    The CSV has columns: cell_id, vertex_x, vertex_y, label_id.
    Each cell_id has multiple rows (polygon vertices).
    """
    from geopandas import GeoDataFrame
    from shapely.geometry import Polygon

    print(f"  Reading {csv_path} ...")
    df = pd.read_csv(csv_path)
    print(f"  {len(df):,} vertex rows")

    # Group vertices by cell_id and build polygons
    grouped = df.groupby("cell_id")
    if max_cells:
        grouped = list(grouped)[:max_cells]
        grouped = {k: v for k, v in grouped}
    else:
        grouped = {k: v for k, v in grouped}

    cell_ids = []
    polygons = []
    skipped = 0

    for cell_id, verts in grouped.items():
        coords = verts[["vertex_x", "vertex_y"]].values
        if len(coords) < 3:
            skipped += 1
            continue
        try:
            poly = Polygon(coords)
            if poly.is_valid and not poly.is_empty:
                cell_ids.append(cell_id)
                polygons.append(poly)
            else:
                skipped += 1
        except Exception:
            skipped += 1

    print(f"  Built {len(polygons):,} polygons ({skipped} skipped)")

    gdf = GeoDataFrame({"geometry": polygons}, index=cell_ids)
    return gdf


def main(data_dir: str, output_path: str):
    data_dir = Path(data_dir)
    output_path = Path(output_path)

    if output_path.exists():
        print(f"Output {output_path} already exists — delete it first to recreate.")
        sys.exit(1)

    # ── 1. Load transcriptome table ──────────────────────────────────
    print("Loading transcriptome h5ad ...")
    st = sc.read_h5ad(data_dir / "transcriptome" / "adata.h5ad")
    print(f"  {st.n_obs:,} cells × {st.n_vars:,} genes")

    # Map high_quality → in_tissue for diffusion analysis
    if "high_quality" in st.obs.columns:
        st.obs["in_tissue"] = st.obs["high_quality"].astype(int)
        print(f"  Mapped high_quality → in_tissue ({st.obs['in_tissue'].sum():,} cells)")

    # Ensure required fields for TableModel
    st.obs["cell_id"] = st.obs_names
    st.obs["region"] = "cell_boundaries"

    table = TableModel.parse(
        st,
        region="cell_boundaries",
        region_key="region",
        instance_key="cell_id",
    )

    # ── 2. Load CODEX proteome table ─────────────────────────────────
    print("Loading CODEX proteome h5ad ...")
    codex = sc.read_h5ad(data_dir / "proteome" / "adata_codex.h5ad")
    print(f"  {codex.n_obs:,} cells × {codex.n_vars:,} proteins")
    print(f"  Proteins: {list(codex.var_names)}")

    codex.obs["cell_id"] = codex.obs_names
    codex.obs["region"] = "codex_cells"

    codex_table = TableModel.parse(
        codex,
        region="codex_cells",
        region_key="region",
        instance_key="cell_id",
    )

    # ── 3. Load cell boundaries ──────────────────────────────────────
    print("Loading cell boundaries ...")
    boundaries_path = data_dir / "segmentation_mask" / "cell_boundaries.csv"
    gdf = load_cell_boundaries(boundaries_path)
    shapes = ShapesModel.parse(
        gdf, transformations={"global": Identity()}
    )

    # ── 4. Build SpatialData object ──────────────────────────────────
    print("Building SpatialData object ...")
    sdata = sd.SpatialData(
        tables={"table": table, "codex_table": codex_table},
        shapes={"cell_boundaries": shapes},
    )

    # ── 5. Write zarr ────────────────────────────────────────────────
    print(f"Writing zarr → {output_path} ...")
    sdata.write(str(output_path))
    print("Done!")

    # Summary
    print(f"\n  Tables:  {list(sdata.tables.keys())}")
    print(f"  Shapes:  {list(sdata.shapes.keys())}")
    print(f"  Images:  {list(sdata.images.keys())}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <data_dir> <output_zarr>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
