"""
Grid Binning Module

Bin spatial data to a uniform resolution grid for cross-platform comparison.
Based on the original SPATCH 2_8um_bin.py analysis.

Supports two input modes:
  1. AnnData table (cell/spot-level expression already aggregated)
  2. Transcript-level DataFrame with x/y coordinates and gene names
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

import spatialdata as sd
from spatialdata.models import TableModel

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class GridBinning(SpatchModule):
    """Bin spatial expression data to a uniform resolution grid.

    This module re-grids spatial transcriptomics data from any platform
    to a common resolution (default 8 µm), enabling apples-to-apples
    cross-platform comparisons.

    Two modes:
      - **table mode**: rebins an existing cell/spot-level AnnData table.
      - **transcript mode**: aggregates raw transcript coordinates directly
        (for Xenium, CosMx, etc.).

    The output is a new AnnData table stored in the SpatialData object.
    """

    name = "grid_binning"
    version = "1.0.0"
    description = "Bin spatial data to a uniform resolution grid"
    category = "analysis"
    requires = ["tables"]
    produces = []  # dynamic key, e.g. tables/binned_8um

    def run(
        self,
        sdata: sd.SpatialData,
        target_resolution_um: float = 8.0,
        mode: str = "table",
        table_key: str = "table",
        coord_key: str = "spatial",
        source_resolution_um: float = 1.0,
        transcript_df: pd.DataFrame | None = None,
        x_col: str = "x_location",
        y_col: str = "y_location",
        gene_col: str = "feature_name",
        output_key: str | None = None,
        **kwargs,
    ) -> ModuleResult:
        """Run grid binning.

        Args:
            sdata: SpatialData object.
            target_resolution_um: Target grid size in µm.
            mode: "table" to rebin an existing AnnData, or "transcript"
                  to aggregate from a transcript DataFrame.
            table_key: Key for the source AnnData in sdata.tables (table mode).
            coord_key: obsm key holding spatial coordinates (table mode).
            source_resolution_um: Original pixel-to-µm factor for the
                coordinates stored in obsm (table mode). E.g. 1.0 if
                coords are already in µm; 0.2125 for raw Xenium pixels.
            transcript_df: DataFrame with per-transcript rows (transcript mode).
                If None and mode=="transcript", looks in sdata.points.
            x_col: Column name for x coordinate (transcript mode).
            y_col: Column name for y coordinate (transcript mode).
            gene_col: Column name for gene/feature name (transcript mode).
            output_key: Key for the binned table in sdata.tables.
                Defaults to ``binned_{target_resolution_um}um``.
            **kwargs: Extra config (ignored).

        Returns:
            ModuleResult with the binned table added to sdata.
        """
        log = []
        output_key = output_key or f"binned_{int(target_resolution_um)}um"

        if mode == "transcript":
            binned = self._bin_transcripts(
                sdata,
                transcript_df,
                x_col,
                y_col,
                gene_col,
                target_resolution_um,
                source_resolution_um,
            )
            log.append(
                f"Binned transcripts to {target_resolution_um} µm grid"
            )
        elif mode == "table":
            if table_key not in sdata.tables:
                raise ValueError(f"Table '{table_key}' not found in sdata")
            adata = sdata.tables[table_key]
            binned = self._bin_table(
                adata,
                coord_key,
                source_resolution_um,
                target_resolution_um,
            )
            log.append(
                f"Rebinned table '{table_key}' from "
                f"{source_resolution_um} µm to {target_resolution_um} µm"
            )
        else:
            raise ValueError(f"Unknown mode '{mode}'. Use 'table' or 'transcript'.")

        log.append(
            f"Result: {binned.n_obs} bins × {binned.n_vars} genes"
        )

        sdata.tables[output_key] = TableModel.parse(binned)

        return ModuleResult(
            sdata=sdata,
            metrics={
                "n_bins": binned.n_obs,
                "n_genes": binned.n_vars,
                "target_resolution_um": target_resolution_um,
            },
            log=log,
        )

    # ── transcript mode ─────────────────────────────────────────────

    def _bin_transcripts(
        self,
        sdata: sd.SpatialData,
        transcript_df: pd.DataFrame | None,
        x_col: str,
        y_col: str,
        gene_col: str,
        target_res: float,
        source_res: float,
    ) -> ad.AnnData:
        """Aggregate raw transcripts into grid bins.

        Mirrors SPATCH ``xenium_8um_bin`` / ``nano_8um_bin``.
        """
        if transcript_df is None:
            # Try to pull from sdata.points
            point_keys = list(sdata.points.keys()) if sdata.points else []
            if not point_keys:
                raise ValueError(
                    "No transcript_df provided and no points found in sdata"
                )
            transcript_df = sdata.points[point_keys[0]].compute()

        df = transcript_df.copy()

        # Compute bin indices
        df["binx"] = (df[x_col] * source_res // target_res).astype(int)
        df["biny"] = (df[y_col] * source_res // target_res).astype(int)

        # Pivot to count matrix
        pivot = (
            df.groupby(["binx", "biny", gene_col])
            .size()
            .unstack(fill_value=0)
        )

        spatial = np.array(list(pivot.index))  # (N, 2)
        gene_names = pivot.columns.tolist()

        adata = ad.AnnData(
            X=csr_matrix(pivot.values),
            var=pd.DataFrame(index=gene_names),
        )
        adata.obsm["spatial"] = spatial * target_res  # back to µm

        return adata

    # ── table mode ──────────────────────────────────────────────────

    def _bin_table(
        self,
        adata: ad.AnnData,
        coord_key: str,
        source_res: float,
        target_res: float,
    ) -> ad.AnnData:
        """Rebin an existing AnnData to a coarser grid.

        Mirrors SPATCH ``bin_8um``.
        """
        coords = adata.obsm[coord_key]

        # Scale coordinates to µm, then quantise to target grid
        grid_indices = (coords * source_res // target_res).astype(int)

        unique_indices, inverse = np.unique(
            grid_indices, axis=0, return_inverse=True
        )

        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

        n_bins = len(unique_indices)
        binned_X = np.zeros((n_bins, X.shape[1]), dtype=np.float32)
        for i in range(n_bins):
            binned_X[i] = X[inverse == i].sum(axis=0)

        binned = ad.AnnData(
            X=csr_matrix(binned_X),
            var=adata.var.copy(),
        )
        binned.obsm["spatial"] = unique_indices.astype(float) * target_res

        return binned
