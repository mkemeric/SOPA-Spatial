"""
QC Preprocessing Module

Filters low-quality cells/genes, normalizes, and computes the PCA/UMAP
embedding used by downstream clustering and spatial-stats modules.

Promotes steps 4-5 ("Preprocessing and QC" / "Dimensionality Reduction")
of notebooks/01_spatch_workflow_example.ipynb into a YAML-configured,
reusable module, so QC thresholds and preprocessing parameters live in
the step's ``config:`` block instead of hardcoded notebook cells.
"""

from __future__ import annotations

import scanpy as sc
import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class QCPreprocessing(SpatchModule):
    """Filter cells/genes by QC thresholds, normalize, and compute PCA/UMAP."""

    name = "qc_preprocessing"
    version = "1.0.0"
    description = "QC filtering, normalization, HVG/PCA/UMAP"
    category = "qc"
    requires = []
    produces = []

    def validate_inputs(self, sdata: sd.SpatialData) -> list[str]:
        table_key = self.config.get("table_key", "table")
        if table_key not in sdata.tables:
            return [f"Missing required table: '{table_key}'"]
        return []

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        min_genes_per_cell: int = 50,
        max_genes_per_cell: int = 10000,
        min_cells_per_gene: int = 5,
        max_pct_mt: float = 25.0,
        mito_prefix: str = "MT-",
        target_sum: float = 10000,
        log1p: bool = True,
        n_top_genes: int = 3000,
        flavor: str = "seurat_v3",
        n_pcs: int = 50,
        n_neighbors: int = 30,
        umap_min_dist: float = 0.3,
        **kwargs,
    ) -> ModuleResult:
        """Run QC filtering + normalization + PCA/UMAP.

        Args:
            sdata: SpatialData object with an aggregated cell x gene table.
            table_key: Key for the cell table in sdata.tables.
            min_genes_per_cell / max_genes_per_cell: Cell filtering bounds.
            min_cells_per_gene: Gene filtering bound.
            max_pct_mt: Max allowed mitochondrial read percentage.
            mito_prefix: Gene name prefix used to flag mitochondrial genes.
            target_sum: Total-count normalization target.
            log1p: Whether to log1p-transform after normalization.
            n_top_genes / flavor: Highly-variable-gene selection.
            n_pcs: Number of PCA components.
            n_neighbors: Neighbors for the KNN graph (feeds Leiden + UMAP).
            umap_min_dist: UMAP `min_dist` parameter.

        Returns:
            ModuleResult with the filtered/normalized table (PCA, neighbor
            graph, and UMAP embedding computed) written back to sdata.
        """
        log = []
        adata = sdata.tables[table_key]

        if hasattr(adata.X, "compute"):
            adata.X = adata.X.compute()
            log.append("Materialized Dask-backed expression matrix")

        cells_before, genes_before = adata.n_obs, adata.n_vars

        adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
        sc.pp.filter_cells(adata, max_genes=max_genes_per_cell)
        adata = adata[adata.obs.pct_counts_mt < max_pct_mt, :].copy()
        sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
        log.append(
            f"Filtered {cells_before}->{adata.n_obs} cells, "
            f"{genes_before}->{adata.n_vars} genes "
            f"(min_genes={min_genes_per_cell}, max_genes={max_genes_per_cell}, "
            f"max_pct_mt={max_pct_mt}, min_cells_per_gene={min_cells_per_gene})"
        )

        sc.pp.normalize_total(adata, target_sum=target_sum)
        if log1p:
            sc.pp.log1p(adata)
        adata.raw = adata.copy()
        log.append("Normalized, log-transformed, and stashed raw counts")

        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)
        sc.pp.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1))
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata, min_dist=umap_min_dist)
        log.append(
            f"Computed {int(adata.var.highly_variable.sum())} HVGs, "
            f"PCA({adata.obsm['X_pca'].shape[1]}), neighbors(k={n_neighbors}), UMAP"
        )

        sdata.tables[table_key] = adata

        return ModuleResult(
            sdata=sdata,
            metrics={
                "cells_before": cells_before,
                "cells_after": adata.n_obs,
                "genes_before": genes_before,
                "genes_after": adata.n_vars,
                "n_hvg": int(adata.var.highly_variable.sum()),
                "n_pcs": int(adata.obsm["X_pca"].shape[1]),
            },
            log=log,
        )
