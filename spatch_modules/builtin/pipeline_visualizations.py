"""
Pipeline Visualizations Module

Generate publication-quality figures from the processed SpatialData object.
Runs lightweight scanpy preprocessing (normalize, PCA, UMAP, Leiden) if not
already present, then produces the plot types listed in the YAML config.

Plot types:
  - spatial_scatter:          Cells on tissue coordinates, colored by cluster
  - umap:                     UMAP embedding colored by cluster
  - cluster_composition:      Bar chart of cluster proportions
  - marker_heatmap:           Mean expression of marker genes per cluster
  - neighborhood_enrichment:  Spatial co-localization matrix (squidpy)
  - cell_shape_distributions: Violin plots of morphology metrics
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless backend — no display needed
import matplotlib.pyplot as plt

import scanpy as sc
import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


# ── helpers ──────────────────────────────────────────────────────────

def _ensure_preprocessed(adata, n_top_genes: int = 3000, n_pcs: int = 50,
                         n_neighbors: int = 30, resolution: float = 0.5):
    """Run scanpy preprocessing if results aren't already present."""
    log = []

    # Stash raw counts before normalizing
    if adata.raw is None:
        adata.raw = adata.copy()
        log.append("Saved raw counts to adata.raw")

    # Normalize + log1p
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        log.append("Normalized and log-transformed")

    # HVG
    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes,
                                    flavor="seurat_v3")
        log.append(f"Selected {n_top_genes} highly variable genes")

    # PCA
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1))
        log.append(f"Computed PCA ({adata.obsm['X_pca'].shape[1]} components)")

    # Neighbors
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        log.append(f"Built neighbor graph (k={n_neighbors})")

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)
        log.append("Computed UMAP embedding")

    # Leiden clustering
    if "leiden" not in adata.obs.columns:
        sc.tl.leiden(adata, resolution=resolution, flavor="igraph",
                     n_iterations=2)
        log.append(f"Leiden clustering (resolution={resolution})")

    return log


def _savefig(fig, path: Path, dpi: int = 300):
    """Save a matplotlib figure and close it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ── module ───────────────────────────────────────────────────────────

@register
class PipelineVisualizations(SpatchModule):
    """Generate figures from the processed SpatialData object."""

    name = "pipeline_visualizations"
    version = "1.0.0"
    description = "Generate pipeline QC and analysis figures"
    category = "analysis"
    requires = ["tables/table"]
    produces = []

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        output_dir: str = "./results/figures",
        figure_format: str = "png",
        figure_dpi: int = 300,
        marker_genes: list[str] | None = None,
        plots: list[str] | None = None,
        leiden_resolution: float = 0.5,
        n_top_genes: int = 3000,
        **kwargs,
    ) -> ModuleResult:
        log = []
        artifacts: dict[str, str] = {}
        out = Path(output_dir)
        fmt = figure_format
        dpi = figure_dpi

        adata = sdata.tables[table_key].copy()

        # Materialize Dask-backed arrays — scanpy requires eager data
        if hasattr(adata.X, 'compute'):
            adata.X = adata.X.compute()
            log.append("Materialized Dask-backed expression matrix")

        # ── preprocessing ────────────────────────────────────────────
        prep_log = _ensure_preprocessed(
            adata,
            n_top_genes=n_top_genes,
            resolution=leiden_resolution,
        )
        log.extend(prep_log)

        all_plots = plots or [
            "spatial_scatter", "umap", "cluster_composition",
            "marker_heatmap", "cell_shape_distributions",
        ]

        # ── spatial scatter ──────────────────────────────────────────
        if "spatial_scatter" in all_plots and "spatial" in adata.obsm:
            path = out / f"spatial_scatter.{fmt}"
            fig, ax = plt.subplots(figsize=(10, 10))
            coords = adata.obsm["spatial"]
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=adata.obs["leiden"].astype("category").cat.codes,
                cmap="tab20", s=0.3, alpha=0.7, rasterized=True,
            )
            ax.set_aspect("equal")
            ax.set_title("Spatial — Leiden clusters")
            ax.set_xlabel("X (pixels)")
            ax.set_ylabel("Y (pixels)")
            ax.invert_yaxis()
            _savefig(fig, path, dpi)
            artifacts["spatial_scatter"] = str(path)
            log.append(f"Saved spatial scatter → {path}")

        # ── UMAP ────────────────────────────────────────────────────
        if "umap" in all_plots and "X_umap" in adata.obsm:
            path = out / f"umap_leiden.{fmt}"
            fig, ax = plt.subplots(figsize=(8, 8))
            sc.pl.umap(adata, color="leiden", ax=ax, show=False,
                       frameon=False, title="UMAP — Leiden clusters")
            _savefig(fig, path, dpi)
            artifacts["umap"] = str(path)
            log.append(f"Saved UMAP → {path}")

        # ── cluster composition ──────────────────────────────────────
        if "cluster_composition" in all_plots and "leiden" in adata.obs.columns:
            path = out / f"cluster_composition.{fmt}"
            counts = adata.obs["leiden"].value_counts().sort_index()
            fig, ax = plt.subplots(figsize=(10, 5))
            counts.plot.bar(ax=ax, color=plt.cm.tab20.colors[:len(counts)])
            ax.set_title("Cells per Leiden cluster")
            ax.set_xlabel("Cluster")
            ax.set_ylabel("Cell count")
            _savefig(fig, path, dpi)
            artifacts["cluster_composition"] = str(path)
            log.append(f"Saved cluster composition → {path}")

        # ── marker heatmap ───────────────────────────────────────────
        if "marker_heatmap" in all_plots and marker_genes:
            available = [g for g in marker_genes if g in adata.var_names]
            if available:
                path = out / f"marker_heatmap.{fmt}"
                sc.pl.heatmap(
                    adata, var_names=available, groupby="leiden",
                    show=False, swap_axes=True, dendrogram=False,
                    figsize=(max(8, len(available) * 0.6), 6),
                )
                _savefig(plt.gcf(), path, dpi)
                artifacts["marker_heatmap"] = str(path)
                log.append(f"Saved marker heatmap ({len(available)} genes) → {path}")
            else:
                log.append("Skipped marker heatmap — no marker genes found in var_names")

        # ── neighborhood enrichment ──────────────────────────────────
        if "neighborhood_enrichment" in all_plots and "spatial" in adata.obsm:
            try:
                import squidpy as sq
                sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=30)
                sq.gr.nhood_enrichment(adata, cluster_key="leiden")
                path = out / f"neighborhood_enrichment.{fmt}"
                fig, ax = plt.subplots(figsize=(8, 7))
                sq.pl.nhood_enrichment(
                    adata, cluster_key="leiden", ax=ax, show=False,
                    title="Neighborhood enrichment (Leiden)",
                )
                _savefig(fig, path, dpi)
                artifacts["neighborhood_enrichment"] = str(path)
                log.append(f"Saved neighborhood enrichment → {path}")
            except Exception as e:
                log.append(f"Skipped neighborhood enrichment — {e}")

        # ── cell shape distributions ─────────────────────────────────
        shape_cols = [c for c in ["area_um2", "circularity", "eccentricity",
                                  "solidity", "aspect_ratio"]
                      if c in adata.obs.columns]
        if "cell_shape_distributions" in all_plots and shape_cols:
            path = out / f"cell_shape_distributions.{fmt}"
            n = len(shape_cols)
            fig, axes = plt.subplots(1, n, figsize=(4 * n, 5))
            if n == 1:
                axes = [axes]
            for ax, col in zip(axes, shape_cols):
                sc.pl.violin(adata, keys=col, groupby="leiden",
                             ax=ax, show=False, rotation=90)
                ax.set_title(col.replace("_", " ").title())
            fig.suptitle("Cell shape metrics by cluster", y=1.02, fontsize=14)
            fig.tight_layout()
            _savefig(fig, path, dpi)
            artifacts["cell_shape_distributions"] = str(path)
            log.append(f"Saved cell shape distributions ({len(shape_cols)} metrics) → {path}")

        return ModuleResult(
            sdata=sdata,
            metrics={
                "n_plots_generated": len(artifacts),
                "n_clusters": int(adata.obs["leiden"].nunique())
                    if "leiden" in adata.obs.columns else 0,
            },
            artifacts=artifacts,
            log=log,
        )
