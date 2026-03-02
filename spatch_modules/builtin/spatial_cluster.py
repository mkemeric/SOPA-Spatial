"""
Spatial Clustering Module

CellCharter-based spatial clustering with latent embedding and
cross-platform cluster proportion correlation.
Based on the original SPATCH 10_spatial_cluster.py analysis.

Requires:
  - scvi-tools     (for scVI / TRVAE latent embeddings)
  - cellcharter    (for spatial neighbor aggregation + GMM clustering)
  - squidpy        (for Delaunay spatial graph)

Install with: ``pip install spatch-modules[spatial_clustering]``
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr

import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


def _ensure_deps():
    """Lazy import of heavy optional deps."""
    deps = {}
    try:
        import squidpy as sq
        deps["sq"] = sq
    except ImportError:
        raise ImportError("squidpy is required for spatial clustering.")
    try:
        import cellcharter as cc
        deps["cc"] = cc
    except ImportError:
        raise ImportError(
            "cellcharter is required. "
            "Install with: pip install spatch-modules[spatial_clustering]"
        )
    return deps


@register
class SpatialCluster(SpatchModule):
    """CellCharter spatial clustering with optional cross-platform correlation.

    Workflow:
      1. Compute latent embedding via scVI (transcriptomics) or TRVAE (protein).
      2. Build Delaunay spatial neighbor graph (Squidpy).
      3. Aggregate neighbor embeddings over *n* layers (CellCharter).
      4. GMM clustering with automatic K selection (or user-specified K).
      5. Optionally compute cluster proportion correlation with a
         reference dataset (e.g. ST vs CODEX).
    """

    name = "spatial_cluster"
    version = "1.0.0"
    description = "CellCharter spatial clustering + cross-platform correlation"
    category = "analysis"
    requires = ["tables/table"]
    produces = []  # adds obs['spatial_cluster']

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        coord_key: str = "spatial",
        data_type: str = "transcriptome",
        n_layers: int = 5,
        n_clusters: tuple[int, int] | int | None = None,
        max_runs: int = 10,
        use_gpu: bool = False,
        cluster_col: str = "spatial_cluster",
        reference_table_key: str | None = None,
        reference_cluster_col: str | None = None,
        correlation_grid_size: float = 100.0,
        random_state: int = 12345,
        **kwargs,
    ) -> ModuleResult:
        """Run spatial clustering.

        Args:
            sdata: SpatialData object.
            table_key: Key for the cell table.
            coord_key: obsm key for spatial coordinates.
            data_type: "transcriptome" (uses scVI) or "protein" (uses TRVAE).
            n_layers: Number of neighbor aggregation layers.
            n_clusters: Either a single int, a (min, max) tuple for auto-K,
                or None for automatic range (5, 10).
            max_runs: Max runs for ClusterAutoK.
            use_gpu: Whether to use GPU for embedding + clustering.
            cluster_col: obs column name for cluster labels.
            reference_table_key: If set, compute cluster proportion
                correlation against this table.
            reference_cluster_col: Cluster column in the reference table.
            correlation_grid_size: Grid cell size (µm) for proportion
                correlation.
            random_state: Seed for reproducibility.
            **kwargs: Extra config (ignored).

        Returns:
            ModuleResult with spatial_cluster column added.
        """
        deps = _ensure_deps()
        sq = deps["sq"]
        cc = deps["cc"]
        log = []

        if table_key not in sdata.tables:
            raise ValueError(f"Table '{table_key}' not found")

        adata = sdata.tables[table_key]

        # ── 1. Latent embedding ──────────────────────────────────────
        if data_type == "transcriptome":
            adata = self._embed_scvi(adata, use_gpu, random_state)
            use_rep = "X_scVI"
            log.append("Computed scVI latent embedding")
        elif data_type == "protein":
            adata = self._embed_trvae(adata, use_gpu, random_state)
            use_rep = "X_trVAE"
            log.append("Computed TRVAE latent embedding")
        else:
            raise ValueError(f"Unknown data_type '{data_type}'")

        # ── 2. Spatial neighbors ─────────────────────────────────────
        sq.gr.spatial_neighbors(
            adata, coord_type="generic", delaunay=True, spatial_key=coord_key
        )
        cc.gr.remove_long_links(adata)
        log.append("Built Delaunay spatial graph")

        # ── 3. Neighbor aggregation ──────────────────────────────────
        cc.gr.aggregate_neighbors(
            adata, n_layers=n_layers, use_rep=use_rep, out_key="X_cellcharter"
        )
        log.append(f"Aggregated neighbors ({n_layers} layers)")

        # ── 4. Clustering ────────────────────────────────────────────
        trainer_params = (
            dict(accelerator="gpu", devices=1) if use_gpu else {}
        )

        if n_clusters is None or isinstance(n_clusters, tuple):
            k_range = n_clusters or (5, 10)
            autok = cc.tl.ClusterAutoK(
                n_clusters=k_range,
                max_runs=max_runs,
                model_params=dict(
                    random_state=random_state,
                    trainer_params=trainer_params,
                ),
            )
            autok.fit(adata, use_rep="X_cellcharter")
            best_k = autok.best_k
            log.append(f"Auto-K selected K={best_k}")
        else:
            best_k = int(n_clusters)
            log.append(f"Using user-specified K={best_k}")

        gmm = cc.tl.Cluster(
            n_clusters=best_k,
            random_state=random_state,
            trainer_params=trainer_params,
        )
        gmm.fit(adata, use_rep="X_cellcharter")
        adata.obs[cluster_col] = gmm.predict(adata, use_rep="X_cellcharter")
        log.append(f"Clustered into {best_k} spatial clusters")

        # Write back
        sdata.tables[table_key] = adata

        # ── 5. Cross-platform correlation ────────────────────────────
        corr_metrics = {}
        if (
            reference_table_key
            and reference_table_key in sdata.tables
            and reference_cluster_col
        ):
            corr_matrix = self._cluster_proportion_correlation(
                adata,
                sdata.tables[reference_table_key],
                cluster_col,
                reference_cluster_col,
                coord_key,
                correlation_grid_size,
            )
            corr_metrics["mean_max_corr"] = float(
                np.nanmax(corr_matrix, axis=1).mean()
            )
            log.append(
                f"Cross-platform correlation: mean max r = "
                f"{corr_metrics['mean_max_corr']:.3f}"
            )

        return ModuleResult(
            sdata=sdata,
            metrics={
                "n_clusters": best_k,
                "n_cells": adata.n_obs,
                **corr_metrics,
            },
            log=log,
        )

    # ── embedding helpers ────────────────────────────────────────────

    @staticmethod
    def _embed_scvi(
        adata: ad.AnnData, use_gpu: bool, seed: int
    ) -> ad.AnnData:
        """Compute scVI latent embedding for transcriptomics data."""
        import scvi

        adata = adata.copy()
        adata.obs["sample"] = "transcriptome"
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)

        scvi.settings.seed = seed
        scvi.model.SCVI.setup_anndata(
            adata, layer="counts", batch_key="sample"
        )
        model = scvi.model.SCVI(adata)
        model.train(early_stopping=True, enable_progress_bar=True)
        adata.obsm["X_scVI"] = model.get_latent_representation(adata).astype(
            np.float32
        )
        return adata

    @staticmethod
    def _embed_trvae(
        adata: ad.AnnData, use_gpu: bool, seed: int
    ) -> ad.AnnData:
        """Compute TRVAE latent embedding for protein data."""
        try:
            import scarches as sca
        except ImportError:
            raise ImportError(
                "scarches is required for TRVAE embedding. "
                "Install with: pip install scarches"
            )

        adata = adata.copy()
        adata.obs["sample"] = "codex"
        adata.X = adata.X.astype(np.float32)

        trvae = sca.models.TRVAE(
            adata=adata,
            condition_key="sample",
            conditions=["codex"],
            hidden_layer_sizes=[128, 128],
        )
        trvae.train(
            n_epochs=500,
            alpha_epoch_anneal=200,
            early_stopping_kwargs={
                "early_stopping_metric": "val_unweighted_loss",
                "patience": 20,
                "threshold": 0,
                "reduce_lr": True,
            },
        )
        adata.obsm["X_trVAE"] = trvae.get_latent().astype(np.float32)
        return adata

    # ── cluster proportion correlation ───────────────────────────────

    @staticmethod
    def _cluster_proportion_correlation(
        adata1: ad.AnnData,
        adata2: ad.AnnData,
        cluster_col1: str,
        cluster_col2: str,
        coord_key: str,
        grid_size: float,
    ) -> np.ndarray:
        """Compute Pearson correlation of cluster proportions on a common grid.

        Direct port of SPATCH ``calculate_cluster_proportion`` +
        ``create_proportion_vectors``.
        """

        def _grid_proportions(adata, cluster_col, coord_key, grid_size):
            coords = adata.obsm[coord_key]
            grid_idx = np.floor(coords / grid_size).astype(int)
            labels = adata.obs[cluster_col].values

            # Build {(gx, gy): {label: count}}
            counts = {}
            for i in range(len(labels)):
                key = (grid_idx[i, 0], grid_idx[i, 1])
                if key not in counts:
                    counts[key] = {}
                lbl = labels[i]
                counts[key][lbl] = counts[key].get(lbl, 0) + 1

            # Normalise to proportions
            proportions = {}
            for key, lbl_counts in counts.items():
                total = sum(lbl_counts.values())
                proportions[key] = {
                    lbl: c / total for lbl, c in lbl_counts.items()
                }
            return proportions

        props1 = _grid_proportions(adata1, cluster_col1, coord_key, grid_size)
        props2 = _grid_proportions(adata2, cluster_col2, coord_key, grid_size)

        common_grids = set(props1.keys()) & set(props2.keys())
        if not common_grids:
            return np.array([[]])

        labels1 = sorted(
            {lbl for p in props1.values() for lbl in p}
        )
        labels2 = sorted(
            {lbl for p in props2.values() for lbl in p}
        )

        def _vectors(proportions, common, labels):
            vecs = {lbl: [] for lbl in labels}
            for grid in common:
                for lbl in labels:
                    vecs[lbl].append(proportions.get(grid, {}).get(lbl, 0))
            return vecs

        vecs1 = _vectors(props1, common_grids, labels1)
        vecs2 = _vectors(props2, common_grids, labels2)

        n1, n2 = len(labels1), len(labels2)
        corr_matrix = np.zeros((n1, n2))
        for i, l1 in enumerate(labels1):
            for j, l2 in enumerate(labels2):
                v1, v2 = vecs1[l1], vecs2[l2]
                if len(v1) > 1 and np.std(v1) > 0 and np.std(v2) > 0:
                    corr_matrix[i, j], _ = pearsonr(v1, v2)

        return corr_matrix
