"""
Diffusion Analysis Module

Compare in-tissue vs out-of-tissue signal to quantify transcript diffusion.
Based on the original SPATCH 4_diffusion.py analysis.
"""

import numpy as np
import pandas as pd
import anndata as ad
from sklearn.neighbors import NearestNeighbors

import spatialdata as sd
from spatialdata.models import TableModel

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class DiffusionAnalysis(SpatchModule):
    """Quantify transcript diffusion from in-tissue to out-of-tissue regions.
    
    This module compares transcript/expression counts between tissue regions
    and surrounding non-tissue areas to estimate signal diffusion effects.
    
    The analysis:
    1. Identifies spots/cells inside vs outside tissue boundaries
    2. Calculates minimum distances from outside spots to tissue edge
    3. Computes per-gene diffusion ratios (outside / total counts)
    4. Stores results as a new table in the SpatialData object
    """
    
    name = "diffusion_analysis"
    version = "1.0.0"
    description = "Compare in-tissue vs out-of-tissue signal to quantify transcript diffusion"
    category = "analysis"
    requires = ["tables/table"]
    produces = ["tables/diffusion_metrics"]

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        tissue_mask_key: str = "tissue_boundary",
        in_tissue_col: str = "in_tissue",
        buffer_distance_um: float = 50.0,
        compute_distances: bool = True,
        batch_size: int = 100000,
        **kwargs
    ) -> ModuleResult:
        """Run diffusion analysis.
        
        Args:
            sdata: SpatialData object with expression table.
            table_key: Key for the AnnData table in sdata.tables.
            tissue_mask_key: Key for tissue boundary shapes (optional).
            in_tissue_col: Column in obs indicating tissue membership.
            buffer_distance_um: Distance buffer for analyzing outside tissue.
            compute_distances: Whether to compute min distances to tissue.
            batch_size: Batch size for distance calculations.
            **kwargs: Additional config overrides.
        
        Returns:
            ModuleResult with diffusion metrics added to sdata.
        """
        log = []
        
        # Get the expression table
        adata = sdata.tables[table_key]
        
        # Determine in-tissue status
        if in_tissue_col in adata.obs.columns:
            in_tissue_mask = adata.obs[in_tissue_col].astype(bool)
        elif tissue_mask_key in sdata.shapes:
            # Use spatial query to determine tissue membership
            in_tissue_mask = self._query_tissue_membership(
                sdata, adata, tissue_mask_key
            )
            adata.obs[in_tissue_col] = in_tissue_mask
        else:
            raise ValueError(
                f"No tissue mask found. Provide either '{in_tissue_col}' column "
                f"in obs or '{tissue_mask_key}' in shapes."
            )
        
        # Split data by tissue status
        adata_in = adata[in_tissue_mask, :].copy()
        adata_out = adata[~in_tissue_mask, :].copy()
        
        log.append(f"In-tissue spots: {adata_in.n_obs}")
        log.append(f"Out-of-tissue spots: {adata_out.n_obs}")
        
        # Global diffusion metrics
        total_in = np.sum(adata_in.X) if hasattr(adata_in.X, 'sum') else np.sum(adata_in.X.toarray())
        total_out = np.sum(adata_out.X) if hasattr(adata_out.X, 'sum') else np.sum(adata_out.X.toarray())
        total = total_in + total_out
        global_diffusion_ratio = total_out / total if total > 0 else 0.0
        
        # Per-gene diffusion metrics
        gene_metrics = self._compute_gene_diffusion(adata_in, adata_out)
        log.append(f"Computed diffusion for {len(gene_metrics)} genes")
        
        # Compute distances from outside spots to nearest in-tissue spot
        if compute_distances and adata_out.n_obs > 0 and adata_in.n_obs > 0:
            distances = self._calculate_min_distances(
                adata_in, adata_out, batch_size
            )
            adata_out.obs["distance_to_tissue"] = distances
            gene_metrics = self._add_distance_correlation(
                gene_metrics, adata_out, distances
            )
            log.append(f"Computed distances to tissue boundary")
        
        # Create diffusion metrics table
        metrics_adata = ad.AnnData(
            obs=gene_metrics.set_index("gene")
        )
        sdata.tables["diffusion_metrics"] = TableModel.parse(metrics_adata)
        
        return ModuleResult(
            sdata=sdata,
            metrics={
                "global_diffusion_ratio": round(global_diffusion_ratio, 4),
                "in_tissue_counts": int(total_in),
                "out_tissue_counts": int(total_out),
                "n_genes_analyzed": len(gene_metrics),
                "n_in_tissue_spots": int(adata_in.n_obs),
                "n_out_tissue_spots": int(adata_out.n_obs)
            },
            log=log
        )

    def _query_tissue_membership(
        self,
        sdata: sd.SpatialData,
        adata: ad.AnnData,
        tissue_mask_key: str
    ) -> np.ndarray:
        """Determine which cells are inside tissue using spatial query."""
        from shapely.geometry import Point
        from shapely.ops import unary_union
        
        tissue_shapes = sdata.shapes[tissue_mask_key]
        tissue_polygon = unary_union(tissue_shapes.geometry)
        
        coords = adata.obsm.get("spatial", None)
        if coords is None:
            raise ValueError("No spatial coordinates found in adata.obsm['spatial']")
        
        in_tissue = np.array([
            tissue_polygon.contains(Point(x, y)) for x, y in coords
        ])
        
        return in_tissue

    def _compute_gene_diffusion(
        self,
        adata_in: ad.AnnData,
        adata_out: ad.AnnData
    ) -> pd.DataFrame:
        """Compute per-gene diffusion metrics."""
        # Get expression matrices
        X_in = adata_in.X.toarray() if hasattr(adata_in.X, 'toarray') else adata_in.X
        X_out = adata_out.X.toarray() if hasattr(adata_out.X, 'toarray') else adata_out.X
        
        # Sum counts per gene
        in_counts = np.array(X_in.sum(axis=0)).flatten()
        out_counts = np.array(X_out.sum(axis=0)).flatten()
        total_counts = in_counts + out_counts
        
        # Compute diffusion ratio
        with np.errstate(divide='ignore', invalid='ignore'):
            diffusion_ratio = np.where(
                total_counts > 0,
                out_counts / total_counts,
                0.0
            )
        
        gene_metrics = pd.DataFrame({
            "gene": adata_in.var_names,
            "in_tissue_counts": in_counts,
            "out_tissue_counts": out_counts,
            "total_counts": total_counts,
            "diffusion_ratio": diffusion_ratio
        })
        
        return gene_metrics

    def _calculate_min_distances(
        self,
        adata_in: ad.AnnData,
        adata_out: ad.AnnData,
        batch_size: int
    ) -> np.ndarray:
        """Calculate minimum distance from each outside spot to nearest in-tissue spot.
        
        Based on the original SPATCH calculate_distance function.
        """
        X_in = adata_in.obsm["spatial"]
        X_out = adata_out.obsm["spatial"]
        
        # Fit nearest neighbors on in-tissue points
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree')
        nbrs.fit(X_in)
        
        # Calculate distances in batches
        min_distances = []
        for i in range(0, len(X_out), batch_size):
            batch = X_out[i:i + batch_size]
            distances, _ = nbrs.kneighbors(batch)
            min_distances.extend(distances.flatten())
        
        return np.array(min_distances)

    def _add_distance_correlation(
        self,
        gene_metrics: pd.DataFrame,
        adata_out: ad.AnnData,
        distances: np.ndarray
    ) -> pd.DataFrame:
        """Add distance-expression correlation to gene metrics."""
        from scipy.stats import spearmanr
        
        X_out = adata_out.X.toarray() if hasattr(adata_out.X, 'toarray') else adata_out.X
        
        distance_correlations = []
        for i in range(X_out.shape[1]):
            gene_expr = X_out[:, i]
            if np.std(gene_expr) > 0 and np.std(distances) > 0:
                corr, _ = spearmanr(distances, gene_expr)
                distance_correlations.append(corr)
            else:
                distance_correlations.append(0.0)
        
        gene_metrics["distance_correlation"] = distance_correlations
        
        return gene_metrics
