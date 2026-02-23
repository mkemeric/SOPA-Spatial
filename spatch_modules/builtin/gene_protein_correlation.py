"""
Gene-Protein Correlation Module

Calculate multi-resolution correlation between spatial transcriptomics
gene expression and CODEX protein levels.
Based on the original SPATCH 5_correlation_with_codex.py analysis.
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

import spatialdata as sd
from spatialdata.models import TableModel

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class GeneProteinCorrelation(SpatchModule):
    """Calculate spatial correlation between gene expression and protein levels.
    
    This module enables cross-modal analysis between spatial transcriptomics
    (e.g., Xenium, Visium HD) and protein imaging (e.g., CODEX) data.
    
    The analysis:
    1. Bins both modalities to common spatial grids at multiple resolutions
    2. Matches gene-protein pairs based on provided mapping
    3. Computes correlation coefficients at each resolution
    4. Stores per-pair and aggregate correlation metrics
    """
    
    name = "gene_protein_correlation"
    version = "1.0.0"
    description = "Multi-resolution ST×CODEX correlation analysis"
    category = "analysis"
    requires = []  # Flexible - accepts multiple table keys
    produces = ["tables/gene_protein_correlation"]

    # Default gene-to-protein mapping (can be overridden in config)
    DEFAULT_GENE_PROTEIN_MAP = {
        "EPCAM": "Pan-Cytokeratin",
        "KRT18": "Pan-Cytokeratin",
        "CD3E": "CD3e",
        "CD3D": "CD3e",
        "CD4": "CD4",
        "CD8A": "CD8",
        "FOXP3": "FOXP3",
        "CD68": "CD68",
        "MS4A1": "CD20",
        "CD20": "CD20",
        "NCAM1": "CD56",
        "CD34": "CD34",
        "ACTA2": "SMA",
    }

    def run(
        self,
        sdata: sd.SpatialData,
        gene_table_key: str = "table",
        protein_table_key: str = "codex_table",
        gene_coord_key: str = "spatial",
        protein_coord_key: str = "spatial",
        resolution_um: list[float] = None,
        gene_protein_map: dict = None,
        method: str = "spearman",
        min_spots_per_bin: int = 1,
        **kwargs
    ) -> ModuleResult:
        """Run gene-protein correlation analysis.
        
        Args:
            sdata: SpatialData object with gene and protein tables.
            gene_table_key: Key for gene expression table.
            protein_table_key: Key for protein expression table.
            gene_coord_key: Key in obsm for gene spatial coordinates.
            protein_coord_key: Key in obsm for protein spatial coordinates.
            resolution_um: List of bin sizes in µm for multi-resolution analysis.
            gene_protein_map: Dict mapping gene names to protein names.
            method: Correlation method ('spearman' or 'pearson').
            min_spots_per_bin: Minimum spots required per bin.
            **kwargs: Additional config overrides.
        
        Returns:
            ModuleResult with correlation table added to sdata.
        """
        log = []
        
        # Set defaults
        resolution_um = resolution_um or self.config.get(
            "resolution_um", [100, 200, 300, 400, 500]
        )
        gene_protein_map = gene_protein_map or self.config.get(
            "gene_protein_map", self.DEFAULT_GENE_PROTEIN_MAP
        )
        method = method or self.config.get("method", "spearman")
        
        # Get tables
        if gene_table_key not in sdata.tables:
            raise ValueError(f"Gene table '{gene_table_key}' not found in sdata")
        if protein_table_key not in sdata.tables:
            raise ValueError(f"Protein table '{protein_table_key}' not found in sdata")
        
        adata_gene = sdata.tables[gene_table_key]
        adata_protein = sdata.tables[protein_table_key]
        
        log.append(f"Gene table: {adata_gene.n_obs} cells, {adata_gene.n_vars} genes")
        log.append(f"Protein table: {adata_protein.n_obs} cells, {adata_protein.n_vars} proteins")
        
        # Find matching gene-protein pairs
        available_genes = set(adata_gene.var_names)
        available_proteins = set(adata_protein.var_names)
        
        valid_pairs = []
        for gene, protein in gene_protein_map.items():
            if gene in available_genes and protein in available_proteins:
                valid_pairs.append((gene, protein))
        
        if not valid_pairs:
            log.append("WARNING: No matching gene-protein pairs found")
            return ModuleResult(
                sdata=sdata,
                metrics={"n_pairs": 0, "warning": "No matching pairs"},
                log=log
            )
        
        log.append(f"Found {len(valid_pairs)} matching gene-protein pairs")
        
        # Run correlation at each resolution
        all_results = []
        
        for bin_size in resolution_um:
            # Bin gene expression data
            gene_binned = self._bin_adata(
                adata_gene, gene_coord_key, bin_size
            )
            
            # Bin protein expression data
            protein_binned = self._bin_adata(
                adata_protein, protein_coord_key, bin_size
            )
            
            # Calculate correlations for each pair
            for gene, protein in valid_pairs:
                result = self._calculate_pair_correlation(
                    gene_binned, protein_binned,
                    gene, protein,
                    bin_size, method, min_spots_per_bin
                )
                if result is not None:
                    all_results.append(result)
        
        # Aggregate results
        results_df = pd.DataFrame(all_results)
        
        if len(results_df) > 0:
            # Summary statistics
            mean_corr = results_df.groupby(["gene", "protein"])["correlation"].mean()
            best_resolution = results_df.loc[
                results_df.groupby(["gene", "protein"])["correlation"].idxmax()
            ][["gene", "protein", "resolution_um", "correlation"]]
        
        # Create results table
        results_adata = ad.AnnData(obs=results_df)
        sdata.tables["gene_protein_correlation"] = TableModel.parse(results_adata)
        
        return ModuleResult(
            sdata=sdata,
            metrics={
                "n_pairs": len(valid_pairs),
                "n_resolutions": len(resolution_um),
                "mean_correlation": float(results_df["correlation"].mean()) if len(results_df) > 0 else 0.0,
                "max_correlation": float(results_df["correlation"].max()) if len(results_df) > 0 else 0.0,
            },
            log=log
        )

    def _bin_adata(
        self,
        adata: ad.AnnData,
        coord_key: str,
        bin_size: float
    ) -> ad.AnnData:
        """Bin spatial data to a regular grid.
        
        Based on the original SPATCH bin_adata function.
        """
        coords = adata.obsm[coord_key]
        
        # Compute grid indices
        grid_indices = (coords // bin_size).astype(int)
        
        # Find unique grid cells
        unique_indices, inverse = np.unique(
            grid_indices, axis=0, return_inverse=True
        )
        
        # Get expression matrix
        X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        
        # Aggregate expression by grid cell (sum)
        n_bins = len(unique_indices)
        n_genes = X.shape[1]
        binned_X = np.zeros((n_bins, n_genes))
        
        for i in range(n_bins):
            mask = inverse == i
            binned_X[i] = X[mask].sum(axis=0)
        
        # Create new AnnData
        binned_adata = ad.AnnData(
            X=binned_X,
            var=adata.var.copy()
        )
        binned_adata.obsm["spatial"] = unique_indices * bin_size
        binned_adata.obs["bin_x"] = unique_indices[:, 0]
        binned_adata.obs["bin_y"] = unique_indices[:, 1]
        binned_adata.obs["n_spots"] = np.bincount(inverse, minlength=n_bins)
        
        return binned_adata

    def _calculate_pair_correlation(
        self,
        gene_binned: ad.AnnData,
        protein_binned: ad.AnnData,
        gene: str,
        protein: str,
        resolution: float,
        method: str,
        min_spots: int
    ) -> dict | None:
        """Calculate correlation for a single gene-protein pair."""
        # Get binned values
        gene_expr = gene_binned[:, gene].X.flatten()
        protein_expr = protein_binned[:, protein].X.flatten()
        
        # Create DataFrames for spatial join
        gene_df = pd.DataFrame({
            "bin_x": gene_binned.obs["bin_x"].values,
            "bin_y": gene_binned.obs["bin_y"].values,
            "gene_expr": gene_expr
        })
        
        protein_df = pd.DataFrame({
            "bin_x": protein_binned.obs["bin_x"].values,
            "bin_y": protein_binned.obs["bin_y"].values,
            "protein_expr": protein_expr
        })
        
        # Join on spatial bins
        merged = pd.merge(gene_df, protein_df, on=["bin_x", "bin_y"], how="inner")
        
        if len(merged) < 3:
            return None
        
        # Replace NaN with 0
        merged = merged.fillna(0)
        
        # Calculate correlation
        if method == "spearman":
            corr, pval = stats.spearmanr(merged["gene_expr"], merged["protein_expr"])
        else:
            corr, pval = stats.pearsonr(merged["gene_expr"], merged["protein_expr"])
        
        return {
            "gene": gene,
            "protein": protein,
            "resolution_um": resolution,
            "correlation": corr,
            "p_value": pval,
            "n_bins": len(merged),
            "method": method
        }
