"""
Spatial Statistics Module

Squidpy-based spatial analyses (neighborhood enrichment, Ripley's
statistics, co-occurrence) keyed on a cluster column produced by
``leiden_clustering``. Promotes step 7 ("Spatial Analysis") of
notebooks/01_spatch_workflow_example.ipynb into a YAML-configured,
reusable module.
"""

from __future__ import annotations

import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class SpatialStats(SpatchModule):
    """Squidpy spatial statistics on a cluster column."""

    name = "spatial_stats"
    version = "1.0.0"
    description = "Spatial neighbor graph + neighborhood enrichment/Ripley/co-occurrence"
    category = "analysis"
    requires = []
    produces = []

    def validate_inputs(self, sdata: sd.SpatialData) -> list[str]:
        table_key = self.config.get("table_key", "table")
        if table_key not in sdata.tables:
            return [f"Missing required table: '{table_key}'"]
        adata = sdata.tables[table_key]
        cluster_key = self.config.get("cluster_key", "leiden_0.5")
        if cluster_key not in adata.obs.columns:
            return [
                f"Missing cluster column '{cluster_key}' — "
                "run leiden_clustering first"
            ]
        return []

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        cluster_key: str = "leiden_0.5",
        coord_type: str = "generic",
        n_neighs: int = 30,
        neighborhood: dict | None = None,
        ripley: dict | None = None,
        cooccurrence: dict | None = None,
        interactions: dict | None = None,
        **kwargs,
    ) -> ModuleResult:
        """Build a spatial neighbor graph and run the requested analyses.

        Args:
            sdata: SpatialData object with a cluster column (see
                leiden_clustering).
            table_key: Key for the cell table in sdata.tables.
            cluster_key: obs column to use as the cluster label for every
                analysis below.
            coord_type / n_neighs: Passed to squidpy's spatial neighbor
                graph construction.
            neighborhood / ripley / cooccurrence / interactions: Each a
                ``{"enabled": bool}`` dict toggling that analysis.
                `neighborhood` defaults to enabled; the others default to
                disabled since they're not exercised by the reference
                workflow and can be memory/time intensive on large panels.

        Returns:
            ModuleResult noting which analyses ran; results are stored in
            adata.uns by the underlying squidpy calls.
        """
        import squidpy as sq

        log = []
        metrics: dict = {}
        adata = sdata.tables[table_key]

        sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)
        log.append(f"Built spatial graph (coord_type={coord_type}, n_neighs={n_neighs})")

        if (neighborhood or {}).get("enabled", True):
            sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
            log.append(f"Neighborhood enrichment on '{cluster_key}'")
            metrics["neighborhood_enrichment"] = True

        if (ripley or {}).get("enabled", False):
            sq.gr.ripley(adata, cluster_key=cluster_key, mode="L")
            log.append(f"Ripley's L statistic on '{cluster_key}'")
            metrics["ripley"] = True

        if (cooccurrence or {}).get("enabled", False):
            sq.gr.co_occurrence(adata, cluster_key=cluster_key)
            log.append(f"Co-occurrence on '{cluster_key}'")
            metrics["cooccurrence"] = True

        if (interactions or {}).get("enabled", False):
            log.append(
                "Skipped ligand-receptor interactions — requires a curated "
                "LR database, not yet implemented in spatial_stats"
            )

        sdata.tables[table_key] = adata

        return ModuleResult(sdata=sdata, metrics=metrics, log=log)
