"""
Leiden Clustering Module

Multi-resolution Leiden clustering on the neighbor graph computed by
``qc_preprocessing``. Promotes step 6 ("Clustering") of
notebooks/01_spatch_workflow_example.ipynb into a YAML-configured,
reusable module — resolutions live in the step's ``config:`` block
instead of hardcoded notebook cells.
"""

from __future__ import annotations

import scanpy as sc
import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class LeidenClustering(SpatchModule):
    """Leiden clustering at one or more resolutions."""

    name = "leiden_clustering"
    version = "1.0.0"
    description = "Multi-resolution Leiden clustering"
    category = "analysis"
    requires = []
    produces = []

    def validate_inputs(self, sdata: sd.SpatialData) -> list[str]:
        table_key = self.config.get("table_key", "table")
        if table_key not in sdata.tables:
            return [f"Missing required table: '{table_key}'"]
        adata = sdata.tables[table_key]
        if "neighbors" not in adata.uns:
            return [
                f"Table '{table_key}' has no neighbor graph — "
                "run qc_preprocessing (or sc.pp.neighbors) first"
            ]
        return []

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        resolutions: list[float] | None = None,
        key_prefix: str = "leiden_",
        flavor: str = "igraph",
        n_iterations: int = 2,
        random_state: int = 0,
        **kwargs,
    ) -> ModuleResult:
        """Run Leiden clustering at each requested resolution.

        Args:
            sdata: SpatialData object with a precomputed neighbor graph
                (see qc_preprocessing).
            table_key: Key for the cell table in sdata.tables.
            resolutions: Leiden resolutions to run. Each produces its own
                obs column named f"{key_prefix}{resolution}"
                (e.g. "leiden_0.5"). Defaults to [0.5].
            key_prefix: Prefix for the obs column name of each resolution.
            flavor / n_iterations: Passed to scanpy's Leiden implementation.
            random_state: Seed for reproducibility.

        Returns:
            ModuleResult with one cluster column per resolution.
        """
        log = []
        metrics: dict = {}
        adata = sdata.tables[table_key]
        resolutions = resolutions or [0.5]

        for res in resolutions:
            key = f"{key_prefix}{res}"
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=key,
                flavor=flavor,
                n_iterations=n_iterations,
                random_state=random_state,
            )
            n_clusters = int(adata.obs[key].nunique())
            metrics[f"n_clusters_{key}"] = n_clusters
            log.append(f"Leiden @ resolution={res} -> {n_clusters} clusters ({key})")

        sdata.tables[table_key] = adata

        return ModuleResult(sdata=sdata, metrics=metrics, log=log)
