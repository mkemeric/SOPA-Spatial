"""
Cell Shape Metrics Module

Compute morphological metrics from cell boundary polygons using Shapely.
Based on the original SPATCH 8_cell_shape.py analysis, but reimplemented
using Shapely/GeoPandas instead of OpenCV for SpatialData compatibility.
"""

import numpy as np
import pandas as pd

import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class CellShapeMetrics(SpatchModule):
    """Compute morphological metrics from cell boundary polygons.
    
    This module calculates cell shape descriptors commonly used in
    spatial biology analysis, including:
    - Area and perimeter
    - Circularity (isoperimetric quotient)
    - Eccentricity (from minimum bounding rectangle)
    - Solidity (area / convex hull area)
    - Aspect ratio
    
    Works directly on SpatialData shapes (GeoPandas GeoDataFrames)
    without requiring OpenCV.
    """
    
    name = "cell_shape_metrics"
    version = "1.0.0"
    description = "Compute Shapely-based cell morphology metrics"
    category = "analysis"
    requires = ["shapes/cell_boundaries", "tables/table"]
    produces = []  # Adds columns to existing table

    def run(
        self,
        sdata: sd.SpatialData,
        boundaries_key: str = "cell_boundaries",
        table_key: str = "table",
        add_to_table: bool = True,
        **kwargs
    ) -> ModuleResult:
        """Compute cell shape metrics.
        
        Args:
            sdata: SpatialData object with cell boundaries.
            boundaries_key: Key for cell boundary shapes in sdata.shapes.
            table_key: Key for cell table in sdata.tables.
            add_to_table: Whether to add metrics to existing table.
            **kwargs: Additional config overrides.
        
        Returns:
            ModuleResult with shape metrics added to sdata.
        """
        log = []
        
        # Get cell boundaries
        if boundaries_key not in sdata.shapes:
            raise ValueError(f"Shapes '{boundaries_key}' not found in sdata")
        
        boundaries = sdata.shapes[boundaries_key]
        log.append(f"Processing {len(boundaries)} cell boundaries")
        
        # Compute metrics for each cell
        metrics = self._compute_shape_metrics(boundaries)
        
        # Add to table if requested
        if add_to_table and table_key in sdata.tables:
            adata = sdata.tables[table_key]
            
            # Match by index if possible
            if len(metrics) == adata.n_obs:
                for col in ["area_um2", "perimeter_um", "circularity", 
                           "eccentricity", "solidity", "aspect_ratio"]:
                    adata.obs[col] = metrics[col].values
                log.append(f"Added shape metrics to {table_key}")
            else:
                log.append(
                    f"WARNING: Shape count ({len(metrics)}) doesn't match "
                    f"table rows ({adata.n_obs}). Metrics not added to table."
                )
        
        # Summary statistics
        summary = {
            "mean_area": float(metrics["area_um2"].mean()),
            "median_area": float(metrics["area_um2"].median()),
            "mean_circularity": float(metrics["circularity"].mean()),
            "mean_eccentricity": float(metrics["eccentricity"].mean()),
            "mean_solidity": float(metrics["solidity"].mean()),
            "n_cells": len(metrics)
        }
        
        return ModuleResult(
            sdata=sdata,
            metrics=summary,
            log=log
        )

    def _compute_shape_metrics(self, boundaries) -> pd.DataFrame:
        """Compute shape metrics for all cell boundaries."""
        areas = []
        perimeters = []
        circularities = []
        eccentricities = []
        solidities = []
        aspect_ratios = []
        
        for geom in boundaries.geometry:
            # Skip invalid or empty geometries
            if geom is None or geom.is_empty:
                areas.append(np.nan)
                perimeters.append(np.nan)
                circularities.append(np.nan)
                eccentricities.append(np.nan)
                solidities.append(np.nan)
                aspect_ratios.append(np.nan)
                continue
            
            # Basic measurements
            area = geom.area
            perimeter = geom.length
            
            areas.append(area)
            perimeters.append(perimeter)
            
            # Circularity: 4π × area / perimeter²
            # Perfect circle = 1.0
            if perimeter > 0:
                circularity = (4 * np.pi * area) / (perimeter ** 2)
            else:
                circularity = 0.0
            circularities.append(circularity)
            
            # Eccentricity from minimum rotated bounding rectangle
            try:
                mbr = geom.minimum_rotated_rectangle
                coords = list(mbr.exterior.coords)
                
                # Calculate edge lengths
                edge1 = np.sqrt(
                    (coords[1][0] - coords[0][0])**2 + 
                    (coords[1][1] - coords[0][1])**2
                )
                edge2 = np.sqrt(
                    (coords[2][0] - coords[1][0])**2 + 
                    (coords[2][1] - coords[1][1])**2
                )
                
                major = max(edge1, edge2)
                minor = min(edge1, edge2)
                
                # Eccentricity: sqrt(1 - (minor/major)²)
                # Circle = 0, line = 1
                if major > 0:
                    ecc = np.sqrt(1 - (minor / major) ** 2)
                    aspect = major / minor if minor > 0 else 0
                else:
                    ecc = 0.0
                    aspect = 1.0
            except Exception:
                ecc = 0.0
                aspect = 1.0
            
            eccentricities.append(ecc)
            aspect_ratios.append(aspect)
            
            # Solidity: area / convex hull area
            try:
                convex_area = geom.convex_hull.area
                if convex_area > 0:
                    solidity = area / convex_area
                else:
                    solidity = 0.0
            except Exception:
                solidity = 0.0
            
            solidities.append(solidity)
        
        return pd.DataFrame({
            "area_um2": areas,
            "perimeter_um": perimeters,
            "circularity": circularities,
            "eccentricity": eccentricities,
            "solidity": solidities,
            "aspect_ratio": aspect_ratios
        })
