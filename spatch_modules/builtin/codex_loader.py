"""
CODEX Loader Module

Loads Akoya CODEX/PhenoCycler data into SpatialData format.
Based on the original SPATCH load_codex() function.
"""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import anndata as ad
from geopandas import GeoDataFrame
from shapely.geometry import Point

import spatialdata as sd
from spatialdata.models import TableModel, ShapesModel, Image2DModel
from spatialdata.transformations import Identity

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class CODEXLoader(SpatchModule):
    """Load Akoya CODEX/PhenoCycler data into SpatialData format.
    
    Supports multiple CODEX data formats:
    - QuPath-processed measurements (TSV files)
    - Standard CODEX output with stitched TIFFs and cell tables
    
    The loader extracts protein expression, cell spatial coordinates,
    and cell type classifications from the CODEX analysis outputs.
    """
    
    name = "codex_loader"
    version = "1.0.0"
    description = "Load Akoya CODEX/PhenoCycler data into SpatialData format"
    category = "loader"
    requires = []
    produces = ["images/codex_image", "tables/codex_table", "shapes/codex_cells"]

    # Cell type mapping from marker classifications to canonical types
    MAJOR_TYPE_MAP = {
        "Pan-Cytokeratin": "Epithelial",
        "CD34": "Endothelial",
        "CD3e": "T",
        "SMA": "Fibroblast",
        "CD68": "Macrophage",
        "CD20": "B",
        "CD56": "NK",
        "CD4": "CD4T",
        "CD8": "CD8T",
        "FOXP3": "Treg"
    }

    def run(
        self,
        sdata: sd.SpatialData = None,
        data_path: str = "",
        load_image: bool = False,
        detection_prob_quantile: float = 0.1,
        coord_system: str = "global",
        **kwargs
    ) -> ModuleResult:
        """Load CODEX data from disk.
        
        Args:
            sdata: Existing SpatialData object (ignored for loader).
            data_path: Path to CODEX output directory.
            load_image: Whether to load the multi-channel TIFF image.
            detection_prob_quantile: Quantile threshold for filtering low-quality cells.
            coord_system: Name for the coordinate system.
            **kwargs: Additional config overrides.
        
        Returns:
            ModuleResult with new SpatialData object containing CODEX data.
        """
        data_path = Path(data_path or self.config.get("data_path", ""))
        load_image = load_image or self.config.get("load_image", False)
        detection_prob_quantile = detection_prob_quantile or self.config.get(
            "detection_prob_quantile", 0.1
        )
        
        log = []
        elements = {}
        
        # Detect CODEX data format and load accordingly
        if (data_path / "measurements_major.tsv").exists():
            # QuPath-processed format
            adata, n_filtered = self._load_qupath_format(
                data_path, detection_prob_quantile
            )
            log.append(f"Loaded QuPath format from {data_path}")
            log.append(f"Filtered {n_filtered} low-quality cells")
        elif (data_path / "FCS" / "compensation" / "cell_table.csv").exists():
            # Standard CODEX format
            adata = self._load_standard_format(data_path)
            log.append(f"Loaded standard CODEX format from {data_path}")
            n_filtered = 0
        else:
            raise FileNotFoundError(
                f"Could not find CODEX data files in {data_path}. "
                "Expected either measurements_major.tsv (QuPath) or "
                "FCS/compensation/cell_table.csv (standard)."
            )
        
        # Create shapes from cell centroids
        geometries = [
            Point(x, y) for x, y in adata.obsm["spatial"]
        ]
        shapes_gdf = GeoDataFrame(
            {"geometry": geometries},
            index=adata.obs_names
        )
        shapes = ShapesModel.parse(
            shapes_gdf,
            transformations={coord_system: Identity()}
        )
        elements["codex_cells"] = shapes
        log.append(f"Created {len(geometries)} cell point geometries")
        
        # Create table
        table = TableModel.parse(
            adata,
            region="codex_cells",
            region_key="region",
            instance_key="cell_id"
        )
        
        # Load image if requested
        if load_image:
            image_element = self._load_image(data_path, coord_system)
            if image_element is not None:
                elements["codex_image"] = image_element
                log.append(f"Loaded multi-channel image")
        
        # Build SpatialData object
        result_sdata = sd.SpatialData(
            images={"codex_image": elements.get("codex_image")} if "codex_image" in elements else {},
            shapes={"codex_cells": shapes},
            tables={"codex_table": table}
        )
        
        return ModuleResult(
            sdata=result_sdata,
            metrics={
                "n_cells": adata.n_obs,
                "n_proteins": adata.n_vars,
                "n_filtered": n_filtered if 'n_filtered' in dir() else 0
            },
            log=log
        )

    def _load_qupath_format(
        self,
        data_path: Path,
        detection_prob_quantile: float
    ) -> tuple[ad.AnnData, int]:
        """Load CODEX data from QuPath-processed TSV files.
        
        This is the format produced by the original SPATCH pipeline
        using StarDist segmentation and QuPath analysis.
        """
        # Load measurement files
        data = pd.read_table(data_path / "measurements_major.tsv")
        
        # Load T cell subtype classifications if available
        t_path = data_path / "measurements_t.tsv"
        foxp3_path = data_path / "measurements_foxp3.tsv"
        
        if t_path.exists():
            t_data = pd.read_table(t_path)
            data["T"] = t_data["Classification"]
        
        if foxp3_path.exists():
            foxp3_data = pd.read_table(foxp3_path)
            data["FOXP3"] = foxp3_data["Classification"]
        
        # Build hierarchical cell type annotations
        data["major"] = data["Classification"]
        data["minor"] = data["Classification"]
        
        # Refine T cell subtypes
        if "T" in data.columns:
            mask_cd4 = (data["major"] == "CD3e") & (data["T"] == "CD4")
            mask_cd8 = (data["major"] == "CD3e") & (data["T"] == "CD8")
            data.loc[mask_cd4, "minor"] = "CD4"
            data.loc[mask_cd8, "minor"] = "CD8"
        
        if "FOXP3" in data.columns:
            mask_foxp3 = (data["minor"] == "CD4") & (data["FOXP3"] == "FOXP3")
            data.loc[mask_foxp3, "minor"] = "FOXP3"
        
        # Remove DAPI channel from intensity data
        data = data.loc[:, ~data.columns.str.contains("DAPI")]
        
        # Extract protein expression matrix (Mean intensities)
        intensity_cols = [c for c in data.columns if "Cell: Mean" in c]
        counts = data[intensity_cols].copy()
        counts.columns = counts.columns.str.replace(": Cell: Mean", "", regex=False)
        
        # Build AnnData object
        adata = ad.AnnData(counts.values.astype(np.float32))
        adata.var_names = counts.columns.tolist()
        adata.obs_names = data["Object ID"].astype(str).tolist()
        
        # Add spatial coordinates
        adata.obsm["spatial"] = data[
            ["Centroid X µm", "Centroid Y µm"]
        ].values.astype(np.float32)
        
        # Add cell type annotations
        adata.obs["major"] = data["major"].map(self.MAJOR_TYPE_MAP).values
        adata.obs["minor"] = data["minor"].map(self.MAJOR_TYPE_MAP).values
        adata.obs["cell_id"] = adata.obs_names
        adata.obs["region"] = "codex_cells"
        
        # Filter low-quality cells based on detection probability
        if "Detection probability" in data.columns:
            prob_threshold = np.quantile(
                data["Detection probability"],
                detection_prob_quantile
            )
            keep_mask = (
                (data["Detection probability"] >= prob_threshold) &
                (~data["major"].isna())
            )
            n_filtered = (~keep_mask).sum()
            adata = adata[keep_mask.values, :].copy()
        else:
            keep_mask = ~data["major"].isna()
            n_filtered = (~keep_mask).sum()
            adata = adata[keep_mask.values, :].copy()
        
        return adata, n_filtered

    def _load_standard_format(self, data_path: Path) -> ad.AnnData:
        """Load CODEX data from standard CSV cell table format."""
        cells_path = data_path / "FCS" / "compensation" / "cell_table.csv"
        cells_df = pd.read_csv(cells_path)
        
        # Identify protein intensity columns (exclude metadata)
        meta_cols = {"x", "y", "cell_id", "area", "size", "Cell ID", "X", "Y"}
        protein_cols = [c for c in cells_df.columns if c not in meta_cols]
        
        # Determine coordinate columns
        x_col = "x" if "x" in cells_df.columns else "X"
        y_col = "y" if "y" in cells_df.columns else "Y"
        id_col = "cell_id" if "cell_id" in cells_df.columns else "Cell ID"
        
        # Build AnnData
        adata = ad.AnnData(
            X=cells_df[protein_cols].values.astype(np.float32),
            obs=pd.DataFrame({
                "cell_id": cells_df[id_col].astype(str),
                "region": "codex_cells"
            }),
            var=pd.DataFrame(index=protein_cols)
        )
        adata.obs_names = adata.obs["cell_id"]
        
        # Add coordinates
        adata.obsm["spatial"] = cells_df[[x_col, y_col]].values.astype(np.float32)
        
        # Add area if available
        if "area" in cells_df.columns:
            adata.obs["area"] = cells_df["area"].values
        
        return adata

    def _load_image(
        self,
        data_path: Path,
        coord_system: str
    ) -> Any:
        """Load multi-channel CODEX TIFF image if available."""
        # Try common image locations
        image_paths = [
            data_path / "stitched" / "reg001_stitched.tif",
            data_path / "stitched.tif",
            data_path / "image.tif",
        ]
        
        for image_path in image_paths:
            if image_path.exists():
                try:
                    from tifffile import imread
                    img = imread(str(image_path))
                    
                    # Ensure CYX dimension order
                    if img.ndim == 3 and img.shape[-1] < img.shape[0]:
                        img = np.moveaxis(img, -1, 0)
                    
                    return Image2DModel.parse(
                        img,
                        transformations={coord_system: Identity()},
                        dims=("c", "y", "x")
                    )
                except Exception as e:
                    print(f"Warning: Failed to load image {image_path}: {e}")
        
        return None
