"""
CODEX Loader Module

Loads Akoya CODEX/PhenoCycler data into SpatialData format.
Based on the original SPATCH load_codex() function.

Supported input layouts:
  1. H5AD multiomics — transcriptome/adata.h5ad + proteome/adata_codex.h5ad
     + segmentation_mask/cell_boundaries.csv  (replaces prepare_codex_sdata.py)
  2. QuPath-processed — measurements_major.tsv (+ optional T / FOXP3 TSVs)
  3. Standard CODEX — FCS/compensation/cell_table.csv + optional stitched TIFF
"""

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from geopandas import GeoDataFrame
from shapely.geometry import Point, Polygon

import spatialdata as sd
from spatialdata.models import TableModel, ShapesModel, Image2DModel
from spatialdata.transformations import Identity

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class CODEXLoader(SpatchModule):
    """Load Akoya CODEX/PhenoCycler data into SpatialData format.
    
    Supports multiple CODEX data formats:
    - H5AD multiomics (transcriptome + proteome h5ads + boundary CSV)
    - QuPath-processed measurements (TSV files)
    - Standard CODEX output with stitched TIFFs and cell tables
    
    The loader extracts protein expression, cell spatial coordinates,
    and cell type classifications from the CODEX analysis outputs.
    """
    
    name = "codex_loader"
    version = "2.0.0"
    description = "Load Akoya CODEX/PhenoCycler data into SpatialData format"
    category = "loader"
    requires = []
    produces = [
        "tables/table", "tables/codex_table",
        "shapes/cell_boundaries", "shapes/codex_cells",
        "images/codex_image",
    ]

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
        
        # ── Detect format and dispatch ────────────────────────────
        # Format 1: H5AD multiomics (transcriptome + proteome + boundaries)
        if (data_path / "transcriptome" / "adata.h5ad").exists():
            return self._run_h5ad_multiomics(
                data_path, coord_system, load_image, **kwargs
            )

        # Format 2: QuPath-processed TSVs
        if (data_path / "measurements_major.tsv").exists():
            adata, n_filtered = self._load_qupath_format(
                data_path, detection_prob_quantile
            )
            log.append(f"Loaded QuPath format from {data_path}")
            log.append(f"Filtered {n_filtered} low-quality cells")

        # Format 3: Standard CODEX CSV + optional TIFF
        elif (data_path / "FCS" / "compensation" / "cell_table.csv").exists():
            adata = self._load_standard_format(data_path)
            log.append(f"Loaded standard CODEX format from {data_path}")
            n_filtered = 0
        else:
            raise FileNotFoundError(
                f"Could not find CODEX data files in {data_path}. "
                "Expected one of: transcriptome/adata.h5ad (h5ad multiomics), "
                "measurements_major.tsv (QuPath), or "
                "FCS/compensation/cell_table.csv (standard)."
            )
        
        # ── Shared path for QuPath / standard formats ─────────────
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
                log.append("Loaded multi-channel image")
        
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

    # ── H5AD multiomics format ────────────────────────────────────

    def _run_h5ad_multiomics(
        self,
        data_path: Path,
        coord_system: str,
        load_image: bool = False,
        **kwargs,
    ) -> ModuleResult:
        """Load pre-processed h5ad multiomics data.

        Expected layout under *data_path*::

            transcriptome/adata.h5ad          # gene expression
            proteome/adata_codex.h5ad         # protein expression
            segmentation_mask/cell_boundaries.csv  # vertex CSV

        This replaces the standalone ``prepare_codex_sdata.py`` script.
        """
        log = []

        # ── 1. Transcriptome table ────────────────────────────────
        st_path = data_path / "transcriptome" / "adata.h5ad"
        print(f"Loading transcriptome h5ad: {st_path}")
        st = sc.read_h5ad(st_path)
        log.append(f"Transcriptome: {st.n_obs:,} cells × {st.n_vars:,} genes")

        self._sanitize_uns_keys(st)

        # Map high_quality → in_tissue (used by diffusion_analysis)
        if "high_quality" in st.obs.columns:
            st.obs["in_tissue"] = st.obs["high_quality"].astype(int)
            log.append(
                f"Mapped high_quality → in_tissue ({st.obs['in_tissue'].sum():,} cells)"
            )

        st.obs["cell_id"] = st.obs_names
        st.obs["region"] = "cell_boundaries"

        table = TableModel.parse(
            st,
            region="cell_boundaries",
            region_key="region",
            instance_key="cell_id",
        )

        # ── 2. Proteome (CODEX) table ─────────────────────────────
        codex_path = data_path / "proteome" / "adata_codex.h5ad"
        print(f"Loading CODEX proteome h5ad: {codex_path}")
        codex = sc.read_h5ad(codex_path)
        log.append(
            f"Proteome: {codex.n_obs:,} cells × {codex.n_vars:,} proteins "
            f"({list(codex.var_names)})"
        )

        self._sanitize_uns_keys(codex)

        codex.obs["cell_id"] = codex.obs_names
        codex.obs["region"] = "codex_cells"

        codex_table = TableModel.parse(
            codex,
            region="codex_cells",
            region_key="region",
            instance_key="cell_id",
        )

        # ── 3. Cell boundary polygons ─────────────────────────────
        boundaries_csv = data_path / "segmentation_mask" / "cell_boundaries.csv"
        print(f"Loading cell boundaries: {boundaries_csv}")
        gdf = self._load_cell_boundaries(boundaries_csv)
        shapes = ShapesModel.parse(
            gdf, transformations={coord_system: Identity()}
        )
        log.append(f"Cell boundaries: {len(gdf):,} polygons")

        # ── 4. Optional image ─────────────────────────────────────
        images = {}
        if load_image:
            img_el = self._load_image(data_path, coord_system)
            if img_el is not None:
                images["codex_image"] = img_el
                log.append("Loaded multi-channel image")

        # ── 5. Assemble SpatialData ───────────────────────────────
        result_sdata = sd.SpatialData(
            tables={"table": table, "codex_table": codex_table},
            shapes={"cell_boundaries": shapes},
            images=images,
        )

        print(f"SpatialData built — tables: {list(result_sdata.tables.keys())}, "
              f"shapes: {list(result_sdata.shapes.keys())}")

        return ModuleResult(
            sdata=result_sdata,
            metrics={
                "n_cells": st.n_obs,
                "n_genes": st.n_vars,
                "n_proteins": codex.n_vars,
                "n_boundaries": len(gdf),
            },
            log=log,
        )

    # ── Shared helpers ──────────────────────────────────────────────

    @staticmethod
    def _sanitize_uns_keys(adata: ad.AnnData) -> None:
        """Rename uns keys containing spaces or special characters.

        Newer spatialdata (>=0.6) requires uns keys to match
        ``[a-zA-Z0-9_.-]+``.
        """
        bad = [
            k for k in list(adata.uns.keys())
            if not re.match(r'^[a-zA-Z0-9_.\-]+$', k)
        ]
        for k in bad:
            new_key = re.sub(r'[^a-zA-Z0-9_.\-]', '_', k)
            adata.uns[new_key] = adata.uns.pop(k)
        if bad:
            print(f"  Sanitized {len(bad)} uns key(s): {bad}")

    @staticmethod
    def _load_cell_boundaries(
        csv_path: Path,
        max_cells: int | None = None,
    ) -> GeoDataFrame:
        """Convert a vertex CSV to a GeoDataFrame of Shapely polygons.

        Expected columns: ``cell_id, vertex_x, vertex_y, label_id``.
        """
        print(f"  Reading {csv_path} ...")
        df = pd.read_csv(csv_path)
        print(f"  {len(df):,} vertex rows")

        grouped = df.groupby("cell_id")
        if max_cells:
            grouped = dict(list(grouped)[:max_cells])
        else:
            grouped = dict(grouped)

        cell_ids: list = []
        polygons: list = []
        skipped = 0

        for cell_id, verts in grouped.items():
            coords = verts[["vertex_x", "vertex_y"]].values
            if len(coords) < 3:
                skipped += 1
                continue
            try:
                poly = Polygon(coords)
                if poly.is_valid and not poly.is_empty:
                    cell_ids.append(cell_id)
                    polygons.append(poly)
                else:
                    skipped += 1
            except Exception:
                skipped += 1

        print(f"  Built {len(polygons):,} polygons ({skipped} skipped)")
        return GeoDataFrame({"geometry": polygons}, index=cell_ids)

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
