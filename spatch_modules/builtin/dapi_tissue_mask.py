"""
DAPI Tissue Mask Module

Generate a binary tissue mask from a DAPI image and classify
spots/cells as in-tissue or out-of-tissue.
Based on the original SPATCH 2_8um_bin.py ``dapi_mask`` and
``tissue_extraction`` functions.

Requires: opencv-python-headless  (install with ``pip install spatch-modules[imaging]``)
"""

import numpy as np
import pandas as pd
import anndata as ad

import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity

from ..base import SpatchModule, ModuleResult
from ..registry import register


def _ensure_cv2():
    """Lazy import of OpenCV with a helpful error message."""
    try:
        import cv2
        return cv2
    except ImportError:
        raise ImportError(
            "opencv-python-headless is required for DAPI tissue masking. "
            "Install it with: pip install spatch-modules[imaging]"
        )


@register
class DapiTissueMask(SpatchModule):
    """Generate a tissue mask from DAPI fluorescence and tag cells.

    Pipeline:
      1. Extract a single DAPI channel from the SpatialData image
      2. Normalise to uint8 and apply Otsu thresholding
      3. Morphological dilation + erosion to fill holes
      4. Select the largest contour as the tissue boundary
      5. Convert the contour to a Shapely polygon (stored in sdata.shapes)
      6. Tag every cell/spot in the table as in_tissue = 0 | 1
    """

    name = "dapi_tissue_mask"
    version = "1.0.0"
    description = "DAPI-based tissue boundary detection and cell tagging"
    category = "qc"
    requires = ["images"]
    produces = ["shapes/tissue_boundary"]

    def run(
        self,
        sdata: sd.SpatialData,
        image_key: str | None = None,
        channel: int | str = 0,
        kernel_width: int = 500,
        dilate_iterations: int = 2,
        erode_iterations: int = 1,
        table_key: str = "table",
        coord_key: str = "spatial",
        coord_system: str = "global",
        output_shapes_key: str = "tissue_boundary",
        in_tissue_col: str = "in_tissue",
        **kwargs,
    ) -> ModuleResult:
        """Run DAPI tissue masking.

        Args:
            sdata: SpatialData object containing at least one image.
            image_key: Key for the image element.  If None the first
                available image is used.
            channel: Channel index or name for the DAPI stain.
            kernel_width: Width of the square morphological kernel.
            dilate_iterations: Dilation iterations.
            erode_iterations: Erosion iterations.
            table_key: Key for the cell/spot table to tag.
            coord_key: obsm key for spatial coordinates in the table.
            coord_system: Coordinate system name for the output shapes.
            output_shapes_key: Key under which the tissue polygon is stored.
            in_tissue_col: Column name added to adata.obs.
            **kwargs: Extra config (ignored).

        Returns:
            ModuleResult with tissue polygon in shapes and in_tissue column.
        """
        cv2 = _ensure_cv2()
        log = []

        # ── 1. Load DAPI channel ────────────────────────────────────
        if image_key is None:
            image_key = next(iter(sdata.images))
        image_element = sdata.images[image_key]

        # Resolve to numpy (handle DataTree / MultiscaleSpatialImage)
        if hasattr(image_element, "values"):
            # DataTree – take highest resolution
            ds = next(iter(image_element.values()))
            img_arr = ds["image"].values if "image" in ds else next(iter(ds.values())).values
        else:
            img_arr = np.asarray(image_element)

        # Select channel (CYX layout)
        if isinstance(channel, str):
            # TODO: look up channel name if available
            channel = 0
        if img_arr.ndim == 3:
            dapi = img_arr[channel]
        else:
            dapi = img_arr

        log.append(f"DAPI image shape: {dapi.shape}")

        # ── 2. Generate binary mask ─────────────────────────────────
        mask = self._dapi_mask(cv2, dapi, kernel_width, dilate_iterations, erode_iterations)
        log.append("Generated binary tissue mask via Otsu + morphology")

        # ── 3. Extract largest contour as polygon ───────────────────
        contours, _ = cv2.findContours(
            mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )
        if not contours:
            log.append("WARNING: No contours found in DAPI mask")
            return ModuleResult(sdata=sdata, metrics={"n_contours": 0}, log=log)

        max_contour = max(contours, key=cv2.contourArea)
        log.append(f"Largest contour area: {cv2.contourArea(max_contour):.0f} px²")

        # Convert contour to Shapely polygon
        from shapely.geometry import Polygon
        from geopandas import GeoDataFrame

        coords = max_contour.squeeze()
        if coords.ndim == 1:
            coords = coords.reshape(-1, 2)
        # OpenCV contours are (x, y)
        polygon = Polygon(coords.tolist())

        gdf = GeoDataFrame(
            {"geometry": [polygon]},
            index=["tissue"],
        )
        shapes = ShapesModel.parse(
            gdf,
            transformations={coord_system: Identity()},
        )
        sdata.shapes[output_shapes_key] = shapes
        log.append(f"Stored tissue polygon in shapes['{output_shapes_key}']")

        # ── 4. Tag cells in table ───────────────────────────────────
        n_tagged = 0
        if table_key in sdata.tables:
            adata = sdata.tables[table_key]
            if coord_key in adata.obsm:
                in_tissue = self._tissue_extraction(
                    adata, coord_key, mask
                )
                adata.obs[in_tissue_col] = in_tissue
                n_tagged = int(in_tissue.sum())
                log.append(
                    f"Tagged {n_tagged}/{adata.n_obs} cells as in-tissue"
                )
            else:
                log.append(
                    f"WARNING: coord_key '{coord_key}' not in table obsm, "
                    "skipping cell tagging"
                )

        return ModuleResult(
            sdata=sdata,
            metrics={
                "mask_area_px": int(mask.sum() / 255),
                "contour_area_px": int(cv2.contourArea(max_contour)),
                "n_in_tissue": n_tagged,
            },
            log=log,
        )

    # ── helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _dapi_mask(
        cv2,
        image: np.ndarray,
        kernel_width: int,
        dilate_iter: int,
        erode_iter: int,
    ) -> np.ndarray:
        """Generate binary mask from DAPI image.

        Direct port of SPATCH ``dapi_mask``.
        """
        # Normalise to uint8
        img = (image / np.max(image) * 255).astype(np.uint8)

        # Otsu threshold
        _, otsu = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

        # Morphological operations
        kernel = np.ones((kernel_width, kernel_width), np.uint8)
        dilated = cv2.dilate(otsu, kernel, iterations=dilate_iter)
        eroded = cv2.erode(dilated, kernel, iterations=erode_iter)

        # Keep only the largest connected component
        contours, _ = cv2.findContours(
            eroded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )
        if not contours:
            return eroded

        max_contour = max(contours, key=cv2.contourArea)
        mask = np.zeros_like(img)
        cv2.drawContours(mask, [max_contour], -1, 255, thickness=cv2.FILLED)
        return mask

    @staticmethod
    def _tissue_extraction(
        adata: ad.AnnData,
        coord_key: str,
        mask: np.ndarray,
    ) -> np.ndarray:
        """Classify each cell as in-tissue (1) or not (0).

        Direct port of SPATCH ``tissue_extraction``.
        """
        loc = adata.obsm[coord_key].astype(int)
        h, w = mask.shape

        in_tissue = np.zeros(adata.n_obs, dtype=int)
        for i in range(adata.n_obs):
            x, y = loc[i, 0], loc[i, 1]
            if 0 <= y < h and 0 <= x < w and mask[y, x] == 255:
                in_tissue[i] = 1

        return in_tissue
