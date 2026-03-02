"""
Image Registration Module

Landmark-based cross-platform image registration using SimpleITK.
Aligns spatial transcriptomics DAPI images to CODEX protein images
(or any pair of images) and applies the resulting transform to both
image channels and AnnData spatial coordinates.

Based on the original SPATCH 3_registeration.py analysis.

Requires: SimpleITK  (install with ``pip install spatch-modules[registration]``)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import anndata as ad

import spatialdata as sd
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

from ..base import SpatchModule, ModuleResult
from ..registry import register


def _ensure_sitk():
    """Lazy import of SimpleITK."""
    try:
        import SimpleITK as sitk
        return sitk
    except ImportError:
        raise ImportError(
            "SimpleITK is required for image registration. "
            "Install it with: pip install spatch-modules[registration]"
        )


@register
class ImageRegistration(SpatchModule):
    """Landmark-based cross-platform image registration.

    Workflow:
      1. Load landmark point-pairs from TSV files (e.g. exported from
         napari or QuPath).
      2. Compute an initial Similarity2D transform from the landmarks.
      3. Refine with multi-resolution gradient descent on Mattes Mutual
         Information.
      4. Apply the transform to:
         a) CODEX multi-channel image → save as registered TIFF / SpatialData image
         b) AnnData spatial coordinates → new ``codex_registered`` obsm key

    The final transform parameters are stored in ``sdata.attrs`` so the
    registration is fully reproducible.
    """

    name = "image_registration"
    version = "1.0.0"
    description = "Landmark-based SimpleITK registration (ST ↔ CODEX)"
    category = "analysis"
    requires = ["images"]
    produces = []  # dynamic: images/codex_registered, tables updated

    def run(
        self,
        sdata: sd.SpatialData,
        fixed_image_key: str = "morphology_focus",
        moving_image_key: str = "codex_image",
        landmark_dir: str = "",
        fixed_landmark_file: str = "nano.tif-points.tsv",
        moving_landmark_file: str = "codex.tif-points.tsv",
        fixed_resolution: float = 0.12,
        codex_resolution: float | None = None,
        n_iterations: int = 500,
        platform: str = "nano",
        apply_to_table: str | None = None,
        coord_key: str = "spatial",
        registered_coord_key: str = "codex_registered",
        save_registered_image: bool = False,
        output_image_key: str = "codex_registered",
        output_tiff_path: str | None = None,
        coord_system: str = "global",
        **kwargs,
    ) -> ModuleResult:
        """Run image registration.

        Args:
            sdata: SpatialData object.
            fixed_image_key: Image key for the fixed (reference) image.
            moving_image_key: Image key for the moving image (e.g. CODEX).
            landmark_dir: Directory containing landmark TSV files.
            fixed_landmark_file: Filename for fixed-image landmarks.
            moving_landmark_file: Filename for moving-image landmarks.
            fixed_resolution: Pixel spacing (µm/px) of the fixed image.
            codex_resolution: Pixel spacing of the moving image.
                If None, will attempt to read from image metadata.
            n_iterations: Max iterations for the optimizer.
            platform: Platform type ("nano" applies a y-flip to landmarks).
            apply_to_table: If set, apply the transform to this table's
                spatial coordinates.
            coord_key: Source coordinate key in adata.obsm.
            registered_coord_key: Destination coordinate key for
                transformed coordinates.
            save_registered_image: Whether to resample & store the
                registered moving image.
            output_image_key: Key for the registered image in sdata.
            output_tiff_path: Optional path to save registered image
                as OME-TIFF.
            coord_system: Coordinate system name.
            **kwargs: Extra config (ignored).

        Returns:
            ModuleResult with transform applied.
        """
        sitk = _ensure_sitk()
        log = []

        # ── 1. Load images as numpy arrays ──────────────────────────
        fix_np = self._image_to_numpy(sdata.images[fixed_image_key])
        mv_np = self._image_to_numpy(sdata.images[moving_image_key])
        log.append(f"Fixed image shape: {fix_np.shape}")
        log.append(f"Moving image shape: {mv_np.shape}")

        # Use first channel if multi-channel
        fix_2d = fix_np[0] if fix_np.ndim == 3 else fix_np
        mv_2d = mv_np[0] if mv_np.ndim == 3 else mv_np

        # ── 2. Load landmarks ───────────────────────────────────────
        lm_dir = Path(landmark_dir)
        fix_lm = pd.read_csv(lm_dir / fixed_landmark_file, sep="\t", header=0)
        mv_lm = pd.read_csv(lm_dir / moving_landmark_file, sep="\t", header=0)

        if codex_resolution is None:
            codex_resolution = self.config.get("codex_resolution", 0.5)

        log.append(
            f"Loaded {len(fix_lm)} fixed / {len(mv_lm)} moving landmarks"
        )

        # ── 3. Run registration ─────────────────────────────────────
        fix_sitk, final_transform = self._register(
            sitk,
            fix_2d,
            mv_2d,
            fix_lm,
            mv_lm,
            fixed_resolution,
            codex_resolution,
            n_iterations,
            platform,
        )
        log.append("Registration completed")

        # Store transform parameters for reproducibility
        transform_params = {
            "type": "Similarity2D",
            "parameters": list(final_transform.GetParameters()),
            "fixed_parameters": list(final_transform.GetFixedParameters()),
            "fixed_resolution": fixed_resolution,
            "codex_resolution": codex_resolution,
        }
        if not hasattr(sdata, "attrs") or sdata.attrs is None:
            sdata.attrs = {}
        sdata.attrs["registration_transform"] = transform_params
        log.append("Stored transform parameters in sdata.attrs")

        # ── 4. Apply transform to coordinates ───────────────────────
        if apply_to_table and apply_to_table in sdata.tables:
            adata = sdata.tables[apply_to_table]
            transformed = self._apply_transform_to_coords(
                sitk, adata, coord_key, codex_resolution, final_transform, mv_lm
            )
            adata.obsm[registered_coord_key] = transformed
            log.append(
                f"Applied transform to table '{apply_to_table}' "
                f"→ obsm['{registered_coord_key}']"
            )

        # ── 5. Optionally register full moving image ────────────────
        if save_registered_image:
            registered_stack = self._apply_transform_to_image(
                sitk, mv_np, fix_sitk, final_transform, codex_resolution, platform
            )
            img_element = Image2DModel.parse(
                registered_stack,
                transformations={coord_system: Identity()},
                dims=("c", "y", "x"),
            )
            sdata.images[output_image_key] = img_element
            log.append(f"Stored registered image as '{output_image_key}'")

            if output_tiff_path:
                self._save_ome_tiff(registered_stack, output_tiff_path, codex_resolution)
                log.append(f"Saved OME-TIFF to {output_tiff_path}")

        return ModuleResult(
            sdata=sdata,
            metrics={
                "n_landmarks": len(fix_lm),
                "n_iterations": n_iterations,
                "transform_type": "Similarity2D",
            },
            log=log,
        )

    # ── image helpers ────────────────────────────────────────────────

    @staticmethod
    def _image_to_numpy(image_element) -> np.ndarray:
        """Extract numpy array from a SpatialData image element."""
        if hasattr(image_element, "values"):
            # DataTree
            ds = next(iter(image_element.values()))
            arr = ds["image"].values if "image" in ds else next(iter(ds.values())).values
        else:
            arr = np.asarray(image_element)
        return arr

    @staticmethod
    def _resample_image(sitk, image: np.ndarray, src_res: float, dst_res: float) -> np.ndarray:
        """Resample an image to a new resolution.

        Direct port of SPATCH ``resample_image``.
        """
        img = sitk.GetImageFromArray(image.astype(np.float32))
        img.SetSpacing([src_res, src_res])

        original_size = img.GetSize()
        original_spacing = img.GetSpacing()
        new_spacing = [dst_res, dst_res]
        new_size = [
            int(round(osz * ospc / nspc))
            for osz, ospc, nspc in zip(original_size, original_spacing, new_spacing)
        ]

        resampler = sitk.ResampleImageFilter()
        resampler.SetOutputSpacing(new_spacing)
        resampler.SetSize(new_size)
        resampler.SetOutputDirection(img.GetDirection())
        resampler.SetOutputOrigin(img.GetOrigin())
        resampler.SetInterpolator(sitk.sitkLinear)

        return sitk.GetArrayFromImage(resampler.Execute(img))

    # ── core registration ────────────────────────────────────────────

    def _register(
        self,
        sitk,
        fix_img: np.ndarray,
        mv_img: np.ndarray,
        fix_lm_df: pd.DataFrame,
        mv_lm_df: pd.DataFrame,
        fix_res: float,
        codex_res: float,
        n_iter: int,
        platform: str,
    ) -> tuple:
        """Landmark-initialised registration with Mattes MI refinement.

        Direct port of SPATCH ``register``.

        Returns:
            (fix_sitk_image, final_transform)
        """
        # Prepare moving image
        if platform == "nano":
            mv_img = np.flip(mv_img, axis=0).copy()

        mv_landmarks = np.array(
            mv_lm_df[["x", "y"]].values * codex_res, dtype=np.float32
        )
        if platform == "nano":
            mv_landmarks[:, 1] = mv_img.shape[0] * codex_res - mv_landmarks[:, 1]
        mv_landmarks_flat = mv_landmarks.reshape(-1).tolist()

        mv_sitk = sitk.GetImageFromArray(mv_img.astype(np.float32))
        mv_sitk.SetSpacing([codex_res, codex_res])

        # Prepare fixed image (resample to codex resolution)
        fix_lm_flat = np.array(
            fix_lm_df[["x", "y"]].values * codex_res, dtype=np.float32
        ).reshape(-1).tolist()

        fix_resampled = self._resample_image(sitk, fix_img, fix_res, codex_res)
        fix_sitk = sitk.GetImageFromArray(fix_resampled.astype(np.float32))
        fix_sitk.SetSpacing([codex_res, codex_res])

        # Initial transform from landmarks
        initial_transform = sitk.LandmarkBasedTransformInitializer(
            sitk.Similarity2DTransform(),
            fix_lm_flat,
            mv_landmarks_flat,
        )

        # Multi-resolution registration
        reg = sitk.ImageRegistrationMethod()
        reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        reg.SetMetricSamplingStrategy(reg.RANDOM)
        reg.SetMetricSamplingPercentage(0.01)
        reg.SetInterpolator(sitk.sitkLinear)

        reg.SetOptimizerAsGradientDescent(
            learningRate=1,
            numberOfIterations=n_iter,
            convergenceMinimumValue=1e-6,
            convergenceWindowSize=10,
        )
        reg.SetOptimizerScalesFromPhysicalShift()

        reg.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
        reg.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
        reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

        reg.SetInitialTransform(initial_transform, inPlace=False)

        final_transform = reg.Execute(
            sitk.Cast(fix_sitk, sitk.sitkFloat32),
            sitk.Cast(mv_sitk, sitk.sitkFloat32),
        )

        return fix_sitk, final_transform

    # ── apply to coordinates ─────────────────────────────────────────

    @staticmethod
    def _apply_transform_to_coords(
        sitk,
        adata: ad.AnnData,
        coord_key: str,
        codex_res: float,
        transform,
        mv_lm_df: pd.DataFrame,
    ) -> np.ndarray:
        """Apply registration transform to AnnData spatial coordinates.

        Direct port of SPATCH ``apply_transform_to_scanpy``.
        """
        coords = adata.obsm[coord_key] / codex_res  # to pixel space

        inner = transform.GetNthTransform(0) if hasattr(transform, "GetNthTransform") else transform
        matrix = np.array(inner.GetMatrix()).reshape(2, 2)
        translation = np.array(inner.GetTranslation())
        center = np.array([
            np.mean(mv_lm_df["x"]),
            np.mean(mv_lm_df["y"]),
        ])

        transformed = (coords - center) @ matrix + center - translation
        return transformed  # in codex pixel coordinates

    # ── apply to image channels ──────────────────────────────────────

    @staticmethod
    def _apply_transform_to_image(
        sitk,
        mv_stack: np.ndarray,
        fix_sitk,
        transform,
        codex_res: float,
        platform: str,
    ) -> np.ndarray:
        """Resample all channels of the moving image.

        Direct port of SPATCH ``apply_transform_to_codex``.
        """
        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(fix_sitk)
        resampler.SetTransform(transform)
        resampler.SetInterpolator(sitk.sitkLinear)

        # If 2D, wrap in a channel dimension
        if mv_stack.ndim == 2:
            mv_stack = mv_stack[np.newaxis]

        registered = []
        for ch_idx in range(mv_stack.shape[0]):
            ch = mv_stack[ch_idx]
            if platform == "nano":
                ch = np.flip(ch, axis=0).copy()
            ch_sitk = sitk.GetImageFromArray(ch.astype(np.float32))
            ch_sitk.SetSpacing([codex_res, codex_res])
            ch_reg = resampler.Execute(sitk.Cast(ch_sitk, sitk.sitkFloat32))
            registered.append(sitk.GetArrayFromImage(ch_reg))

        stack = np.stack(registered, axis=0)
        return np.round(stack).astype(np.uint16)

    # ── OME-TIFF output ──────────────────────────────────────────────

    @staticmethod
    def _save_ome_tiff(
        stack: np.ndarray,
        path: str,
        resolution: float,
    ) -> None:
        """Save a CYX stack as OME-TIFF."""
        import tifffile
        tifffile.imwrite(
            path,
            stack,
            metadata={
                "Pixels": {
                    "PhysicalSizeX": resolution,
                    "PhysicalSizeXUnit": "µm",
                    "PhysicalSizeY": resolution,
                    "PhysicalSizeYUnit": "µm",
                },
                "axes": "CYX",
            },
            ome=True,
            compression="lzw",
        )
