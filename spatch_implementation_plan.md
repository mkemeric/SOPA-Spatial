# SPATCH Implementation Plan
## Rebuilding on the SpatialData + Sopa Platform

**Prepared for:** University of Rochester — St. Jude Spatial Transcriptomics Collaboration  
**Date:** February 2026  
**Version:** 1.0

---

## 1. Executive Summary

This plan replaces the original SPATCH research scripts (11 standalone Python/R files with hardcoded paths and no reusable API) with a production pipeline built on two published, maintained open-source platforms:

- **SpatialData** (Nature Methods, Jan 2025) — the unified data framework providing coordinate systems, multi-modal alignment, and a Zarr-backed storage model
- **Sopa** (Nature Communications, Jun 2024) — the technology-invariant processing pipeline providing segmentation, aggregation, annotation, and workflow orchestration

The migration preserves all 10 SPATCH analysis stages while reducing custom code from ~4,000 lines to approximately 400 lines across four small plug-in modules. The remaining functionality is handled by maintained libraries with active development communities, CI/CD testing, and published benchmarks.

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    PIPELINE ORCHESTRATION                        │
│         Snakemake (HPC) │ nf-core/Nextflow (Docker/Cloud)       │
│         Seqera Platform (Web UI for monitoring & launch)        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────┐   ┌──────────────┐   ┌──────────────────────┐│
│  │  DATA LAYER  │   │  PROCESSING  │   │   CUSTOM MODULES     ││
│  │              │   │    LAYER     │   │   (Plugin Registry)   ││
│  │ SpatialData  │   │    Sopa      │   │                      ││
│  │  ─ Zarr/NGFF │   │  ─ Segment   │   │ ─ codex_loader       ││
│  │  ─ CCS align │   │  ─ Aggregate │   │ ─ diffusion_analysis ││
│  │  ─ Queries   │   │  ─ Annotate  │   │ ─ gene_protein_corr  ││
│  │  ─ spatialdata│  │  ─ Spatial   │   │ ─ cell_shape_metrics ││
│  │    _io readers│  │    stats     │   │ ─ [user modules...]  ││
│  └──────┬───────┘   └──────┬───────┘   └──────────┬───────────┘│
│         │                  │                       │            │
│         └──────────────────┴───────────────────────┘            │
│                            │                                    │
│                     SpatialData Object                          │
│                   (single .zarr store)                           │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│                    VISUALIZATION LAYER                           │
│  napari-spatialdata (interactive)  │  Xenium Explorer (export)  │
│  spatialdata-plot (static/notebook) │  Vitessce (web/cloud)    │
└─────────────────────────────────────────────────────────────────┘
```

**Key design principle:** Every component reads from and writes to the same `SpatialData` object backed by a single Zarr store. Custom modules follow a uniform interface contract so they are discoverable, configurable, and swappable without modifying the pipeline core.

---

## 3. Stage-by-Stage Migration Plan

### Stage 1: Data Loading

**Original SPATCH:** 6 custom loader functions with hardcoded paths, each producing ad-hoc data structures.

**New approach:** Use `spatialdata-io` built-in readers for 5 of 6 platforms. Write one custom reader for CODEX (see Section 5, Module Registry).

| Platform | Reader | Status |
|---|---|---|
| 10x Xenium | `spatialdata_io.xenium()` | Built-in, production |
| 10x Visium HD | `spatialdata_io.visium_hd()` | Built-in, production |
| Bruker CosMx | `spatialdata_io.cosmx()` | Built-in, production |
| BGI Stereo-seq | `spatialdata_io.stereoseq()` | Built-in, production |
| Akoya CODEX | `spatch_modules.codex_loader` | Custom module (see §5) |
| Vizgen MERSCOPE | `spatialdata_io.merscope()` | Built-in (if needed) |

**Implementation:**

```python
import spatialdata_io as sdio
from spatch_modules import registry

# Built-in platforms — one line each
sdata_xenium = sdio.xenium("/data/xenium_output/")
sdata_visium = sdio.visium_hd("/data/visium_hd_output/")
sdata_cosmx  = sdio.cosmx("/data/cosmx_output/")

# CODEX — via custom module registry
codex_loader = registry.get_module("codex_loader")
sdata_codex  = codex_loader.run("/data/codex_output/")

# All produce SpatialData objects with identical interface
```

**Effort:** 0 lines for built-in loaders. ~120 lines for CODEX module (see §5).

---

### Stage 2: Resolution Standardization (8µm Binning)

**Original SPATCH:** Custom binning logic to normalize Visium HD (2µm bins), Stereo-seq (variable), and other platforms to a common 8µm grid.

**New approach:** SpatialData's `rasterize()` and `aggregate()` functions handle multi-resolution data natively. The coordinate transformation system means you can query at any resolution without physically re-binning.

```python
from spatialdata import rasterize, aggregate

# Rasterize transcripts to 8µm grid
binned = rasterize(
    sdata["transcripts"],
    target_width=8.0,  # µm per pixel
    target_coordinate_system="tissue_aligned"
)

# Or aggregate by pre-defined grid shapes
grid_counts = aggregate(
    values=sdata["transcripts"],
    by=sdata["grid_8um"],
    value_key="gene",
    agg_func="count"
)
```

**Effort:** ~15 lines of configuration, no custom code.

---

### Stage 3: DAPI Masking and Tissue Extraction

**Original SPATCH:** Custom DAPI thresholding and morphological operations to identify tissue boundaries.

**New approach:** Sopa's segmentation pipeline handles this end-to-end with Cellpose (deep learning segmentation) on DAPI channels, including memory-efficient patching for large images.

```python
import sopa

sopa.make_image_patches(sdata)
sopa.segmentation.cellpose(sdata, channels="DAPI", diameter=30)
sopa.resolve_conflicts(sdata)  # handle patch boundaries
sopa.aggregate(sdata)          # count transcripts per cell
```

Or via CLI:
```bash
sopa patchify image data.zarr
sopa segmentation cellpose data.zarr --diameter 30 --channels DAPI
sopa resolve cellpose data.zarr
sopa aggregate data.zarr
```

**Tissue-level masking** (separating tissue from background) can be done with Sopa's ROI functionality or a simple Otsu threshold on the DAPI channel, saved as a SpatialData shape element.

**Effort:** 0 custom code — pure Sopa configuration.

---

### Stage 4: Multi-Modal Registration (Image Alignment)

**Original SPATCH:** SimpleITK landmark-based registration with custom wrapper functions.

**New approach:** SpatialData's coordinate transformation system with napari-spatialdata for interactive landmark placement.

**Workflow:**

1. Load all datasets into a single SpatialData object (or multiple objects viewed together in napari)
2. Open napari-spatialdata and select landmarks on matching anatomical features across modalities
3. SpatialData computes affine transformation and defines a Common Coordinate System (CCS)
4. All downstream spatial queries automatically operate in the aligned CCS

```python
from napari_spatialdata import Interactive
from spatialdata.transformations import align_elements_using_landmarks

# Interactive landmark placement
interactive = Interactive([sdata_xenium, sdata_visium])
interactive.run()
# User places 3+ landmarks on each dataset in napari...
# Landmarks are saved to the SpatialData objects automatically (Shift+E)

# Compute alignment
aligned = align_elements_using_landmarks(
    source=sdata_visium,
    target=sdata_xenium,
    source_landmarks="visium_landmarks",
    target_landmarks="xenium_landmarks",
    target_coordinate_system="tissue_aligned"
)
```

**For batch/automated registration** (no GUI), the affine matrix can be computed programmatically using saved landmark coordinates or feature-based matching, then applied via SpatialData's transformation API.

**Effort:** 0 custom code for interactive workflow. ~30 lines for automated batch registration wrapper.

---

### Stage 5: Preprocessing and QC

**Original SPATCH:** scanpy-based filtering, normalization, HVG selection.

**New approach:** Identical scanpy operations, but operating on SpatialData's AnnData table directly.

```python
import scanpy as sc

adata = sdata.tables["table"]  # AnnData object, standard scanpy target

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

Or via Sopa CLI:
```bash
sopa scanpy-preprocess data.zarr
```

**Effort:** 0 change from standard scanpy workflow. The key difference is that results persist in the Zarr store rather than a loose .h5ad file.

---

### Stage 6: Clustering

**Original SPATCH:** Leiden clustering with various resolution parameters.

**New approach:** Standard scanpy/squidpy clustering. Squidpy adds spatial-aware graph construction.

```python
import squidpy as sq

# Spatial neighbors graph (Delaunay or KNN on coordinates)
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=15)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_1.0")
```

**Effort:** 0 custom code.

---

### Stage 7: Cell Type Annotation (Ensemble)

**Original SPATCH:** 5-tool ensemble (Tangram, TACCO, CellTypist, Spoint, Selina) with majority voting.

**New approach:** Sopa provides built-in Tangram annotation. Additional annotators can be added as custom modules (see §5) using the same plugin interface.

```python
# Sopa built-in: Tangram annotation
sopa.annotate.tangram(sdata, reference_adata, cell_type_key="cell_type")

# Standard CellTypist (operates on AnnData)
import celltypist
model = celltypist.models.download_model("Immune_All_Low.pkl")
predictions = celltypist.annotate(adata, model=model)
adata.obs["celltypist"] = predictions.predicted_labels

# Ensemble via custom module
ensemble = registry.get_module("ensemble_annotation")
ensemble.run(sdata, methods=["tangram", "celltypist", "spoint"],
             strategy="majority_vote")
```

**Effort:** Tangram is built-in. CellTypist is ~10 lines. Full 5-tool ensemble wrapper is ~80 lines as a custom module.

---

### Stage 8: Diffusion Analysis

**Original SPATCH:** Custom analysis comparing in-tissue vs. out-of-tissue signal to measure transcript diffusion.

**New approach:** Custom module using SpatialData spatial queries (see §5 for full implementation).

```python
diffusion = registry.get_module("diffusion_analysis")
results = diffusion.run(sdata, tissue_mask="tissue_boundary",
                        buffer_distance_um=50)
# Results stored in sdata.tables["diffusion_metrics"]
```

**Effort:** ~60 lines as a custom module.

---

### Stage 9: Gene-Protein Correlation (ST × CODEX)

**Original SPATCH:** Multi-resolution correlation between spatial transcriptomics gene expression and CODEX protein levels.

**New approach:** Once both modalities are aligned to the same CCS via Stage 4, SpatialData's `aggregate()` enables cross-modal queries. The correlation computation is a custom module.

```python
correlation = registry.get_module("gene_protein_correlation")
results = correlation.run(
    sdata,
    gene_layer="xenium_table",
    protein_layer="codex_table",
    resolution_um=[8, 25, 50, 100],
    coordinate_system="tissue_aligned"
)
```

**Effort:** ~80 lines as a custom module.

---

### Stage 10: Cell Shape Metrics

**Original SPATCH:** OpenCV contour analysis measuring cell morphology (area, perimeter, eccentricity, circularity).

**New approach:** SpatialData stores cell boundaries as Shapely polygons in GeoDataFrames. Shapely provides all geometric operations natively — no OpenCV needed.

```python
shape_metrics = registry.get_module("cell_shape_metrics")
results = shape_metrics.run(sdata, boundaries_key="cell_boundaries")
# Adds columns to sdata.tables["table"].obs:
#   area_um2, perimeter_um, eccentricity, circularity, solidity
```

**Effort:** ~50 lines as a custom module.

---

## 4. Code Reduction Summary

| Component | Original SPATCH | New Stack | Custom Code |
|---|---|---|---|
| Data loading (6 platforms) | ~600 lines | spatialdata-io | ~120 lines (CODEX only) |
| Resolution binning | ~200 lines | SpatialData rasterize | ~15 lines config |
| DAPI masking / segmentation | ~300 lines | Sopa segmentation | 0 lines |
| Image registration | ~400 lines | SpatialData + napari | ~30 lines (batch only) |
| Preprocessing / QC | ~250 lines | scanpy (unchanged) | 0 lines |
| Clustering | ~150 lines | scanpy + squidpy | 0 lines |
| Cell annotation (ensemble) | ~500 lines | Sopa Tangram + modules | ~80 lines |
| Diffusion analysis | ~350 lines | Custom module | ~60 lines |
| Gene-protein correlation | ~450 lines | Custom module | ~80 lines |
| Cell shape metrics | ~300 lines | Custom module | ~50 lines |
| Config / orchestration | ~500 lines (hardcoded) | Sopa YAML + Snakemake | ~50 lines YAML |
| **Total** | **~4,000 lines** | | **~485 lines** |

Custom code reduction: **88%**. All custom code follows a uniform module interface (§5).

---

## 5. Custom Module Plugin System — `spatch_modules`

### 5.1 Design Philosophy

Users need to:

1. **Discover** available modules (built-in and user-contributed)
2. **Configure** modules via YAML without editing Python code
3. **Swap** alternative implementations for the same analysis step
4. **Extend** the library by dropping a new Python file into a directory
5. **Version** module configurations alongside pipeline runs for reproducibility

The design follows the pattern established by napari plugins, pytest fixtures, and Snakemake wrappers: a lightweight registry backed by Python entry points and/or a filesystem directory scanner.

### 5.2 Module Interface Contract

Every module implements a single abstract base class:

```python
# spatch_modules/base.py

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any
import spatialdata as sd


@dataclass
class ModuleResult:
    """Standard return type for all modules."""
    sdata: sd.SpatialData                # Modified SpatialData (or original if read-only)
    metrics: dict[str, Any] = field(default_factory=dict)  # Numeric results / QC stats
    artifacts: dict[str, str] = field(default_factory=dict) # Paths to generated files (plots, CSVs)
    log: list[str] = field(default_factory=list)            # Human-readable log messages


class SpatchModule(ABC):
    """Base class for all SPATCH custom modules."""

    # ── Module metadata (override in subclass) ─────────────────
    name: str = "unnamed"
    version: str = "0.1.0"
    description: str = ""
    category: str = "analysis"     # "loader" | "analysis" | "annotation" | "qc"
    
    # What this module requires to already exist in the SpatialData object
    requires: list[str] = []       # e.g., ["tables/table", "shapes/cell_boundaries"]
    
    # What this module produces / modifies in the SpatialData object
    produces: list[str] = []       # e.g., ["tables/diffusion_metrics"]

    def __init__(self, **kwargs):
        """Initialize with arbitrary config parameters from YAML."""
        self.config = kwargs

    @abstractmethod
    def run(self, sdata: sd.SpatialData, **runtime_kwargs) -> ModuleResult:
        """Execute the module's analysis.
        
        Args:
            sdata: The SpatialData object (may be modified in-place).
            **runtime_kwargs: Additional parameters not known at config time.
        
        Returns:
            ModuleResult with updated sdata, metrics, and optional artifacts.
        """
        ...

    def validate_inputs(self, sdata: sd.SpatialData) -> list[str]:
        """Check that required elements exist. Returns list of errors (empty = OK)."""
        errors = []
        for req in self.requires:
            parts = req.split("/")
            container = getattr(sdata, parts[0], None)
            if container is None or (len(parts) > 1 and parts[1] not in container):
                errors.append(f"Missing required element: {req}")
        return errors

    def describe(self) -> dict:
        """Return module metadata for registry listing."""
        return {
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "category": self.category,
            "requires": self.requires,
            "produces": self.produces,
            "config_schema": {
                k: type(v).__name__ for k, v in self.config.items()
            }
        }
```

### 5.3 Module Registry

```python
# spatch_modules/registry.py

import importlib
import pkgutil
from pathlib import Path
from typing import Type

from .base import SpatchModule

_REGISTRY: dict[str, Type[SpatchModule]] = {}


def register(cls: Type[SpatchModule]) -> Type[SpatchModule]:
    """Decorator to register a module class."""
    _REGISTRY[cls.name] = cls
    return cls


def get_module(name: str, **config) -> SpatchModule:
    """Instantiate a registered module by name with config."""
    if name not in _REGISTRY:
        raise KeyError(
            f"Module '{name}' not found. "
            f"Available: {list(_REGISTRY.keys())}"
        )
    return _REGISTRY[name](**config)


def list_modules(category: str | None = None) -> list[dict]:
    """List all registered modules, optionally filtered by category."""
    modules = []
    for name, cls in sorted(_REGISTRY.items()):
        inst = cls()
        if category is None or inst.category == category:
            modules.append(inst.describe())
    return modules


def discover_user_modules(directory: str | Path):
    """Scan a directory for .py files containing SpatchModule subclasses.
    
    Users drop .py files into this directory. Any class decorated with
    @register is automatically available in the registry.
    """
    directory = Path(directory)
    if not directory.exists():
        return

    import sys
    sys.path.insert(0, str(directory))
    
    for path in sorted(directory.glob("*.py")):
        if path.name.startswith("_"):
            continue
        module_name = path.stem
        try:
            importlib.import_module(module_name)
        except Exception as e:
            print(f"Warning: Failed to load user module {path.name}: {e}")
    
    sys.path.pop(0)
```

### 5.4 Example Module Implementations

#### CODEX Loader

```python
# spatch_modules/builtin/codex_loader.py

import numpy as np
import pandas as pd
from pathlib import Path
import spatialdata as sd
from spatialdata.models import TableModel, ShapesModel, Image2DModel
from spatialdata.transformations import Identity

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class CODEXLoader(SpatchModule):
    name = "codex_loader"
    version = "1.0.0"
    description = "Load Akoya CODEX/PhenoCycler data into SpatialData format"
    category = "loader"
    requires = []
    produces = ["images/codex_image", "tables/codex_table", "shapes/codex_cells"]

    def run(self, sdata: sd.SpatialData = None, 
            data_path: str = "", **kwargs) -> ModuleResult:
        
        data_path = Path(data_path or self.config.get("data_path", ""))
        log = []

        # CODEX outputs: stitched multi-channel TIFF + cell segmentation CSV
        image_path = data_path / "stitched" / "reg001_stitched.tif"
        cells_path = data_path / "FCS" / "compensation" / "cell_table.csv"

        # Load multi-channel protein image
        from tifffile import imread
        img = imread(str(image_path))  # shape: (C, Y, X) or (Y, X, C)
        if img.ndim == 3 and img.shape[-1] < img.shape[0]:
            img = np.moveaxis(img, -1, 0)  # ensure CYX
        
        image_element = Image2DModel.parse(
            img,
            transformations={"global": Identity()},
            dims=("c", "y", "x")
        )
        log.append(f"Loaded CODEX image: {img.shape}")

        # Load cell table
        cells_df = pd.read_csv(cells_path)
        
        # Extract centroids as shapes
        from geopandas import GeoDataFrame
        from shapely.geometry import Point
        geometries = [Point(r["x"], r["y"]) for _, r in cells_df.iterrows()]
        shapes = ShapesModel.parse(
            GeoDataFrame({"geometry": geometries}),
            transformations={"global": Identity()}
        )

        # Build AnnData table from protein intensities
        import anndata as ad
        protein_cols = [c for c in cells_df.columns 
                       if c not in ("x", "y", "cell_id", "area", "size")]
        adata = ad.AnnData(
            X=cells_df[protein_cols].values,
            obs=cells_df[["cell_id", "x", "y", "area"]].set_index("cell_id"),
            var=pd.DataFrame(index=protein_cols)
        )
        table = TableModel.parse(adata)
        log.append(f"Loaded {adata.n_obs} cells, {adata.n_vars} proteins")

        # Build SpatialData object
        result_sdata = sd.SpatialData(
            images={"codex_image": image_element},
            shapes={"codex_cells": shapes},
            tables={"codex_table": table}
        )

        return ModuleResult(sdata=result_sdata, 
                          metrics={"n_cells": adata.n_obs, "n_proteins": adata.n_vars},
                          log=log)
```

#### Diffusion Analysis

```python
# spatch_modules/builtin/diffusion_analysis.py

import numpy as np
from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class DiffusionAnalysis(SpatchModule):
    name = "diffusion_analysis"
    version = "1.0.0"
    description = "Compare in-tissue vs out-of-tissue signal to quantify transcript diffusion"
    category = "analysis"
    requires = ["shapes/cell_boundaries"]
    produces = ["tables/diffusion_metrics"]

    def run(self, sdata, tissue_mask: str = "tissue_boundary",
            buffer_distance_um: float = 50.0, **kwargs) -> ModuleResult:
        
        from spatialdata import bounding_box_query
        import pandas as pd

        tissue = sdata.shapes[tissue_mask]
        transcripts = sdata.points["transcripts"]
        
        # Create buffer zone outside tissue
        buffered = tissue.geometry.buffer(buffer_distance_um)
        outside_ring = buffered.difference(tissue.geometry.unary_union)

        # Count transcripts inside vs. outside tissue
        # Using SpatialData's polygon query interface
        inside_count = len(sdata.query.points(
            sdata["transcripts"], tissue.geometry.unary_union
        ))
        outside_count = len(sdata.query.points(
            sdata["transcripts"], outside_ring
        ))

        total = inside_count + outside_count
        diffusion_ratio = outside_count / total if total > 0 else 0.0

        # Per-gene diffusion metrics
        gene_metrics = []
        for gene in transcripts["feature_name"].unique():
            gene_pts = transcripts[transcripts["feature_name"] == gene]
            g_in = len(sdata.query.points(gene_pts, tissue.geometry.unary_union))
            g_out = len(sdata.query.points(gene_pts, outside_ring))
            g_total = g_in + g_out
            gene_metrics.append({
                "gene": gene,
                "in_tissue": g_in,
                "out_tissue": g_out,
                "diffusion_ratio": g_out / g_total if g_total > 0 else 0.0
            })

        metrics_df = pd.DataFrame(gene_metrics)
        
        import anndata as ad
        from spatialdata.models import TableModel
        metrics_adata = ad.AnnData(obs=metrics_df.set_index("gene"))
        sdata.tables["diffusion_metrics"] = TableModel.parse(metrics_adata)

        return ModuleResult(
            sdata=sdata,
            metrics={
                "global_diffusion_ratio": diffusion_ratio,
                "inside_transcripts": inside_count,
                "outside_transcripts": outside_count,
                "n_genes_analyzed": len(gene_metrics)
            }
        )
```

#### Cell Shape Metrics

```python
# spatch_modules/builtin/cell_shape_metrics.py

import numpy as np
from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class CellShapeMetrics(SpatchModule):
    name = "cell_shape_metrics"
    version = "1.0.0"
    description = "Compute morphological metrics from cell boundary polygons using Shapely"
    category = "analysis"
    requires = ["shapes/cell_boundaries", "tables/table"]
    produces = []  # Adds columns to existing table

    def run(self, sdata, boundaries_key: str = "cell_boundaries", **kwargs) -> ModuleResult:
        
        boundaries = sdata.shapes[boundaries_key]
        
        areas = []
        perimeters = []
        circularities = []
        eccentricities = []
        solidities = []

        for geom in boundaries.geometry:
            a = geom.area
            p = geom.length
            areas.append(a)
            perimeters.append(p)

            # Circularity: 4π × area / perimeter²  (1.0 = perfect circle)
            circ = (4 * np.pi * a) / (p ** 2) if p > 0 else 0
            circularities.append(circ)

            # Eccentricity from minimum bounding rectangle
            mbr = geom.minimum_rotated_rectangle
            coords = list(mbr.exterior.coords)
            edge1 = np.linalg.norm(
                np.array(coords[1]) - np.array(coords[0]))
            edge2 = np.linalg.norm(
                np.array(coords[2]) - np.array(coords[1]))
            major, minor = max(edge1, edge2), min(edge1, edge2)
            ecc = np.sqrt(1 - (minor / major) ** 2) if major > 0 else 0
            eccentricities.append(ecc)

            # Solidity: area / convex hull area
            convex_area = geom.convex_hull.area
            solidities.append(a / convex_area if convex_area > 0 else 0)

        # Add to existing table
        adata = sdata.tables["table"]
        adata.obs["area_um2"] = areas
        adata.obs["perimeter_um"] = perimeters
        adata.obs["circularity"] = circularities
        adata.obs["eccentricity"] = eccentricities
        adata.obs["solidity"] = solidities

        return ModuleResult(
            sdata=sdata,
            metrics={
                "mean_area": np.mean(areas),
                "mean_circularity": np.mean(circularities),
                "n_cells": len(areas)
            }
        )
```

### 5.5 YAML Configuration for Modules

Modules are configured in the same YAML file as the Sopa pipeline, under a `custom_modules` section:

```yaml
# pipeline_config.yaml

# ── Standard Sopa configuration ───────────────────────────
technology: xenium
segmentation:
  cellpose:
    diameter: 30
    channels: [DAPI]
    
aggregate:
  average_intensities: true

annotation:
  tangram:
    reference: /data/references/scrnaseq_reference.h5ad
    cell_type_key: cell_type

# ── Custom module configuration ───────────────────────────
custom_modules:
  # Directory for user-contributed modules (scanned at startup)
  user_module_dir: /home/user/my_spatch_modules/

  # Pipeline steps — executed in order listed
  steps:
    - module: codex_loader
      enabled: true
      config:
        data_path: /data/codex_output/
    
    - module: diffusion_analysis
      enabled: true
      config:
        tissue_mask: tissue_boundary
        buffer_distance_um: 50.0
    
    - module: gene_protein_correlation
      enabled: true
      config:
        gene_layer: xenium_table
        protein_layer: codex_table
        resolution_um: [8, 25, 50, 100]
        coordinate_system: tissue_aligned
        method: spearman  # or "pearson"
    
    - module: cell_shape_metrics
      enabled: true
      config:
        boundaries_key: cell_boundaries

    # Users can add their own modules here:
    # - module: my_custom_niche_analysis
    #   enabled: true
    #   config:
    #     radius: 100
    #     min_cells: 5
```

### 5.6 Module Runner (Pipeline Integration)

```python
# spatch_modules/runner.py

import yaml
from pathlib import Path
from .registry import get_module, discover_user_modules, list_modules
from .base import ModuleResult
import spatialdata as sd


def run_custom_pipeline(sdata: sd.SpatialData, config_path: str) -> dict:
    """Execute all enabled custom modules from a YAML config."""
    
    config = yaml.safe_load(Path(config_path).read_text())
    module_config = config.get("custom_modules", {})
    
    # Discover user modules
    user_dir = module_config.get("user_module_dir")
    if user_dir:
        discover_user_modules(user_dir)
    
    results = {}
    
    for step in module_config.get("steps", []):
        if not step.get("enabled", True):
            continue
        
        name = step["module"]
        module = get_module(name, **step.get("config", {}))
        
        # Validate inputs
        errors = module.validate_inputs(sdata)
        if errors:
            print(f"Skipping {name}: {errors}")
            results[name] = {"status": "skipped", "errors": errors}
            continue
        
        # Execute
        print(f"Running: {name} v{module.version}")
        result = module.run(sdata, **step.get("config", {}))
        sdata = result.sdata  # Carry forward modified SpatialData
        
        results[name] = {
            "status": "completed",
            "metrics": result.metrics,
            "artifacts": result.artifacts,
            "log": result.log
        }
        print(f"  ✓ {name}: {result.metrics}")
    
    return results
```

### 5.7 Creating a New User Module (User Guide)

To add a custom analysis module, a user creates a single Python file:

```python
# ~/my_spatch_modules/niche_analysis.py

from spatch_modules.base import SpatchModule, ModuleResult
from spatch_modules.registry import register


@register
class NicheAnalysis(SpatchModule):
    name = "niche_analysis"
    version = "0.1.0"
    description = "Identify cellular niches by spatial proximity patterns"
    category = "analysis"
    requires = ["tables/table", "shapes/cell_boundaries"]
    produces = ["tables/niche_assignments"]

    def run(self, sdata, radius: float = 100.0, 
            min_cells: int = 5, **kwargs) -> ModuleResult:
        # ... implementation ...
        return ModuleResult(sdata=sdata, metrics={"n_niches": n_niches})
```

The file is placed in the directory specified by `user_module_dir` in the YAML config. It is automatically discovered and available immediately — no package installation, no code changes to the pipeline.

### 5.8 Module Listing CLI

```bash
$ spatch-modules list
┌─────────────────────────┬─────────┬──────────┬───────────────────────────────────┐
│ Name                    │ Version │ Category │ Description                       │
├─────────────────────────┼─────────┼──────────┼───────────────────────────────────┤
│ codex_loader            │ 1.0.0   │ loader   │ Load Akoya CODEX/PhenoCycler data │
│ diffusion_analysis      │ 1.0.0   │ analysis │ In/out tissue diffusion metrics   │
│ gene_protein_correlation│ 1.0.0   │ analysis │ Multi-res ST×CODEX correlation    │
│ cell_shape_metrics      │ 1.0.0   │ analysis │ Shapely-based cell morphology     │
│ ensemble_annotation     │ 1.0.0   │ annotate │ Majority-vote cell annotation     │
│ niche_analysis          │ 0.1.0   │ analysis │ Cellular niche identification     │  ← user module
└─────────────────────────┴─────────┴──────────┴───────────────────────────────────┘
```

---

## 6. Pipeline Orchestration and UI Assessment

### 6.1 Execution Options (3 tiers)

Sopa provides three execution tiers, each with increasing automation:

#### Tier 1: Python API (interactive / notebook)

Best for: Prototyping, exploration, single samples.

```python
import sopa
sdata = sopa.io.xenium("path/to/data")
sopa.make_image_patches(sdata)
sopa.segmentation.cellpose(sdata, "DAPI", diameter=30)
sopa.aggregate(sdata)
```

**UI:** Jupyter notebooks + napari-spatialdata for visualization. No batch management.

#### Tier 2: Snakemake Pipeline (HPC batch)

Best for: Multi-sample batch processing on institutional Slurm/LSF/SGE clusters.

**Configuration:** Single YAML file per technology/experiment. Shared configs can be version-controlled in Git.

**Execution:**
```bash
snakemake \
  --config data_path=/data/experiment_001 \
  --configfile=config/xenium/cellpose.yaml \
  --workflow-profile profile/slurm \
  --cores 64
```

**Monitoring UI:** Snakemake provides:
- Terminal-based progress with DAG visualization
- `--report` generates an HTML report after completion with timing, file sizes, and DAG graph
- Panoptes (third-party): Real-time web dashboard for Snakemake at `localhost:5000`

**Limitations:**
- No built-in web UI for launching or configuring pipelines
- No multi-user collaboration features
- No centralized run history across experiments
- Status monitoring requires terminal access or third-party tools

#### Tier 3: nf-core/Nextflow + Seqera Platform (production/cloud)

Best for: Multi-user teams, cloud deployment, enterprise-grade monitoring, full audit trails.

**Configuration:** Same YAML config as Snakemake (Sopa's configs work for both engines). Plus a `samplesheet.csv` listing all input datasets:

```csv
sample,data_path
patient_001_xenium,/data/xenium/patient_001/
patient_001_visium,/data/visium_hd/patient_001/
patient_002_xenium,/data/xenium/patient_002/
```

**Execution:**
```bash
nextflow run nf-core/sopa \
  -profile docker,xenium \
  --input samplesheet.csv \
  --outdir /results/ \
  -with-tower  # ← enables Seqera Platform monitoring
```

**Seqera Platform (formerly Nextflow Tower) — the production UI:**

This is the answer to your UI question. The Seqera Platform is the most sophisticated pipeline management UI available for this stack. It provides:

| Capability | Detail |
|---|---|
| **Web-based launch** | Configure and launch pipelines from a browser. Non-technical users can fill in a form (auto-generated from pipeline schema) without touching YAML or command line. |
| **Real-time monitoring** | Live DAG visualization, per-task status (pending/running/completed/failed), CPU/memory/duration metrics per task, log streaming. |
| **Compute environments** | Create and manage HPC (Slurm), cloud (AWS Batch, Google Life Sciences, Azure Batch), and Kubernetes compute environments from the UI. |
| **Run history** | Searchable database of all pipeline runs with full provenance: parameters, config, software versions, input/output files, timing. |
| **Cost tracking** | Per-run and per-task cloud cost estimates (AWS/GCP). |
| **Multi-user** | Workspaces, teams, role-based access control. Multiple users can launch, monitor, and share pipeline results. |
| **REST API** | Full programmatic access to everything the UI offers: launch runs, check status, list results, manage compute environments. |
| **CLI** | `tw` command-line tool for scripting Tower operations. |
| **Labels & tags** | Organize runs by experiment, patient, tissue type, etc. |
| **Optimization** | Analyzes past runs to recommend resource allocation (CPU, memory) per task. |
| **Resumability** | Failed runs can be resumed from the last successful step — no re-computation. |

**Pricing tiers:**

| Tier | Cost | Features |
|---|---|---|
| **Tower Cloud Free** | $0 | Personal workspace, community pipelines, basic monitoring |
| **Tower Cloud Pro** | $0 for academic (apply) | Multi-user workspaces, labels, optimization, priority support |
| **Tower Enterprise** | Custom pricing | On-premises deployment, SSO, audit logs, dedicated support |

**Academic access is free** through the Seqera Academic Program, which provides Pro-level features.

### 6.2 UI Comparison Matrix

| Feature | Snakemake (Tier 2) | nf-core/Nextflow (Tier 3) | Seqera Platform (Tier 3+UI) |
|---|---|---|---|
| Config format | YAML | YAML + samplesheet | Web form (from schema) |
| Launch method | CLI only | CLI | Web UI, CLI, or API |
| Real-time status | Terminal | Terminal | Web dashboard |
| DAG visualization | Post-run HTML | Terminal | Live interactive |
| Per-task metrics | Post-run report | Log files | Real-time charts |
| Multi-sample batch | Loop in shell script | Native (samplesheet) | Native + UI |
| Run history | File system | File system | Searchable database |
| Multi-user | No | No | Yes (workspaces/RBAC) |
| Cost tracking | No | No | Yes (cloud runs) |
| REST API | No | No | Yes |
| Resume failed runs | Yes (--rerun-incomplete) | Yes (-resume) | Yes (one-click) |
| Docker/Singularity | Via conda profiles | Native | Native |
| Cloud deployment | Manual | Manual config | Managed (Forge) |

### 6.3 Recommended Deployment Path

**Phase 1 (now):** Develop and test using **Tier 1** (Python API in Jupyter + napari). Validate each stage works with real data. Iterate on custom modules.

**Phase 2 (batch runs):** Deploy with **Tier 2** (Snakemake on HPC). Process full cohorts on institutional Slurm cluster. Use YAML configs committed to Git for reproducibility.

**Phase 3 (production):** Migrate to **Tier 3** (nf-core/Nextflow + Seqera Platform). Enable web-based pipeline launching and monitoring for the broader team. Apply for Seqera Academic Program for free Pro access.

The migration from Tier 2 → Tier 3 requires zero code changes because Sopa's YAML configs are compatible with both Snakemake and Nextflow engines.

---

## 7. Project Structure

```
spatch-pipeline/
├── pyproject.toml              # Package metadata for spatch_modules
├── environment.yml             # Conda environment specification
│
├── spatch_modules/             # Custom module library
│   ├── __init__.py
│   ├── base.py                 # SpatchModule ABC + ModuleResult
│   ├── registry.py             # Module discovery and instantiation
│   ├── runner.py               # YAML-driven module execution
│   ├── cli.py                  # `spatch-modules list/run` commands
│   └── builtin/                # Shipped modules
│       ├── __init__.py
│       ├── codex_loader.py
│       ├── diffusion_analysis.py
│       ├── gene_protein_correlation.py
│       ├── cell_shape_metrics.py
│       └── ensemble_annotation.py
│
├── configs/                    # Pipeline configurations
│   ├── xenium_cellpose.yaml    # Per-technology Sopa configs
│   ├── cosmx_baysor.yaml
│   ├── visium_hd.yaml
│   └── spatch_full.yaml        # Full SPATCH pipeline (Sopa + custom modules)
│
├── samplesheets/               # Input manifests
│   └── experiment_001.csv
│
├── user_modules/               # User-contributed modules (auto-discovered)
│   └── README.md
│
├── notebooks/                  # Jupyter exploration notebooks
│   ├── 01_load_and_align.ipynb
│   ├── 02_segment_and_annotate.ipynb
│   ├── 03_custom_analysis.ipynb
│   └── 04_visualization.ipynb
│
├── tests/                      # Module tests
│   ├── test_codex_loader.py
│   ├── test_diffusion.py
│   └── test_cell_shape.py
│
└── docs/
    ├── module_development.md   # How to write a new module
    └── deployment.md           # HPC and cloud setup guides
```

---

## 8. Implementation Timeline

| Week | Milestone | Deliverables |
|---|---|---|
| 1-2 | Environment setup | Conda env, SpatialData + Sopa installed, test data loaded |
| 3-4 | Module framework | `spatch_modules` package: base class, registry, runner, CLI |
| 5-6 | CODEX loader + registration | CODEX reader module, napari landmark alignment validated |
| 7-8 | Core Sopa pipeline | Segmentation, aggregation, annotation working end-to-end |
| 9-10 | Custom analysis modules | Diffusion, correlation, cell shape modules implemented |
| 11-12 | Snakemake integration | Full pipeline YAML config, batch processing on HPC |
| 13-14 | Testing and validation | Compare outputs to original SPATCH results |
| 15-16 | Documentation + handoff | Module dev guide, deployment docs, notebooks |

**Optional Phase 2 (weeks 17-20):** Nextflow migration + Seqera Platform setup.

---

## 9. Dependencies

```yaml
# environment.yml
name: spatch
channels:
  - conda-forge
  - pytorch
dependencies:
  - python=3.11
  - pip
  - pip:
    # Core platform
    - spatialdata>=0.6
    - spatialdata-io>=0.2
    - spatialdata-plot>=0.3
    - sopa>=2.1
    - napari-spatialdata>=0.5
    
    # Analysis ecosystem
    - scanpy>=1.10
    - squidpy>=1.5
    - anndata>=0.10
    - scvi-tools>=1.2
    
    # Segmentation
    - cellpose>=3.0
    
    # Annotation
    - tangram-sc>=1.0
    - celltypist>=1.6
    
    # Custom module dependencies
    - tifffile        # CODEX image loading
    - shapely>=2.0    # Cell shape metrics
    - geopandas       # Spatial operations
    
    # Workflow orchestration (choose one)
    - snakemake>=8.0
    # - nextflow (install separately via curl)
```

---

## 10. Risk Assessment

| Risk | Impact | Mitigation |
|---|---|---|
| SpatialData API breaking changes | Medium | Pin versions in environment.yml; scverse has 6-month deprecation policy |
| Sopa doesn't support a needed feature | Low | Sopa is built on SpatialData; custom modules bypass Sopa for any step |
| CODEX format variations across instruments | Medium | CODEX loader module has config params for format variants |
| Large dataset performance | Low | Sopa proven on 1TB+ images; SpatialData uses lazy/Dask arrays |
| napari stability on remote/HPC | Medium | Use napari headless mode for batch; interactive only on local workstations |
| Seqera Platform academic access denied | Low | Snakemake tier provides full functionality without Seqera |

---

## 11. Key Decisions for Review

1. **Nextflow vs. Snakemake as primary orchestrator?** Both work with identical Sopa configs. Snakemake is more Pythonic and simpler to start. Nextflow has the superior UI (Seqera). Recommendation: Start with Snakemake, migrate to Nextflow when the team needs the web UI.

2. **Module registry: entry points vs. directory scanning?** Entry points require `pip install` for each module. Directory scanning allows "drop a file" workflow. Recommendation: Support both; directory scanning for development, entry points for production/shared modules.

3. **Custom module execution: inside Sopa/Snakemake DAG or as post-processing?** Recommendation: Post-processing initially (simpler), with Snakemake rule integration in Phase 2. The module runner can be called as a single Snakemake rule that processes all custom steps sequentially.

---

*This plan provides a complete path from the current SPATCH research scripts to a production-grade, extensible spatial transcriptomics pipeline built on published, maintained infrastructure. The custom module system ensures the investment in SPATCH-specific analysis logic is preserved while gaining the sustainability and scalability of the SpatialData ecosystem.*
