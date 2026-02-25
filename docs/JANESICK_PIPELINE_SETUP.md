# Janesick Breast Cancer Pipeline Setup

This document describes how the SPATCH pipeline was configured to analyze the Janesick et al. (Nature Communications, 2023) FFPE human breast cancer dataset using the SpatialData and Sopa frameworks.

## Table of Contents
1. [Dataset Overview](#dataset-overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Configuration Details](#configuration-details)
4. [Running the Pipeline](#running-the-pipeline)
5. [Expected Outputs](#expected-outputs)
6. [Troubleshooting](#troubleshooting)

## Dataset Overview

### Sample Information
- **Publication**: Janesick A, et al. (2023) "High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis." Nature Communications 14:8353
- **Technology**: 10x Genomics Xenium In Situ
- **Sample**: FFPE Human Breast Cancer Replicate 1
- **Pathology**: Invasive Ductal Carcinoma with Ductal Carcinoma In Situ (DCIS)
- **Data Location**: `../data/outs/`
- **Source**: https://github.com/10XGenomics/janesick_nature_comms_2023_companion

### Dataset Contents
```
../data/outs/
├── cells.parquet                    # Cell segmentation results (~70K cells)
├── cell_boundaries.parquet          # Cell boundary polygons
├── cell_feature_matrix.h5           # Gene x Cell expression matrix (313 genes)
├── transcripts.parquet              # Individual transcript locations (~100M transcripts)
├── morphology_mip.ome.tif          # DAPI maximum intensity projection
├── gene_panel.json                  # Gene panel metadata
└── experiment.xenium                # Experiment metadata
```

### Data Statistics
- **Cells**: ~70,000 cells
- **Genes**: 313 targeted genes (Xenium Human Breast Cancer Panel)
- **Transcripts**: ~100 million individual transcripts
- **Resolution**: Subcellular resolution (~0.2 µm/pixel)
- **Tissue Area**: ~11 mm²

## Pipeline Architecture

The SPATCH pipeline integrates three main frameworks:

### 1. Sopa + SpatialData (stock from PyPI)
- **Sopa** (Blampey et al., Nature Communications 2024): Config-driven pipeline for segmentation, aggregation, annotation, preprocessing, and Xenium Explorer export. Includes a Snakemake workflow and CLI.
- **SpatialData** (Marconato et al., Nature Methods 2025): Unified Zarr-backed data model for spatial omics — coordinates, images, tables, and shapes in one object.
- **SpatialData-IO**: Built-in readers for Xenium, MERSCOPE, CosMx, Visium HD, etc.

### 2. SPATCH Custom Modules (`spatch_modules`)
- **Purpose**: Domain-specific analyses not covered by stock sopa
- **Janesick modules**: `diffusion_analysis` (transcript diffusion quantification) and `cell_shape_metrics` (morphological descriptors)
- **Architecture**: Lightweight plugin system — modules are auto-discovered and driven by YAML config

## Configuration Details

The pipeline configuration is defined in `configs/janesick_breast_cancer.yaml`.

### Key Configuration Decisions

#### 1. Segmentation Strategy
```yaml
segmentation:
  use_existing: true
```

**Rationale**: 
- Xenium provides high-quality segmentation from Xenium Ranger
- Re-segmentation not typically needed unless testing alternative methods
- Saves ~2-3 hours of compute time

**Alternative** (if re-segmentation desired):
```yaml
segmentation:
  method: cellpose
  cellpose:
    diameter: 35  # Breast cancer cells ~30-40 pixels
    channels: [DAPI]
    model_type: cyto2  # General cell model
    use_gpu: true
```

#### 2. Quality Control Parameters
```yaml
preprocessing:
  min_genes_per_cell: 50      # Relaxed for tumor heterogeneity
  max_genes_per_cell: 10000   # Higher for transcriptionally active cells
  max_pct_mt: 25              # Permissive for FFPE tissue
```

**Rationale**:
- **Lower gene threshold**: Tumor cells can be transcriptionally suppressed
- **Higher upper limit**: Some tumor cells are highly proliferative
- **Mitochondrial %**: FFPE tissue degradation can increase MT reads

#### 3. Dimensionality Reduction
```yaml
preprocessing:
  n_top_genes: 3000           # More genes for TME complexity
  n_pcs: 50
  n_neighbors: 30             # Dense Xenium data benefits from more neighbors
  umap_min_dist: 0.3
```

**Rationale**:
- **3000 HVGs**: Captures tumor microenvironment (TME) heterogeneity
- **30 neighbors**: Xenium's high cell density supports larger neighborhoods
- **UMAP min_dist**: Balanced local vs global structure preservation

#### 4. Multi-Resolution Clustering
```yaml
clustering:
  leiden:
    resolutions: [0.3, 0.5, 0.8, 1.0, 1.5]
```

**Rationale**:
- **Low resolution (0.3)**: Major cell type compartments (tumor, immune, stroma)
- **Medium resolution (0.8)**: Cell type refinement (invasive vs DCIS, T cell subtypes)
- **High resolution (1.5)**: Fine-grained tumor heterogeneity

#### 5. Cell Type Markers
The configuration includes breast cancer-specific marker genes:

| Cell Type | Markers | Biological Significance |
|-----------|---------|------------------------|
| Invasive Tumor | FASN, CEACAM6, MKI67 | Lipid metabolism, proliferation |
| DCIS | CEACAM6, ESR1 | Precursor lesion markers |
| Myoepithelial | ACTA2, KRT15, KRT14 | Tumor boundary cells |
| Stromal | POSTN, COL1A1, DCN | Cancer-associated fibroblasts |
| T cells | CD3E, CD3D, IL7R | Adaptive immunity |
| Macrophages | ITGAX, CD68, LYZ | Tumor-associated macrophages |
| B cells | MS4A1, CD79A | Humoral immunity |
| Endothelial | VWF, PECAM1 | Vasculature |

#### 6. Spatial Analysis Parameters
```yaml
spatial_analysis:
  neighborhood:
    enabled: true  # Cell neighborhood enrichment
  interactions:
    enabled: true  # Ligand-receptor interactions
  ripley:
    enabled: true  # Spatial point pattern analysis
```

**Rationale**:
- **Neighborhood enrichment**: Identifies which cell types co-localize
- **Interactions**: Predicts cell-cell communication pathways
- **Ripley's statistics**: Tests for clustering vs random distribution

#### 7. SPATCH Custom Modules

##### Cell Shape Metrics
```yaml
- module: cell_shape_metrics
  enabled: true
  config:
    metrics:
      - area
      - perimeter
      - circularity
      - eccentricity
      - solidity
      - aspect_ratio
```

**Biological Rationale**:
- Tumor cells often have irregular morphology vs normal epithelium
- Invasive cells may show different shape characteristics than DCIS
- Morphology correlates with cell state and function

##### Diffusion Analysis
```yaml
- module: diffusion_analysis
  enabled: true
  config:
    buffer_distance_um: 50.0
    compute_distances: true
```

**Biological Rationale**:
- Quantifies transcript leakage from cells
- Important for understanding DCIS boundary dynamics
- Validates quality of transcript assignment

## Running the Pipeline

The Janesick pipeline runs in two phases:
1. **Phase 1 (sopa)**: Data loading, segmentation, aggregation, preprocessing — handled entirely by stock sopa.
2. **Phase 2 (SPATCH modules)**: Custom analysis extensions (diffusion analysis, cell shape metrics).

### Phase 1: Sopa Standard Pipeline

#### Option A: Sopa CLI (step by step)

```bash
# 1. Convert Xenium output to SpatialData
sopa convert ../data/outs/ --sdata-path results/janesick.zarr --technology xenium

# 2. Patchify for parallel segmentation
sopa patchify image results/janesick.zarr --patch-width-pixel 2048 --patch-overlap-pixel 100

# 3. Segment with Cellpose
sopa segmentation cellpose results/janesick.zarr --diameter 35 --channels DAPI --model-type cyto2

# 4. Resolve boundary conflicts
sopa resolve cellpose results/janesick.zarr

# 5. Aggregate transcripts per cell
sopa aggregate results/janesick.zarr --min-transcripts 10

# 6. Export to Xenium Explorer
sopa explorer write results/janesick.zarr --output-path results/janesick_explorer/
```

#### Option B: Sopa Snakemake (fully automated)

Uses `configs/janesick_sopa.yaml` which defines all parameters in sopa's native config format:

```bash
snakemake --snakefile sopa/workflow/Snakefile \
  --configfile configs/janesick_sopa.yaml \
  --config data_path=../data/outs/ sdata_path=results/janesick.zarr \
  --cores 4
```

This runs the entire standard pipeline end-to-end and produces:
- SpatialData Zarr store at `results/janesick.zarr`
- Xenium Explorer output for interactive visualization
- QC report

### Phase 2: SPATCH Custom Modules

After sopa completes, run the SPATCH-specific analyses. These are configured in `configs/janesick_breast_cancer.yaml` under the `custom_modules` section.

#### From Python / Notebook

```python
import spatialdata as sd
from spatch_modules.runner import run_custom_pipeline

# Load the sopa-processed data
sdata = sd.read_zarr("results/janesick.zarr")

# Run SPATCH custom modules
results = run_custom_pipeline(sdata, "configs/janesick_breast_cancer.yaml", verbose=True)

# Save updated data
sdata.write("results/janesick_breast_cancer/processed.zarr")
```

#### What the custom modules do

**diffusion_analysis** — Compares in-tissue vs out-of-tissue transcript signal:
- Splits cells by tissue membership (`in_tissue` column)
- Computes per-gene diffusion ratios (outside / total counts)
- Calculates minimum distances from out-of-tissue spots to tissue edge
- Stores results as `diffusion_metrics` table in sdata
- Key for understanding DCIS boundary dynamics in the Janesick sample

**cell_shape_metrics** — Computes morphological descriptors from cell boundaries:
- Area, perimeter, circularity (isoperimetric quotient)
- Eccentricity (from minimum bounding rectangle)
- Solidity (area / convex hull area), aspect ratio
- Adds metrics as columns to the cell table
- Useful for distinguishing tumor from normal epithelium based on shape

#### From CLI

```bash
spatch-modules run \
  -i results/janesick.zarr \
  -c configs/janesick_breast_cancer.yaml \
  -o results/janesick_breast_cancer/processed.zarr
```

### Interactive Exploration (Notebook)

Open `notebooks/02_janesick_breast_cancer_analysis.ipynb` for interactive analysis including:
- Detailed scanpy preprocessing (multi-resolution clustering, marker genes)
- Spatial visualization of cell types and markers
- SPATCH custom module results exploration
- Publication-quality figure generation

## Expected Outputs

After running the pipeline, the following outputs will be generated:

### Directory Structure
```
results/janesick_breast_cancer/
├── processed.zarr/              # Unified SpatialData object
│   ├── images/                  # Morphology images
│   ├── shapes/                  # Cell boundaries
│   ├── points/                  # Transcript coordinates
│   └── tables/                  # AnnData tables
├── processed_table.h5ad         # Standalone AnnData for scanpy/squidpy
├── marker_genes.csv             # Differential expression results
└── figures/
    ├── qc_metrics.png
    ├── spatial_clusters.png
    ├── spatial_markers.png
    ├── cell_shape_distributions.png
    ├── neighborhood_enrichment.png
    ├── cooccurrence.png
    └── marker_heatmap.png
```

### File Descriptions

#### 1. `processed.zarr/`
- **Format**: Zarr store (SpatialData)
- **Size**: ~5-10 GB
- **Contents**: All spatial data in unified coordinate system
- **Usage**: Load with `sdata = sd.read_zarr("processed.zarr")`

#### 2. `processed_table.h5ad`
- **Format**: AnnData H5AD
- **Size**: ~500 MB
- **Contents**: 
  - `.X`: Normalized gene expression (log1p)
  - `.raw`: Raw normalized counts
  - `.obs`: Cell metadata (QC, clusters, spatial coords, shape metrics)
  - `.var`: Gene metadata
  - `.obsm`: PCA, UMAP embeddings
  - `.obsp`: Neighbor graphs (transcriptional, spatial)
  - `.uns`: Analysis parameters, marker genes
- **Usage**: Load with `adata = sc.read_h5ad("processed_table.h5ad")`

#### 3. `marker_genes.csv`
- **Format**: CSV
- **Contents**: Differential expression results for all clusters
- **Columns**: gene, cluster, logfoldchange, pval, pval_adj

#### 4. Figures
All figures saved at 300 DPI for publication quality.

## Troubleshooting

### Import errors

**`ModuleNotFoundError: No module named 'spatch_modules'`**: Run `pip install -e . --no-deps` from the project root.

**`ModuleNotFoundError: No module named 'sopa'`**: Run `./setup_environment.sh`.

### Segmentation boundaries not found

`KeyError: 'cell_boundaries'` — Check available shapes: `print(sdata.shapes.keys())`. The Xenium reader may use a different key; check with `sdata` to see the full structure.

### Missing marker genes

The Xenium Human Breast Cancer Panel has 313 genes. Verify availability:
```python
print(adata.var_names)
```

### Memory errors

- Sopa's patching (`patch_width_pixel: 2048`) handles large images automatically
- For downstream analysis, reduce `n_top_genes` or `n_neighbors` in the config
- Use `sc.pp.subsample(adata, n_obs=10000)` for quick exploration

### GPU acceleration

Cellpose supports GPU: add `gpu: true` to the cellpose section of `janesick_sopa.yaml`.

## Next Steps

After completing the pipeline setup:
1. Review [Results Exploration Guide](JANESICK_RESULTS_EXPLORATION.md)
2. Customize analysis for specific biological questions
3. Integrate with other modalities (H&E, IF imaging)
4. Compare with original Janesick publication analyses

## References

1. **Janesick A, et al.** (2023) High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. Nature Communications 14:8353

2. **Marconato L, et al.** (2025) SpatialData: an open and universal framework for processing spatial omics data. Nature Methods

3. **Blampey Q, et al.** (2024) Sopa: a technology-invariant pipeline for analyses of image-based spatial omics. Nature Communications

4. **Wolf FA, et al.** (2018) SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19:15

5. **Palla G, et al.** (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19:171-178
