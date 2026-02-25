# Janesick Breast Cancer Results Exploration Guide

This document describes how to explore and interpret the results from the SPATCH analysis pipeline applied to the Janesick breast cancer dataset.

## Table of Contents
1. [Loading Results](#loading-results)
2. [Quality Assessment](#quality-assessment)
3. [Cell Type Analysis](#cell-type-analysis)
4. [Spatial Analysis](#spatial-analysis)
5. [Morphological Analysis](#morphological-analysis)
6. [Advanced Exploration](#advanced-exploration)
7. [Interpretation Guidelines](#interpretation-guidelines)
8. [Export and Sharing](#export-and-sharing)

## Loading Results

### Python/Jupyter Environment

```python
import spatialdata as sd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns

# Load complete SpatialData object
sdata = sd.read_zarr("results/janesick_breast_cancer/processed.zarr")

# Load AnnData table for analysis
adata = sc.read_h5ad("results/janesick_breast_cancer/processed_table.h5ad")

print(f"Loaded {adata.n_obs:,} cells and {adata.n_vars:,} genes")
```

### Quick Data Inspection

```python
# View SpatialData structure
print(sdata)

# Check available data layers
print("Images:", list(sdata.images.keys()))
print("Shapes:", list(sdata.shapes.keys()))
print("Points:", list(sdata.points.keys()))
print("Tables:", list(sdata.tables.keys()))

# Check AnnData structure
print("\nAnnData structure:")
print("  .X: Expression matrix")
print("  .obs columns:", adata.obs.columns.tolist())
print("  .obsm keys:", adata.obsm.keys())
print("  .obsp keys:", adata.obsp.keys())
```

## Quality Assessment

### 1. Cell Quality Metrics

```python
# Distribution of key QC metrics
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Total counts
adata.obs['total_counts'].hist(bins=100, ax=axes[0])
axes[0].set_xlabel('Total Counts')
axes[0].set_title(f'Median: {adata.obs["total_counts"].median():.0f}')

# Genes detected
adata.obs['n_genes_by_counts'].hist(bins=100, ax=axes[1])
axes[1].set_xlabel('Genes Detected')
axes[1].set_title(f'Median: {adata.obs["n_genes_by_counts"].median():.0f}')

# Mitochondrial percentage
adata.obs['pct_counts_mt'].hist(bins=100, ax=axes[2])
axes[2].set_xlabel('MT %')
axes[2].set_title(f'Median: {adata.obs["pct_counts_mt"].median():.1f}%')

plt.tight_layout()
```

**Interpretation**:
- **Total counts**: Typical range 50-500 transcripts/cell for FFPE Xenium
- **Genes detected**: Expect 20-150 genes with 313-gene panel
- **MT%**: FFPE should be <25%; high values indicate RNA degradation

### 2. Spatial QC

```python
# Visualize cell density across tissue
sq.pl.spatial_scatter(
    adata, 
    color='total_counts',
    size=2,
    figsize=(12, 12),
    cmap='viridis',
    title='Transcript Counts (Spatial)'
)
```

**Look for**:
- Uniform distribution of high-quality cells
- Areas of low signal (necrosis, processing artifacts)
- Edge effects (drop in quality at tissue boundaries)

## Cell Type Analysis

### 1. Clustering Overview

```python
# Cluster sizes at different resolutions
for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    col = f'leiden_{res}'
    if col in adata.obs.columns:
        n_clusters = adata.obs[col].nunique()
        print(f"Resolution {res}: {n_clusters} clusters")
        print(adata.obs[col].value_counts().head())
        print()
```

### 2. UMAP Visualization

```python
# Multi-panel UMAP
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Cluster identity
sc.pl.umap(adata, color='leiden_0.8', ax=axes[0, 0], show=False, 
           legend_loc='on data', legend_fontsize=8)

# Key markers
markers = ['FASN', 'ACTA2', 'POSTN', 'CD3E', 'ITGAX']
for i, gene in enumerate(markers):
    if gene in adata.var_names:
        sc.pl.umap(adata, color=gene, ax=axes.flatten()[i+1], 
                   show=False, cmap='viridis')

plt.tight_layout()
```

### 3. Spatial Cell Type Distribution

```python
# Major cell types on tissue
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()

cell_type_scores = [col for col in adata.obs.columns if col.endswith('_score')]

for i, score_col in enumerate(cell_type_scores[:8]):
    sq.pl.spatial_scatter(
        adata,
        color=score_col,
        ax=axes[i],
        size=1,
        cmap='viridis',
        title=score_col.replace('_score', '')
    )

plt.tight_layout()
```

**Interpretation**:
- **Tumor regions**: High FASN, CEACAM6 scores
- **Myoepithelial boundary**: ACTA2+ cells surrounding DCIS
- **Stromal regions**: POSTN+ cancer-associated fibroblasts
- **Immune infiltrates**: CD3E (T cells), ITGAX (macrophages)

### 4. Marker Gene Analysis

```python
# Top markers per cluster
import pandas as pd

marker_df = pd.read_csv("results/janesick_breast_cancer/marker_genes.csv")

# Top 5 genes per cluster
for cluster in adata.obs['leiden_0.8'].unique()[:5]:
    print(f"\n=== Cluster {cluster} ===")
    cluster_markers = marker_df[marker_df['cluster'] == cluster].head(5)
    for _, row in cluster_markers.iterrows():
        print(f"  {row['gene']}: log2FC={row['logfoldchange']:.2f}, "
              f"p={row['pval_adj']:.2e}")
```

### 5. Cell Type Annotation

**Manual annotation based on marker expression**:

```python
# Example annotation logic
def annotate_cell_type(row):
    # Tumor cells
    if row['Tumor_Invasive_score'] > 0.2:
        return 'Tumor_Invasive'
    elif row['Tumor_DCIS_score'] > 0.2:
        return 'Tumor_DCIS'
    # Myoepithelial
    elif row['Myoepithelial_score'] > 0.15:
        return 'Myoepithelial'
    # Immune cells
    elif row['T_cells_score'] > 0.15:
        return 'T_cells'
    elif row['Macrophages_score'] > 0.15:
        return 'Macrophages'
    elif row['B_cells_score'] > 0.1:
        return 'B_cells'
    # Stromal/other
    elif row['Stromal_score'] > 0.15:
        return 'Stromal'
    elif row['Endothelial_score'] > 0.1:
        return 'Endothelial'
    else:
        return 'Unknown'

adata.obs['cell_type'] = adata.obs.apply(annotate_cell_type, axis=1)

# Visualize
sc.pl.umap(adata, color='cell_type', legend_loc='on data')
sq.pl.spatial_scatter(adata, color='cell_type', size=2, figsize=(15, 15))
```

## Spatial Analysis

### 1. Neighborhood Enrichment

```python
# Compute and visualize neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')
sq.pl.nhood_enrichment(adata, cluster_key='cell_type', figsize=(10, 10))
```

**Interpretation**:
- **Red (enriched)**: Cell types that co-localize more than expected
- **Blue (depleted)**: Cell types that avoid each other
- **Expected patterns**:
  - Tumor cells enrich with tumor cells
  - T cells enrich near tumor boundaries
  - Myoepithelial cells enrich with DCIS

### 2. Co-occurrence Analysis

```python
# Test spatial co-occurrence
sq.gr.co_occurrence(adata, cluster_key='cell_type')

# Visualize for specific cell types
sq.pl.co_occurrence(
    adata,
    cluster_key='cell_type',
    clusters=['Tumor_Invasive', 'Tumor_DCIS', 'Myoepithelial', 'T_cells'],
    figsize=(12, 5)
)
```

**Interpretation**:
- Identifies cell types enriched/depleted at specific distances
- Useful for finding immune-tumor interactions
- Can reveal distinct tumor microenvironment zones

### 3. Spatial Patterns

```python
# Ripley's statistics for cell types
cell_types_to_test = ['Tumor_Invasive', 'T_cells', 'Macrophages']

for ct in cell_types_to_test:
    sq.gr.ripley(adata, cluster_key='cell_type', mode='L', 
                 spatial_key='spatial')
    
fig, axes = plt.subplots(1, len(cell_types_to_test), figsize=(15, 5))
for i, ct in enumerate(cell_types_to_test):
    sq.pl.ripley(adata, cluster_key='cell_type', mode='L', 
                 ax=axes[i], title=ct)
```

**Interpretation**:
- Above CSR line: Clustering (cells aggregate)
- Below CSR line: Dispersion (cells avoid each other)
- On CSR line: Random distribution

### 4. Ligand-Receptor Interactions

```python
# Cell-cell interaction analysis
sq.gr.ligrec(
    adata,
    n_perms=100,
    cluster_key='cell_type',
    copy=False
)

# Visualize top interactions
sq.pl.ligrec(
    adata,
    source_groups=['Tumor_Invasive', 'Tumor_DCIS'],
    target_groups=['T_cells', 'Macrophages', 'Stromal'],
    means_range=(0.3, np.inf),
    alpha=0.001,
    swap_axes=True
)
```

**Biological Insights**:
- Immune checkpoint interactions (PD-L1/PD-1)
- Growth factor signaling (VEGF, TGF-β)
- Cell adhesion molecules

## Morphological Analysis

### 1. Cell Shape Distributions

```python
# Overall shape metrics
shape_cols = ['area_um2', 'circularity', 'eccentricity', 'solidity']

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for i, col in enumerate(shape_cols):
    if col in adata.obs.columns:
        adata.obs[col].hist(bins=100, ax=axes[i], edgecolor='black', alpha=0.7)
        axes[i].set_xlabel(col)
        axes[i].set_ylabel('Frequency')
        axes[i].axvline(adata.obs[col].median(), color='red', 
                       linestyle='--', label='Median')
        axes[i].legend()

plt.tight_layout()
```

**Expected values**:
- **Area**: 50-300 µm² (tumor cells often larger)
- **Circularity**: 0.6-0.9 (1.0 = perfect circle)
- **Eccentricity**: 0-0.8 (0 = circle, 1 = line)
- **Solidity**: 0.85-1.0 (convexity measure)

### 2. Shape by Cell Type

```python
# Compare morphology across cell types
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for i, metric in enumerate(shape_cols):
    if metric in adata.obs.columns:
        ax = axes.flatten()[i]
        
        # Boxplot
        adata.obs.boxplot(
            column=metric,
            by='cell_type',
            ax=ax,
            rot=45
        )
        ax.set_title(f'{metric} by Cell Type')
        ax.set_xlabel('')

plt.tight_layout()
```

**Biological Interpretation**:
- **Large, irregular cells**: Likely malignant (invasive tumor)
- **Elongated cells**: Myoepithelial, fibroblasts
- **Small, circular**: Lymphocytes (T/B cells)
- **Medium, irregular**: Macrophages

### 3. Shape-Phenotype Correlations

```python
# Correlate morphology with expression
marker_genes = ['FASN', 'ACTA2', 'CD3E', 'POSTN']

for gene in marker_genes:
    if gene in adata.var_names:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Area vs expression
        axes[0].scatter(adata.obs['area_um2'], adata[:, gene].X.toarray().flatten(),
                       alpha=0.1, s=1)
        axes[0].set_xlabel('Cell Area (µm²)')
        axes[0].set_ylabel(f'{gene} Expression')
        axes[0].set_title(f'Area vs {gene}')
        
        # Circularity vs expression
        axes[1].scatter(adata.obs['circularity'], adata[:, gene].X.toarray().flatten(),
                       alpha=0.1, s=1)
        axes[1].set_xlabel('Circularity')
        axes[1].set_ylabel(f'{gene} Expression')
        axes[1].set_title(f'Circularity vs {gene}')
        
        plt.tight_layout()
        plt.show()
```

## Advanced Exploration

### 1. Spatial Domains/Niches

```python
# Identify cellular neighborhoods using clustering
from sklearn.cluster import KMeans

# Get spatial neighborhood composition
nhood_composition = pd.DataFrame(
    sq.gr.spatial_neighbors(adata, n_neighs=30, coord_type='generic', 
                            key_added='spatial_neighbors')
)

# Cluster based on neighborhood composition
# (This is a simplified example - see squidpy docs for full implementation)
```

### 2. Gene Expression Gradients

```python
# Identify genes with spatial gradients
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100)

# Plot top spatially variable genes
spatially_variable = adata.uns['moranI'].sort_values('I', ascending=False).head(20)
print(spatially_variable)

# Visualize
top_sv_genes = spatially_variable.index[:6].tolist()
sq.pl.spatial_scatter(adata, color=top_sv_genes, ncols=3, size=2, figsize=(15, 10))
```

### 3. Trajectory Analysis

For understanding tumor progression (DCIS → Invasive):

```python
import scanpy as sc

# Subset to epithelial compartment
epithelial = adata[adata.obs['cell_type'].isin(['Tumor_DCIS', 'Tumor_Invasive'])].copy()

# Run diffusion pseudotime
sc.tl.diffmap(epithelial)
sc.tl.dpt(epithelial, n_dcs=10)

# Visualize
sc.pl.diffmap(epithelial, color=['cell_type', 'dpt_pseudotime'])
sq.pl.spatial_scatter(epithelial, color='dpt_pseudotime', size=3, cmap='viridis')
```

## Interpretation Guidelines

### Key Biological Questions

#### 1. Where is the DCIS-to-Invasive Transition?

**Approach**:
- Plot DCIS (CEACAM6+) and invasive (FASN+) markers spatially
- Identify myoepithelial boundary (ACTA2+, KRT15+)
- Look for regions where myoepithelial is disrupted

```python
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for i, gene in enumerate(['CEACAM6', 'FASN', 'ACTA2']):
    sq.pl.spatial_scatter(adata, color=gene, ax=axes[i], size=2, 
                         cmap='viridis', title=gene)
```

**Biological Insight**:
- Intact myoepithelial barrier → DCIS
- Disrupted barrier with FASN+ cells → Invasion

#### 2. How Does the Immune Infiltrate Distribute?

**Approach**:
- Identify T cell (CD3E+) and macrophage (ITGAX+) locations
- Test enrichment near tumor vs stroma
- Measure distances to tumor

```python
# Calculate distance of immune cells to nearest tumor cell
from scipy.spatial.distance import cdist

tumor_coords = adata[adata.obs['cell_type'] == 'Tumor_Invasive'].obsm['spatial']
tcell_coords = adata[adata.obs['cell_type'] == 'T_cells'].obsm['spatial']

distances = cdist(tcell_coords, tumor_coords).min(axis=1)

# Plot distribution
plt.hist(distances, bins=50)
plt.xlabel('Distance to Nearest Tumor Cell (µm)')
plt.ylabel('Number of T Cells')
plt.title('T Cell Distribution Relative to Tumor')
```

#### 3. What is the Stromal Architecture?

**Approach**:
- Identify fibroblast (POSTN+) patterns
- Measure density and organization
- Correlate with tumor progression markers

```python
# Stromal cell density in tumor vs normal regions
from scipy.spatial import Voronoi

# Compute local density
sq.gr.spatial_neighbors(adata, n_neighs=30, coord_type='generic')
adata.obs['neighbor_density'] = adata.obsp['spatial_connectivities'].sum(axis=1)

# Compare by region
stromal_cells = adata[adata.obs['cell_type'] == 'Stromal']
sq.pl.spatial_scatter(stromal_cells, color='neighbor_density', 
                      size=3, cmap='viridis')
```

## Export and Sharing

### 1. Export for Xenium Explorer

```python
# Export back to Xenium Explorer format
# (Requires xenium-exporter package)
from xenium_exporter import export_to_xenium

export_to_xenium(
    sdata,
    output_dir="results/janesick_breast_cancer/xenium_export",
    include_clusters=['leiden_0.8', 'cell_type'],
    include_scores=['Tumor_Invasive_score', 'T_cells_score']
)
```

### 2. Export Figures for Publication

```python
# High-resolution figure export
fig = plt.figure(figsize=(16, 16))

# Main spatial plot
sq.pl.spatial_scatter(
    adata,
    color='cell_type',
    size=3,
    legend_loc='right margin',
    title='Cell Type Annotation - Janesick Breast Cancer'
)

plt.savefig('results/figures/main_figure.pdf', dpi=600, bbox_inches='tight')
plt.savefig('results/figures/main_figure.png', dpi=300, bbox_inches='tight')
```

### 3. Create Interactive Visualizations

```python
# Napari visualization (if installed)
import napari

viewer = napari.Viewer()
sdata.pl.render(viewer)  # Render SpatialData in Napari
```

### 4. Share Results

```python
# Create summary report
summary = {
    'n_cells': adata.n_obs,
    'n_genes': adata.n_vars,
    'cell_types': adata.obs['cell_type'].value_counts().to_dict(),
    'median_counts': adata.obs['total_counts'].median(),
    'median_genes': adata.obs['n_genes_by_counts'].median()
}

import json
with open('results/janesick_breast_cancer/summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
```

## Further Resources

- **SpatialData**: https://spatialdata.scverse.org/
- **Scanpy tutorials**: https://scanpy.readthedocs.io/
- **Squidpy tutorials**: https://squidpy.readthedocs.io/
- **Xenium Explorer**: https://www.10xgenomics.com/support/software/xenium-explorer
- **Original Publication**: Janesick et al., Nat Commun (2023) DOI: 10.1038/s41467-023-43458-x
