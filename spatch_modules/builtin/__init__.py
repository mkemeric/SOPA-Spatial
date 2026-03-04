"""
Built-in SPATCH modules.

These modules are automatically registered when spatch_modules is imported.
"""

from . import codex_loader
from . import diffusion_analysis
from . import gene_protein_correlation
from . import cell_shape_metrics
from . import grid_binning
from . import dapi_tissue_mask
from . import image_registration
from . import annotation_consensus
from . import spatial_cluster
from . import pipeline_visualizations

__all__ = [
    "codex_loader",
    "diffusion_analysis",
    "gene_protein_correlation",
    "cell_shape_metrics",
    "grid_binning",
    "dapi_tissue_mask",
    "image_registration",
    "annotation_consensus",
    "spatial_cluster",
    "pipeline_visualizations",
]
