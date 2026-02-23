"""
Built-in SPATCH modules.

These modules are automatically registered when spatch_modules is imported.
"""

from . import codex_loader
from . import diffusion_analysis
from . import gene_protein_correlation
from . import cell_shape_metrics

__all__ = [
    "codex_loader",
    "diffusion_analysis",
    "gene_protein_correlation",
    "cell_shape_metrics",
]
