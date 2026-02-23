"""
SPATCH Custom Modules

A plugin system for custom spatial transcriptomics analysis modules
built on the SpatialData + Sopa platform.
"""

from .base import SpatchModule, ModuleResult
from .registry import register, get_module, list_modules, discover_user_modules
from .runner import run_custom_pipeline

__version__ = "1.0.0"

__all__ = [
    "SpatchModule",
    "ModuleResult",
    "register",
    "get_module",
    "list_modules",
    "discover_user_modules",
    "run_custom_pipeline",
]

# Auto-discover builtin modules on import
from . import builtin  # noqa: F401
