"""
Base classes for SPATCH custom modules.

Every module implements the SpatchModule abstract base class,
providing a uniform interface for discovery, configuration, and execution.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

import spatialdata as sd


@dataclass
class ModuleResult:
    """Standard return type for all modules.
    
    Attributes:
        sdata: Modified SpatialData object (or original if read-only)
        metrics: Numeric results / QC statistics
        artifacts: Paths to generated files (plots, CSVs)
        log: Human-readable log messages
    """
    sdata: sd.SpatialData
    metrics: dict[str, Any] = field(default_factory=dict)
    artifacts: dict[str, str] = field(default_factory=dict)
    log: list[str] = field(default_factory=list)


class SpatchModule(ABC):
    """Base class for all SPATCH custom modules.
    
    Subclasses must implement the `run` method and define module metadata
    via class attributes.
    
    Example:
        @register
        class MyModule(SpatchModule):
            name = "my_module"
            version = "1.0.0"
            description = "My custom analysis"
            category = "analysis"
            requires = ["tables/table"]
            produces = ["tables/my_results"]
            
            def run(self, sdata, **kwargs):
                # ... analysis logic ...
                return ModuleResult(sdata=sdata, metrics={"result": 42})
    """

    # ── Module metadata (override in subclass) ─────────────────
    name: str = "unnamed"
    version: str = "0.1.0"
    description: str = ""
    category: str = "analysis"  # "loader" | "analysis" | "annotation" | "qc"
    
    # What this module requires to already exist in the SpatialData object
    requires: list[str] = []  # e.g., ["tables/table", "shapes/cell_boundaries"]
    
    # What this module produces / modifies in the SpatialData object
    produces: list[str] = []  # e.g., ["tables/diffusion_metrics"]

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
        """Check that required elements exist.
        
        Returns:
            List of error messages (empty list = OK).
        """
        errors = []
        for req in self.requires:
            parts = req.split("/")
            if len(parts) == 2:
                container_name, element_name = parts
                container = getattr(sdata, container_name, None)
                if container is None:
                    errors.append(f"Missing required container: {container_name}")
                elif element_name not in container:
                    errors.append(f"Missing required element: {req}")
            else:
                # Single-level requirement (e.g., just check container exists)
                container = getattr(sdata, parts[0], None)
                if container is None:
                    errors.append(f"Missing required container: {parts[0]}")
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
