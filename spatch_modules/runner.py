"""
Pipeline runner for SPATCH custom modules.

Executes modules from YAML configuration in sequence,
managing state and collecting results.
"""

from pathlib import Path
from typing import Any

import yaml
import spatialdata as sd

from .registry import get_module, discover_user_modules


def run_custom_pipeline(
    sdata: sd.SpatialData,
    config_path: str | Path,
    verbose: bool = True
) -> dict[str, Any]:
    """Execute all enabled custom modules from a YAML config.
    
    Args:
        sdata: Input SpatialData object (will be modified in-place).
        config_path: Path to YAML configuration file.
        verbose: Whether to print progress messages.
    
    Returns:
        Dictionary mapping module names to their results, including:
        - status: "completed", "skipped", or "failed"
        - metrics: Module-specific metrics
        - artifacts: Paths to generated files
        - log: Module log messages
        - errors: Error messages (if any)
    """
    config = yaml.safe_load(Path(config_path).read_text())
    module_config = config.get("custom_modules", {})
    
    # Discover user modules if directory specified
    user_dir = module_config.get("user_module_dir")
    if user_dir:
        loaded = discover_user_modules(user_dir)
        if verbose and loaded:
            print(f"Loaded user modules: {loaded}")
    
    results = {}
    
    for step in module_config.get("steps", []):
        if not step.get("enabled", True):
            continue
        
        name = step["module"]
        step_config = step.get("config", {})
        
        try:
            module = get_module(name, **step_config)
        except KeyError as e:
            if verbose:
                print(f"✗ {name}: Module not found")
            results[name] = {"status": "failed", "errors": [str(e)]}
            continue
        
        # Validate inputs
        errors = module.validate_inputs(sdata)
        if errors:
            if verbose:
                print(f"⊘ {name}: Skipped - {errors}")
            results[name] = {"status": "skipped", "errors": errors}
            continue
        
        # Execute
        if verbose:
            print(f"▶ Running: {name} v{module.version}")
        
        try:
            result = module.run(sdata, **step_config)
            sdata = result.sdata  # Carry forward modified SpatialData
            
            results[name] = {
                "status": "completed",
                "metrics": result.metrics,
                "artifacts": result.artifacts,
                "log": result.log
            }
            
            if verbose:
                metrics_str = ", ".join(f"{k}={v}" for k, v in result.metrics.items())
                print(f"  ✓ {name}: {metrics_str}")
                
        except Exception as e:
            if verbose:
                print(f"  ✗ {name}: Failed - {e}")
            results[name] = {
                "status": "failed",
                "errors": [str(e)]
            }
    
    return results


def run_single_module(
    sdata: sd.SpatialData,
    module_name: str,
    **config
) -> sd.SpatialData:
    """Run a single module by name.
    
    Convenience function for interactive/notebook use.
    
    Args:
        sdata: Input SpatialData object.
        module_name: Name of the registered module.
        **config: Configuration parameters for the module.
    
    Returns:
        Modified SpatialData object.
    """
    module = get_module(module_name, **config)
    
    errors = module.validate_inputs(sdata)
    if errors:
        raise ValueError(f"Input validation failed: {errors}")
    
    result = module.run(sdata, **config)
    
    print(f"✓ {module_name}: {result.metrics}")
    
    return result.sdata
