# User Modules Directory

Drop custom `.py` module files in this directory. Any class that:

1. Subclasses `SpatchModule`
2. Is decorated with `@register`

Will be automatically discovered and available in the pipeline.

## Template

```python
# my_custom_analysis.py

from spatch_modules.base import SpatchModule, ModuleResult
from spatch_modules.registry import register


@register
class MyCustomAnalysis(SpatchModule):
    """Short description of what this module does."""
    
    name = "my_custom_analysis"  # Used to reference in YAML config
    version = "1.0.0"
    description = "Longer description for module listing"
    category = "analysis"  # One of: loader, analysis, annotation, qc
    
    # SpatialData elements this module requires
    requires = ["tables/table"]
    
    # SpatialData elements this module produces
    produces = ["tables/my_results"]

    def run(self, sdata, param1=10, param2="default", **kwargs):
        """Execute the analysis.
        
        Args:
            sdata: SpatialData object to analyze
            param1: Example numeric parameter
            param2: Example string parameter
            **kwargs: Additional parameters from YAML config
        
        Returns:
            ModuleResult with updated sdata and metrics
        """
        log = []
        
        # Your analysis logic here
        # ...
        
        # Example: add results to sdata
        # sdata.tables["my_results"] = ...
        
        return ModuleResult(
            sdata=sdata,
            metrics={"example_metric": 42},
            artifacts={},
            log=log
        )
```

## Using in Configuration

```yaml
custom_modules:
  user_module_dir: ./user_modules/
  
  steps:
    - module: my_custom_analysis
      enabled: true
      config:
        param1: 20
        param2: "custom_value"
```
