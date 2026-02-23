"""
Module registry for SPATCH custom modules.

Provides discovery, registration, and instantiation of modules.
Supports both decorator-based registration and filesystem scanning
for user-contributed modules.
"""

import importlib
from pathlib import Path
from typing import Type

from .base import SpatchModule

_REGISTRY: dict[str, Type[SpatchModule]] = {}


def register(cls: Type[SpatchModule]) -> Type[SpatchModule]:
    """Decorator to register a module class.
    
    Example:
        @register
        class MyModule(SpatchModule):
            name = "my_module"
            ...
    """
    _REGISTRY[cls.name] = cls
    return cls


def get_module(name: str, **config) -> SpatchModule:
    """Instantiate a registered module by name with config.
    
    Args:
        name: The registered module name.
        **config: Configuration parameters passed to module constructor.
    
    Returns:
        Instantiated SpatchModule.
    
    Raises:
        KeyError: If module name is not found in registry.
    """
    if name not in _REGISTRY:
        raise KeyError(
            f"Module '{name}' not found. "
            f"Available: {list(_REGISTRY.keys())}"
        )
    return _REGISTRY[name](**config)


def list_modules(category: str | None = None) -> list[dict]:
    """List all registered modules, optionally filtered by category.
    
    Args:
        category: Optional filter by module category
                  ("loader", "analysis", "annotation", "qc").
    
    Returns:
        List of module metadata dictionaries.
    """
    modules = []
    for name, cls in sorted(_REGISTRY.items()):
        inst = cls()
        if category is None or inst.category == category:
            modules.append(inst.describe())
    return modules


def discover_user_modules(directory: str | Path) -> list[str]:
    """Scan a directory for .py files containing SpatchModule subclasses.
    
    Users drop .py files into this directory. Any class decorated with
    @register is automatically available in the registry.
    
    Args:
        directory: Path to directory containing user module files.
    
    Returns:
        List of successfully loaded module names.
    """
    directory = Path(directory)
    loaded = []
    
    if not directory.exists():
        return loaded

    import sys
    sys.path.insert(0, str(directory))
    
    for path in sorted(directory.glob("*.py")):
        if path.name.startswith("_"):
            continue
        module_name = path.stem
        try:
            importlib.import_module(module_name)
            loaded.append(module_name)
        except Exception as e:
            print(f"Warning: Failed to load user module {path.name}: {e}")
    
    sys.path.pop(0)
    return loaded


def get_registry() -> dict[str, Type[SpatchModule]]:
    """Return a copy of the current registry."""
    return _REGISTRY.copy()
