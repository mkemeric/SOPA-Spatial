"""
SPATCH Environment Helpers

Verifies that stock sopa/spatialdata/spatch_modules are importable
and configures matplotlib/scanpy for notebook use.

Usage in notebooks:
    from setup_local_imports import setup_notebook_environment
    setup_notebook_environment()
"""

import importlib
import importlib.metadata


def verify_imports(verbose=True):
    """
    Verify that all required packages can be imported.

    Returns
    -------
    dict
        Import status for each package.
    """
    packages = [
        'spatialdata',
        'spatialdata_io',
        'sopa',
        'spatch_modules',
        'scanpy',
        'squidpy',
        'numpy',
        'pandas',
    ]

    status = {}

    if verbose:
        print("=" * 70)
        print("Verifying package imports")
        print("=" * 70)

    for pkg_name in packages:
        try:
            mod = importlib.import_module(pkg_name)
            version = getattr(mod, '__version__', '?')
            status[pkg_name] = {'importable': True, 'version': version}
            if verbose:
                print(f"  {pkg_name:20s} \u2713  v{version}")
        except ImportError as e:
            status[pkg_name] = {'importable': False, 'error': str(e)}
            if verbose:
                print(f"  {pkg_name:20s} \u2717  {e}")

    if verbose:
        failed = [k for k, v in status.items() if not v['importable']]
        print("=" * 70)
        if failed:
            print(f"\n\u26a0\ufe0f  {len(failed)} package(s) failed to import.")
            print("   Run: ./setup_environment.sh")
        else:
            print("\nAll packages OK.")
        print()

    return status


def get_package_versions(verbose=True):
    """
    Get version information for key packages.

    Returns
    -------
    dict
        Version information.
    """
    packages = [
        'spatialdata', 'sopa', 'scanpy', 'squidpy', 'numpy', 'pandas',
    ]

    versions = {}

    if verbose:
        print("=" * 70)
        print("Package Versions")
        print("=" * 70)

    for pkg in packages:
        try:
            version = importlib.metadata.version(pkg)
            versions[pkg] = version
            if verbose:
                print(f"  {pkg:20s} {version}")
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = None
            if verbose:
                print(f"  {pkg:20s} not installed")

    if verbose:
        print("=" * 70)
        print()

    return versions


def setup_notebook_environment(verbose=True):
    """
    Complete setup for Jupyter notebooks.

    1. Verifies all imports work
    2. Configures matplotlib/plotting
    3. Configures scanpy defaults

    Usage at top of notebook:
        from setup_local_imports import setup_notebook_environment
        setup_notebook_environment()
    """
    # Verify imports
    status = verify_imports(verbose=verbose)
    failed = [k for k, v in status.items() if not v['importable']]
    if failed:
        if verbose:
            print("\u26a0\ufe0f  Some packages failed to import!")
            print("   Run: ./setup_environment.sh")
        return False

    # Configure matplotlib for notebooks
    try:
        import matplotlib.pyplot as plt
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['figure.facecolor'] = 'white'
        if verbose:
            print("\u2713 Matplotlib configured")
    except ImportError:
        pass

    # Configure scanpy
    try:
        import scanpy as sc
        sc.settings.verbosity = 2
        sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
        if verbose:
            print("\u2713 Scanpy configured")
    except ImportError:
        pass

    if verbose:
        print("\nNotebook environment ready!\n")

    return True


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("SPATCH Environment Diagnostic")
    print("=" * 70 + "\n")

    verify_imports(verbose=True)
    get_package_versions(verbose=True)
