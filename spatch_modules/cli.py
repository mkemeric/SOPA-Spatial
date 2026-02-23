"""
Command-line interface for SPATCH modules.

Provides commands for listing, running, and managing modules.
"""

import argparse
import sys
from pathlib import Path

from .registry import list_modules, discover_user_modules, get_module


def cmd_list(args):
    """List all registered modules."""
    if args.user_dir:
        discover_user_modules(args.user_dir)
    
    modules = list_modules(category=args.category)
    
    if not modules:
        print("No modules registered.")
        return
    
    # Calculate column widths
    name_width = max(len(m["name"]) for m in modules) + 2
    ver_width = max(len(m["version"]) for m in modules) + 2
    cat_width = max(len(m["category"]) for m in modules) + 2
    
    # Header
    header = (
        f"{'Name':<{name_width}}"
        f"{'Version':<{ver_width}}"
        f"{'Category':<{cat_width}}"
        f"Description"
    )
    print(header)
    print("-" * len(header))
    
    # Rows
    for m in modules:
        print(
            f"{m['name']:<{name_width}}"
            f"{m['version']:<{ver_width}}"
            f"{m['category']:<{cat_width}}"
            f"{m['description']}"
        )


def cmd_describe(args):
    """Show detailed info about a module."""
    if args.user_dir:
        discover_user_modules(args.user_dir)
    
    try:
        module = get_module(args.module)
    except KeyError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    info = module.describe()
    
    print(f"Module: {info['name']} v{info['version']}")
    print(f"Category: {info['category']}")
    print(f"Description: {info['description']}")
    print()
    print(f"Requires: {info['requires'] or '(none)'}")
    print(f"Produces: {info['produces'] or '(none)'}")
    
    if info["config_schema"]:
        print()
        print("Config parameters:")
        for param, ptype in info["config_schema"].items():
            print(f"  - {param}: {ptype}")


def cmd_run(args):
    """Run the custom module pipeline."""
    import spatialdata as sd
    from .runner import run_custom_pipeline
    
    if args.user_dir:
        discover_user_modules(args.user_dir)
    
    # Load input data
    print(f"Loading: {args.input}")
    sdata = sd.read_zarr(args.input)
    
    # Run pipeline
    results = run_custom_pipeline(sdata, args.config, verbose=True)
    
    # Save output
    if args.output:
        print(f"Saving: {args.output}")
        sdata.write(args.output)
    
    # Summary
    completed = sum(1 for r in results.values() if r["status"] == "completed")
    skipped = sum(1 for r in results.values() if r["status"] == "skipped")
    failed = sum(1 for r in results.values() if r["status"] == "failed")
    
    print()
    print(f"Summary: {completed} completed, {skipped} skipped, {failed} failed")


def main():
    """Main entry point for CLI."""
    parser = argparse.ArgumentParser(
        prog="spatch-modules",
        description="SPATCH Custom Module CLI"
    )
    parser.add_argument(
        "--user-dir",
        help="Directory containing user module files"
    )
    
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # List command
    list_parser = subparsers.add_parser("list", help="List registered modules")
    list_parser.add_argument(
        "--category",
        choices=["loader", "analysis", "annotation", "qc"],
        help="Filter by category"
    )
    list_parser.set_defaults(func=cmd_list)
    
    # Describe command
    desc_parser = subparsers.add_parser("describe", help="Show module details")
    desc_parser.add_argument("module", help="Module name")
    desc_parser.set_defaults(func=cmd_describe)
    
    # Run command
    run_parser = subparsers.add_parser("run", help="Run module pipeline")
    run_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input SpatialData zarr path"
    )
    run_parser.add_argument(
        "-c", "--config",
        required=True,
        help="Pipeline YAML config file"
    )
    run_parser.add_argument(
        "-o", "--output",
        help="Output SpatialData zarr path"
    )
    run_parser.set_defaults(func=cmd_run)
    
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
