"""
Typer-based CLI for SPATCH custom modules.

Provides the ``spatch`` command with subcommands that mirror sopa's
conventions (positional sdata_path, YAML-driven config).  Also usable
as a sopa subcommand group via ``spatch_modules.sopa_plugin``.

Usage (standalone)::

    spatch run results.zarr --config configs/janesick_breast_cancer.yaml
    spatch run results.zarr -c configs/janesick_breast_cancer.yaml -m dapi_tissue_mask --save
    spatch list
    spatch describe dapi_tissue_mask

Usage (via sopa, after spatch-modules is installed)::

    sopa spatch run results.zarr --config configs/janesick_breast_cancer.yaml
"""

from __future__ import annotations

import typer

# Typer app that can be mounted standalone *or* onto sopa's app.
app_spatch = typer.Typer(
    name="spatch",
    help="SPATCH custom analysis modules for spatial transcriptomics.",
    no_args_is_help=True,
)

# Standalone entry-point app (mirrors app_spatch, adds --version callback).
app = typer.Typer(no_args_is_help=True)


# ── helpers ──────────────────────────────────────────────────────────

def _ensure_modules():
    """Import builtin modules so the registry is populated."""
    import spatch_modules.builtin  # noqa: F401


def _format_table(rows: list[dict], columns: list[str]) -> str:
    """Simple column-aligned text table."""
    widths = {
        c: max(len(c), *(len(str(r.get(c, ""))) for r in rows)) + 2
        for c in columns
    }
    header = "".join(f"{c:<{widths[c]}}" for c in columns)
    sep = "-" * len(header)
    lines = [header, sep]
    for r in rows:
        lines.append(
            "".join(f"{str(r.get(c, '')):<{widths[c]}}" for c in columns)
        )
    return "\n".join(lines)


# ── run (shared implementation) ──────────────────────────────────────

def _cmd_run(
    sdata_path: str,
    config: str,
    module: str | None,
    save: bool,
    user_dir: str | None,
):
    """Backend for both ``spatch run`` and ``sopa spatch run``."""
    import spatialdata as sd
    from spatch_modules.runner import run_custom_pipeline
    from spatch_modules.registry import discover_user_modules

    _ensure_modules()

    if user_dir:
        discover_user_modules(user_dir)

    typer.echo(f"Loading: {sdata_path}")
    sdata = sd.read_zarr(sdata_path)

    results = run_custom_pipeline(
        sdata,
        config,
        module_name=module,
        verbose=True,
    )

    if save:
        typer.echo(f"Saving: {sdata_path}")
        sdata.write(sdata_path, overwrite=True)

    completed = sum(1 for r in results.values() if r["status"] == "completed")
    skipped = sum(1 for r in results.values() if r["status"] == "skipped")
    failed = sum(1 for r in results.values() if r["status"] == "failed")
    typer.echo(f"\nSummary: {completed} completed, {skipped} skipped, {failed} failed")

    if failed:
        raise typer.Exit(code=1)


# ── subcommands (on app_spatch — shared with sopa plugin) ───────────

@app_spatch.command()
def run(
    sdata_path: str = typer.Argument(
        help="Path to the SpatialData .zarr directory",
    ),
    config: str = typer.Option(
        ..., "--config", "-c", help="YAML config file",
    ),
    module: str = typer.Option(
        None, "--module", "-m",
        help="Run only this module (must match a step name in the config). "
             "Omit to run all enabled modules.",
    ),
    save: bool = typer.Option(
        False, "--save", "-s",
        help="Write the modified SpatialData back to sdata_path after running. "
             "Required for chained container execution.",
    ),
    user_dir: str = typer.Option(
        None, "--user-dir",
        help="Directory containing user-contributed module .py files",
    ),
):
    """Run SPATCH custom modules on a SpatialData object."""
    _cmd_run(sdata_path, config, module, save, user_dir)


@app_spatch.command("list")
def list_modules(
    category: str = typer.Option(
        None,
        help="Filter by category (loader, analysis, annotation, qc)",
    ),
):
    """List all registered SPATCH modules."""
    from spatch_modules.registry import list_modules as _list

    _ensure_modules()
    modules = _list(category=category)

    if not modules:
        typer.echo("No modules registered.")
        raise typer.Exit()

    typer.echo(
        _format_table(modules, ["name", "version", "category", "description"])
    )


@app_spatch.command()
def describe(
    name: str = typer.Argument(help="Module name"),
    user_dir: str = typer.Option(None, "--user-dir"),
):
    """Show detailed information about a SPATCH module."""
    from spatch_modules.registry import get_module, discover_user_modules

    _ensure_modules()
    if user_dir:
        discover_user_modules(user_dir)

    try:
        mod = get_module(name)
    except KeyError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)

    info = mod.describe()
    typer.echo(f"Module:      {info['name']} v{info['version']}")
    typer.echo(f"Category:    {info['category']}")
    typer.echo(f"Description: {info['description']}")
    typer.echo(f"Requires:    {info['requires'] or '(none)'}")
    typer.echo(f"Produces:    {info['produces'] or '(none)'}")
    if info["config_schema"]:
        typer.echo("\nConfig parameters:")
        for param, ptype in info["config_schema"].items():
            typer.echo(f"  - {param}: {ptype}")


# ── standalone entry-point mirrors ───────────────────────────────────
# Register the same commands on ``app`` so ``spatch run`` works.

app.command()(run)
app.command("list")(list_modules)
app.command()(describe)


def _version_callback(value: bool):
    if value:
        from spatch_modules import __version__

        typer.echo(f"spatch-modules {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: bool = typer.Option(
        None, "--version", callback=_version_callback,
        help="Show the spatch-modules version and exit.",
    ),
):
    """SPATCH — custom spatial transcriptomics analysis modules."""
    pass
