"""
Sopa CLI plugin — extends sopa with SPATCH subcommands at runtime.

This module imports sopa's Typer app, attaches the SPATCH command group,
and re-exports the combined app.  **No sopa source files are modified.**

When registered as an entry point::

    [project.scripts]
    sopa = "spatch_modules.sopa_plugin:app"

it shadows sopa's own ``sopa`` entry point so that ``sopa spatch run …``
works alongside all the original sopa commands (convert, patchify, etc.).

If spatch-modules is uninstalled, running ``pip install sopa`` restores
the original entry point.
"""

from sopa.cli.app import app  # sopa's existing Typer app (all built-in cmds)
from spatch_modules.cli import app_spatch  # our SPATCH subcommand group

app.add_typer(
    app_spatch,
    name="spatch",
    help="Run SPATCH custom analysis modules (tissue mask, morphology, "
         "diffusion, visualizations, and more).",
)

__all__ = ["app"]
