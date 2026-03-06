"""
Pipeline runner for SPATCH custom modules.

Executes modules from YAML configuration in sequence,
managing state and collecting results.

Each module's tabular output is persisted as a separate parquet file
under ``{output_dir}/spatch/``.  This avoids the spatialdata Dask-lock
issue that blocks zarr write-back, and enables true per-module
parallelism / container execution.
"""

from pathlib import Path
from typing import Any

import yaml
import spatialdata as sd

from .registry import get_module, discover_user_modules
from .io import (
    save_obs_columns,
    save_table,
    load_obs_parquets,
    load_table_parquet,
    list_saved_outputs,
    update_manifest,
)


# ── helpers ──────────────────────────────────────────────────────────

def _resolve_output_dir(config: dict, fallback: str = ".") -> str:
    """Pull output_dir from the YAML config's ``output`` section."""
    return config.get("output", {}).get("output_dir", fallback)


def merge_prior_outputs(
    sdata: sd.SpatialData,
    output_dir: str | Path,
    table_key: str = "table",
    verbose: bool = False,
) -> list[str]:
    """Load saved parquet outputs and merge into the in-memory *sdata*.

    Only columns / tables that are **not** already present are added,
    so re-running a module that already executed in the same session
    won't overwrite its in-memory results.

    Returns:
        List of log messages describing what was merged.
    """
    import anndata as ad
    from spatialdata.models import TableModel

    log: list[str] = []
    saved = list_saved_outputs(output_dir)

    # ── obs-level columns ────────────────────────────────────────
    obs_extra = load_obs_parquets(output_dir)
    if obs_extra is not None and table_key in sdata.tables:
        adata = sdata.tables[table_key]
        added: list[str] = []
        for col in obs_extra.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = obs_extra[col].reindex(adata.obs.index)
                added.append(col)
        if added:
            log.append(
                f"Merged {len(added)} prior obs columns: {', '.join(added)}"
            )

    # ── standalone tables ────────────────────────────────────────
    for tname in saved.get("tables", []):
        if tname in sdata.tables:
            continue
        df = load_table_parquet(tname, output_dir)
        if df is not None:
            sdata.tables[tname] = TableModel.parse(ad.AnnData(obs=df))
            log.append(f"Loaded standalone table '{tname}' from parquet")

    if verbose:
        for msg in log:
            print(f"  ↻ {msg}")

    return log


def _persist_module_outputs(
    sdata: sd.SpatialData,
    name: str,
    obs_before: set[str],
    tables_before: set[str],
    output_dir: str | Path,
    table_key: str = "table",
    verbose: bool = False,
) -> list[str]:
    """Diff obs columns / tables and write new outputs as parquet."""
    saved_files: list[str] = []

    # ── new obs columns ──────────────────────────────────────────
    if table_key in sdata.tables:
        obs_after = set(sdata.tables[table_key].obs.columns)
        new_cols = sorted(obs_after - obs_before)
        if new_cols:
            path = save_obs_columns(
                sdata.tables[table_key].obs, new_cols, name, output_dir,
            )
            if path:
                saved_files.append(str(path))
                if verbose:
                    print(f"    💾 Saved obs columns → {path.name}")

    # ── new standalone tables ────────────────────────────────────
    tables_after = set(sdata.tables.keys())
    for tname in sorted(tables_after - tables_before):
        tbl = sdata.tables[tname]
        path = save_table(tbl.obs, tname, output_dir)
        saved_files.append(str(path))
        if verbose:
            print(f"    💾 Saved table → {path.name}")

    if saved_files:
        update_manifest(name, output_dir, saved_files)

    return saved_files


# ── main pipeline entry point ────────────────────────────────────────

def run_custom_pipeline(
    sdata: sd.SpatialData,
    config_path: str | Path,
    module_names: list[str] | str | None = None,
    output_dir: str | None = None,
    verbose: bool = True,
    # Legacy alias
    module_name: str | None = None,
) -> dict[str, Any]:
    """Execute custom modules from a YAML config.

    Args:
        sdata: Input SpatialData object (will be modified in-place).
        config_path: Path to YAML configuration file.
        module_names: If set, run **only** these steps (by name).
            Accepts a single string or list of strings.
            Omit to run all enabled steps.
        output_dir: Directory for parquet outputs.  Defaults to the
            ``output.output_dir`` value in the YAML config.
        verbose: Whether to print progress messages.
        module_name: Deprecated — use module_names instead.

    Returns:
        Dictionary mapping module names to their results, including:
        - status: "completed", "skipped", or "failed"
        - metrics: Module-specific metrics
        - artifacts: Paths to generated files
        - log: Module log messages
        - errors: Error messages (if any)
        - saved_files: Parquet files written by this step
    """
    # Normalize module_names
    if module_names is None and module_name is not None:
        module_names = [module_name]
    elif isinstance(module_names, str):
        module_names = [module_names]
    run_set: set[str] | None = set(module_names) if module_names else None

    config = yaml.safe_load(Path(config_path).read_text())
    module_config = config.get("custom_modules", {})

    # Resolve output directory
    if output_dir is None:
        output_dir = _resolve_output_dir(config)

    # Discover user modules if directory specified
    user_dir = module_config.get("user_module_dir")
    if user_dir:
        loaded = discover_user_modules(user_dir)
        if verbose and loaded:
            print(f"Loaded user modules: {loaded}")

    # ── merge prior parquet outputs into sdata ───────────────────
    merge_log = merge_prior_outputs(sdata, output_dir, verbose=verbose)

    results: dict[str, Any] = {}

    for step in module_config.get("steps", []):
        if not step.get("enabled", True):
            continue

        name = step["module"]

        # If specific modules were requested, skip everything else
        if run_set is not None and name not in run_set:
            continue

        step_config = step.get("config", {})
        table_key = step_config.get("table_key", "table")

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

        # Snapshot before execution
        obs_before: set[str] = set()
        if table_key in sdata.tables:
            obs_before = set(sdata.tables[table_key].obs.columns)
        tables_before = set(sdata.tables.keys())

        # Execute
        if verbose:
            print(f"▶ Running: {name} v{module.version}")

        try:
            result = module.run(sdata, **step_config)
            sdata = result.sdata  # Carry forward modified SpatialData

            # Persist new tabular outputs as parquet
            saved = _persist_module_outputs(
                sdata, name, obs_before, tables_before,
                output_dir, table_key, verbose,
            )

            results[name] = {
                "status": "completed",
                "metrics": result.metrics,
                "artifacts": result.artifacts,
                "log": result.log,
                "saved_files": saved,
            }

            if verbose:
                metrics_str = ", ".join(
                    f"{k}={v}" for k, v in result.metrics.items()
                )
                print(f"  ✓ {name}: {metrics_str}")

        except Exception as e:
            if verbose:
                print(f"  ✗ {name}: Failed - {e}")
            results[name] = {
                "status": "failed",
                "errors": [str(e)],
            }

    return results


def run_single_module(
    sdata: sd.SpatialData,
    module_name: str,
    output_dir: str | None = None,
    **config,
) -> sd.SpatialData:
    """Run a single module by name.

    Convenience function for interactive/notebook use.  If *output_dir*
    is provided, prior parquet outputs are merged first and new outputs
    are saved after execution.

    Args:
        sdata: Input SpatialData object.
        module_name: Name of the registered module.
        output_dir: Optional directory for parquet I/O.
        **config: Configuration parameters for the module.

    Returns:
        Modified SpatialData object.
    """
    if output_dir:
        merge_prior_outputs(sdata, output_dir)

    module = get_module(module_name, **config)

    errors = module.validate_inputs(sdata)
    if errors:
        raise ValueError(f"Input validation failed: {errors}")

    table_key = config.get("table_key", "table")
    obs_before = set(
        sdata.tables[table_key].obs.columns
    ) if table_key in sdata.tables else set()
    tables_before = set(sdata.tables.keys())

    result = module.run(sdata, **config)

    if output_dir:
        _persist_module_outputs(
            result.sdata, module_name, obs_before, tables_before,
            output_dir, table_key,
        )

    print(f"✓ {module_name}: {result.metrics}")

    return result.sdata
