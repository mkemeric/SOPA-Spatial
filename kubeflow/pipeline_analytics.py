"""Kubeflow Pipeline: SPATCH analytics on an existing SpatialData zarr.

Runs the YAML-configured `spatch run` custom-module pipeline — QC/
preprocessing, clustering, spatial stats, shape metrics, diffusion
analysis, visualizations (see spatch_modules/builtin/) — against a zarr
produced by e.g. kubeflow/pipeline_raw_to_zarr.py, or any existing
SpatialData store with an aggregated cell x gene table.

Unlike pipeline_raw_to_zarr.py, this pipeline exposes NO per-module
parameters (QC thresholds, Leiden resolutions, marker genes, ...) — every
one of those already lives in the spatch config YAML (see configs/*.yaml),
which is the single source of truth `spatch run` reads from. The config's
full text is passed into the pipeline as a string parameter at submit
time (see kubeflow/submit_analytics_run.py), so changing analysis
parameters, or which zarr to run against, never requires touching this
file or recompiling the pipeline — only editing the config YAML and/or
the small run-params YAML that points at it.

If `mlflow_tracking_uri` is set, every module's metrics/artifacts (the
same dict `run_custom_pipeline` returns — see spatch_modules/runner.py)
are logged to that MLflow server under one run, alongside the config
text and the full parquet output directory. Leave it blank to skip
MLflow entirely.

Usage:
    python kubeflow/pipeline_analytics.py
    # writes kubeflow/pipeline_analytics.yaml, compiled and ready to
    # submit via kubeflow/submit_analytics_run.py
"""

from kfp import compiler, dsl
from kfp.dsl import Dataset, Input, Output

# Build this image from kubeflow/Dockerfile.analytics (spatch-modules +
# its dependency stack). Pin the tag to a specific build, never `:latest`,
# for reproducibility.
SPATCH_IMAGE = "ghcr.io/your-org/spatch-analytics:1.0.0"


@dsl.component(base_image=SPATCH_IMAGE)
def spatch_run_op(
    sdata_in: Input[Dataset],
    spatch_config_yaml: str,
    results: Output[Dataset],
    modules: list = [],  # noqa: B006 — empty = run every enabled step
    mlflow_tracking_uri: str = "",
    mlflow_experiment_name: str = "spatch-analytics",
):
    """Runs `run_custom_pipeline` in-process (rather than shelling out to
    the `spatch` CLI) so every module's metrics/artifacts dict is
    available to log to MLflow — the CLI only prints them to stdout.

    `spatch_config_yaml` is the *content* of a spatch config file (e.g.
    configs/janesick_breast_cancer.yaml), passed as a string rather than
    baked into the image — the pipeline never needs to be recompiled when
    the config changes.

    The input zarr is never mutated in place: only parquet outputs (and
    figures, if `pipeline_visualizations` is enabled) are written under
    `results.path`.
    """
    from pathlib import Path

    import spatialdata as sd
    import yaml

    from spatch_modules.runner import run_custom_pipeline

    config_path = Path("/tmp/spatch_config.yaml")
    config_path.write_text(spatch_config_yaml)
    Path(results.path).mkdir(parents=True, exist_ok=True)

    sdata = sd.read_zarr(sdata_in.path)
    module_results = run_custom_pipeline(
        sdata,
        config_path,
        module_names=modules or None,
        output_dir=results.path,
        verbose=True,
    )

    failed = {
        name: r["errors"] for name, r in module_results.items()
        if r["status"] == "failed"
    }

    if mlflow_tracking_uri:
        import mlflow

        def _flatten(d: dict, prefix: str = "") -> dict:
            """Dotted-key scalars only — nested lists (e.g. custom_modules.steps)
            aren't flattened, the full YAML is logged separately as an artifact.
            """
            flat = {}
            for k, v in d.items():
                key = f"{prefix}{k}"
                if isinstance(v, dict):
                    flat.update(_flatten(v, prefix=f"{key}."))
                elif isinstance(v, (str, int, float, bool)) or v is None:
                    flat[key] = v
            return flat

        mlflow.set_tracking_uri(mlflow_tracking_uri)
        mlflow.set_experiment(mlflow_experiment_name)
        with mlflow.start_run():
            mlflow.log_text(spatch_config_yaml, "spatch_config.yaml")
            mlflow.log_params(_flatten(yaml.safe_load(spatch_config_yaml)))
            mlflow.log_param("modules_requested", ",".join(modules) or "all")

            for name, r in module_results.items():
                mlflow.log_param(f"{name}.status", r["status"])
                for metric_name, value in r.get("metrics", {}).items():
                    if isinstance(value, (int, float)) and not isinstance(value, bool):
                        mlflow.log_metric(f"{name}.{metric_name}", value)
                    else:
                        mlflow.log_param(f"{name}.{metric_name}", value)

            mlflow.log_artifacts(results.path, artifact_path="spatch_outputs")

    if failed:
        raise RuntimeError(f"SPATCH modules failed: {failed}")


@dsl.pipeline(
    name="spatch-analytics",
    description="Run SPATCH's YAML-configured custom-module analysis on an existing SpatialData zarr.",
)
def spatch_analytics_pipeline(
    sdata_uri: str,
    spatch_config_yaml: str,
    modules: list = [],  # noqa: B006
    mlflow_tracking_uri: str = "",
    mlflow_experiment_name: str = "spatch-analytics",
) -> Dataset:
    """
    Args:
        sdata_uri: Location of the input SpatialData zarr — an S3/GCS URI
            or a path on cluster-visible storage (e.g. a mounted PVC).
            Same convention as `raw_data_uri` in pipeline_raw_to_zarr.py;
            brought in via `dsl.importer` without copying it ahead of time.
        spatch_config_yaml: Full text of a spatch config YAML (any file
            under configs/). Every module's parameters live here — this
            pipeline has no opinion on QC thresholds, resolutions, marker
            genes, etc.
        modules: Optional list of step names to run (matches `-m/--module`
            on `spatch run`, repeatable). Empty = run every step in the
            config with `enabled: true`.
        mlflow_tracking_uri: MLflow tracking server URL (e.g.
            `http://mlflow.mlflow.svc.cluster.local:5000` from inside the
            cluster). Leave blank to skip MLflow logging entirely.
        mlflow_experiment_name: MLflow experiment to log the run under.

    Returns:
        The `results` output directory: parquet outputs for every step
        (under `spatch/`), plus figures if `pipeline_visualizations` ran.
    """
    sdata = dsl.importer(
        artifact_uri=sdata_uri,
        artifact_class=Dataset,
        reimport=False,
    )

    run_task = spatch_run_op(
        sdata_in=sdata.output,
        spatch_config_yaml=spatch_config_yaml,
        modules=modules,
        mlflow_tracking_uri=mlflow_tracking_uri,
        mlflow_experiment_name=mlflow_experiment_name,
    )
    # Reliability: retry transient failures (e.g. a shared-cluster OOM or
    # a flaky read of sdata_uri) instead of failing the whole run outright.
    run_task.set_retry(num_retries=2, backoff_duration="30s", backoff_factor=2)
    run_task.set_cpu_request("4").set_cpu_limit("8")
    run_task.set_memory_request("16Gi").set_memory_limit("32Gi")

    return run_task.outputs["results"]


if __name__ == "__main__":
    output_path = __file__.replace(".py", ".yaml")
    compiler.Compiler().compile(
        pipeline_func=spatch_analytics_pipeline,
        package_path=output_path,
    )
    print(f"Compiled pipeline -> {output_path}")
