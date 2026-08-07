"""Submit a spatch-analytics Kubeflow run from a small YAML params file.

This is the file customers actually touch day-to-day: point a run-params
YAML (see kubeflow/runs/*.run.yaml) at a new zarr or a different spatch
config and resubmit. kubeflow/pipeline_analytics.py itself, and its
compiled kubeflow/pipeline_analytics.yaml, never change for any of that —
only for changes to the pipeline's DAG structure.

Usage:
    python kubeflow/pipeline_analytics.py   # once, or after a DAG change
    python kubeflow/submit_analytics_run.py kubeflow/runs/janesick_breast_cancer.run.yaml
"""

import sys
from pathlib import Path

import yaml
import kfp

PIPELINE_PACKAGE = Path(__file__).parent / "pipeline_analytics.yaml"


def main(params_path: str) -> None:
    params = yaml.safe_load(Path(params_path).read_text(encoding="utf-8"))

    config_text = Path(params["spatch_config_path"]).read_text(encoding="utf-8")

    client = kfp.Client(host=params["host"])
    run = client.create_run_from_pipeline_package(
        pipeline_file=str(PIPELINE_PACKAGE),
        arguments={
            "sdata_uri": params["sdata_uri"],
            "spatch_config_yaml": config_text,
            "modules": params.get("modules") or [],
            "mlflow_tracking_uri": params.get("mlflow_tracking_uri") or "",
            "mlflow_experiment_name": params.get("mlflow_experiment_name") or "spatch-analytics",
        },
        run_name=params.get("run_name"),
        experiment_name=params.get("experiment_name"),
    )
    print(f"Submitted run: {run.run_id}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {Path(__file__).name} <params.yaml>", file=sys.stderr)
        raise SystemExit(1)
    main(sys.argv[1])
