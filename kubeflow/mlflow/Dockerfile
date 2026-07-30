# MLflow tracking server image for the SPATCH kubeflow pipelines.
# Backend store: SQLite file on a PVC. Artifact store: local dir on the
# same PVC. Fine for single-node/dev use (see kubeflow/mlflow/deployment.yaml);
# swap both for Postgres + S3/GCS for a real multi-node deployment.
FROM python:3.12-slim

RUN pip install --no-cache-dir mlflow==2.19.0

EXPOSE 5000

CMD ["mlflow", "server", \
     "--backend-store-uri", "sqlite:////mlflow/db/mlflow.db", \
     "--default-artifact-root", "/mlflow/artifacts", \
     "--host", "0.0.0.0", \
     "--port", "5000"]
