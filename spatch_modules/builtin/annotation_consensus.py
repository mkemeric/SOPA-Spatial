"""
Annotation Consensus Module

Run multiple cell-type annotation tools and produce a majority-vote
consensus label.  Based on the original SPATCH 9_st_annotation.py.

Each external tool is imported lazily and gated behind a try/except,
so users only need the tools they actually enable.

Supported tools:
  - Tangram       (``tangram-sc``)     — transcript-based mapping
  - CellTypist    (``celltypist``)     — logistic-regression classifier
  - TACCO         (``tacco``)          — transfer annotation
  - Spoint        (``cell2location``/custom)  — spatial deconvolution
  - Selina        (user-provided path) — autoencoder classifier
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import spatialdata as sd

from ..base import SpatchModule, ModuleResult
from ..registry import register


@register
class AnnotationConsensus(SpatchModule):
    """Multi-tool cell-type annotation with majority-vote consensus.

    Workflow:
      1. For each enabled tool, run annotation and store per-cell predictions.
      2. Apply row-wise majority vote across all predictions.
      3. For cells where the vote is "unknown" / tied, fall back to the
         single tool that has the highest overall agreement with the
         consensus (best-correlated method).
      4. Store individual tool columns *and* the final consensus in
         ``adata.obs``.
    """

    name = "annotation_consensus"
    version = "1.0.0"
    description = "Multi-tool cell-type annotation with majority-vote consensus"
    category = "annotation"
    requires = ["tables/table"]
    produces = []  # adds obs columns

    def run(
        self,
        sdata: sd.SpatialData,
        table_key: str = "table",
        reference_path: str = "",
        reference_celltype_key: str = "major",
        tools: dict[str, dict] | None = None,
        consensus_col: str = "cell_type_consensus",
        unknown_label: str = "unknown",
        **kwargs,
    ) -> ModuleResult:
        """Run annotation consensus.

        Args:
            sdata: SpatialData object.
            table_key: Key for the cell expression table.
            reference_path: Path to the scRNA-seq reference h5ad.
            reference_celltype_key: obs column in reference with cell types.
            tools: Dict of tool_name → tool-specific config.  Example::

                    tools:
                      tangram:
                        enabled: true
                      celltypist:
                        enabled: true
                        model: Immune_All_Low.pkl
                      tacco:
                        enabled: true
                      spoint:
                        enabled: false
                      selina:
                        enabled: false
                        selina_path: /path/to/Selina

            consensus_col: Name of the final consensus obs column.
            unknown_label: Label assigned when no majority exists.
            **kwargs: Extra config (ignored).

        Returns:
            ModuleResult with annotation columns added.
        """
        log = []

        if table_key not in sdata.tables:
            raise ValueError(f"Table '{table_key}' not found")
        adata = sdata.tables[table_key]

        tools = tools or self.config.get("tools", {})
        if not tools:
            raise ValueError("No annotation tools configured")

        # Load reference once (shared across tools that need it)
        sc_ref = None
        if reference_path:
            sc_ref = sc.read_h5ad(reference_path)
            log.append(f"Loaded reference: {sc_ref.n_obs} cells, "
                       f"{sc_ref.n_vars} genes")

        # ── Run each tool ────────────────────────────────────────────
        predictions: dict[str, np.ndarray] = {}

        for tool_name, tool_cfg in tools.items():
            if not tool_cfg.get("enabled", True):
                continue
            try:
                preds = self._run_tool(
                    tool_name, adata, sc_ref, reference_celltype_key, tool_cfg
                )
                col = f"ct_{tool_name}"
                adata.obs[col] = preds
                predictions[col] = preds
                log.append(f"✓ {tool_name}: {len(set(preds))} types")
            except Exception as e:
                log.append(f"✗ {tool_name}: {e}")

        if not predictions:
            log.append("WARNING: No tools produced predictions")
            return ModuleResult(sdata=sdata, metrics={"n_tools": 0}, log=log)

        # ── Majority vote ────────────────────────────────────────────
        pred_df = pd.DataFrame(predictions, index=adata.obs_names)
        consensus = pred_df.apply(
            lambda row: self._majority_vote(row, unknown_label), axis=1
        )

        # Fallback: for unknowns, use the best-correlated single tool
        n_unknown = (consensus == unknown_label).sum()
        if n_unknown > 0 and len(predictions) > 1:
            best_tool = self._best_correlated_tool(pred_df, consensus)
            fallback = pred_df[best_tool]
            consensus = consensus.where(
                consensus != unknown_label, fallback
            )
            log.append(f"Fallback tool for {n_unknown} unknowns: {best_tool}")

        adata.obs[consensus_col] = consensus.values
        log.append(
            f"Consensus: {len(consensus.unique())} types, "
            f"{(consensus == unknown_label).sum()} remaining unknowns"
        )

        return ModuleResult(
            sdata=sdata,
            metrics={
                "n_tools": len(predictions),
                "n_cell_types": len(consensus.unique()),
                "n_unknown": int((consensus == unknown_label).sum()),
            },
            log=log,
        )

    # ── voting helpers ───────────────────────────────────────────────

    @staticmethod
    def _majority_vote(row: pd.Series, unknown_label: str) -> str:
        """Row-wise majority vote across tool columns."""
        counts = Counter(row.values)
        if not counts:
            return unknown_label
        most_common = counts.most_common(2)
        # Require strict majority (top > second)
        if len(most_common) == 1 or most_common[0][1] > most_common[1][1]:
            return most_common[0][0]
        return unknown_label

    @staticmethod
    def _best_correlated_tool(
        pred_df: pd.DataFrame, consensus: pd.Series
    ) -> str:
        """Find the tool column most consistent with the consensus."""
        agreements = {}
        for col in pred_df.columns:
            agreements[col] = (pred_df[col] == consensus).sum()
        return max(agreements, key=agreements.get)

    # ── tool runners ─────────────────────────────────────────────────

    def _run_tool(
        self,
        name: str,
        adata: ad.AnnData,
        sc_ref: ad.AnnData | None,
        celltype_key: str,
        cfg: dict,
    ) -> np.ndarray:
        """Dispatch to the correct tool runner."""
        runners = {
            "tangram": self._run_tangram,
            "celltypist": self._run_celltypist,
            "tacco": self._run_tacco,
            "spoint": self._run_spoint,
            "selina": self._run_selina,
        }
        if name not in runners:
            raise ValueError(
                f"Unknown tool '{name}'. Available: {list(runners)}"
            )
        return runners[name](adata, sc_ref, celltype_key, cfg)

    # ── Tangram ──────────────────────────────────────────────────────

    @staticmethod
    def _run_tangram(
        adata: ad.AnnData,
        sc_ref: ad.AnnData,
        celltype_key: str,
        cfg: dict,
    ) -> np.ndarray:
        import tangram as tg

        markers = cfg.get("markers", None)
        if markers is None:
            # Compute DEGs as markers
            sc_ref_copy = sc_ref.copy()
            sc.tl.rank_genes_groups(sc_ref_copy, groupby=celltype_key)
            df = pd.DataFrame(sc_ref_copy.uns["rank_genes_groups"]["names"])
            markers = list(np.unique(df.iloc[:50].values.flatten()))

        tg.pp_adatas(sc_ref, adata, genes=markers)
        ad_map = tg.map_cells_to_space(
            adata_sc=sc_ref,
            adata_sp=adata,
            device=cfg.get("device", "cpu"),
            mode="clusters",
            cluster_label=celltype_key,
        )
        tg.project_cell_annotations(ad_map, adata, annotation=celltype_key)
        pred_key = f"tangram_ct_pred"
        if pred_key in adata.obsm:
            return adata.obsm[pred_key].idxmax(axis=1).values
        return np.full(adata.n_obs, "unknown")

    # ── CellTypist ───────────────────────────────────────────────────

    @staticmethod
    def _run_celltypist(
        adata: ad.AnnData,
        sc_ref: ad.AnnData,
        celltype_key: str,
        cfg: dict,
    ) -> np.ndarray:
        import celltypist

        st_ad = adata.copy()
        sc.pp.normalize_total(st_ad, target_sum=1e4)
        sc.pp.log1p(st_ad)

        model_path = cfg.get("model", None)
        if model_path:
            model = celltypist.models.Model.load(model=model_path)
        elif sc_ref is not None:
            ref_copy = sc_ref.copy()
            sc.pp.normalize_total(ref_copy, target_sum=1e4)
            sc.pp.log1p(ref_copy)
            model = celltypist.train(
                ref_copy,
                labels=celltype_key,
                n_jobs=cfg.get("n_jobs", 10),
                feature_selection=True,
            )
        else:
            raise ValueError("CellTypist requires a model path or reference")

        predictions = celltypist.annotate(
            st_ad, model=model, majority_voting=False
        )
        return predictions.predicted_labels["predicted_labels"].values

    # ── TACCO ────────────────────────────────────────────────────────

    @staticmethod
    def _run_tacco(
        adata: ad.AnnData,
        sc_ref: ad.AnnData,
        celltype_key: str,
        cfg: dict,
    ) -> np.ndarray:
        import tacco as tc

        st_ad = adata.copy()
        st_ad.X = st_ad.X.astype(np.float32)
        ref = sc_ref.copy()
        ref.X = ref.X.astype(np.float32)

        tc.tl.annotate(st_ad, ref, celltype_key, result_key=celltype_key)
        if celltype_key in st_ad.obs:
            return st_ad.obs[celltype_key].values
        # TACCO stores proportions in obsm
        if celltype_key in st_ad.obsm:
            return st_ad.obsm[celltype_key].idxmax(axis=1).values
        return np.full(adata.n_obs, "unknown")

    # ── Spoint ───────────────────────────────────────────────────────

    @staticmethod
    def _run_spoint(
        adata: ad.AnnData,
        sc_ref: ad.AnnData,
        celltype_key: str,
        cfg: dict,
    ) -> np.ndarray:
        """Spoint spatial deconvolution."""
        try:
            from Spoint import Spoint
        except ImportError:
            raise ImportError("Spoint is required. Install it manually.")

        import torch

        spoint_model = Spoint.init_model(
            sc_ref, adata,
            celltype_key=celltype_key,
            deg_method="t-test",
            sm_size=cfg.get("sm_size", 100000),
            use_gpu=cfg.get("use_gpu", False),
        )
        spoint_model.train(
            max_steps=cfg.get("max_steps", 5000),
            batch_size=cfg.get("batch_size", 1024),
            save_mode="best",
        )
        st_data = torch.tensor(spoint_model.st_data).to("cpu")
        model = spoint_model.model.to("cpu")
        model.eval()
        _, pre, _ = model(st_data)
        pre = pre.detach().numpy()
        pre_df = pd.DataFrame(
            pre, columns=spoint_model.clusters, index=adata.obs_names
        )
        return pre_df.idxmax(axis=1).values

    # ── Selina ───────────────────────────────────────────────────────

    @staticmethod
    def _run_selina(
        adata: ad.AnnData,
        sc_ref: ad.AnnData,
        celltype_key: str,
        cfg: dict,
    ) -> np.ndarray:
        """Selina autoencoder classifier."""
        import sys
        selina_path = cfg.get("selina_path", "")
        if selina_path:
            sys.path.insert(0, selina_path)

        try:
            from selina.selina import (
                preprocessing, label2dic, train, Autoencoder,
                Normal_Classifier, test,
            )
        except ImportError:
            raise ImportError(
                "Selina is required. Provide selina_path in config."
            )

        import torch
        device = torch.device(cfg.get("device", "cpu"))
        params_train = cfg.get("params_train", {})

        common_genes = np.intersect1d(
            sc_ref.var_names.tolist(), adata.var_names.tolist()
        )
        train_data = sc_ref[:, common_genes].to_df().T
        st_data = adata[:, common_genes.tolist()]

        train_data_proc, celltypes, platforms, genes = preprocessing(
            train_data, sc_ref.obs
        )
        ct_dic = label2dic(celltypes)
        plat_dic = label2dic(platforms)
        nfeatures = train_data_proc.shape[0]
        nct = len(ct_dic)
        nplat = len(plat_dic)

        network = train(
            train_data_proc, params_train, celltypes, platforms,
            nfeatures, nct, nplat, ct_dic, plat_dic, device,
        )
        network = Autoencoder(network, nfeatures, nct)
        network = Normal_Classifier(network).to(device)

        pred_labels, _ = test(st_data, network, ct_dic, device)
        return np.array(pred_labels)
