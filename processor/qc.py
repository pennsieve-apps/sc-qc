"""
The per-sample method after ambient removal: QC metrics, outlier flags by median absolute
deviations (Heumos et al. 2023, sc-best-practices "Quality control"), doublets with Scrublet
(scanpy.pp.scrublet, Wolock et al. 2019), normalization (normalize_total + log1p).

Observations are never reordered: removal is a boolean mask in file order.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

log = logging.getLogger("sc-qc")

# In the order they are listed in qc.csv's `reason` column.
REASONS = ("ambient_empty", "low_counts", "high_counts", "low_genes", "high_genes", "low_top20", "high_top20", "high_mito",
           "doublet")
MIN_CELLS_FOR_SCRUBLET = 100       # Scrublet's kNN and 30 PCs need more cells than this
CELL_PROBABILITY_CUTOFF = 0.5      # CellBender's own cutoff for calling a droplet a cell


@dataclass
class Thresholds:
    """One row per metric for qc_summary.csv and thresholds.csv."""
    rows: list[dict] = field(default_factory=list)


def mito_ribo(adata: AnnData) -> None:
    """var["mt"], var["ribo"] (kept from sc-prep, else MT-/mt- and RPS/RPL prefixes, case-insensitive)."""
    upper = adata.var_names.str.upper()
    if "mt" not in adata.var:
        adata.var["mt"] = upper.str.startswith("MT-")
    if "ribo" not in adata.var:
        adata.var["ribo"] = upper.str.startswith(("RPS", "RPL"))
    adata.var["mt"] = adata.var["mt"].astype(bool)
    adata.var["ribo"] = adata.var["ribo"].astype(bool)


def qc_metrics(counts: sp.csr_matrix, var: pd.DataFrame, obs_names) -> pd.DataFrame:
    """total_counts, n_genes_by_counts, pct_counts_mt, pct_counts_ribo, pct_counts_in_top_20_genes
    of a count matrix (scanpy.pp.calculate_qc_metrics)."""
    a = AnnData(X=counts, var=var[["mt", "ribo"]].copy(), obs=pd.DataFrame(index=obs_names))
    top = [20] if a.n_vars > 20 else None
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=top, log1p=False, inplace=True)
    out = a.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]].astype(np.float64)
    out["pct_counts_in_top_20_genes"] = a.obs["pct_counts_in_top_20_genes"] if top else 100.0
    # A cell without counts has 0/0 % mt; call it 0 so medians stay defined (it fails the count cut anyway).
    return out.fillna(0.0)


def _median_mad(values: np.ndarray) -> tuple[float, float]:
    med = float(np.median(values))
    return med, float(np.median(np.abs(values - med)))


def flag(metrics: pd.DataFrame, params: dict, basis: np.ndarray | None = None) -> tuple[pd.DataFrame, Thresholds]:
    """One boolean column per QC reason (not doublets or empty droplets), plus the thresholds.

    As in sc-best-practices: log1p(total_counts), log1p(n_genes_by_counts) and the top-20-gene share
    outside median ± MADS × MAD, both directions. Mitochondrial: pct_counts_mt above MAX_PCT_MT (a
    fixed cap; see README). `basis` (a boolean mask) restricts the cells the medians and MADs are
    computed from; every cell is flagged against them."""
    k = params["MADS"]
    thr = Thresholds()
    flags = pd.DataFrame(False, index=metrics.index, columns=list(REASONS))
    for metric, low, high, log1p in (("total_counts", "low_counts", "high_counts", True),
                                     ("n_genes_by_counts", "low_genes", "high_genes", True),
                                     ("pct_counts_in_top_20_genes", "low_top20", "high_top20", False)):
        raw = metrics[metric].to_numpy()
        values = np.log1p(raw) if log1p else raw
        med, mad = _median_mad(values if basis is None or not basis.any() else values[basis])
        back = np.expm1 if log1p else (lambda v: v)
        lo, hi = back(med - k * mad), back(med + k * mad)
        if mad == 0:
            log.warning("MAD of %s is 0; no cut on it", metric)
        else:
            flags[low] = raw < lo
            flags[high] = raw > hi
        thr.rows.append({"metric": metric, "median": float(back(med)), "mad": mad, "scale": "log1p" if log1p else "linear",
                         "lower": float(lo) if mad else None, "upper": float(hi) if mad else None,
                         "rule": f"median ± {k:g} MAD"})
    flags["high_mito"] = metrics["pct_counts_mt"].to_numpy() > params["MAX_PCT_MT"]
    thr.rows.append({"metric": "pct_counts_mt", "median": float(metrics["pct_counts_mt"].median()), "mad": None,
                     "scale": "linear", "lower": None, "upper": params["MAX_PCT_MT"], "rule": "> MAX_PCT_MT"})
    return flags, thr


def scrublet(counts: sp.csr_matrix, obs_names, params: dict) -> tuple[pd.DataFrame, dict]:
    """scanpy.pp.scrublet on the counts of the cells that passed QC. Returns obs (doublet_score,
    predicted_doublet) and a summary (threshold, simulated scores, rates)."""
    a = AnnData(X=counts.astype(np.float32), obs=pd.DataFrame(index=obs_names))
    sc.pp.scrublet(a, expected_doublet_rate=params["DOUBLET_RATE"], random_state=params["SEED"], verbose=False)
    info = a.uns["scrublet"]
    threshold = info.get("threshold")
    if threshold is None:
        log.warning("Scrublet found no threshold (the simulated doublet scores are not bimodal); nothing called a doublet")
    obs = pd.DataFrame({"doublet_score": a.obs["doublet_score"].to_numpy(dtype=np.float32),
                        "predicted_doublet": a.obs["predicted_doublet"].to_numpy(dtype=bool)}, index=obs_names)
    summary = {"threshold": None if threshold is None else float(threshold),
               "simulated_scores": np.asarray(info["doublet_scores_sim"], dtype=np.float32),
               "n_predicted": int(obs["predicted_doublet"].sum())}
    return obs, summary


def reasons_of(flags: pd.DataFrame) -> pd.Series:
    """'low_counts;high_mito' per flagged cell, '' when kept."""
    out = pd.Series("", index=flags.index, dtype=object)
    for r in REASONS:
        hit = flags[r].to_numpy()
        out[hit] = np.where(out[hit] == "", r, out[hit] + ";" + r)
    return out


def normalize(adata: AnnData, target_sum: float) -> None:
    """X = log1p(normalize_total(layers["counts"], target_sum)), float32 CSR."""
    adata.X = adata.layers["counts"].astype(np.float32)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.X = sp.csr_matrix(adata.X, dtype=np.float32)


def as_counts(x) -> sp.csr_matrix:
    """int32 CSR of an integer-valued count matrix; ValueError when it is not counts."""
    x = sp.csr_matrix(x)
    x.sum_duplicates()
    if x.nnz and (x.data.min() < 0 or not np.all(np.mod(x.data, 1) == 0)):
        raise ValueError("X is not raw counts (negative or non-integer values); sc-qc needs sc-prep's output")
    return sp.csr_matrix(x, dtype=np.int32)
