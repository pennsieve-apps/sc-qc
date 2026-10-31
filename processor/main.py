"""
sc-qc: per-sample quality control of single-cell / single-nucleus RNA-seq, the step after
sc-prep. For each <sample_id>.h5ad (port `cells`), one at a time:

  1. ambient RNA removal with CellBender when the sample's raw droplet matrix
     <sample_id>.raw.h5 (port `droplets`) is present (AMBIENT),
  2. QC metrics and outlier flags by MADs (counts, genes, top-20 share) and a mitochondrial cap,
  3. doublets with Scrublet on the cells that passed 2,
  4. removal of the flagged cells, normalize_total + log1p, counts in layers.

Writes cells/<sample_id>.h5ad and qc/ (qc_summary.csv, thresholds.csv, <sample_id>_qc.csv,
figures, summary.json); provenance.json in both.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import platform
import shutil
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("sc-qc")
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "h5py"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "sc-qc"
REPO_ROOT = Path(__file__).resolve().parent.parent

DEFAULTS = {
    "AMBIENT": "auto",
    "AMBIENT_EPOCHS": "150",
    "AMBIENT_FPR": "0.01",
    "DOUBLETS": "remove",
    "DOUBLET_RATE": "0.05",
    "MADS": "5",
    "MAX_PCT_MT": "5",
    "TARGET_SUM": "10000",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)
AMBIENT_MODES = ("auto", "cellbender", "none")
DOUBLET_MODES = ("remove", "flag")

# Checked on Crossref and PubMed (see README "Citations").
CITATIONS = {
    "best_practices": {"doi": "10.1038/s41576-023-00586-w", "pmid": "37002403",
                       "text": "Heumos et al. Nature Reviews Genetics 2023"},
    "scanpy": {"doi": "10.1186/s13059-017-1382-0", "pmid": "29409532", "text": "Wolf et al. Genome Biology 2018"},
    "anndata": {"doi": "10.21105/joss.04371", "text": "Virshup et al. Journal of Open Source Software 2024"},
    "scrublet": {"doi": "10.1016/j.cels.2018.11.005", "pmid": "30954476", "text": "Wolock et al. Cell Systems 2019"},
    "cellbender": {"doi": "10.1038/s41592-023-01943-7", "pmid": "37550580", "text": "Fleming et al. Nature Methods 2023"},
}


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {"AMBIENT": raw["AMBIENT"].strip().lower(), "DOUBLETS": raw["DOUBLETS"].strip().lower()}
    if p["AMBIENT"] not in AMBIENT_MODES:
        raise ParameterError(f"AMBIENT must be one of {', '.join(AMBIENT_MODES)}, got {raw['AMBIENT']!r}")
    if p["DOUBLETS"] not in DOUBLET_MODES:
        raise ParameterError(f"DOUBLETS must be one of {', '.join(DOUBLET_MODES)}, got {raw['DOUBLETS']!r}")
    try:
        for k in ("AMBIENT_FPR", "DOUBLET_RATE", "MADS", "MAX_PCT_MT", "TARGET_SUM"):
            p[k] = float(raw[k])
        for k in ("AMBIENT_EPOCHS", "SEED", "THREADS"):
            p[k] = int(raw[k])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if not 1 <= p["AMBIENT_EPOCHS"] <= 1000:
        raise ParameterError("AMBIENT_EPOCHS must be 1-1000")
    if not 0 < p["AMBIENT_FPR"] < 1:
        raise ParameterError("AMBIENT_FPR must be between 0 and 1")
    if not 0 < p["DOUBLET_RATE"] < 0.5:
        raise ParameterError("DOUBLET_RATE must be between 0 and 0.5")
    if p["MADS"] <= 0:
        raise ParameterError("MADS must be > 0")
    if not 0 <= p["MAX_PCT_MT"] <= 100:
        raise ParameterError("MAX_PCT_MT is a percentage, 0-100")
    if p["TARGET_SUM"] <= 0:
        raise ParameterError("TARGET_SUM must be > 0")
    if p["THREADS"] < 0:
        raise ParameterError("THREADS must be >= 0")
    return p


def _fmt(v):
    if isinstance(v, (bool,)) or type(v).__name__ == "bool_":
        return "true" if v else "false"
    if v is None:
        return ""
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.6g}"
    return v


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r.get(c)) for c in columns])


def _clean_json(obj):
    import numpy as np

    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return None if math.isnan(v) or math.isinf(v) else v
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, np.ndarray)):
        return [_clean_json(v) for v in obj]
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    out = {"python": platform.python_version()}
    for name in ("scanpy", "anndata", "cellbender", "torch", "pyro-ppl", "scikit-image", "numpy", "scipy", "pandas",
                 "h5py", "scikit-learn", "matplotlib"):
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


def peak_rss_mb() -> float | None:
    import resource

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return round(max(peak, child) / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)


def cgroup_peak_mb() -> float | None:
    try:
        return round(int(open("/sys/fs/cgroup/memory.peak").read()) / 1024 / 1024, 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def cuda_available() -> bool:
    import torch

    return bool(torch.cuda.is_available())


def ambient_step(sample, adata, raw_counts, params: dict, cuda: bool, threads: int, out_qc: Path) -> dict:
    """CellBender for one sample when it applies. Returns {"status", "note", "counts", "probability", ...}."""
    from processor import ambient as A

    mode = params["AMBIENT"]
    if mode == "none":
        return {"status": "off", "note": "AMBIENT=none"}
    if sample.raw is None:
        log.warning("%s: no raw droplet matrix (<sample_id>.raw.h5); ambient RNA is not removed", sample.sample_id)
        return {"status": "skipped", "note": "no raw droplet matrix"}
    if mode == "auto" and not cuda:
        log.warning("%s: no GPU; ambient RNA is not removed (CellBender takes about an hour per sample on CPU: "
                    "use the gpu variant, or AMBIENT=cellbender to run it on CPU anyway)", sample.sample_id)
        return {"status": "skipped", "note": "no GPU (AMBIENT=auto runs CellBender on a GPU only)"}
    workdir = Path(tempfile.mkdtemp(prefix=f"cellbender-{sample.sample_id}-", dir=os.environ.get("TMPDIR", "/tmp")))
    try:
        run = A.run_cellbender(sample.raw, workdir, epochs=params["AMBIENT_EPOCHS"],
                               fpr=params["AMBIENT_FPR"], cuda=cuda, threads=threads)
        counts, probability = A.corrected_counts(adata, raw_counts, run.output)
        if run.metrics is not None:
            shutil.copy(run.metrics, out_qc / f"{sample.sample_id}_cellbender_metrics.csv")
        return {"status": "cellbender", "note": f"CellBender on {run.device}", "counts": counts,
                "probability": probability, "seconds": round(run.seconds, 1),
                "seconds_per_epoch": run.seconds_per_epoch, "device": run.device}
    except Exception as exc:   # auto: go on without it; cellbender: the user asked for it, the sample fails
        if mode == "cellbender":
            raise
        log.warning("%s: CellBender failed (%s); ambient RNA is not removed", sample.sample_id, exc)
        return {"status": "failed", "note": f"CellBender failed: {exc}"}
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def process_sample(sample, params: dict, cuda: bool, threads: int, out_cells: Path, out_qc: Path) -> dict:
    import anndata as ad
    import numpy as np
    import pandas as pd

    from processor import figures as F
    from processor import qc as Q

    sid = sample.sample_id
    t0 = time.monotonic()
    adata = ad.read_h5ad(sample.path)
    meta = dict(adata.uns["pennsieve"]) if isinstance(adata.uns.get("pennsieve"), dict) else {}
    if "qc" in meta:
        raise ValueError("already through sc-qc (uns['pennsieve']['qc'] is set); give sc-qc sc-prep's output")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"empty matrix ({adata.n_obs} cells x {adata.n_vars} genes)")
    if "gene_ids" not in adata.var:
        adata.var["gene_ids"] = adata.var_names
    raw_counts = Q.as_counts(adata.layers["counts"] if "counts" in adata.layers else adata.X)
    Q.mito_ribo(adata)
    n_before = adata.n_obs

    t1 = time.monotonic()
    amb = ambient_step(sample, adata, raw_counts, params, cuda, threads, out_qc)
    t_ambient = time.monotonic() - t1
    counts = amb.get("counts", raw_counts)
    ambient_ran = amb["status"] == "cellbender"

    metrics = Q.qc_metrics(counts, adata.var, adata.obs_names)
    flags, thr = Q.flag(metrics, params)
    per_cell = metrics.copy()
    if ambient_ran:
        raw_total = np.asarray(raw_counts.sum(axis=1)).ravel().astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            per_cell["ambient_fraction"] = np.where(raw_total > 0, 1.0 - metrics["total_counts"] / raw_total, np.nan)
        per_cell["total_counts_raw"] = raw_total
        per_cell["cellbender_cell_probability"] = amb["probability"]
        flags["ambient_empty"] = np.nan_to_num(amb["probability"], nan=1.0) < Q.CELL_PROBABILITY_CUTOFF

    t2 = time.monotonic()
    passed = ~flags.any(axis=1).to_numpy()
    per_cell["doublet_score"] = np.nan
    per_cell["predicted_doublet"] = False
    dbl = {"threshold": None, "simulated_scores": np.array([]), "n_predicted": 0, "status": "skipped"}
    if passed.sum() >= Q.MIN_CELLS_FOR_SCRUBLET:
        scored, dbl = Q.scrublet(counts[np.flatnonzero(passed)], adata.obs_names[passed], params)
        dbl["status"] = "ok"
        per_cell.loc[scored.index, "doublet_score"] = scored["doublet_score"]
        per_cell.loc[scored.index, "predicted_doublet"] = scored["predicted_doublet"]
        if params["DOUBLETS"] == "remove":
            flags["doublet"] = per_cell["predicted_doublet"].to_numpy(dtype=bool)
    else:
        log.warning("%s: %d cells passed QC, fewer than %d; no doublet detection", sid, int(passed.sum()),
                    Q.MIN_CELLS_FOR_SCRUBLET)
    t_doublets = time.monotonic() - t2

    reasons = Q.reasons_of(flags)
    kept = (reasons == "").to_numpy()
    n_after = int(kept.sum())
    if n_after < 2:
        raise ValueError(f"{n_after} of {n_before} cells would be kept; check the thresholds")

    # The output h5ad: kept cells in input order, every gene, counts and metrics after ambient removal.
    t3 = time.monotonic()
    idx = np.flatnonzero(kept)
    out = adata[idx].copy()
    out.layers["counts"] = counts[idx]
    if ambient_ran:
        out.layers["counts_raw"] = raw_counts[idx]
        out.obs["ambient_fraction"] = per_cell["ambient_fraction"].to_numpy(dtype=np.float32)[idx]
        out.obs["cellbender_cell_probability"] = per_cell["cellbender_cell_probability"].to_numpy(dtype=np.float32)[idx]
    for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"):
        out.obs[c] = per_cell[c].to_numpy()[idx].astype(np.float32 if c.startswith("pct") else np.float64)
    out.obs["doublet_score"] = per_cell["doublet_score"].to_numpy(dtype=np.float32)[idx]
    out.obs["predicted_doublet"] = per_cell["predicted_doublet"].to_numpy(dtype=bool)[idx]
    Q.normalize(out, params["TARGET_SUM"])
    meta["qc"] = {"version": app_version(), "parameters": params, "ambient": amb["status"],
                  "ambient_note": amb["note"], "doublet_threshold": dbl["threshold"], "n_cells_before": n_before,
                  "n_cells_after": n_after}
    meta["normalized"] = True
    meta["target_sum"] = params["TARGET_SUM"]
    out.uns["pennsieve"] = meta
    out.write_h5ad(out_cells / f"{sid}.h5ad", compression="gzip")
    t_write = time.monotonic() - t3

    table = per_cell.copy()
    table.insert(0, "obs_id", table.index)
    table.insert(0, "sample_id", sid)
    table["kept"] = kept
    table["reason"] = reasons.to_numpy()
    columns = [c for c in ("sample_id", "obs_id", "kept", "reason", "total_counts", "n_genes_by_counts", "pct_counts_mt",
                           "pct_counts_ribo", "pct_counts_in_top_20_genes", "total_counts_raw", "ambient_fraction",
                           "cellbender_cell_probability", "doublet_score", "predicted_doublet") if c in table]
    write_csv(out_qc / f"{sid}_qc.csv", table[columns].to_dict("records"), columns)

    t4 = time.monotonic()
    figures = []
    try:
        lines = {r["metric"]: (r["lower"], r["upper"]) for r in thr.rows}
        figures.append(F.qc_overview(sid, metrics, kept, lines, out_qc / f"{sid}_qc.png"))
        if dbl["status"] == "ok":
            figures.append(F.doublets(sid, per_cell["doublet_score"].dropna().to_numpy(), dbl["simulated_scores"],
                                      dbl["threshold"], out_qc / f"{sid}_doublets.png"))
        if ambient_ran:
            figures.append(F.ambient(sid, per_cell["ambient_fraction"].to_numpy(), out_qc / f"{sid}_ambient.png"))
    except Exception as exc:   # figures are a convenience; the tables and the h5ad stand without them
        log.exception("%s: figures failed: %s", sid, exc)
    t_fig = time.monotonic() - t4

    n_by = {r: int(flags[r].sum()) for r in Q.REASONS}
    row = {"sample_id": sid, "n_cells_before": n_before, "n_cells_after": n_after,
           "pct_removed": 100.0 * (n_before - n_after) / n_before, "n_genes": adata.n_vars,
           "ambient_method": "cellbender" if ambient_ran else "none", "ambient_status": amb["status"],
           "ambient_note": amb["note"], "median_ambient_fraction": float(np.nanmedian(per_cell["ambient_fraction"])) if ambient_ran else None,
           "median_counts": float(np.median(per_cell["total_counts"][kept])),
           "median_genes": float(np.median(per_cell["n_genes_by_counts"][kept])),
           "median_pct_mt": float(np.median(per_cell["pct_counts_mt"][kept])),
           "doublet_threshold": dbl["threshold"], "n_predicted_doublets": int(per_cell["predicted_doublet"].sum()),
           **{f"n_removed_{r}": n_by[r] for r in Q.REASONS}, "status": "ok"}
    thr_rows = [{"sample_id": sid, **r} for r in thr.rows]
    timings = {"ambient": t_ambient, "doublets": t_doublets, "normalize_write": t_write, "figures": t_fig,
               "total": time.monotonic() - t0}
    log.info("%s: %d -> %d cells (%.1f %% removed: %s), %d genes, ambient %s%s, doublets %d (threshold %s), "
             "%.1fs (ambient %.1f, doublets %.1f, write %.1f, figures %.1f)",
             sid, n_before, n_after, row["pct_removed"], ", ".join(f"{k} {v}" for k, v in n_by.items() if v) or "none",
             adata.n_vars, amb["status"],
             f" (median {100 * row['median_ambient_fraction']:.1f} %)" if ambient_ran else "",
             row["n_predicted_doublets"], "none" if dbl["threshold"] is None else f"{dbl['threshold']:.2f}",
             timings["total"], t_ambient, t_doublets, t_write, t_fig)
    return {"row": row, "thresholds": thr_rows, "figures": [f.name for f in figures], "ambient": amb,
            "timings": {k: round(v, 2) for k, v in timings.items()}}


SUMMARY_COLUMNS = ["sample_id", "n_cells_before", "n_cells_after", "pct_removed",
                   *(f"n_removed_{r}" for r in ("ambient_empty", "low_counts", "high_counts", "low_genes", "high_genes",
                                                "high_top20", "high_mito", "doublet")),
                   "ambient_method", "ambient_status", "ambient_note", "median_ambient_fraction", "n_genes",
                   "median_counts", "median_genes", "median_pct_mt", "doublet_threshold", "n_predicted_doublets", "status"]
THRESHOLD_COLUMNS = ["sample_id", "metric", "scale", "median", "mad", "lower", "upper", "rule"]


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    threads = int(os.environ.get("SC_QC_THREADS") or params["THREADS"] or 2)

    from processor.inputs import discover

    found = discover(input_dir)
    for rel, reason in found.skipped.items():
        log.warning("%s: skipped, %s", rel, reason)
    if not found.samples:
        log.error("no <sample_id>.h5ad found under %s", input_dir)
        return 1
    cuda = params["AMBIENT"] != "none" and any(s.raw for s in found.samples) and cuda_available()
    log.info("parameters: %s; threads %d (%s); device %s; %d samples, %d with a raw droplet matrix", params, threads,
             os.environ.get("SC_QC_THREADS_SOURCE", "environment"), "cuda" if cuda else "cpu", len(found.samples),
             sum(1 for s in found.samples if s.raw))
    out_cells, out_qc = output_dir / "cells", output_dir / "qc"
    out_cells.mkdir(parents=True, exist_ok=True)
    out_qc.mkdir(parents=True, exist_ok=True)

    rows, thr_rows, statuses, figures = [], [], [], []
    for s in found.samples:
        t = time.monotonic()
        status = {"sample_id": s.sample_id, "file": str(s.path.relative_to(input_dir)),
                  "raw": str(s.raw.relative_to(input_dir)) if s.raw else None}
        try:
            result = process_sample(s, params, cuda, threads, out_cells, out_qc)
            rows.append(result["row"])
            thr_rows.extend(result["thresholds"])
            figures.extend(result["figures"])
            amb = {k: v for k, v in result["ambient"].items() if k not in ("counts", "probability")}
            status.update(status="ok", n_cells_before=result["row"]["n_cells_before"],
                          n_cells_after=result["row"]["n_cells_after"], ambient=amb,
                          timings_seconds=result["timings"], peak_rss_mb=peak_rss_mb())
        except Exception as exc:   # one bad sample is skipped, not fatal
            status.update(status="failed", error=f"{type(exc).__name__}: {exc}", seconds=round(time.monotonic() - t, 2))
            log.exception("%s: failed: %s", s.sample_id, exc)
            rows.append({"sample_id": s.sample_id, "status": f"failed: {type(exc).__name__}: {exc}"})
        statuses.append(status)
    ok = [r for r in rows if r.get("status") == "ok"]

    write_csv(out_qc / "qc_summary.csv", rows, SUMMARY_COLUMNS)
    write_csv(out_qc / "thresholds.csv", thr_rows, THRESHOLD_COLUMNS)
    from processor.qc import REASONS

    ambient_ran = any(r["ambient_method"] == "cellbender" for r in ok)
    summary = {
        "n_samples": len(statuses), "n_ok": len(ok), "n_failed": len(statuses) - len(ok),
        "n_cells_before": sum(r["n_cells_before"] for r in ok), "n_cells_after": sum(r["n_cells_after"] for r in ok),
        "n_removed_by_reason": {k: sum(r[f"n_removed_{k}"] for r in ok) for k in REASONS},
        "n_samples_ambient_removed": sum(1 for r in ok if r["ambient_method"] == "cellbender"),
        "device": "cuda" if cuda else "cpu",
        "samples": [{k: r.get(k) for k in SUMMARY_COLUMNS} for r in rows],
        "parameters": params, "figures": sorted(figures),
    }
    write_json(out_qc / "summary.json", summary)
    citations = [CITATIONS[k] for k in ("best_practices", "scanpy", "anndata", "scrublet")]
    if ambient_ran:
        citations.append(CITATIONS["cellbender"])
    provenance = {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": params,
        "threads": {"n": threads, "source": os.environ.get("SC_QC_THREADS_SOURCE", "environment")},
        "device": "cuda" if cuda else "cpu",
        "method": {
            "ambient": "cellbender remove-background (its own cell-count and droplet heuristics, --epochs AMBIENT_EPOCHS, "
                       "--fpr AMBIENT_FPR) on <sample_id>.raw.h5; counts of the sc-prep cells replaced by CellBender's; "
                       "cells with CellBender cell probability < 0.5 removed (ambient_empty)",
            "metrics": "scanpy.pp.calculate_qc_metrics(qc_vars=['mt','ribo'], percent_top=[20]) on the counts after ambient removal",
            "outliers": "log1p(total_counts), log1p(n_genes_by_counts) outside median ± MADS*MAD; pct_counts_in_top_20_genes > "
                        "median + MADS*MAD; pct_counts_mt > MAX_PCT_MT",
            "doublets": "scanpy.pp.scrublet(expected_doublet_rate=DOUBLET_RATE, random_state=SEED) on the counts of the cells "
                        "that passed the outlier flags; automatic threshold",
            "normalization": "layers['counts'] = counts (after ambient removal), layers['counts_raw'] = sc-prep's counts when "
                             "CellBender ran; X = log1p(normalize_total(target_sum=TARGET_SUM)), float32 CSR; every gene kept",
        },
        "citations": citations,
        "inputs": found.seen, "skipped": found.skipped,
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": peak_rss_mb(), "cgroup_memory_peak_mb": cgroup_peak_mb(),
        "samples": statuses, "n_ok": len(ok), "n_failed": len(statuses) - len(ok),
    }
    write_json(out_qc / "provenance.json", provenance)
    write_json(out_cells / "provenance.json", provenance)
    log.info("done: %d samples (%d failed), %d -> %d cells, %.1fs, peak RSS %s MB, cgroup peak %s MB -> %s",
             len(ok), len(statuses) - len(ok), summary["n_cells_before"], summary["n_cells_after"],
             time.monotonic() - t0, provenance["peak_rss_mb"], provenance["cgroup_memory_peak_mb"], output_dir)
    if not ok:
        log.error("every sample failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
