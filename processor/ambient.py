"""
Ambient RNA removal with CellBender remove-background (Fleming et al. 2023), run per sample as
its own process on the sample's raw droplet matrix (<sample_id>.raw.h5).

CellBender learns the ambient ("soup") profile from the empty droplets and removes the
expected background counts from every cell. Its corrected counts replace the counts of the
sc-prep cells (matched by barcode and Ensembl id); the per-cell ambient fraction is
1 - corrected / original total counts, and cellbender_cell_probability is CellBender's own
posterior that the droplet holds a cell.
"""
from __future__ import annotations

import logging
import re
import subprocess
import sys
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp

log = logging.getLogger("sc-qc")

EPOCH_TIME = re.compile(r"\(([0-9.]+) seconds per epoch\)")


@dataclass
class CellBenderRun:
    output: Path              # CellBender's output.h5
    metrics: Path | None      # its output_metrics.csv
    seconds: float
    seconds_per_epoch: float | None
    device: str


def run_cellbender(raw: Path, workdir: Path, *, epochs: int, fpr: float, cuda: bool,
                   threads: int) -> CellBenderRun:
    """cellbender remove-background in `workdir`, which is also its working directory: its
    checkpoint (~1 GB, written once after training) is written and read back there, a fresh folder
    per sample, so no checkpoint from anywhere else is ever loaded.
    Raises RuntimeError with the end of its log when it fails."""
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "output.h5"
    cmd = [sys.executable, "-m", "cellbender.base_cli", "remove-background", "--input", str(raw.resolve()), "--output", str(out),
           "--epochs", str(epochs), "--fpr", str(fpr),
           "--cpu-threads", str(threads), "--checkpoint-mins", "100000"]   # no mid-run checkpoints: never resumed
    if cuda:
        cmd.append("--cuda")
    t0 = time.monotonic()
    tail: deque[str] = deque(maxlen=3)
    per_epoch = None
    with open(workdir / "cellbender.log", "w") as log_file, subprocess.Popen(
            cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
        for line in proc.stdout:
            log_file.write(line)
            tail.append(line.rstrip())
            m = EPOCH_TIME.search(line)
            if m and per_epoch is None:
                per_epoch = float(m.group(1))
                log.info("%s: CellBender %.1f s per epoch on %s, %d epochs: ~%.0f min of training (then the posterior)",
                         raw.name, per_epoch, "cuda" if cuda else "cpu", epochs, per_epoch * epochs / 60.0)
        code = proc.wait()
    if code != 0 or not out.exists():
        raise RuntimeError(f"cellbender exited with {code}: " + " | ".join(tail))
    metrics = workdir / "output_metrics.csv"
    return CellBenderRun(out, metrics if metrics.exists() else None, time.monotonic() - t0, per_epoch,
                         "cuda" if cuda else "cpu")


def corrected_counts(adata, raw_counts: sp.csr_matrix, output: Path) -> tuple[sp.csr_matrix, np.ndarray]:
    """CellBender's counts for the cells and genes of `adata` (obs_names <sample_id>:<barcode>,
    var["gene_ids"]), as int32 CSR in adata's order, and each cell's CellBender cell probability.
    A cell CellBender did not analyse (fewer counts than its low-count cutoff) keeps its counts and
    gets probability 0: CellBender saw nothing in it that looked like a cell."""
    from cellbender.remove_background.downstream import anndata_from_h5

    cb = anndata_from_h5(str(output), analyzed_barcodes_only=True)
    barcodes = np.array([name.split(":", 1)[-1] for name in adata.obs_names])
    rows = cb.obs_names.get_indexer(barcodes)
    first = {gid: i for i, gid in reversed(list(enumerate(cb.var["gene_id"].astype(str))))}
    cols = np.array([first.get(g, -1) for g in adata.var["gene_ids"].astype(str)])
    if (cols < 0).any():
        raise ValueError(f"{int((cols < 0).sum())} of {adata.n_vars} genes are not in the raw matrix")
    analysed = rows >= 0
    corrected = sp.csr_matrix(cb.X)[rows[analysed]][:, cols]
    x = sp.vstack([corrected, raw_counts[~analysed]]).tocsr()
    order = np.argsort(np.r_[np.flatnonzero(analysed), np.flatnonzero(~analysed)], kind="stable")
    x = sp.csr_matrix(x[order], dtype=np.float64)
    x.data = np.rint(x.data)          # integer counts (CellBender's estimate is integer; rounding makes it certain)
    x = sp.csr_matrix(x, dtype=np.int32)
    x.eliminate_zeros()
    probability = np.zeros(adata.n_obs)
    probability[analysed] = cb.obs["cell_probability"].to_numpy(dtype=float)[rows[analysed]]
    return x, probability
