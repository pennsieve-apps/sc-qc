"""
PNG figures per sample (150 dpi), in the colour system of the other workflow apps: a blue for
kept cells, the reserved critical red for removed cells and thresholds, a second categorical
slot (orange) for simulated doublets, muted inks for chrome.

- <sample_id>_qc.png: total counts, genes and % mitochondrial counts before (every input cell)
  and after QC (kept cells) as violins with the thresholds, and counts vs genes per cell
  coloured kept / removed.
- <sample_id>_doublets.png: Scrublet doublet scores of the observed cells and of the simulated
  doublets, with the threshold.
- <sample_id>_ambient.png: per-cell ambient fraction removed by CellBender (when it ran).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
CRITICAL = "#d03b3b"
KEPT = "#2a78d6"
SIMULATED = "#eb6834"


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": INK, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
    })
    return plt


def _title(fig, text: str) -> None:
    fig.suptitle(text, x=0.01, ha="left", fontsize=11, fontweight="bold")


def _violin(ax, before: np.ndarray, after: np.ndarray, label: str, lines=(), log=False):
    data = [np.log10(np.maximum(v, 1)) if log else v for v in (before, after)]
    data = [d if d.size > 1 else np.r_[d, d, 0.0][:2] for d in data]
    parts = ax.violinplot(data, positions=[0, 1], widths=0.8, showextrema=False, showmedians=True)
    for body, colour in zip(parts["bodies"], (MUTED, KEPT)):
        body.set_facecolor(colour)
        body.set_edgecolor(SURFACE)
        body.set_alpha(0.75)
    parts["cmedians"].set_color(INK)
    for value in lines:
        if value is not None and np.isfinite(value):
            ax.axhline(np.log10(max(value, 1)) if log else value, color=CRITICAL, linewidth=1.1, linestyle=(0, (4, 2)))
    ax.set_xticks([0, 1], [f"before\n{before.size:,}", f"after\n{after.size:,}"])
    ax.set_ylabel(f"log10 {label}" if log else label)
    ax.set_title(label, loc="left")
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def qc_overview(sample_id: str, metrics, kept: np.ndarray, thresholds: dict, path: Path) -> Path:
    """metrics: the per-cell QC table of every input cell; thresholds: metric -> (lower, upper)."""
    plt = _plt()
    fig, axes = plt.subplots(1, 4, figsize=(14.5, 3.6), gridspec_kw={"width_ratios": [1, 1, 1, 1.4]})
    counts = metrics["total_counts"].to_numpy(dtype=float)
    genes = metrics["n_genes_by_counts"].to_numpy(dtype=float)
    mt = metrics["pct_counts_mt"].to_numpy(dtype=float)
    _violin(axes[0], counts, counts[kept], "total counts", thresholds.get("total_counts", ()), log=True)
    _violin(axes[1], genes, genes[kept], "genes detected", thresholds.get("n_genes_by_counts", ()), log=True)
    _violin(axes[2], mt, mt[kept], "% mitochondrial counts", thresholds.get("pct_counts_mt", ()))
    ax = axes[3]
    size = float(np.clip(3000.0 / max(len(counts), 1), 1.0, 12.0))
    ax.scatter(np.maximum(counts[kept], 1), np.maximum(genes[kept], 1), s=size, c=KEPT, linewidths=0,
               rasterized=True, label=f"kept ({kept.sum():,})")
    ax.scatter(np.maximum(counts[~kept], 1), np.maximum(genes[~kept], 1), s=size, c=CRITICAL, linewidths=0,
               rasterized=True, label=f"removed ({(~kept).sum():,})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total counts")
    ax.set_ylabel("genes detected")
    ax.set_title("counts vs genes per cell", loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper left", markerscale=max(1.0, 4.0 / size))
    pct = 100.0 * (~kept).sum() / max(len(kept), 1)
    _title(fig, f"{sample_id} · QC: {len(kept):,} → {kept.sum():,} cells ({pct:.1f} % removed)")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def doublets(sample_id: str, observed: np.ndarray, simulated: np.ndarray, threshold, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    bins = np.linspace(0, max(float(np.max(observed, initial=0)), float(np.max(simulated, initial=0)), 0.1), 50)
    ax.hist(simulated, bins=bins, color=SIMULATED, alpha=0.6, density=True, linewidth=0,
            label=f"simulated doublets ({simulated.size:,})")
    ax.hist(observed, bins=bins, color=KEPT, alpha=0.75, density=True, linewidth=0, label=f"cells ({observed.size:,})")
    if threshold is not None:
        ax.axvline(threshold, color=CRITICAL, linewidth=1.2, linestyle=(0, (4, 2)))
        ax.set_title(f"threshold {threshold:.2f}: {int((observed > threshold).sum()):,} cells called doublets", loc="left")
    else:
        ax.set_title("no threshold found (simulated scores not bimodal): nothing called a doublet", loc="left")
    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("density")
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    _title(fig, f"{sample_id} · doublets")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def ambient(sample_id: str, fraction: np.ndarray, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    values = fraction[np.isfinite(fraction)]
    ax.hist(values, bins=np.linspace(0, max(float(np.max(values, initial=0)), 0.05), 50), color=KEPT, linewidth=0)
    med = float(np.median(values)) if values.size else float("nan")
    if np.isfinite(med):
        ax.axvline(med, color=INK_2, linewidth=1.0)
        ax.text(med, ax.get_ylim()[1] * 0.97, f" median {100 * med:.1f} %", ha="left", va="top", fontsize=8, color=INK_2)
    ax.set_xlabel("fraction of the cell's counts removed as ambient RNA")
    ax.set_ylabel("cells")
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    _title(fig, f"{sample_id} · ambient RNA (CellBender)")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
