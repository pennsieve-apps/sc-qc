"""
A synthetic sc-prep output that follows the workflow contract (cells/<sample_id>.h5ad and
droplets/<sample_id>.raw.h5), for the unit tests and `make run`:

  python -m tests.make_fixture data/fixture [--layout linked|flat|both] [--seed 7]

Two single-nucleus samples, 300 genes (13 MT-, 10 RPS/RPL), three cell types with their own
marker genes:

- `nuc_a`: 400 nuclei with a raw droplet matrix (the nuclei plus 70,000 empty droplets that
  hold the ambient profile, a mix of all three types; every nucleus carries ~10 % ambient
  counts). Planted: 12 nuclei with ~4 % of the normal counts, 10 with ~30 % mitochondrial
  counts, 20 doublets (the sum of two nuclei of different types).
- `nuc_b`: 300 nuclei, no raw matrix (ambient removal is skipped). Planted: 8 nuclei with ~30 %
  mitochondrial counts.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_GENES = 300
PLANTED = {"nuc_a": {"low_counts": 12, "high_mito": 10, "doublet": 20}, "nuc_b": {"high_mito": 8}}
GENES = [f"MT-G{i}" for i in range(13)] + [f"RPL{i}" for i in range(10)] + [f"GENE{i}" for i in range(N_GENES - 23)]
GENE_IDS = [f"ENSG{i:011d}" for i in range(N_GENES)]


def _profiles(rng: np.random.Generator) -> np.ndarray:
    """Three cell types: a shared base and 30 markers each, normalized to sum 1; ~1 % mitochondrial."""
    base = rng.gamma(1.0, 1.0, N_GENES)
    base[:13] = 0.04
    prof = np.tile(base, (3, 1))
    for t in range(3):
        prof[t, 23 + 30 * t: 53 + 30 * t] *= 12.0
    return prof / prof.sum(1, keepdims=True)


def _draw(rng, prof: np.ndarray, totals: np.ndarray) -> np.ndarray:
    lam = prof * totals[:, None]
    return rng.poisson(rng.gamma(4.0, lam / 4.0)).astype(np.int64)


def make_sample(rng: np.random.Generator, sid: str, n: int, planted: dict, with_raw: bool):
    prof = _profiles(rng)
    types = rng.integers(0, 3, n)
    totals = rng.lognormal(np.log(1500), 0.25, n)
    X = _draw(rng, prof[types], totals)
    label = np.array([""] * n, dtype=object)
    free = list(rng.permutation(n))
    if planted.get("low_counts"):
        idx = [free.pop() for _ in range(planted["low_counts"])]
        X[idx] = rng.binomial(X[idx], 0.04)
        label[idx] = "low_counts"
    if planted.get("high_mito"):
        idx = [free.pop() for _ in range(planted["high_mito"])]
        X[np.ix_(idx, range(13))] = rng.poisson(X[idx, 13:].sum(1, keepdims=True) * 0.43 / 13)
        label[idx] = "high_mito"
    if planted.get("doublet"):
        for _ in range(planted["doublet"]):
            i = free.pop()
            other = rng.choice(np.flatnonzero(types != types[i]))
            X[i] = X[i] + _draw(rng, prof[types[other]][None, :], np.array([totals[other]]))[0]
            label[i] = "doublet"
    ambient = prof.mean(0)
    if with_raw:
        X = X + _draw(rng, ambient[None, :], X.sum(1) * 0.11)
    barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(n)]
    obs = pd.DataFrame({"sample_id": pd.Categorical([sid] * n), "cell_type_truth": pd.Categorical(types.astype(str)),
                        "planted": pd.Categorical(label.astype(str))}, index=[f"{sid}:{b}" for b in barcodes])
    var = pd.DataFrame({"gene_ids": GENE_IDS, "feature_types": "Gene Expression",
                        "mt": [g.startswith("MT-") for g in GENES], "ribo": [g.startswith("RPL") for g in GENES]},
                       index=GENES)
    adata = ad.AnnData(X=sp.csr_matrix(X.astype(np.int32)), obs=obs, var=var)
    a = adata.copy()
    import scanpy as sc

    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
    for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"):
        adata.obs[c] = a.obs[c].to_numpy()
    adata.uns["pennsieve"] = {"source": f"{sid}/", "chemistry": "3p", "organism": "human", "prep_version": "0.1.0"}
    raw = None
    if with_raw:
        n_empty = 70000   # CellBender takes ~20,000 droplets beyond the cells, then learns from the emptier ones
        empty_tot = np.r_[rng.lognormal(np.log(25), 0.4, 8000), rng.lognormal(np.log(9), 0.4, n_empty - 8000)]
        empties = _draw(rng, ambient[None, :].repeat(n_empty, 0), empty_tot)
        raw_x = np.vstack([X, empties])
        raw_bc = barcodes + [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(n_empty)]
        raw = (sp.csr_matrix(raw_x.astype(np.int32)), raw_bc)
    return adata, raw


def write_10x_h5(path: Path, x: sp.csr_matrix, barcodes: list[str]) -> None:
    """A Cell Ranger v3 feature-barcode matrix (HDF5), barcodes x genes in, genes x barcodes CSC on disk."""
    csc = sp.csr_matrix(x).T.tocsc()   # genes x barcodes: CSC column pointers run over barcodes
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        m.create_dataset("data", data=csc.data.astype(np.int32))
        m.create_dataset("indices", data=csc.indices.astype(np.int64))
        m.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array(csc.shape, dtype=np.int32))
        feat = m.create_group("features")
        feat.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
        feat.create_dataset("feature_type", data=np.array([b"Gene Expression"] * N_GENES))
        feat.create_dataset("genome", data=np.array([b"GRCh38"] * N_GENES))
        feat.create_dataset("id", data=np.array(GENE_IDS, dtype="S"))
        feat.create_dataset("name", data=np.array(GENES, dtype="S"))


def make_fixture(root: Path, layout: str = "both", seed: int = 7) -> dict[str, Path]:
    rng = np.random.default_rng(seed)
    linked = root / "linked"
    (linked / "cells").mkdir(parents=True, exist_ok=True)
    (linked / "droplets").mkdir(parents=True, exist_ok=True)
    for sid, n, with_raw in (("nuc_a", 400, True), ("nuc_b", 300, False)):
        adata, raw = make_sample(rng, sid, n, PLANTED[sid], with_raw)
        adata.write_h5ad(linked / "cells" / f"{sid}.h5ad")
        if raw is not None:
            write_10x_h5(linked / "droplets" / f"{sid}.raw.h5", *raw)
    out = {"linked": linked}
    if layout in ("flat", "both"):
        flat = root / "flat"
        flat.mkdir(parents=True, exist_ok=True)
        for p in [*(linked / "cells").iterdir(), *(linked / "droplets").iterdir()]:
            shutil.copy(p, flat / p.name)
        out["flat"] = flat
    if layout == "flat":
        shutil.rmtree(linked)
        out.pop("linked")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    for name, path in make_fixture(args.out, args.layout, args.seed).items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
