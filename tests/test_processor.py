"""The synthetic fixture through the processor (both input layouts, with and without CellBender),
the contract on the written h5ad, the modes and the parameter checks."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from processor import main as app
from processor import qc as Q
from tests.make_fixture import PLANTED, make_fixture

QC_FILES = {"qc_summary.csv", "thresholds.csv", "summary.json", "provenance.json"}


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


@pytest.fixture(scope="session")
def cellbender_run(fixture_root, tmp_path_factory) -> Path:
    """One CellBender run (5 epochs, a few seconds on the fixture) shared by the tests that need it."""
    out = tmp_path_factory.mktemp("cellbender") / "out"
    with pytest.MonkeyPatch.context() as mp:
        assert run(fixture_root["linked"], out, mp, AMBIENT="cellbender", AMBIENT_EPOCHS="5") == 0
    return out


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, str(value))
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def source(root: Path, sid: str) -> ad.AnnData:
    linked = root / "cells" / f"{sid}.h5ad"
    return ad.read_h5ad(linked if linked.exists() else root / f"{sid}.h5ad")


def check_contract(src: ad.AnnData, out: ad.AnnData) -> None:
    # Cells: a subset of the input in input order; every gene kept.
    assert list(out.obs_names) == list(src.obs_names[src.obs_names.isin(out.obs_names)])
    assert list(out.var_names) == list(src.var_names)
    counts = out.layers["counts"]
    assert sp.isspmatrix_csr(counts) and np.issubdtype(counts.dtype, np.integer)
    if "counts_raw" in out.layers:
        raw = out.layers["counts_raw"]
        assert (raw != src[out.obs_names].X).nnz == 0
        assert (counts > raw).nnz == 0          # ambient removal only takes counts away
    else:
        assert (counts != src[out.obs_names].X).nnz == 0
    assert sp.isspmatrix_csr(out.X) and out.X.dtype == np.float32
    totals = np.asarray(np.expm1(out.X).sum(axis=1)).ravel()
    assert np.allclose(totals, 10000, rtol=1e-3)
    assert np.allclose(out.obs["total_counts"], np.asarray(counts.sum(axis=1)).ravel())
    for c in ("sample_id", "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "doublet_score",
              "predicted_doublet"):
        assert c in out.obs
    assert out.uns["pennsieve"]["normalized"] and out.uns["pennsieve"]["target_sum"] == 10000
    assert out.uns["pennsieve"]["qc"]["n_cells_after"] == out.n_obs


def removed_by(qc_csv: Path) -> dict[str, list[str]]:
    return {r["obs_id"]: r["reason"].split(";") for r in rows(qc_csv) if r["kept"] == "false"}


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_cells_are_removed(fixture_root, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixture_root[layout], out, monkeypatch) == 0     # auto, no GPU here: no CellBender
    assert sorted(p.name for p in out.iterdir()) == ["cells", "qc"]
    assert QC_FILES <= {p.name for p in (out / "qc").iterdir()}
    summary = {r["sample_id"]: r for r in rows(out / "qc" / "qc_summary.csv")}
    assert set(summary) == set(PLANTED) and all(r["status"] == "ok" and r["ambient_method"] == "none" for r in summary.values())
    for sid, planted in PLANTED.items():
        src = source(fixture_root[layout], sid)
        res = ad.read_h5ad(out / "cells" / f"{sid}.h5ad")
        check_contract(src, res)
        got = removed_by(out / "qc" / f"{sid}_qc.csv")
        assert len(rows(out / "qc" / f"{sid}_qc.csv")) == src.n_obs
        for why, n in planted.items():
            ids = set(src.obs_names[src.obs["planted"] == why])
            assert len(ids) == n
            hit = sum(1 for oid in ids if oid in got and why in got[oid])
            assert hit >= (0.8 * n if why == "doublet" else n), f"{sid}: {hit} of {n} {why} removed for that reason"
        planted_ids = set(src.obs_names[src.obs["planted"] != ""])
        assert len(set(got) - planted_ids) <= 0.03 * src.n_obs, f"{sid}: too many unplanted cells removed"
        assert res.n_obs == src.n_obs - len(got) == int(summary[sid]["n_cells_after"])
        assert not res.obs["predicted_doublet"].any()
        assert (out / "qc" / f"{sid}_qc.png").exists() and (out / "qc" / f"{sid}_doublets.png").exists()
    js = json.loads((out / "qc" / "summary.json").read_text())
    assert js["n_ok"] == 2 and js["n_samples_ambient_removed"] == 0
    prov = json.loads((out / "qc" / "provenance.json").read_text())
    assert prov["app"] == "sc-qc" and prov["version"] == "0.1.0" and prov["device"] == "cpu"
    dois = {c["doi"] for c in prov["citations"]}
    assert "10.1016/j.cels.2018.11.005" in dois and "10.1038/s41592-023-01943-7" not in dois
    assert all(set(c) <= {"doi", "pmid", "text"} and c["text"] for c in prov["citations"])
    assert {"app", "version", "tools", "parameters", "inputs", "started", "finished", "samples", "peak_rss_mb"} <= set(prov)
    assert (out / "cells" / "provenance.json").exists()


def test_cellbender_removes_ambient(fixture_root, cellbender_run):
    out = cellbender_run
    summary = {r["sample_id"]: r for r in rows(out / "qc" / "qc_summary.csv")}
    assert summary["nuc_a"]["ambient_method"] == "cellbender" and summary["nuc_b"]["ambient_method"] == "none"
    assert 0 < float(summary["nuc_a"]["median_ambient_fraction"]) < 0.5
    res = ad.read_h5ad(out / "cells" / "nuc_a.h5ad")
    src = source(fixture_root["linked"], "nuc_a")
    check_contract(src, res)
    assert "counts_raw" in res.layers and (res.obs["ambient_fraction"] > 0).mean() > 0.9
    assert res.layers["counts"].sum() < res.layers["counts_raw"].sum()
    # The near-empty planted nuclei: CellBender calls them empty droplets.
    got = removed_by(out / "qc" / "nuc_a_qc.csv")
    low = set(src.obs_names[src.obs["planted"] == "low_counts"])
    assert all("ambient_empty" in got[oid] for oid in low)
    assert {"nuc_a_ambient.png", "nuc_a_cellbender_metrics.csv"} <= {p.name for p in (out / "qc").iterdir()}
    assert "counts_raw" not in ad.read_h5ad(out / "cells" / "nuc_b.h5ad").layers
    prov = json.loads((out / "qc" / "provenance.json").read_text())
    assert "10.1038/s41592-023-01943-7" in {c["doi"] for c in prov["citations"]}
    assert prov["samples"][0]["ambient"]["status"] == "cellbender"


def test_doublets_flagged_not_removed(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, DOUBLETS="flag", AMBIENT="none") == 0
    res = ad.read_h5ad(out / "cells" / "nuc_a.h5ad")
    src = source(fixture_root["linked"], "nuc_a")
    doublets = set(src.obs_names[src.obs["planted"] == "doublet"])
    kept = doublets & set(res.obs_names)        # a doublet can still go for its high counts
    assert len(kept) >= 0.8 * len(doublets)
    assert res.obs.loc[list(kept), "predicted_doublet"].mean() >= 0.8
    assert not any("doublet" in why for why in removed_by(out / "qc" / "nuc_a_qc.csv").values())


def test_mito_cap(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, MAX_PCT_MT="50", AMBIENT="none") == 0
    assert all(int(r["n_removed_high_mito"]) == 0 for r in rows(out / "qc" / "qc_summary.csv"))


@pytest.mark.parametrize("params", [{"AMBIENT": "soupx"}, {"DOUBLETS": "drop"}, {"MADS": "0"}, {"MAX_PCT_MT": "150"},
                                    {"TARGET_SUM": "0"}, {"AMBIENT_EPOCHS": "x"}, {"DOUBLET_RATE": "0.9"}])
def test_bad_parameters(fixture_root, tmp_path, monkeypatch, params):
    assert run(fixture_root["linked"], tmp_path / "out", monkeypatch, **params) == 2


def test_a_bad_sample_is_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixture_root["linked"], src)
    (src / "cells" / "nuc_b.h5ad").write_bytes(b"not an h5ad")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch, AMBIENT="none") == 0
    prov = json.loads((out / "qc" / "provenance.json").read_text())
    assert prov["n_failed"] == 1 and prov["n_ok"] == 1
    assert not (out / "cells" / "nuc_b.h5ad").exists()


def test_qc_output_is_refused(cellbender_run, tmp_path, monkeypatch):
    """sc-qc's own output is not raw counts in X: the sample fails rather than being filtered twice."""
    assert run(cellbender_run / "cells", tmp_path / "out", monkeypatch, AMBIENT="none") == 1


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_links_outside_input_are_skipped(fixture_root, tmp_path):
    from processor.inputs import discover

    outside = tmp_path / "elsewhere"
    outside.mkdir()
    shutil.copy(fixture_root["flat"] / "nuc_b.h5ad", outside / "nuc_c.h5ad")
    src = tmp_path / "in"
    (src / "cells" / "sub").mkdir(parents=True)
    shutil.copy(fixture_root["flat"] / "nuc_a.h5ad", src / "cells" / "nuc_a.h5ad")
    (src / "cells" / "sub" / "nuc_c.h5ad").symlink_to(outside / "nuc_c.h5ad")
    found = discover(src)
    assert [s.sample_id for s in found.samples] == ["nuc_a"]
    assert found.skipped == {"cells/sub/nuc_c.h5ad": "links outside INPUT_DIR"}


def test_mad_flags():
    rng = np.random.default_rng(0)
    n = 1000
    counts = np.exp(rng.normal(8, 0.2, n))
    counts[:5] = 30
    counts[5:8] = 1e6
    metrics = pd.DataFrame({"total_counts": counts, "n_genes_by_counts": counts / 3,
                            "pct_counts_in_top_20_genes": rng.normal(30, 2, n),
                            "pct_counts_mt": np.clip(rng.normal(1, 0.5, n), 0, None)}, index=[str(i) for i in range(n)])
    metrics.loc["10", "pct_counts_mt"] = 12
    flags, thr = Q.flag(metrics, {"MADS": 5.0, "MAX_PCT_MT": 5.0})
    assert flags["low_counts"][:5].all() and flags["high_counts"][5:8].all() and flags.loc["10", "high_mito"]
    med = np.median(np.log1p(counts))
    mad = np.median(np.abs(np.log1p(counts) - med))
    assert thr.rows[0]["lower"] == pytest.approx(np.expm1(med - 5 * mad))
    assert flags.any(axis=1).sum() <= 15
    assert Q.reasons_of(flags)["0"] == "low_counts;low_genes"
