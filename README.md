# sc-qc

Quality control of single-cell and single-nucleus RNA-seq samples, the step after `sc-prep` in the
single-cell workflows:

- **A. Single-nucleus Atlas Mapping & Condition Comparison:**
  `sc-prep → sc-qc → sc-integrate → sc-annotate → {sc-markers, st-compare} → sc-report`
- **B. Multiome Gene Regulation:** the RNA half goes through A's chain, which gives the cell types.

Not every barcode Cell Ranger calls a cell is one good cell. Droplets also pick up **ambient RNA**,
transcripts floating in the suspension from lysed cells; in nuclei preparations of nervous tissue
(dorsal root ganglia, spinal cord) the ambient pool is rich in neuronal transcripts, so glial nuclei
appear to express neuronal genes. Some droplets hold **two cells** (doublets), which look like
intermediate cell types. Some hold **damaged cells or debris**: few genes, or a transcriptome
dominated by a few genes or by mitochondrial reads. Left in, all three create false cell types and
false "differential expression". sc-qc removes them sample by sample and leaves every sample
normalized the same way, so `sc-integrate` and everything after it starts from clean cells.

For the computational analyst: per sample, streaming one `.h5ad` at a time, (1) CellBender
`remove-background` on the raw droplet matrix (when linked, and on a GPU unless asked otherwise);
(2) `scanpy.pp.calculate_qc_metrics` and median-absolute-deviation outlier flags as in the
single-cell best-practices guide (Heumos et al. 2023), with a fixed mitochondrial cap suited to
nuclei; (3) Scrublet (`scanpy.pp.scrublet`) on the cells that passed (2); (4) removal of the
flagged cells and `normalize_total` + `log1p`, with the counts in `layers`. Cells are only removed,
never reordered; every gene is kept.

## Method

1. **Ambient RNA** (`AMBIENT`). When the sample's raw (unfiltered) droplet matrix is linked on the
   `droplets` port, CellBender `remove-background` (Fleming et al. 2023) learns the ambient profile
   from the empty droplets and each droplet's probability of holding a cell, and removes the
   expected background counts. It runs as its own process in a fresh temporary folder, with its own
   cell-count and droplet heuristics (no `--expected-cells`: those heuristics are what its authors
   recommend, and sc-prep may have subsampled the cells), `--epochs AMBIENT_EPOCHS` (150) and
   `--fpr AMBIENT_FPR` (0.01). The corrected counts of sc-prep's cells (matched by barcode and
   Ensembl id) replace their counts; the originals go to `layers["counts_raw"]`. Cells CellBender
   gives a cell probability below 0.5 (its own cutoff) are removed as `ambient_empty`; a cell it did
   not analyse at all (below its low-count cutoff) gets probability 0.
   - `auto` (default): CellBender when the step has a GPU (the `gpu` variant). Without one it is
     skipped with a warning: on CPU CellBender took about an hour for the smallest demo sample even
     at 10 epochs (see Runtime), so a CPU budget would almost never be met for real 10x data and
     would only waste the time it took to find out. `cellbender` runs it anyway (also on CPU) and
     fails the sample when CellBender fails; `none` never runs it.
   - **Why no fallback correction on CPU.** SoupX (Young & Behjati 2020) and DecontX are R packages;
     there is no maintained, licence-clean Python port, and re-implementing SoupX's contamination
     estimate here would be an unreviewed method. Samples without correction keep their counts and
     say so (`ambient_method = none` in `qc_summary.csv`, a warning in the log and provenance).
2. **QC metrics** on the counts that will be written (after ambient removal when it ran):
   `total_counts`, `n_genes_by_counts`, `pct_counts_mt`, `pct_counts_ribo`,
   `pct_counts_in_top_20_genes` (mito = `var["mt"]` from sc-prep, else `MT-`/`mt-`).
3. **Outlier flags** (sc-best-practices, "Quality Control" chapter of Heumos et al. 2023): a cell
   is an outlier when `log1p(total_counts)` or `log1p(n_genes_by_counts)` is more than `MADS` (5)
   median absolute deviations from the sample's median, on either side (`low_counts`,
   `high_counts`, `low_genes`, `high_genes`), or when its share of counts in the 20 most expressed
   genes is more than `MADS` MADs above the median (`high_top20`). **Deviation from the guide:**
   mitochondrial content is a fixed cap, `pct_counts_mt > MAX_PCT_MT` (5) → `high_mito`, instead of
   the guide's "3 MADs or > 8 %". In nuclei the mitochondrial reads are cytoplasmic carry-over and
   near zero for most nuclei, so a MAD cut on that narrow distribution would remove good nuclei;
   5 % is the usual nuclei cap (whole cells: 10-20 %, more for muscle and heart).
4. **Doublets** (`DOUBLETS`): Scrublet (Wolock et al. 2019, scanpy's port) on the counts of the
   cells that passed 1-3, `expected_doublet_rate = DOUBLET_RATE`, `random_state = SEED`, its
   automatic threshold (the minimum between the two modes of the simulated doublets' scores). With
   `remove` (default) predicted doublets are removed (`doublet`); with `flag` they are kept with
   `obs["predicted_doublet"] = true`. Samples with fewer than 100 cells passing 1-3 are not scored.
   Scrublet rather than scDblFinder (R) or SOLO (scvi-tools, a second model to train per sample):
   it is BSD-licensed inside scanpy, needs no GPU, and the guide lists it among the methods to use.
5. **Normalization.** `layers["counts"]` = the counts (int32, integers; ambient-corrected when
   CellBender ran), `X` = `log1p(normalize_total(target_sum = TARGET_SUM))` (10,000 for every
   sample, so samples share a scale), float32 CSR. **No highly variable genes and no gene
   filtering here**: `sc-integrate` selects HVGs across all samples, and dropping genes per sample
   would shrink the gene set the samples share.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `cells` | input (subfolder), required | `<sample_id>.h5ad` | sc-prep output: raw counts in `X`. Found in `INPUT_DIR/cells/`, else anywhere under `INPUT_DIR`. |
| `droplets` | input (subfolder), optional | `<sample_id>.raw.h5` | sc-prep's copy of the raw 10x matrix. Found in `INPUT_DIR/droplets/`, else anywhere. |
| `cells` | output | `cells/<sample_id>.h5ad`, `cells/provenance.json` | The kept cells; contract below. |
| `qc` | output | `qc/*.csv`, `qc/*.png`, `qc/summary.json`, `qc/provenance.json` | Tables and figures, below. |

Symlinks are followed; a file whose real path leaves `INPUT_DIR` (and the platform's own links at
its top level) is skipped. A `.h5ad` that already went through sc-qc is refused (its `X` is not counts).

### `cells/<sample_id>.h5ad`

Everything sc-prep wrote for the kept cells, every gene, plus:

| Slot | Content |
|---|---|
| `X` | log-normalized, float32 CSR |
| `layers["counts"]` | counts, int32 CSR (after ambient removal when it ran) |
| `layers["counts_raw"]` | sc-prep's counts, only when CellBender ran |
| `obs["total_counts"]`, `["n_genes_by_counts"]`, `["pct_counts_mt"]`, `["pct_counts_ribo"]` | recomputed on `layers["counts"]` |
| `obs["doublet_score"]`, `obs["predicted_doublet"]` | Scrublet (NaN / false when not scored) |
| `obs["ambient_fraction"]`, `obs["cellbender_cell_probability"]` | only when CellBender ran |
| `uns["pennsieve"]` | adds `qc` (version, parameters, ambient status and note, doublet threshold, cells before/after), `normalized: true`, `target_sum` |

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `AMBIENT` | `auto` \| `cellbender` \| `none` | `auto` | CellBender on a GPU (`auto`), always (`cellbender`) or never. |
| `AMBIENT_EPOCHS` | integer 1-1000 | `150` | CellBender training epochs (its recommended default). |
| `AMBIENT_FPR` | number | `0.01` | CellBender `--fpr`: higher removes more background and more signal. |
| `DOUBLETS` | `remove` \| `flag` | `remove` | Remove Scrublet's doublets or only flag them. |
| `DOUBLET_RATE` | number | `0.05` | Scrublet's expected doublet rate (~0.8 % per 1,000 cells loaded on 10x Chromium). |
| `MADS` | number | `5` | MADs for the count, gene and top-20 outlier flags (5 permissive, 3 strict). |
| `MAX_PCT_MT` | number 0-100 | `5` | Mitochondrial cap (nuclei 5; whole cells 10-20). |
| `TARGET_SUM` | number | `10000` | Counts per cell before `log1p`. |
| `SEED` | integer | `42` | Scrublet's random state. |
| `THREADS` | integer | `0` | BLAS / numba / CellBender CPU threads; 0 = the container's CPU quota. |

Bad values exit with code 2 before anything is read.

## Outputs

### `qc/qc_summary.csv`

One row per sample: `sample_id`, `n_cells_before`, `n_cells_after`, `pct_removed`,
`n_removed_<reason>` for `ambient_empty`, `low_counts`, `high_counts`, `low_genes`, `high_genes`,
`high_top20`, `high_mito`, `doublet` (a cell can fail several), `ambient_method` (`cellbender` or
`none`), `ambient_status` (`cellbender`, `skipped`, `off`, `failed`), `ambient_note` (why),
`median_ambient_fraction`, `n_genes`, `median_counts`, `median_genes`, `median_pct_mt` (of the kept
cells), `doublet_threshold`, `n_predicted_doublets`, `status` (`ok` or `failed: <error>`).

### `qc/<sample_id>_qc.csv`

Every input cell: `sample_id, obs_id, kept, reason` (`;`-joined reasons, empty when kept), the five
QC metrics, and when CellBender ran `total_counts_raw, ambient_fraction,
cellbender_cell_probability`, then `doublet_score, predicted_doublet`.

### `qc/thresholds.csv`

Per sample × metric: `scale` (log1p or linear), `median`, `mad`, `lower`, `upper` (linear scale),
`rule`.

### Figures

- `qc/<sample_id>_qc.png`: total counts, genes detected and % mitochondrial counts as violins of all
  input cells (grey) and the kept cells (blue) with the thresholds (red dashes), and counts vs genes
  per cell, kept and removed.
- `qc/<sample_id>_doublets.png`: Scrublet scores of the cells and of the simulated doublets, with
  the threshold. A good call shows a bimodal simulated distribution with the threshold between the
  modes; a threshold far in the tail means Scrublet found little to call.
- `qc/<sample_id>_ambient.png` (CellBender ran): per-cell fraction of counts removed as ambient.
- `qc/<sample_id>_cellbender_metrics.csv` (CellBender ran): CellBender's own run metrics.

### `qc/summary.json`, `provenance.json`

`summary.json`: totals, removals by reason, samples with ambient removal, the `qc_summary.csv` rows.
`provenance.json` (both outputs): `app`, `version`, `tools`, `parameters`, `inputs`, `skipped`,
`started`, `finished`, `samples` (status, ambient details, timings), `peak_rss_mb`,
`cgroup_memory_peak_mb`, `device`, `method`, and `citations` as `{doi, pmid, text}`.

## Results on the demo data

drg-c2 (Discover 475 v1, 12 C2 dorsal root ganglia from neck-fusion surgery, 10x Flex on nuclei)
through sc-prep 0.1.0 (seeded subsample to ≤ 3,000 nuclei per sample, 18,082-18,099 genes), then
sc-qc with the defaults on CPU (so no ambient removal: `ambient_method = none` for every sample).
Outputs in `workflow-single-cell/data/work/drg-c2/qc/{cells,qc}/` for sc-integrate.

| Sample | Nuclei before → after (% removed) | Median counts / genes / % mt (kept) | Removed: high counts / low genes / high genes / top-20 / mito / doublet | Scrublet threshold |
|---|---|---|---|---|
| CC0009-DRGL | 1,265 → 1,100 (13.0 %) | 832 / 698 / 1.38 | 33 / 0 / 24 / 52 / 98 / 13 | 0.28 |
| CC0021-DRGL | 1,593 → 1,509 (5.3 %) | 1,173 / 992 / 0.32 | 7 / 0 / 0 / 67 / 6 / 5 | 0.30 |
| CC0024-DRGR | 186 → 165 (11.3 %) | 760 / 646 / 0.65 | 5 / 0 / 2 / 9 / 3 / 5 | 0.13 |
| CC0026-DRGL | 235 → 209 (11.1 %) | 1,391 / 1,065 / 0.38 | 0 / 0 / 0 / 10 / 4 / 12 | 0.13 |
| CC0027-DRGL | 293 → 260 (11.3 %) | 1,005 / 810 / 1.56 | 7 / 0 / 3 / 8 / 10 / 9 | 0.15 |
| CC0029-DRGR | 3,000 → 2,724 (9.2 %) | 654 / 570 / 2.56 | 35 / 1 / 18 / 69 / 176 / 3 | 0.31 |
| CC0030-DRGR | 1,472 → 1,243 (15.6 %) | 1,098 / 855 / 2.20 | 26 / 0 / 18 / 8 / 189 / 4 | 0.29 |
| CC0032-DRGR | 3,000 → 2,710 (9.7 %) | 2,329 / 1,625 / 1.27 | 45 / 0 / 1 / 120 / 186 / 2 | 0.31 |
| CC0033-DRGL | 371 → 328 (11.6 %) | 804 / 687 / 1.17 | 14 / 0 / 7 / 17 / 11 / 3 | 0.16 |
| CC0037-DRGR | 3,000 → 2,862 (4.6 %) | 1,574 / 1,167 / 0.16 | 5 / 0 / 1 / 100 / 44 / 0 | 0.46 |
| CC0048-DRGR | 3,000 → 2,520 (16.0 %) | 1,236 / 979 / 1.85 | 32 / 0 / 6 / 49 / 417 / 3 | 0.36 |
| CC0050-DRGL | 3,000 → 2,809 (6.4 %) | 1,534 / 1,202 / 1.77 | 44 / 0 / 1 / 84 / 108 / 0 | 0.38 |

20,415 → 18,439 nuclei (9.7 % removed). (Sample ids are `sam-UTD-<id>`.)

- **Mitochondrial content** is the largest cause (1,252 nuclei): most samples have a median of 1-2.5 %
  and a tail above 5 % (CC0048: 14 % of nuclei above 5 %). Raise `MAX_PCT_MT` to keep more.
- **Counts and genes**: nothing falls below the lower MAD cuts (Cell Ranger's cell calling already
  removed near-empty barcodes; the lower count cut is ~90 counts); 253 nuclei are above the upper
  count cut (more than ~13× the median), likely multiplets or debris clumps.
- **Top-20 share** removes 593 nuclei (2.9 %) whose counts are dominated by a few genes.
- **Doublets**: Scrublet calls few (59, 0.3 %), and on the larger samples its automatic threshold sits
  in the tail of the simulated scores (0.36-0.46), i.e. the simulated doublets are not bimodal.
  On a full, not subsampled, sample (CC0050, 19,372 nuclei) it called none. Flex can pool up to 16
  samples per run, and doublets between samples cannot be seen within one; within-sample
  doublets are rare at a few hundred to 3,000 nuclei per sample.
- **Ambient RNA, CellBender on CPU** (CC0024, 186 nuclei, `AMBIENT=cellbender`, 10 epochs, a check of
  the path, not a result): CellBender analysed 9,529 droplets plus 32,783 empties, estimated 1,889
  cells where Cell Ranger called 186, and after 10 epochs (its convergence indicator needs > 20)
  removed 5 % of counts overall but nothing from the 186 Cell Ranger nuclei (cell probability 1.0,
  ambient fraction 0). A real run needs the default 150 epochs, i.e. the gpu variant.

## Runtime

`computeTypes: [standard]` for the cpu variant (4 vCPU / 16 GiB, 8 h timeout); `[gpu]` with 14 GiB
for the gpu variant. Samples are processed one at a time; memory follows the largest sample.

Measured in the linux/amd64 image under emulation on an M-series Mac (Docker Desktop, `--cpus 2`,
`--network none`, read-only root, other builders running in parallel); peak = the cgroup's
`memory.peak`:

| Input | Cells × genes | Wall time | Peak memory |
|---|---|---|---|
| drg-c2, 12 samples (≤ 3,000 nuclei each), no CellBender | 20,415 nuclei in all | 139 s | 2.76 GB |
| CC0050 full (not subsampled), no CellBender | 19,372 × 18,099 (32 M non-zeros) | 164 s | 5.60 GB |
| fixture, CellBender 5 epochs | 400 + 70,000 empty droplets × 300 | 51 s | 1.70 GB |

Memory grows with the cells of the largest sample (Scrublet simulates twice as many doublets and
works on dense PCA input): about 5.6 GB at 19 k nuclei, so 16 GiB covers samples up to roughly
50 k cells.

**CellBender on CPU** (native arm64, 4 threads, sample CC0024: 9.5 k candidate droplets + 33 k
empties × 15.6 k features): 12.5-14.8 s per training epoch, then the posterior took 40-50 min;
55-84 min wall in all for 10 epochs (17-18 min of CPU time; other builds competed for the CPU), peak
3.4-4.1 GB. At the default 150 epochs that is ~35 min of training plus the posterior, i.e. 1-2 h per
sample on 4 vCPU, longer under amd64 emulation (not measured there; hours). **On a GPU** nothing was
measured: the gpu image was built and imports CUDA PyTorch, but there is no NVIDIA GPU here. Expect
tens of minutes at most per 10x sample on a T4 at 150 epochs (an estimate, to be confirmed on the
platform's first gpu run; the timeout is 8 h). CellBender writes a ~1 GB checkpoint
into the run's temporary folder under `/tmp`, deleted after each sample.

Images (`docker image inspect`): cpu 492 MB, gpu 3.28 GB (CUDA 13.2 PyTorch 2.15.0+cu132; builds and imports, `torch.cuda.is_available()` false here).

## Variants

| Variant | Image | Compute | Ambient RNA by default |
|---|---|---|---|
| `cpu` (default) | CPU-only PyTorch | `standard` | no (`AMBIENT=auto` skips it; set `cellbender` to force) |
| `gpu` | CUDA 13.2 PyTorch (NVIDIA driver ≥ 580, T4 or newer) | `gpu`, 14 GiB | yes |

Both are built from the same `Dockerfile` (`TORCH_FLAVOUR=cpu|cuda`); parameters and ports are
identical. The code picks the device with `torch.cuda.is_available()` and logs it.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt \
  && .venv/bin/pip install --no-deps cellbender==0.4.0
make test PYTHON=.venv/bin/python       # pytest on a synthetic fixture (no Docker), incl. a 5-epoch CellBender run
make run PYTHON=.venv/bin/python        # build the image, run it over data/fixture/linked -> data/output
make validate PYTHON=.venv/bin/python   # app.yml against the manifest schema, CITATION.cff
```

The fixture (`tests/make_fixture.py`) has a 400-nucleus sample with a raw droplet matrix (planted
near-empty, high-mitochondrial and doublet nuclei, ~10 % ambient counts) and a 300-nucleus sample
without one; the tests check that the planted cells are removed for their reason, the AnnData
contract of the output, CellBender's correction, both input layouts, parameter errors and path
checks. Real input: `sample-input/README.md`.

## Limitations

- **No ambient correction without a GPU** (by default) and none at all without the raw matrix. The
  uncorrected counts of neuron-rich tissue carry neuronal transcripts in glial nuclei; check marker
  genes in unexpected cell types before reading too much into them.
- **CellBender** needs a real raw matrix with tens of thousands of empty droplets; it may disagree
  with Cell Ranger's cell calls (on the Flex demo it estimated ten times more cells). Cells it calls
  empty are removed (`ambient_empty`); look at that count in `qc_summary.csv`.
- **Thresholds are per sample.** MADs follow each sample's distribution, so a sample that is poor
  throughout keeps most of its cells; compare the medians in `qc_summary.csv` across samples.
- **The mitochondrial cap is fixed.** 5 % suits nuclei; whole-cell data of muscle, heart or kidney
  need a much higher value.
- **Scrublet** misses doublets of similar cells (homotypic) and cannot see doublets between
  multiplexed samples (Flex, hashing); its automatic threshold fails on unimodal score distributions
  (logged; nothing is called then).
- **High top-20 share** can also be biology (cells dominated by a few very highly expressed genes);
  with 5 MADs it removed 3 % of the demo nuclei.

## Licences

Apache-2.0 for this app (see `LICENSE`). Bundled components and their licences: `THIRD_PARTY.md`
(CellBender BSD-3-Clause, PyTorch, Pyro Apache-2.0, scanpy BSD-3-Clause, …; NVIDIA CUDA libraries
under the NVIDIA EULA in the gpu image only; no GPL code).

## Credits

- **CellBender** by Stephen J. Fleming, Mehrtash Babadi and colleagues at the Broad Institute removes
  the ambient RNA and gives each droplet's cell probability.
- **Scrublet** by Samuel L. Wolock, Romain Lopez and Allon M. Klein is the doublet detector, run
  through scanpy's port.
- **Scanpy** by F. Alexander Wolf, Philipp Angerer, Fabian J. Theis and the scverse community
  computes the QC metrics, runs Scrublet and normalizes.
- **anndata** by Isaac Virshup and colleagues (scverse) reads and writes every `.h5ad`.
- **The single-cell best-practices guide** by Lukas Heumos, Fabian J. Theis and the Single-cell Best
  Practices Consortium sets the MAD-based outlier rules this app follows.
- **PyTorch** (the PyTorch Foundation) and **Pyro** (Uber AI Labs and the Pyro team) run CellBender's
  model; **scikit-image** provides Scrublet's automatic threshold.
- **Demo data**: Asta Arendt-Tranholm, Ishwarya Sankaranarayanan and Theodore J. Price (University of
  Texas at Dallas) shared the C2 dorsal root ganglia nuclei on Pennsieve Discover (dataset 475).

## Citations

> Fleming SJ, Chaffin MD, Arduini A, Akkad AD, Banks E, Marioni JC, Philippakis AA, Ellinor PT,
> Babadi M. **Unsupervised removal of systematic background noise from droplet-based single-cell
> experiments using CellBender.** Nat Methods 20(9):1323-1335 (2023).
> https://doi.org/10.1038/s41592-023-01943-7, PMID 37550580.

> Wolock SL, Lopez R, Klein AM. **Scrublet: Computational Identification of Cell Doublets in
> Single-Cell Transcriptomic Data.** Cell Syst 8(4):281-291.e9 (2019).
> https://doi.org/10.1016/j.cels.2018.11.005, PMID 30954476.

> Heumos L, Schaar AC, Lance C, Litinetskaya A, Drost F, Zappia L, Lücken MD, Strobl DC, Henao J,
> Curion F, Single-cell Best Practices Consortium, Aliee H, Ansari M, Badia-i-Mompel P, Büttner M,
> Dann E, Dimitrov D, Dony L, Frishberg A, He D, Hediyeh-zadeh S, Hetzel L, Ibarra IL, Jones MG,
> Lotfollahi M, Martens LD, Müller CL, Nitzan M, Ostner J, Palla G, Patro R, Piran Z,
> Ramírez-Suástegui C, Saez-Rodriguez J, Sarkar H, Schubert B, Sikkema L, Srivastava A, Tanevski J,
> Virshup I, Weiler P, Schiller HB, Theis FJ. **Best practices for single-cell analysis across
> modalities.** Nat Rev Genet 24(8):550-572 (2023). https://doi.org/10.1038/s41576-023-00586-w,
> PMID 37002403.

> Wolf FA, Angerer P, Theis FJ. **SCANPY: large-scale single-cell gene expression data analysis.**
> Genome Biol 19(1):15 (2018). https://doi.org/10.1186/s13059-017-1382-0, PMID 29409532.

> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. **anndata: Access and store annotated data
> matrices.** J Open Source Softw 9(101):4371 (2024). https://doi.org/10.21105/joss.04371 (not in
> PubMed).

Demo data:

> Arendt-Tranholm A, Sankaranarayanan I, Price TJ. **Single nucleus sequencing of human C2 dorsal
> root ganglia tissue from fusion surgeries** (Version 1). SPARC Portal / Pennsieve Discover dataset
> 475 (2025). https://doi.org/10.26275/cxuf-ftav (CC BY 4.0).

> Arendt-Tranholm A, Sankaranarayanan I, Payne C, Mancilla Moreno M, Mazhar K, Yap N, Chiu AP,
> Barry A, Patel PJ, Inturi NN, Tavares-Ferreira D, Amin A, Karandikar M, Jarvik JG, Turner JA,
> Hofstetter CP, Curatolo M, Price TJ. **Transcriptome of the human C2 dorsal root ganglia in C1–2
> arthrodesis surgery: insight for neck pain.** Brain 149(4):1349-1364 (2026; online 2025).
> https://doi.org/10.1093/brain/awaf370, PMID 41032656.

Every paper was checked on Crossref (title, full author list and order, journal, volume, issue,
pages) and PubMed (E-utilities `esearch` by DOI, then `esummary`: PMID, title, pages, DOI); the
dataset on DataCite (creators with ORCIDs, title, year, CC BY 4.0) and the Discover API (dataset 475,
version 1, contributors, licence). Software licences were read from each package at its pinned
version (CellBender's `LICENSE` in the 0.4.0 wheel). All are in `CITATION.cff`.
