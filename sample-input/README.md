# sample-input

This app reads the output of `sc-prep`: `cells/<sample_id>.h5ad` (raw counts in `X`) and, optionally,
`droplets/<sample_id>.raw.h5` (the raw 10x droplet matrix, for ambient RNA removal). The contract is
in `workflow-single-cell/docs/conventions.md`.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/{cells,qc}/
```

Two single-nucleus samples of 300 genes and three cell types: `nuc_a`, 400 nuclei with a raw droplet
matrix (70,000 empty droplets holding the ambient profile; ~10 % ambient counts in every nucleus), 12
near-empty nuclei, 10 with ~30 % mitochondrial counts and 20 doublets; `nuc_b`, 300 nuclei without a
raw matrix, 8 with ~30 % mitochondrial counts. `make run` uses `AMBIENT=cellbender` with 5 epochs
(`dev.env`), a few seconds on the fixture.

## Real data

The drg-c2 demo (Discover 475 v1, *Single nucleus sequencing of human C2 dorsal root ganglia tissue
from fusion surgeries*, doi:10.26275/cxuf-ftav, CC BY 4.0) through `sc-prep`, then:

```sh
docker run --rm --platform linux/amd64 --memory 4g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "/path/to/sc-prep/output:/data/input:ro" -v "$PWD/data/output:/data/output" sc-qc:cpu
```
