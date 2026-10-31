# Third-party software

This app is Apache-2.0 (see `LICENSE`). The images it builds contain the following third-party
components, installed from PyPI (PyTorch for the cpu variant from the PyTorch CPU index) at the pinned
versions in `requirements.txt`, `requirements-cpu.txt` and `requirements-cuda.txt`; their licences apply
to the images. Licences were read from each package's metadata at the pinned version (CellBender's also
from the `LICENSE` file in its 0.4.0 wheel). No GPL Python packages are installed (sc-qc does no
clustering, so neither igraph nor leidenalg is); the Debian base layer contains GPL system tools.

| Component | What | Licence | Source |
|---|---|---|---|
| CellBender 0.4.0 | `remove-background`: ambient RNA removal and cell probabilities; its `downstream.anndata_from_h5` reads the result. Installed without its Jupyter report extras. | BSD-3-Clause (Broad Institute) | https://github.com/broadinstitute/CellBender |
| PyTorch 2.14.1 (`+cpu` in the cpu image, `+cu126` in the gpu image) | CellBender's model. | PyTorch licence (package metadata: Apache-2.0 AND Apache-2.0 WITH LLVM-exception AND …, BSD-style parts) | https://pytorch.org |
| NVIDIA CUDA 12.6 runtime wheels (`nvidia-*-cu12`, `cuda-toolkit`, `cuda-bindings`, `cuda-pathfinder`), Triton 3.8.0 | gpu image only: the CUDA libraries PyTorch loads. | NVIDIA CUDA Toolkit EULA (redistributable runtime components); Triton MIT | https://pypi.org/project/nvidia-cublas/ |
| Pyro 1.9.2, pyro-api 0.1.2, opt_einsum 3.4.0 | Probabilistic programming used by CellBender. | Apache-2.0 (opt_einsum MIT) | https://pyro.ai |
| scanpy 1.12.5 | `calculate_qc_metrics`, `scrublet` (its port of Scrublet), `normalize_total`, `log1p`. | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata 0.13.4 | Reads and writes `.h5ad`. | BSD-3-Clause | https://github.com/scverse/anndata |
| scikit-image 0.26.0 (+ imageio 2.38.1, tifffile, lazy-loader) | `threshold_minimum`, Scrublet's automatic doublet threshold. | BSD-3-Clause (imageio BSD-2-Clause) | https://scikit-image.org |
| PyTables 3.11.1, blosc2 4.14.1, ndindex, numexpr, msgpack, py-cpuinfo | HDF5 reader CellBender uses for 10x files. | BSD-3-Clause / MIT / Apache-2.0 (msgpack) | https://www.pytables.org |
| dill 0.4.1, psutil 7.2.2 | CellBender's checkpoint serialization (its own, in the run's temporary folder) and memory reporting. | BSD-3-Clause | PyPI |
| IPython 9.18.0 and its dependencies (traitlets, jedi, parso, prompt_toolkit, pexpect, ptyprocess, wcwidth, stack-data, executing, asttokens, pure-eval, matplotlib-inline, ipython-pygments-lexers, Pygments) | Imported by CellBender's report module (the report itself is not written). | BSD-3-Clause / MIT / ISC (pexpect, ptyprocess) / Apache-2.0 (asttokens) / BSD-2-Clause (Pygments) | PyPI |
| NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.0, h5py 3.16.0 | Arrays, CSR matrices, tables, HDF5. | BSD-3-Clause | numpy.org, scipy.org, pandas.pydata.org, h5py.org |
| scikit-learn 1.9.1, numba 0.68.1, llvmlite 0.50.0 | Scrublet's PCA and nearest neighbours (through scanpy); compiled kernels. | BSD-3-Clause; BSD-2-Clause (llvmlite also Apache-2.0 with LLVM exception) | PyPI |
| Matplotlib 3.11.2, Pillow 13.0.0 | PNG figures. | Matplotlib licence (PSF-based); MIT-CMU | matplotlib.org |
| umap-learn, pynndescent, statsmodels, patsy, seaborn, networkx, joblib, threadpoolctl, natsort, session-info2, legacy-api-wrap, fast-array-utils, scverse-misc, zarr, numcodecs, donfig, google-crc32c, formulaic, narwhals, pydantic, rich, httpx and the remaining pins in `requirements.txt` | Dependencies of the above. | BSD / MIT / Apache-2.0 / MPL-2.0 (scverse helpers; file-level copyleft that leaves our code untouched) / PSF | PyPI |
| sympy 1.14.0, mpmath, Jinja2, MarkupSafe, filelock, fsspec, setuptools, typing_extensions | Dependencies of PyTorch. | BSD-3-Clause / MIT | PyPI |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Not bundled

- **Demo data** (not in the image or the repository): Arendt-Tranholm A, Sankaranarayanan I, Price TJ,
  *Single nucleus sequencing of human C2 dorsal root ganglia tissue from fusion surgeries*, Pennsieve
  Discover dataset 475 v1, doi:10.26275/cxuf-ftav (CC BY 4.0).
