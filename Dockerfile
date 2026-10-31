# Base pinned by its multi-arch index digest (built for linux/amd64).
FROM python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=120

WORKDIR /processor

# PyTorch build: cpu or cuda. app.yml's variants set it (cpu variant: CPU-only torch from the PyTorch CPU
# index, no NVIDIA libraries; gpu variant: torch 2.14.1+cu126 for the platform's driver-550 GPU nodes, which
# CellBender uses with --cuda; the app exits 1 when this build finds no GPU).
# The default is cpu for a builder that passes no build arguments.
ARG TORCH_FLAVOUR=cpu

# The complete pinned dependency set, binary wheels only, no resolver (see requirements.txt).
COPY requirements.txt requirements-${TORCH_FLAVOUR}.txt /processor/
RUN pip install --no-deps --only-binary=:all: \
        -r /processor/requirements.txt -r /processor/requirements-${TORCH_FLAVOUR}.txt \
    && find /usr/local/lib/python3.12/site-packages -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12/site-packages -path '*/torch/*' -prune -o -path '*/tables/*' -prune -o -type d -name tests -prune -exec rm -rf {} + \
    && rm -rf /usr/local/lib/python3.12/site-packages/torch/include

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    HOME=/tmp \
    TMPDIR=/tmp \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    XDG_CACHE_HOME=/tmp/cache \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility \
    TORCH_FLAVOUR=${TORCH_FLAVOUR}

# Check that the stack imports; nothing is cached in the image (matplotlib's font cache and numba's cache
# are built under /tmp at run time, the only writable place besides OUTPUT_DIR).
RUN MPLCONFIGDIR=/tmp/build-mpl NUMBA_CACHE_DIR=/tmp/build-numba \
        python -c "import matplotlib.pyplot, scanpy, skimage.filters, cellbender.remove_background.run" \
    && rm -rf /tmp/build-mpl /tmp/build-numba

# Runs as uid 1000 on the platform; nothing here needs root. PyTorch looks the user up (getpass), so
# uid 1000 gets a passwd entry.
RUN useradd --uid 1000 --no-create-home --home-dir /tmp --shell /usr/sbin/nologin processor
USER 1000

ENTRYPOINT ["python3", "-m", "processor.entry"]
