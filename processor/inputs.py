"""
Finding the inputs: sc-prep's <sample_id>.h5ad files (port `cells`) and, optionally, the raw
droplet matrices <sample_id>.raw.h5 (port `droplets`, 10x HDF5).

Linked layout (placement: subfolder): INPUT_DIR/cells/ and INPUT_DIR/droplets/. Without links the
same files are found by name anywhere under INPUT_DIR. Which files may be read is decided by the
shared safe_inputs.InputDir (the platform's symlink rule); every HDF5 file is checked with
check_h5 (no external links, external storage or virtual datasets) before it is used.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from processor.safe_inputs import InputDir, UnsafeInput

SAMPLE_ID = re.compile(r"^[A-Za-z0-9_-]+$")
RAW_SUFFIX = ".raw.h5"


@dataclass
class Sample:
    sample_id: str
    path: str                 # resolved, checked path of the .h5ad
    rel: str                  # relative to INPUT_DIR, for logs and provenance
    raw: str | None = None    # resolved, checked path of the .raw.h5
    raw_rel: str | None = None


@dataclass
class Inputs:
    samples: list[Sample]
    seen: list[str] = field(default_factory=list)              # relative paths of every file considered
    skipped: dict[str, str] = field(default_factory=dict)       # relative path -> reason


def _port_files(files: list[Path], idir: InputDir, port: str) -> list[Path]:
    """A port's files: those under INPUT_DIR/<port>/ when that folder exists, else every file."""
    if not (idir.root / port).is_dir():
        return files
    return [p for p in files if idir.rel(p).split("/", 1)[0] == port]


def _collect(files: list[Path], idir: InputDir, suffix: str, inputs: Inputs) -> dict[str, tuple[str, str]]:
    by_id: dict[str, tuple[str, str]] = {}
    for p in files:
        name = p.name.lower()
        if not name.endswith(suffix) or (suffix == ".h5ad" and name.endswith(RAW_SUFFIX)):
            continue
        rel = idir.rel(p)
        inputs.seen.append(rel)
        sid = p.name[: -len(suffix)]
        if not SAMPLE_ID.match(sid):
            inputs.skipped[rel] = f"file name is not <sample_id>{suffix} with sample_id in [A-Za-z0-9_-]"
            continue
        if sid in by_id:
            inputs.skipped[rel] = f"a second file for sample {sid}"
            continue
        try:
            idir.check_h5(p)
            by_id[sid] = (idir.check(p), rel)
        except UnsafeInput as exc:
            inputs.skipped[rel] = str(exc)
    return by_id


def discover(input_dir: Path) -> Inputs:
    idir = InputDir(input_dir)
    files = idir.walk()
    inputs = Inputs(samples=[])
    cells = _collect(_port_files(files, idir, "cells"), idir, ".h5ad", inputs)
    raws = _collect(_port_files(files, idir, "droplets"), idir, RAW_SUFFIX, inputs)
    for sid, (_, rel) in sorted(raws.items()):
        if sid not in cells:
            inputs.skipped[rel] = "raw matrix without a <sample_id>.h5ad"
    for sid, (path, rel) in sorted(cells.items()):
        raw, raw_rel = raws.get(sid, (None, None))
        inputs.samples.append(Sample(sid, path, rel, raw, raw_rel))
    return inputs
