"""
Finding the inputs: sc-prep's <sample_id>.h5ad files (port `cells`) and, optionally, the raw
droplet matrices <sample_id>.raw.h5 (port `droplets`, 10x HDF5).

Linked layout (placement: subfolder): INPUT_DIR/cells/ and INPUT_DIR/droplets/. Without links
the same files are found by name anywhere under INPUT_DIR. Symlinks are followed: the platform's
own links (the entries of INPUT_DIR) may point anywhere, but a file or folder further down whose
real path leaves INPUT_DIR and those link targets is skipped.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

SAMPLE_ID = re.compile(r"^[A-Za-z0-9_-]+$")
RAW_SUFFIX = ".raw.h5"


@dataclass
class Sample:
    sample_id: str
    path: Path
    raw: Path | None = None


@dataclass
class Inputs:
    samples: list[Sample]
    seen: list[str] = field(default_factory=list)              # relative paths of every file considered
    skipped: dict[str, str] = field(default_factory=dict)       # relative path -> reason


def allowed_roots(input_dir: Path) -> list[Path]:
    """INPUT_DIR and the real paths of its top-level entries (the platform's links)."""
    return [input_dir.resolve()] + [entry.resolve() for entry in input_dir.iterdir()]


def inside(path: Path, roots: list[Path]) -> bool:
    real = path.resolve()
    return any(real == r or real.is_relative_to(r) for r in roots)


def walk(root: Path, roots: list[Path]) -> list[Path]:
    """Every file under root, following symlinks that stay inside `roots`, each real folder once."""
    found: list[Path] = []
    visited: set[Path] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = Path(dirpath).resolve()
        if real in visited or not inside(Path(dirpath), roots):
            dirnames[:] = []
            continue
        visited.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return found


def _port_files(input_dir: Path, port: str, roots: list[Path]) -> list[Path]:
    folder = input_dir / port
    return walk(folder if folder.is_dir() else input_dir, roots)


def _collect(files: list[Path], input_dir: Path, roots: list[Path], suffix: str, inputs: Inputs) -> dict[str, Path]:
    by_id: dict[str, Path] = {}
    for p in files:
        if not p.name.lower().endswith(suffix) or (suffix == ".h5ad" and p.name.lower().endswith(RAW_SUFFIX)):
            continue
        rel = str(p.relative_to(input_dir))
        inputs.seen.append(rel)
        sid = p.name[: -len(suffix)]
        if not inside(p, roots):
            inputs.skipped[rel] = "links outside INPUT_DIR"
        elif not SAMPLE_ID.match(sid):
            inputs.skipped[rel] = f"file name is not <sample_id>{suffix} with sample_id in [A-Za-z0-9_-]"
        elif sid in by_id:
            inputs.skipped[rel] = f"a second file for sample {sid}"
        else:
            by_id[sid] = p
    return by_id


def discover(input_dir: Path) -> Inputs:
    roots = allowed_roots(input_dir)
    inputs = Inputs(samples=[])
    cells = _collect(_port_files(input_dir, "cells", roots), input_dir, roots, ".h5ad", inputs)
    raws = _collect(_port_files(input_dir, "droplets", roots), input_dir, roots, RAW_SUFFIX, inputs)
    for sid, path in sorted(raws.items()):
        if sid not in cells:
            inputs.skipped[str(path.relative_to(input_dir))] = "raw matrix without a <sample_id>.h5ad"
    inputs.samples = [Sample(sid, path, raws.get(sid)) for sid, path in sorted(cells.items())]
    return inputs
