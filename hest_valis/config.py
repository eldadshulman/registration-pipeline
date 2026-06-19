"""Tiny config + samples loader (stdlib only, no yaml dependency).

config.json   : shared parameters (paths, pixel size). See examples/config.json.
samples.csv   : one row per slide. Required columns:
                  sample_id, he_path, dapi_path, xenium_cells
                he_path     : H&E whole-slide image (moving)
                dapi_path   : Xenium DAPI morphology_focus ch0 (fixed reference)
                xenium_cells: Xenium cells.parquet (centroids in microns) for QC
"""
import csv
import json
import os

DEFAULTS = {
    "pixel_um": 0.2125,          # Xenium DAPI um/pixel
    "output_dir": "./output",    # per-sample work + results land here
    "valis_env": "",             # path to the valis_hest_env (for register/warp jobs)
    "centroids_in_um": True,     # Xenium cells.parquet centroids are microns
}

REQUIRED_COLS = ["sample_id", "he_path", "dapi_path", "xenium_cells"]


def load_config(path):
    cfg = dict(DEFAULTS)
    if path and os.path.exists(path):
        with open(path) as _f:
            cfg.update(json.load(_f))
    return cfg


def load_samples(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            missing = [c for c in REQUIRED_COLS if not r.get(c)]
            if missing:
                raise ValueError(f"samples row {r.get('sample_id','?')} missing: {missing}")
            rows.append(r)
    return rows


def get_sample(samples, sample_id):
    for r in samples:
        if r["sample_id"] == sample_id:
            return r
    raise KeyError(f"sample_id {sample_id} not in samples file")
