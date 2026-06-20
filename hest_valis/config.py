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

    # --- QC / rescue thresholds: the SINGLE source of truth ---------------------------------
    # The report's cohort triage lines (report.py) AND the provenance acceptance gate
    # (provenance.py) both read THESE, so "what looks good" and "what's safe to move" can
    # never drift apart. Override per-cohort in config.json; do not hard-code elsewhere.
    "density_r_accept": 0.70,    # accept gate / cohort "good" line: density_r must be >= this
    "median_um_accept": 5.0,     # accept gate / cohort "good" line: median offset (um) <= this
    "rescue_trigger_r": 0.10,    # selected density_r below this -> attempt the coarse/orient rescue
    #                              (run_qc COARSE_TRIGGER default; the TNBC cohort used 0.60 because
    #                               its 90-deg-rotated slides scored ~0.2-0.3, above 0.10)
    "rescue_delta_min": 0.20,    # rescued slide with (post_r - pre_r) < this -> flag for eyes
}

REQUIRED_COLS = ["sample_id", "he_path", "dapi_path", "xenium_cells"]
# samples.csv may also carry an optional `batch` column (scanner-run key); provenance.py uses it
# for per-batch orientation-consistency checks. Absent -> all slides share batch "UNKNOWN".


def thresholds(cfg):
    """The single QC/rescue threshold set (config overrides, else DEFAULTS).

    Returns a dict: density_r_accept, median_um_accept, rescue_trigger_r, rescue_delta_min.
    report.py and provenance.py both pull from here so their lines cannot diverge.
    """
    return {k: cfg.get(k, DEFAULTS[k]) for k in
            ("density_r_accept", "median_um_accept", "rescue_trigger_r", "rescue_delta_min")}


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
