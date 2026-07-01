"""Tests for run_cohort_qc.py -- the cohort QC overlay PDF (selected protocol per slide).

Covers the pure logic (disposition from the real gate), the standard-output slide assembly
(qc.json + he_nuclei_<chosen>.npy + xenium parquet), and that the renderer writes a non-empty PDF
with 1 summary page + one page per slide.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_cohort_qc as R
from hest_valis import config as hvconfig

TH = hvconfig.thresholds(hvconfig.DEFAULTS)


def test_disposition_uses_real_gate():
    assert R.disposition("nomicro", 0.92, 3.0, TH)[0] == "accepted"
    assert R.disposition("micro", 0.50, 3.0, TH)[0] == "manual-review"   # r < 0.70
    assert R.disposition("micro", 0.92, 9.0, TH)[0] == "manual-review"   # um > 5.0
    assert R.disposition("coarse", 0.92, 3.0, TH)[0].startswith("rescued")
    assert R.disposition(None, None, None, TH)[0] == "QC-FAILED"


def _seed_sample(root, sid, chosen, r, u, n=500):
    sdir = os.path.join(root, sid); os.makedirs(sdir, exist_ok=True)
    json.dump({"decision": {"chosen": chosen, "rule": "lower_um", "sel_density_r": r, "sel_median_um": u},
               "metrics": {}}, open(os.path.join(sdir, "qc.json"), "w"))
    np.save(os.path.join(sdir, f"he_nuclei_{chosen}.npy"),
            np.random.default_rng(0).random((n, 2)) * 100.0)
    cells = os.path.join(sdir, "cells.parquet")
    rng = np.random.default_rng(1)
    pd.DataFrame({"x_centroid": rng.random(n) * 21.25, "y_centroid": rng.random(n) * 21.25}
                 ).to_parquet(cells)
    return sdir, cells


def test_load_slide_from_standard_outputs(tmp_path):
    root = tmp_path / "out"
    sdir, cells = _seed_sample(str(root), "S1", "nomicro", 0.88, 3.4)
    s = R._load_slide("S1", sdir, cells, pixel_um=0.2125)
    assert s["protocol"] == "nomicro"
    assert s["he"].shape[1] == 2 and len(s["he"]) == 500
    assert s["xen"].shape[1] == 2 and len(s["xen"]) == 500     # loaded + converted um->px
    assert s["density_r"] == 0.88 and s["median_um"] == 3.4


def test_missing_qc_or_nuclei_is_handled(tmp_path):
    root = tmp_path / "out"; os.makedirs(root)
    assert R._load_slide("NONE", str(root / "NONE"), "", 0.2125) is None    # no qc.json -> skip
    sdir = root / "S2"; sdir.mkdir()
    json.dump({"decision": {"chosen": "micro", "sel_density_r": 0.8, "sel_median_um": 3.0}},
              open(sdir / "qc.json", "w"))
    s = R._load_slide("S2", str(sdir), "", 0.2125)              # qc.json but no he_nuclei_micro.npy
    assert s["he"] is None and "not found" in s["note"]


def test_render_writes_pdf_with_page_per_slide(tmp_path):
    slides = [
        {"sample_id": "S1", "protocol": "nomicro", "density_r": 0.90, "median_um": 3.1,
         "xen": np.random.default_rng(0).random((300, 2)) * 100,
         "he":  np.random.default_rng(1).random((280, 2)) * 100},
        {"sample_id": "S2", "protocol": "micro", "density_r": 0.55, "median_um": 4.0,
         "xen": np.random.default_rng(2).random((300, 2)) * 100,
         "he":  np.random.default_rng(3).random((280, 2)) * 100},
    ]
    out = tmp_path / "cohort_qc.pdf"
    n = R.render_cohort_pdf(slides, str(out), TH)
    assert n == 3                                   # 1 summary + 2 slides
    assert out.exists() and out.stat().st_size > 1000
