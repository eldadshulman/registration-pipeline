"""End-to-end tests for hest_valis.report: a synthetic <out>/<sample> must render a non-empty
PDF, and the no_nuclei case must render without crashing. Also pins the concordance refactor
(mutual_nn_pairs / density_grids) that the report reuses.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from hest_valis import concordance, report

UM = 0.2125


def _write_qc(sample_dir, metrics, decision):
    os.makedirs(sample_dir, exist_ok=True)
    with open(os.path.join(sample_dir, "qc.json"), "w") as f:
        json.dump({"sample_id": os.path.basename(sample_dir),
                   "metrics": metrics, "decision": decision}, f)


def _aligned_points(seed=0, n=1000, span=6000.0, noise_px=2.0):
    rng = np.random.default_rng(seed)
    xen = rng.uniform(0, span, (n, 2))
    he = xen + rng.normal(0, noise_px, (n, 2))   # well-aligned warped H&E nuclei
    return he, xen


# --------------------------------------------------------------------------------------
# concordance refactor: report reuses these, so they must stay consistent with the metrics
# --------------------------------------------------------------------------------------
def test_mutual_nn_pairs_matches_distance_view():
    he, xen = _aligned_points()
    cutoff_px = concordance.CUTOFF_UM / UM
    ia, ib, dpx = concordance.mutual_nn_pairs(xen, he, cutoff_px)
    d_only = concordance._mutual_nn(xen, he, cutoff_px)
    assert len(ia) == len(ib) == len(dpx) == len(d_only)
    assert np.allclose(np.sort(dpx), np.sort(d_only))
    assert (dpx <= cutoff_px).all()


def test_density_grids_shape_and_footprint():
    he, xen = _aligned_points()
    bin_px = concordance.BIN_UM / UM
    gh, gx, foot, nbx, nby = concordance.density_grids(he, xen, bin_px)
    assert gh.shape == gx.shape == foot.shape == (nby, nbx)
    assert foot.sum() >= 10            # enough footprint bins for a valid density_r
    # histogram2d silently clips points outside [0, n*bin_px] (e.g. he = xen + noise can dip
    # below 0 near the origin), so the binned count is <= n and should hold nearly all of them
    assert 0 < gh.sum() <= len(he)
    assert gh.sum() >= 0.98 * len(he)


# --------------------------------------------------------------------------------------
# report end-to-end
# --------------------------------------------------------------------------------------
def test_render_sample_end_to_end(tmp_path):
    sample_dir = str(tmp_path / "S1")
    os.makedirs(sample_dir)
    he, xen = _aligned_points(seed=1)
    np.save(os.path.join(sample_dir, "he_nuclei_nomicro.npy"), he)

    m = concordance.compute_qc(he, xen, UM)          # real schema, real numbers
    assert isinstance(m.get("density_r"), float)     # synthetic data is non-degenerate
    decision = {"chosen": "nomicro", "rule": "lower_um",
                "sel_median_um": m["nucleus_coincidence"]["median_um"],
                "sel_density_r": m["density_r"]}
    _write_qc(sample_dir, {"nomicro": m}, decision)

    pdf = report.render_sample(sample_dir, xen, UM, sample_id="S1", save_png=True)
    assert os.path.exists(pdf) and os.path.getsize(pdf) > 1000
    assert os.path.exists(os.path.join(sample_dir, "report.png"))


def test_render_sample_no_nuclei(tmp_path):
    sample_dir = str(tmp_path / "S2")
    os.makedirs(sample_dir)
    # no he_nuclei_*.npy on disk; metric flags no_nuclei
    metric = {"status": "no_nuclei", "pixel_um": UM, "n_he": 0, "n_xenium": 5}
    decision = {"chosen": "nomicro", "rule": "nomicro_unavailable",
                "sel_median_um": None, "sel_density_r": None}
    _write_qc(sample_dir, {"nomicro": metric}, decision)

    xen = np.random.default_rng(0).uniform(0, 100, (5, 2))
    pdf = report.render_sample(sample_dir, xen, UM, sample_id="S2")
    assert os.path.exists(pdf) and os.path.getsize(pdf) > 500    # banner-only page still valid


def test_render_sample_coarse_no_raster(tmp_path):
    """A coarse-rescued slide (no warped raster) renders point panels only, no crash."""
    sample_dir = str(tmp_path / "S3")
    os.makedirs(sample_dir)
    he, xen = _aligned_points(seed=3)
    np.save(os.path.join(sample_dir, "he_nuclei_coarse.npy"), he)
    m = concordance.compute_qc(he, xen, UM)
    decision = {"chosen": "coarse", "rule": "coarse_rescue_negative_density_r",
                "sel_median_um": m["nucleus_coincidence"]["median_um"],
                "sel_density_r": m["density_r"]}
    _write_qc(sample_dir, {"nomicro": m, "coarse": m}, decision)
    pdf = report.render_sample(sample_dir, xen, UM, sample_id="S3")
    assert os.path.exists(pdf) and os.path.getsize(pdf) > 1000


def test_render_cohort(tmp_path):
    out = str(tmp_path / "output")
    for i in range(3):
        sd = os.path.join(out, f"P{i}")
        os.makedirs(sd)
        he, xen = _aligned_points(seed=i)
        m = concordance.compute_qc(he, xen, UM)
        _write_qc(sd, {"nomicro": m},
                  {"chosen": "nomicro", "rule": "lower_um",
                   "sel_median_um": m["nucleus_coincidence"]["median_um"],
                   "sel_density_r": m["density_r"]})
    pdf = report.render_cohort(out)
    assert os.path.exists(pdf) and os.path.getsize(pdf) > 1000
