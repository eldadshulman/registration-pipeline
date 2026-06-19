"""Unit tests for hest_valis.concordance._mutual_nn and compute_qc edge cases."""
import math
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest
from hest_valis.concordance import _mutual_nn, compute_qc


def test_mutual_nn_perfect_match():
    A = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0]])
    B = A.copy()
    d = _mutual_nn(A, B, cutoff_px=1.0)
    assert len(d) == 3
    assert np.allclose(d, 0.0)


def test_mutual_nn_no_match_outside_cutoff():
    A = np.array([[0.0, 0.0]])
    B = np.array([[100.0, 0.0]])
    d = _mutual_nn(A, B, cutoff_px=1.0)
    assert len(d) == 0


def test_mutual_nn_asymmetric():
    # A has an extra point far from B -> only the close pair matches mutually
    A = np.array([[0.0, 0.0], [50.0, 50.0]])
    B = np.array([[0.5, 0.0]])
    d = _mutual_nn(A, B, cutoff_px=2.0)
    assert len(d) == 1
    assert d[0] < 1.0


def test_compute_qc_empty_he():
    he = np.empty((0, 2), dtype=float)
    xen = np.array([[0.0, 0.0], [10.0, 0.0]])
    result = compute_qc(he, xen, pixel_um=0.2125)
    assert result.get("status") == "no_nuclei"
    assert result["n_he"] == 0


def test_compute_qc_empty_xen():
    he = np.array([[0.0, 0.0], [10.0, 0.0]])
    xen = np.empty((0, 2), dtype=float)
    result = compute_qc(he, xen, pixel_um=0.2125)
    assert result.get("status") == "no_nuclei"
    assert result["n_xenium"] == 0


def test_compute_qc_runs_with_valid_data():
    rng = np.random.default_rng(42)
    pts = rng.uniform(0, 1000, (200, 2))
    result = compute_qc(pts, pts + rng.uniform(-1, 1, (200, 2)), pixel_um=0.2125)
    assert "density_r" in result
    assert "nucleus_coincidence" in result
    assert isinstance(result["density_r"], float)
