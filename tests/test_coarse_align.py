"""Unit tests for hest_valis.coarse_align (rotation convention + rescue recovery).

Pure numpy; no valis / stardist / GPU. Guards the cardinal-rotation convention because it is on
the run_rescue critical path: a wrong convention would silently mis-orient a rescued slide.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest
from hest_valis import coarse_align

W, H = 3000, 2000
SCALE = 0.262774 / 0.2125
# inverse cardinal rotation that undoes a clockwise rotation by r0
INV = {0: 0, 90: 270, 180: 180, 270: 90}


def _tissue(n=6000, seed=0):
    """Asymmetric, chiral tissue so rotations and flips are distinguishable."""
    rng = np.random.RandomState(seed)
    blob = rng.normal([2200, 650], [170, 110], (n, 2))               # off-centre dense blob
    arm = np.c_[rng.uniform(400, 2200, n // 3), rng.normal(1450, 35, n // 3)]  # horizontal arm
    return np.r_[blob, arm]


def test_rotate_points_roundtrip():
    """Four 90-degree rotations (with the dimension swap) return to the start."""
    p0 = _tissue()
    p1 = coarse_align.rotate_points(p0, 90, W, H)
    p2 = coarse_align.rotate_points(p1, 90, H, W)
    p3 = coarse_align.rotate_points(p2, 90, W, H)
    p4 = coarse_align.rotate_points(p3, 90, H, W)
    assert np.allclose(p4, p0)


@pytest.mark.parametrize("r0", [0, 90, 180, 270])
def test_cardinal_rotation_recovers(r0):
    """An H&E rotated clockwise by r0 must be recovered by the inverse rotation, with high r."""
    xen = _tissue()
    he_src = coarse_align.rotate_points(xen, r0, W, H) / SCALE          # H&E px, rotated by r0
    he_wh = (H, W) if r0 in (90, 270) else (W, H)                       # dims after rotation
    rot, r = coarse_align.cardinal_rotation(he_src, xen, SCALE, he_wh)
    assert rot == INV[r0]
    assert r > 0.9


def test_cardinal_rotation_rejects_wrong():
    """A 90-degree-off candidate should score clearly worse than the correct one."""
    xen = _tissue()
    he_src = xen / SCALE                                                # already aligned -> 0 deg
    rot, r = coarse_align.cardinal_rotation(he_src, xen, SCALE, (W, H))
    assert rot == 0 and r > 0.9


def test_coarse_align_recovers_flip():
    """coarse_align searches flip x angle, so a mirrored + rotated H&E is recovered."""
    xen = _tissue()
    he = xen.copy()
    he[:, 0] = xen[:, 0].max() - he[:, 0]                               # mirror (x-flip)
    th = np.deg2rad(30)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    he = (he @ R.T) / SCALE
    _, params, r = coarse_align.coarse_align(he, xen, SCALE, angle_step=6)
    assert r > 0.8
    assert params["flip"] is True
