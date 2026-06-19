"""Unit tests for hest_valis.select.choose()."""
import math
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from hest_valis import select


def _m(median_um, density_r):
    """Build a minimal metric dict for choose()."""
    return {
        "nucleus_coincidence": {"median_um": median_um},
        "density_r": density_r,
    }


def test_both_none():
    r = select.choose(None, None)
    assert r["chosen"] is None
    assert r["rule"] == "no_data"


def test_micro_none():
    r = select.choose(None, _m(3.0, 0.7))
    assert r["chosen"] == "nomicro"
    assert r["rule"] == "micro_unavailable"


def test_nomicro_none():
    r = select.choose(_m(3.0, 0.7), None)
    assert r["chosen"] == "micro"
    assert r["rule"] == "nomicro_unavailable"


def test_lower_um_wins():
    # micro is clearly better on median_um and density_r difference is not alarming
    r = select.choose(_m(2.0, 0.8), _m(4.0, 0.75))
    assert r["chosen"] == "micro"
    assert r["rule"] == "lower_um"


def test_tie_picks_higher_density_r():
    # medians within TIE_UM=0.15; nomicro has higher density_r
    r = select.choose(_m(3.00, 0.6), _m(3.10, 0.9))
    assert r["chosen"] == "nomicro"
    assert "tie" in r["rule"]


def test_overfit_guard_flips():
    # micro wins on um by 0.3 (< GUARD_UM=0.5) but density_r is 0.20 below nomicro (> GUARD_R=0.10)
    r = select.choose(_m(2.0, 0.50), _m(2.3, 0.75))
    assert r["chosen"] == "nomicro"
    assert r["rule"] == "overfit_guard_flip"


def test_nan_density_r_does_not_win():
    # A slide with NaN density_r should not beat a good slide
    r = select.choose(_m(1.0, float("nan")), _m(3.0, 0.8))
    # micro has lower um but NaN density_r; guard should not kick in; it wins on um
    # (NaN -> -1; overfit guard: lR-wR = 0.8 - (-1) = 1.8 > GUARD_R, um diff = 2.0 > GUARD_UM=0.5 -> no flip)
    assert r["chosen"] == "micro"
    assert r["rule"] == "lower_um"


def test_one_none_returns_other():
    r = select.choose(None, _m(5.0, 0.5))
    assert r["chosen"] == "nomicro"
    r2 = select.choose(_m(5.0, 0.5), None)
    assert r2["chosen"] == "micro"
