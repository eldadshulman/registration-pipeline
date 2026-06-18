"""Concordance QC: does the registered H&E agree with the Xenium molecular data?

Four checks, all computed from two point sets in the Xenium/DAPI pixel frame:
  he_px  : H&E nuclei (StarDist) warped into the Xenium frame
  xen_px : Xenium nuclei (cell centroids from the Xenium output)

  1. nucleus coincidence  -> median offset (um) between mutually-nearest nuclei. Lower better.
  2. density correlation  -> Pearson r of H&E vs Xenium nucleus-density maps (50 um bins),
                             restricted to the tissue footprint.
  3. tissue occupancy     -> fraction of Xenium cells that fall on tissue (needs a tissue mask).
  4. negative control     -> shift the H&E nuclei by +/-100 um; a trustworthy density-r and
                             match-rate MUST collapse. The 'collapse' is density_r minus the
                             worst shifted density_r.

compute_qc() returns a JSON-serialisable dict.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from scipy.ndimage import binary_dilation

CUTOFF_UM = 20.0     # max distance to count two nuclei as a coincidence match
BIN_UM = 50.0        # density-map bin size
NEG_SHIFT_UM = 100.0 # negative-control displacement


def _mutual_nn(A, B, cutoff_px):
    ta, tb = cKDTree(A), cKDTree(B)
    dab, nab = tb.query(A, k=1)
    _, nba = ta.query(B, k=1)
    ia = np.arange(len(A))
    keep = (nba[nab] == ia) & (dab <= cutoff_px)
    return dab[keep]


def _density(P, nby, nbx, bin_px):
    H, _, _ = np.histogram2d(P[:, 1], P[:, 0], bins=[nby, nbx],
                             range=[[0, nby * bin_px], [0, nbx * bin_px]])
    return H


def compute_qc(he_px, xen_px, pixel_um, tissue_mask=None, mask_pixel_um=None):
    """he_px, xen_px : (N, 2) nuclei in Xenium-frame pixels. pixel_um : um per pixel.

    tissue_mask    : optional 2D bool array (tissue vs background) for occupancy.
    mask_pixel_um  : um per pixel of tissue_mask (defaults to pixel_um).
    """
    um = pixel_um
    cutoff_px = CUTOFF_UM / um
    bin_px = BIN_UM / um
    R = {"pixel_um": um, "n_he": int(len(he_px)), "n_xenium": int(len(xen_px))}

    # 1. nucleus coincidence
    d = _mutual_nn(xen_px, he_px, cutoff_px) * um
    R["nucleus_coincidence"] = {
        "n_matched": int(len(d)),
        "frac_matched": float(len(d) / max(1, len(xen_px))),
        "median_um": float(np.median(d)) if len(d) else None,
        "pct_within_5um": float((d <= 5).mean() * 100) if len(d) else None,
        "pct_within_10um": float((d <= 10).mean() * 100) if len(d) else None,
    }

    # 2. density correlation (footprint-masked)
    xmax = max(he_px[:, 0].max(), xen_px[:, 0].max())
    ymax = max(he_px[:, 1].max(), xen_px[:, 1].max())
    nbx, nby = int(xmax / bin_px) + 1, int(ymax / bin_px) + 1
    gh, gx = _density(he_px, nby, nbx, bin_px), _density(xen_px, nby, nbx, bin_px)
    foot = binary_dilation(gh > 0, iterations=2)

    def dens_r(P):
        g = _density(P, nby, nbx, bin_px)
        m = foot & ((gx + g) > 0)
        return float(pearsonr(gx[m], g[m])[0]) if m.sum() >= 10 else float("nan")

    r = dens_r(he_px)
    R["density_r"] = round(r, 3)

    # 3. occupancy
    if tissue_mask is not None:
        mum = mask_pixel_um or um
        scale = um / mum
        def occ(P):
            xi = np.clip((P[:, 0] * scale).astype(int), 0, tissue_mask.shape[1] - 1)
            yi = np.clip((P[:, 1] * scale).astype(int), 0, tissue_mask.shape[0] - 1)
            return float(tissue_mask[yi, xi].mean())
        R["occupancy"] = round(occ(xen_px), 3)

    # 4. negative control (bidirectional)
    shifts = {"+x": (NEG_SHIFT_UM, 0), "-x": (-NEG_SHIFT_UM, 0),
              "+y": (0, NEG_SHIFT_UM), "-y": (0, -NEG_SHIFT_UM)}
    neg = {}
    for k, (dx, dy) in shifts.items():
        sh = np.array([dx / um, dy / um])
        Bs = he_px + sh
        dr = _mutual_nn(xen_px, Bs, cutoff_px) * um
        neg[k] = {"density_r": round(dens_r(Bs), 3),
                  "median_um": round(float(np.median(dr)), 2) if len(dr) else None,
                  "frac_matched": round(float(len(dr) / max(1, len(xen_px))), 3)}
    worst = max(v["density_r"] for v in neg.values())
    R["negative_control_100um"] = {"by_dir": neg, "worst_density_r": round(worst, 3),
                                   "density_collapse": round(r - worst, 3)}
    return R
