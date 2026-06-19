"""Coarse global-alignment fallback (self-healing).

VALIS feature matching can lock onto a wrong solution when the H&E is grossly mis-oriented
relative to the DAPI (e.g. a 90/180/270 degree rotation or a mirror). The slide then comes
back with a NEGATIVE density-r: locally nothing coincides even though the footprints roughly
overlap. No amount of rigid/non-rigid/reflection inside VALIS fixes this, because the starting
orientation is wrong.

This recovers the gross orientation directly from the nuclei: it searches rotation x flip and,
for each, finds the best translation by FFT phase correlation, scoring by the agreement of the
two nuclei-density maps. The best transform is applied to the H&E nuclei. It is meant as a
fallback that runs only when the normal registration returns a negative (or near-zero) density-r.

coarse_align() returns (aligned_he_px, params, density_r) with aligned_he_px in the Xenium/DAPI
pixel frame, so it plugs straight into the same concordance QC.
"""
import numpy as np
from numpy.fft import fft2, ifft2


def _grid(P, ext, nb):
    H, _, _ = np.histogram2d(P[:, 1], P[:, 0], bins=[nb, nb], range=[[-ext, ext], [-ext, ext]])
    return H


def coarse_align(he_src_px, xen_px, scale, dapi_um=0.2125, bin_um=40.0, angle_step=6):
    """he_src_px : ORIGINAL H&E nuclei in H&E pixels (pre-registration, from segment).
    xen_px      : Xenium nuclei in DAPI pixels.
    scale       : he_pixel_um / dapi_pixel_um (puts H&E nuclei on the DAPI pixel scale).
    Returns (aligned_he_px, {angle, flip, dx, dy}, density_r).
    """
    he = he_src_px * scale
    B = bin_um / dapi_um
    allp = np.vstack([he - he.mean(0), xen_px - xen_px.mean(0)])
    ext = np.abs(allp).max() * 1.15
    nb = int(2 * ext / B) + 1
    Xc = xen_px - xen_px.mean(0)
    Bx = _grid(Xc, ext, nb)
    Fx = fft2(Bx)

    best = (-2.0, None)
    for flip in (False, True):
        for ang in range(0, 360, angle_step):
            h = he - he.mean(0)
            if flip:
                h = h * np.array([-1, 1])
            th = np.deg2rad(ang)
            R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
            h = h @ R.T
            Fh = fft2(_grid(h, ext, nb))
            cps = Fx * np.conj(Fh)
            cps /= np.abs(cps) + 1e-9
            cc = np.abs(ifft2(cps))
            pk = np.unravel_index(np.argmax(cc), cc.shape)
            sy = pk[0] - (nb if pk[0] > nb // 2 else 0)
            sx = pk[1] - (nb if pk[1] > nb // 2 else 0)
            hp = h + np.array([sx * B, sy * B])
            gh = _grid(hp, ext, nb)
            m = (gh > 0) & ((gh + Bx) > 0)
            r = float(np.corrcoef(gh[m], Bx[m])[0, 1]) if m.sum() > 10 else -1.0
            if r > best[0]:
                best = (r, (ang, flip, sx * B, sy * B))

    r, (ang, flip, dx, dy) = best
    h = he - he.mean(0)
    if flip:
        h = h * np.array([-1, 1])
    th = np.deg2rad(ang)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    aligned = (h @ R.T) + np.array([dx, dy]) + xen_px.mean(0)
    return aligned, {"angle": ang, "flip": flip, "dx": float(dx), "dy": float(dy)}, r
