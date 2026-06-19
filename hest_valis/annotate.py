"""Per-cell annotation transfer: tag each Xenium cell with its H&E region.

After registration, build a morphological region map from the registered H&E
(high_density / low_density / background) and assign every Xenium cell the region it falls in.
This mirrors annotation-transfer pipelines that overlap each cell with an aligned mask, but the
mask here is derived from H&E morphology rather than a hand-drawn pathology annotation.

Region map: per `bin_um` bin, the registered H&E nuclear density. A 2-component GMM splits tissue
bins into the higher-density cluster ('high_density', often but not always tumour) and the lower
('low_density', often stroma); empty bins are 'background'. These are density-based, NOT
pathologist-validated (see the REGIONS note below). When real pathologist masks become available,
replace region_map() with a lookup into that mask; assign_cells() stays the same.
"""
import numpy as np
import pandas as pd
from scipy.ndimage import binary_closing
from sklearn.mixture import GaussianMixture

# NOTE: labels are density-based (H&E nuclear density GMM), not pathologist-validated.
# High-density tissue bins are labelled "high_density" (often but not always tumour);
# low-density bins are "low_density" (often stroma). Replace with a real mask lookup
# when pathologist annotations are available.
REGIONS = ("high_density", "low_density", "background")


def region_map(he_px, pixel_um, nbx, nby, bin_um=50.0):
    """Return (region grid [nby, nbx] of strings, density grid)."""
    bpx = bin_um / pixel_um
    dens, _, _ = np.histogram2d(he_px[:, 1], he_px[:, 0], bins=[nby, nbx],
                                range=[[0, nby * bpx], [0, nbx * bpx]])
    tissue = binary_closing(dens >= 1, iterations=2)
    region = np.full((nby, nbx), "background", dtype=object)
    if tissue.sum() >= 4:
        v = np.log1p(dens[tissue]).reshape(-1, 1)
        gm = GaussianMixture(2, random_state=0).fit(v)
        lab = gm.predict(v)
        tumor = int(np.argmax(gm.means_.ravel()))
        ti = np.argwhere(tissue)
        region[ti[:, 0], ti[:, 1]] = np.where(lab == tumor, "high_density", "low_density")
    return region, dens


def assign_cells(he_px, xen_cells, pixel_um, bin_um=50.0):
    """he_px    : registered H&E nuclei in DAPI pixels.
    xen_cells : DataFrame with cell_id, x_um, y_um, x_px, y_px (Xenium cells in the DAPI frame).
    Returns DataFrame: cell_id, x_um, y_um, he_region.
    """
    bpx = bin_um / pixel_um
    nbx = int(max(he_px[:, 0].max(), xen_cells.x_px.max()) / bpx) + 2
    nby = int(max(he_px[:, 1].max(), xen_cells.y_px.max()) / bpx) + 2
    region, _ = region_map(he_px, pixel_um, nbx, nby, bin_um)
    cx = np.clip((xen_cells.x_px.to_numpy() / bpx).astype(int), 0, nbx - 1)
    cy = np.clip((xen_cells.y_px.to_numpy() / bpx).astype(int), 0, nby - 1)
    return pd.DataFrame({"cell_id": xen_cells.cell_id.to_numpy(),
                         "x_um": xen_cells.x_um.to_numpy(),
                         "y_um": xen_cells.y_um.to_numpy(),
                         "he_region": region[cy, cx]})
