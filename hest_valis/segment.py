"""StarDist nuclei detection on an H&E whole-slide image (tiled, level 0).

Returns nucleus centroids in H&E pixel space. These are then warped into the Xenium frame
(registration.warp_points) and compared to the Xenium nuclei in the concordance QC.

Reads the WSI with tifffile + zarr (no openslide dependency). Runs the StarDist
2D_versatile_he model; a GPU makes it much faster but it works on CPU too.
"""
import numpy as np
import tifffile
import zarr


def _level0_zarr(wsi_path):
    tf = tifffile.TiffFile(wsi_path)
    za = zarr.open(tf.series[0].aszarr(), mode="r")
    # multiscale series open as a group keyed '0'..'n'; a plain array has .shape
    if hasattr(za, "shape"):
        return za
    return za["0"]


def segment_he(wsi_path, tile=2048, overlap=128, prob_blank=235, model=None):
    """Detect nuclei on an H&E WSI. Returns (N, 2) array of (x, y) centroids in level-0 px.

    tile/overlap : tiling for memory; nuclei whose centroid falls in the tile core are kept.
    prob_blank   : skip tiles whose mean intensity is above this (near-white background).
    model        : a loaded StarDist2D model; if None, loads '2D_versatile_he'.
    """
    from csbdeep.utils import normalize
    from skimage.measure import regionprops_table
    if model is None:
        from stardist.models import StarDist2D
        model = StarDist2D.from_pretrained("2D_versatile_he")

    z = _level0_zarr(wsi_path)
    H, W = z.shape[0], z.shape[1]
    rows = []
    for y0 in range(0, H, tile):
        for x0 in range(0, W, tile):
            ya, xa = max(0, y0 - overlap), max(0, x0 - overlap)
            yb, xb = min(H, y0 + tile + overlap), min(W, x0 + tile + overlap)
            t = np.asarray(z[ya:yb, xa:xb])[..., :3]
            if t.mean() > prob_blank or t.shape[0] < 16 or t.shape[1] < 16:
                continue
            lab, _ = model.predict_instances(normalize(t, 1, 99.8, axis=(0, 1)),
                                             n_tiles=(1, 1, 1), show_tile_progress=False, verbose=False)
            if lab.max() == 0:
                continue
            p = regionprops_table(lab, properties=("centroid",))
            cy, cx = p["centroid-0"] + ya, p["centroid-1"] + xa
            core = (cy >= y0) & (cy < min(y0 + tile, H)) & (cx >= x0) & (cx < min(x0 + tile, W))
            rows.extend(zip(cx[core], cy[core]))
    return np.asarray(rows, dtype=float)
