"""Load Xenium nuclei/cells into the DAPI pixel frame, and build a tissue mask.

The Xenium DAPI image and the Xenium cell table share the same micron coordinate system
(origin at the DAPI top-left), so cell-centroid pixels = centroid_um / pixel_um.
"""
import numpy as np
import pandas as pd
import tifffile
import zarr


# common column names for cell centroids across Xenium / HEST exports
_X_COLS = ["x_centroid", "x_um", "he_x", "x"]
_Y_COLS = ["y_centroid", "y_um", "he_y", "y"]


def load_xenium_nuclei(cells_path, pixel_um, in_um=True):
    """Load Xenium cell centroids -> (N, 2) array in DAPI pixels.

    cells_path : a parquet with centroid columns (Xenium cells.parquet uses x_centroid/y_centroid
                 in microns; HEST nucleus-centroid parquets use x_um/y_um or he_x/he_y).
    pixel_um   : DAPI um per pixel (e.g. 0.2125 for Xenium).
    in_um      : True if the centroid columns are microns (divide by pixel_um to get pixels);
                 False if they are already pixels.
    """
    df = pd.read_parquet(cells_path)
    xc = next(c for c in _X_COLS if c in df.columns)
    yc = next(c for c in _Y_COLS if c in df.columns)
    xy = np.c_[df[xc].to_numpy(float), df[yc].to_numpy(float)]
    if in_um:
        xy = xy / pixel_um
    return xy


def tissue_mask_from_dapi(dapi_path, level=4, thresh_frac=0.02):
    """Coarse tissue mask from the DAPI image (signal above background).

    Returns (mask, mask_pixel_um_scale_factor) where the mask is at a downsampled level.
    The caller passes mask_pixel_um = dapi_pixel_um * (full_dim / mask_dim).
    """
    tf = tifffile.TiffFile(dapi_path)
    s = tf.series[0]
    za = zarr.open(s.aszarr(), mode="r")
    arr = za if hasattr(za, "shape") else za["0"]
    # DAPI may be (C, Y, X); take channel 0
    full_h = arr.shape[-2]
    step = max(1, full_h // 4000)
    img = np.asarray(arr[0, ::step, ::step]) if arr.ndim == 3 else np.asarray(arr[::step, ::step])
    thr = np.percentile(img, 100 * (1 - thresh_frac))
    mask = img > max(1, thr * 0.15)
    return mask, step
