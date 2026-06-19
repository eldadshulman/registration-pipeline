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


# Target maximum image dimension for the tissue mask (pixels); coarser = faster.
_MASK_MAX_DIM = 4000

def tissue_mask_from_dapi(dapi_path, thresh_frac=0.02):
    """Coarse tissue mask from the DAPI image (signal above background).

    Downsamples so the longer axis is at most _MASK_MAX_DIM px.
    Returns (mask, step) where step = full_dim / mask_dim (integer downsampling factor).
    The caller passes mask_pixel_um = dapi_pixel_um * step.
    """
    tf = tifffile.TiffFile(dapi_path)
    s = tf.series[0]
    za = zarr.open(s.aszarr(), mode="r")
    arr = za if hasattr(za, "shape") else za["0"]
    # DAPI may be (C, Y, X); take channel 0
    full_h = arr.shape[-2]
    step = max(1, full_h // _MASK_MAX_DIM)
    img = np.asarray(arr[0, ::step, ::step]) if arr.ndim == 3 else np.asarray(arr[::step, ::step])
    thr = np.percentile(img, 100 * (1 - thresh_frac))
    mask = img > max(1, thr * 0.15)
    return mask, step


def load_xenium_cells(cells_path, pixel_um, in_um=True):
    """Load Xenium cells as a DataFrame for annotation transfer.

    Returns columns: cell_id, x_um, y_um, x_px, y_px (DAPI pixel frame).
    """
    df = pd.read_parquet(cells_path)
    xc = next(c for c in _X_COLS if c in df.columns)
    yc = next(c for c in _Y_COLS if c in df.columns)
    cid = next((c for c in ("cell_id", "cell", "id") if c in df.columns), None)
    x = df[xc].to_numpy(float)
    y = df[yc].to_numpy(float)
    x_um, y_um = (x, y) if in_um else (x * pixel_um, y * pixel_um)
    x_px, y_px = (x / pixel_um, y / pixel_um) if in_um else (x, y)
    out = pd.DataFrame({"x_um": x_um, "y_um": y_um, "x_px": x_px, "y_px": y_px})
    out.insert(0, "cell_id", df[cid].to_numpy() if cid else np.arange(len(df)))
    return out


# Fallback H&E um/pixel (Aperio 20x standard; override via config["he_pixel_um"]).
HE_FALLBACK_MPP = 0.2628


def he_pixel_um(he_path, fallback=HE_FALLBACK_MPP):
    """Read the H&E microns-per-pixel (Aperio 'MPP' tag); fallback if not present."""
    import re
    try:
        with tifffile.TiffFile(he_path) as tf:
            desc = tf.pages[0].description or ""
        m = re.search(r"MPP\s*=\s*([0-9.]+)", desc)
        if m:
            return float(m.group(1))
    except Exception:
        pass
    return fallback
