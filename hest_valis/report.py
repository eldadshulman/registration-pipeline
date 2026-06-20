"""Non-interactive per-slide alignment QC report (PDF).

Renders the same diagnostics a human would flip through in a notebook
(cf. Tyler Lam's check_alignment.ipynb, Yang's align_tissue_block.ipynb), but
batch-friendly: Agg backend, no plt.show / input(), runnable inside a SLURM array.

It READS what the pipeline already produced -- it never re-registers and only reloads a
whole-slide image when the warped raster already exists on disk:

  <out>/<sample>/qc.json                  metrics + chosen protocol (run_qc.py)
  <out>/<sample>/he_nuclei_<chosen>.npy   warped H&E nuclei for the chosen variant
  Xenium/DAPI centroids                    passed in as xen_px (driver loads via xenium.py)
  <out>/<sample>/registered/*.ome.tif*     warped H&E raster (run_wsi.py) -- optional
  <out>/<sample>/cell_labels.parquet       per-cell region labels (run_annotate.py) -- optional

Panels (one page):
  POINT-BASED (always, from the nuclei arrays + qc.json):
    - DAPI vs warped-H&E centroid scatter, matched vs unmatched
    - displacement quiver of matched pairs
    - histogram of matched NN distances (um), median line = nucleus_coincidence.median_um
    - the two binned nucleus-density maps (DAPI, warped H&E), same bins density_r uses
    - text banner (density_r, median_um, occupancy, density_collapse, status, chosen + rule)
  RASTER (only if a warped raster exists): magenta/green DAPI-vs-H&E overlay + sitk checkerboard
  ANNOTATION (only if cell_labels.parquet exists): cells by region label + background-cell panel

All matching/binning is reused from hest_valis.concordance (mutual_nn_pairs, density_grids),
not reimplemented, so the report and the QC metrics agree by construction.
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")                       # headless; never plt.show()
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42           # keep PDF text editable (TrueType, not paths)
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from hest_valis import concordance

# region colours mirror run_annotate.COL so the report and the overlay agree
REGION_COL = {"high_density": "#c0392b", "low_density": "#2c7fb8", "background": "#dddddd"}
# chosen-protocol colours for the cohort triage page
PROTOCOL_COL = {"micro": "#1b9e77", "nomicro": "#7570b3", "coarse": "#d95f02",
                "rescued": "#e7298a", None: "#999999"}
QC_DENSITY_R = 0.5      # soft "good alignment" line on the cohort density-r plot (red flag below 0)
QC_MEDIAN_UM = 10.0     # one-cell target on the cohort median-um plot

_RASTER_MAX_DIM = 1500  # downsample rasters so a panel is cheap to draw


# --------------------------------------------------------------------------------------
# qc.json / array readers
# --------------------------------------------------------------------------------------
def read_qc(sample_dir):
    """Load <sample_dir>/qc.json (returns {} if absent)."""
    p = os.path.join(sample_dir, "qc.json")
    if not os.path.exists(p):
        return {}
    with open(p) as f:
        return json.load(f)


def chosen_protocol(qc):
    """Return (chosen, decision_dict). chosen is one of micro/nomicro/coarse/rescued or None."""
    dec = qc.get("decision", {}) or {}
    return dec.get("chosen"), dec


_NUCLEI_NAME = {"micro": "he_nuclei_micro.npy", "nomicro": "he_nuclei_nomicro.npy",
                "coarse": "he_nuclei_coarse.npy", "rescued": "he_nuclei_rescued.npy"}


def chosen_metric(qc, chosen):
    """The metric dict for the chosen variant (or None)."""
    if not chosen:
        return None
    return (qc.get("metrics", {}) or {}).get(chosen)


def load_chosen_nuclei(sample_dir, chosen):
    """Warped H&E nuclei for the chosen variant; empty (0,2) array if missing/unknown."""
    name = _NUCLEI_NAME.get(chosen)
    if name:
        p = os.path.join(sample_dir, name)
        if os.path.exists(p):
            return np.load(p).astype(float)
    return np.empty((0, 2), dtype=float)


def find_raster(sample_dir):
    """Path to the warped H&E raster from run_wsi.py, or None."""
    hits = sorted(glob.glob(os.path.join(sample_dir, "registered", "*.ome.tif*")))
    return hits[0] if hits else None


def find_annotation(sample_dir):
    """Path to the per-cell region labels from run_annotate.py, or None."""
    p = os.path.join(sample_dir, "cell_labels.parquet")
    return p if os.path.exists(p) else None


# --------------------------------------------------------------------------------------
# small array helpers
# --------------------------------------------------------------------------------------
def _norm01(img):
    img = np.asarray(img, dtype=float)
    lo, hi = np.percentile(img, 1), np.percentile(img, 99)
    if hi <= lo:
        hi = img.max() if img.max() > lo else lo + 1.0
    return np.clip((img - lo) / (hi - lo), 0, 1)


def _resize_nn(img, out_hw):
    """Nearest-neighbour resize (pure numpy, no skimage) to out_hw=(H,W)."""
    h, w = img.shape[:2]
    oh, ow = out_hw
    yi = np.clip((np.arange(oh) * h / oh).astype(int), 0, h - 1)
    xi = np.clip((np.arange(ow) * w / ow).astype(int), 0, w - 1)
    return img[yi][:, xi]


def _load_raster_2d(path, max_dim=_RASTER_MAX_DIM):
    """Downsampled 2D intensity from a (possibly pyramidal) TIFF/OME-TIFF.

    Handles RGB H&E (Y,X,3) and multi-channel DAPI (C,Y,X) -- channel 0 for DAPI, mean for RGB.
    Imports tifffile/zarr lazily so the unit test does not need them.
    """
    import tifffile
    import zarr
    with tifffile.TiffFile(path) as tf:
        za = zarr.open(tf.series[0].aszarr(), mode="r")
        arr = za if hasattr(za, "shape") else za[0]
        if arr.ndim == 3 and arr.shape[-1] in (3, 4):          # (Y, X, RGB[A])
            step = max(1, arr.shape[0] // max_dim)
            img = np.asarray(arr[::step, ::step, :3]).mean(-1)
        elif arr.ndim == 3:                                     # (C, Y, X) -> channel 0
            step = max(1, arr.shape[1] // max_dim)
            img = np.asarray(arr[0, ::step, ::step])
        else:                                                   # (Y, X)
            step = max(1, arr.shape[0] // max_dim)
            img = np.asarray(arr[::step, ::step])
    return img.astype(float)


# --------------------------------------------------------------------------------------
# point-based panels (always rendered)
# --------------------------------------------------------------------------------------
def panel_match_scatter(ax, xen_px, he_px, ib, um):
    """DAPI nuclei (context) + warped-H&E nuclei coloured matched vs unmatched."""
    matched = np.zeros(len(he_px), dtype=bool)
    matched[ib] = True
    ax.scatter(xen_px[:, 0], xen_px[:, 1], s=1, c="#cccccc", alpha=0.5, lw=0, label="DAPI")
    ax.scatter(he_px[~matched, 0], he_px[~matched, 1], s=2, c="#d81b60", alpha=0.6, lw=0,
               label=f"H&E unmatched ({(~matched).sum():,})")
    ax.scatter(he_px[matched, 0], he_px[matched, 1], s=2, c="#1b9e3f", alpha=0.6, lw=0,
               label=f"H&E matched ({matched.sum():,})")
    ax.set_title("nucleus matching (mutual-NN, <=20um)", fontsize=9)
    ax.legend(fontsize=6, markerscale=3, loc="upper right", framealpha=0.85)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])


def panel_quiver(ax, xen_px, he_px, ia, ib, um, max_arrows=2000):
    """Residual displacement DAPI->warped-H&E for matched pairs (shift/rotation at a glance)."""
    if len(ia) == 0:
        ax.text(0.5, 0.5, "no matched pairs", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("displacement", fontsize=9); ax.set_xticks([]); ax.set_yticks([]); return
    src = xen_px[ia]
    dxy = (he_px[ib] - src)
    if len(ia) > max_arrows:                                   # thin for legibility only
        sel = np.linspace(0, len(ia) - 1, max_arrows).astype(int)
        src, dxy = src[sel], dxy[sel]
    mag_um = np.hypot(dxy[:, 0], dxy[:, 1]) * um
    # angles/scale_units="xy", scale=1 draws src -> src+(U,V) in DATA coords; invert_yaxis()
    # below flips the display, so pass the true dxy (do NOT negate V) for correct directions.
    q = ax.quiver(src[:, 0], src[:, 1], dxy[:, 0], dxy[:, 1], mag_um,
                  cmap="viridis", angles="xy", scale_units="xy", scale=1,
                  width=0.003, clim=(0, 20))
    plt.colorbar(q, ax=ax, fraction=0.046, pad=0.02, label="um")
    ax.set_title("matched-pair displacement", fontsize=9)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])


def panel_dist_hist(ax, dist_um, median_um):
    """Histogram of matched NN distances (um) with the median line from qc.json."""
    if len(dist_um):
        ax.hist(dist_um, bins=40, range=(0, 20), color="#4c72b0", alpha=0.85)
    med = median_um if median_um is not None else (float(np.median(dist_um)) if len(dist_um) else None)
    if med is not None:
        ax.axvline(med, color="#c0392b", lw=1.5, label=f"median {med:.2f} um")
        ax.legend(fontsize=7)
    ax.set_xlabel("matched NN distance (um)", fontsize=8)
    ax.set_ylabel("pairs", fontsize=8)
    ax.set_title("offset distribution", fontsize=9)


def panel_density(ax, grid, title):
    """One binned nucleus-density map (same 50um bins density_r uses)."""
    ax.imshow(np.log1p(grid), origin="upper", cmap="magma", aspect="equal")
    ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])


def _fmt(v, nd=3):
    return "n/a" if v is None else (f"{v:.{nd}f}" if isinstance(v, float) else str(v))


def panel_banner(ax, sample_id, chosen, decision, metric, status):
    """Text banner: the headline numbers + the decision that was made."""
    ax.axis("off")
    nc = (metric or {}).get("nucleus_coincidence", {}) or {}
    neg = (metric or {}).get("negative_control_100um", {}) or {}
    lines = [
        f"{sample_id}",
        f"chosen: {chosen}   ({decision.get('rule', 'n/a')})",
        f"status: {status or 'ok'}",
        "",
        f"density_r:        {_fmt((metric or {}).get('density_r'))}",
        f"median offset:    {_fmt(nc.get('median_um'))} um",
        f"occupancy:        {_fmt((metric or {}).get('occupancy'))}",
        f"density_collapse: {_fmt(neg.get('density_collapse'))}   (neg-control, want large +)",
        "",
        f"n H&E nuclei:     {_fmt((metric or {}).get('n_he'), 0)}",
        f"n Xenium cells:   {_fmt((metric or {}).get('n_xenium'), 0)}",
        f"matched:          {_fmt(nc.get('n_matched'), 0)}"
        f"  ({_fmt(nc.get('frac_matched'))} of cells)",
    ]
    if chosen in ("coarse", "rescued"):
        lines += ["", "(rescued alignment; raster overlay may be absent)"]
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace",
            fontsize=9, transform=ax.transAxes)


# --------------------------------------------------------------------------------------
# raster panels (only if a warped raster exists)
# --------------------------------------------------------------------------------------
def panel_raster_overlay(ax, dapi2d, he2d):
    """Magenta (DAPI) / green (warped H&E tissue) overlay -- NOT red/blue."""
    h = min(dapi2d.shape[0], he2d.shape[0]); w = min(dapi2d.shape[1], he2d.shape[1])
    d = _norm01(_resize_nn(dapi2d, (h, w)))
    e = 1.0 - _norm01(_resize_nn(he2d, (h, w)))               # invert: H&E tissue -> high
    rgb = np.zeros((h, w, 3))
    rgb[..., 0] = d; rgb[..., 2] = d                          # magenta = DAPI
    rgb[..., 1] = e                                           # green   = H&E
    ax.imshow(rgb, origin="upper")
    ax.set_title("DAPI (magenta) vs warped H&E (green)", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])


def panel_checkerboard(ax, dapi2d, he2d, tiles=8):
    """sitk checkerboard mosaic; misalignment shows as breaks at tile seams."""
    h = min(dapi2d.shape[0], he2d.shape[0]); w = min(dapi2d.shape[1], he2d.shape[1])
    a = _norm01(_resize_nn(dapi2d, (h, w)))
    b = 1.0 - _norm01(_resize_nn(he2d, (h, w)))
    try:
        import SimpleITK as sitk
        ia, ib = sitk.GetImageFromArray(a.astype(np.float32)), sitk.GetImageFromArray(b.astype(np.float32))
        cb = sitk.CheckerBoard(ia, ib, [tiles, tiles])
        board = sitk.GetArrayFromImage(cb)
    except Exception:                                          # numpy fallback if sitk absent
        ty, tx = max(1, h // tiles), max(1, w // tiles)
        m = ((np.arange(h)[:, None] // ty) + (np.arange(w)[None, :] // tx)) % 2 == 0
        board = np.where(m, a, b)
    ax.imshow(board, origin="upper", cmap="gray")
    ax.set_title("checkerboard (DAPI / H&E)", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])


# --------------------------------------------------------------------------------------
# annotation panels (only if cell_labels.parquet exists)
# --------------------------------------------------------------------------------------
def panel_region_cells(ax, df, um):
    """Xenium cells coloured by transferred H&E region label."""
    xp, yp = df["x_um"].to_numpy() / um, df["y_um"].to_numpy() / um
    reg = df["he_region"].to_numpy()
    for k, c in REGION_COL.items():
        m = reg == k
        if m.any():
            ax.scatter(xp[m], yp[m], s=1, c=c, lw=0, alpha=0.6, label=f"{k} ({m.sum():,})")
    ax.legend(fontsize=6, markerscale=4, loc="upper right", framealpha=0.85)
    ax.set_title("cells by transferred region", fontsize=9)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])


def panel_background_cells(ax, df, um):
    """Background cells highlighted (cells landing off tissue = a registration smell)."""
    xp, yp = df["x_um"].to_numpy() / um, df["y_um"].to_numpy() / um
    bg = df["he_region"].to_numpy() == "background"
    ax.scatter(xp[~bg], yp[~bg], s=1, c="#dddddd", lw=0, alpha=0.5, label="on tissue")
    ax.scatter(xp[bg], yp[bg], s=2, c="#000000", lw=0, alpha=0.7,
               label=f"background ({bg.sum():,})")
    ax.legend(fontsize=6, markerscale=4, loc="upper right", framealpha=0.85)
    ax.set_title("background (off-tissue) cells", fontsize=9)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])


# --------------------------------------------------------------------------------------
# top-level: per-slide figure + PDF
# --------------------------------------------------------------------------------------
def build_sample_figure(sample_dir, xen_px, pixel_um, sample_id=None, dapi_path=None):
    """Build (and return) the one-page report Figure for a slide. Pure rendering; reads only
    files already in sample_dir (+ a raster reload only when registered/*.ome.tif* exists)."""
    sample_id = sample_id or os.path.basename(os.path.normpath(sample_dir))
    um = float(pixel_um)
    qc = read_qc(sample_dir)
    chosen, decision = chosen_protocol(qc)
    metric = chosen_metric(qc, chosen)
    status = (metric or {}).get("status") if metric else ("no_data" if not chosen else None)
    he_px = load_chosen_nuclei(sample_dir, chosen)
    xen_px = np.asarray(xen_px, dtype=float).reshape(-1, 2)

    degenerate = (chosen is None or metric is None or status == "no_nuclei"
                  or len(he_px) == 0 or len(xen_px) == 0)

    # --- degenerate (no_nuclei / no_data): single banner page, never crash ---
    if degenerate:
        fig = plt.figure(figsize=(7, 4))
        ax = fig.add_subplot(111)
        status = status or "no_nuclei"
        panel_banner(ax, sample_id, chosen, decision, metric, status)
        ax.text(0.0, 0.02, "point/raster panels skipped (no usable nuclei)",
                transform=ax.transAxes, fontsize=8, color="#888888")
        fig.suptitle(f"{sample_id} -- alignment QC report", fontsize=11)
        fig.tight_layout()
        return fig

    # --- matching + binning, reused from concordance (single source of truth) ---
    cutoff_px = concordance.CUTOFF_UM / um
    ia, ib, dpx = concordance.mutual_nn_pairs(xen_px, he_px, cutoff_px)
    dist_um = dpx * um
    gh, gx, _foot, _nbx, _nby = concordance.density_grids(he_px, xen_px, concordance.BIN_UM / um)
    median_um = (metric.get("nucleus_coincidence", {}) or {}).get("median_um")

    # --- optional inputs ---
    raster = find_raster(sample_dir)
    have_raster = bool(raster and dapi_path and os.path.exists(dapi_path))
    dapi2d = he2d = None
    if have_raster:
        try:
            dapi2d, he2d = _load_raster_2d(dapi_path), _load_raster_2d(raster)
        except Exception as e:
            print(f"[{sample_id}] raster panels skipped ({e})", flush=True)
            have_raster = False

    ann_path = find_annotation(sample_dir)
    ann_df = None
    if ann_path:
        try:
            import pandas as pd
            ann_df = pd.read_parquet(ann_path)
        except Exception as e:
            print(f"[{sample_id}] annotation panel skipped ({e})", flush=True)

    # --- compose grid: 2 base rows (+raster, +annotation) ---
    nrows = 2 + (1 if have_raster else 0) + (1 if ann_df is not None else 0)
    fig = plt.figure(figsize=(15, 4.7 * nrows))
    gs = GridSpec(nrows, 3, figure=fig, hspace=0.28, wspace=0.18)

    panel_match_scatter(fig.add_subplot(gs[0, 0]), xen_px, he_px, ib, um)
    panel_quiver(fig.add_subplot(gs[0, 1]), xen_px, he_px, ia, ib, um)
    panel_dist_hist(fig.add_subplot(gs[0, 2]), dist_um, median_um)
    panel_density(fig.add_subplot(gs[1, 0]), gx, "Xenium/DAPI nucleus density (50um bins)")
    panel_density(fig.add_subplot(gs[1, 1]), gh, "warped H&E nucleus density (50um bins)")
    panel_banner(fig.add_subplot(gs[1, 2]), sample_id, chosen, decision, metric, status)

    row = 2
    if have_raster:
        panel_raster_overlay(fig.add_subplot(gs[row, 0]), dapi2d, he2d)
        panel_checkerboard(fig.add_subplot(gs[row, 1]), dapi2d, he2d)
        fig.add_subplot(gs[row, 2]).axis("off")
        row += 1
    if ann_df is not None:
        panel_region_cells(fig.add_subplot(gs[row, 0]), ann_df, um)
        panel_background_cells(fig.add_subplot(gs[row, 1]), ann_df, um)
        fig.add_subplot(gs[row, 2]).axis("off")

    note = "" if have_raster else "   (no warped raster yet -> no image overlay)"
    fig.suptitle(f"{sample_id} -- alignment QC report  [chosen={chosen}, density_r="
                 f"{_fmt(metric.get('density_r'))}]{note}", fontsize=12)
    return fig


def render_sample(sample_dir, xen_px, pixel_um, sample_id=None, dapi_path=None,
                  out_pdf=None, save_png=False):
    """Build the slide report and write <sample_dir>/report.pdf. Returns the PDF path.

    save_png=True also writes report.png (a thin viewer / the review notebook displays it).
    """
    out_pdf = out_pdf or os.path.join(sample_dir, "report.pdf")
    fig = build_sample_figure(sample_dir, xen_px, pixel_um, sample_id=sample_id, dapi_path=dapi_path)
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    if save_png:
        fig.savefig(os.path.splitext(out_pdf)[0] + ".png", dpi=110,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_pdf


# --------------------------------------------------------------------------------------
# top-level: cohort triage page
# --------------------------------------------------------------------------------------
def collect_cohort(output_dir):
    """Read every <output_dir>/*/qc.json into a list of per-slide rows for the triage page."""
    rows = []
    for qp in sorted(glob.glob(os.path.join(output_dir, "*", "qc.json"))):
        sample = os.path.basename(os.path.dirname(qp))
        try:
            with open(qp) as f:
                qc = json.load(f)
        except Exception:
            continue
        chosen, dec = chosen_protocol(qc)
        metric = chosen_metric(qc, chosen) or {}
        neg = metric.get("negative_control_100um", {}) or {}
        rows.append({
            "sample": sample,
            "chosen": chosen,
            "density_r": dec.get("sel_density_r", metric.get("density_r")),
            "median_um": dec.get("sel_median_um",
                                 (metric.get("nucleus_coincidence", {}) or {}).get("median_um")),
            "density_collapse": neg.get("density_collapse"),
            "status": metric.get("status", "ok"),
        })
    return rows


def _bar(ax, rows, key, title, threshold=None, lower_better=False, xlabel=""):
    vals = [(r["sample"], r[key], r["chosen"]) for r in rows if isinstance(r.get(key), (int, float))]
    vals.sort(key=lambda t: (t[1] if t[1] is not None else 0), reverse=not lower_better)
    if not vals:
        ax.text(0.5, 0.5, f"no {key}", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title, fontsize=10); return
    names = [v[0] for v in vals]
    y = np.arange(len(vals))
    ax.barh(y, [v[1] for v in vals], color=[PROTOCOL_COL.get(v[2], "#999999") for v in vals])
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=6); ax.invert_yaxis()
    if threshold is not None:
        ax.axvline(threshold, color="#c0392b", ls="--", lw=1, label=f"threshold {threshold:g}")
        ax.legend(fontsize=7)
    ax.set_title(title, fontsize=10); ax.set_xlabel(xlabel, fontsize=8)
    ax.grid(axis="x", alpha=0.25)


def render_cohort(output_dir, out_pdf=None):
    """Cohort triage page -> <output_dir>/cohort_report.pdf. Returns the path."""
    rows = collect_cohort(output_dir)
    out_pdf = out_pdf or os.path.join(output_dir, "cohort_report.pdf")
    n = max(1, len(rows))
    fig = plt.figure(figsize=(16, max(4.0, 0.22 * n + 1.5)))
    gs = GridSpec(1, 3, figure=fig, wspace=0.5)
    _bar(fig.add_subplot(gs[0, 0]), rows, "density_r", "density_r (ranked)",
         threshold=QC_DENSITY_R, xlabel="density_r")
    _bar(fig.add_subplot(gs[0, 1]), rows, "median_um", "median offset (ranked, lower better)",
         threshold=QC_MEDIAN_UM, lower_better=True, xlabel="um")
    _bar(fig.add_subplot(gs[0, 2]), rows, "density_collapse",
         "neg-control density collapse", xlabel="density_r drop")
    handles = [plt.Line2D([0], [0], marker="s", ls="", color=c, label=str(k))
               for k, c in PROTOCOL_COL.items() if k is not None]
    fig.legend(handles=handles, fontsize=8, loc="lower center", ncol=4, title="chosen protocol")
    fig.suptitle(f"cohort alignment QC -- {len(rows)} slides", fontsize=13)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.92, bottom=0.12, wspace=0.5)
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_pdf
