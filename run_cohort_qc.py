#!/usr/bin/env python
"""Cohort QC overlay PDF -- one figure per slide showing the SELECTED (final) registration only.

Runs on RESULTS (no re-registration): for each sample it reads the authoritative `qc.json`
decision, loads the chosen protocol's warped H&E nuclei (`he_nuclei_<chosen>.npy`) and the Xenium
cell centroids, and renders:

  * page 1  -- a cohort summary: the selected/final protocol per slide + density_r / median_um +
               disposition (accepted / manual-review), from the real gate (provenance.gate).
  * page k  -- one page per slide: the H&E nuclei (red) overlaid on the Xenium DAPI nuclei (grey),
               for the SELECTED protocol only. The losing protocol and any failed/rescued
               intermediate are NOT drawn -- this is the final alignment, per slide.

Usage:
    python run_cohort_qc.py --samples samples.csv --config config.json [--out cohort_qc.pdf]

The render core (`render_cohort_pdf`) is importable so other drivers (e.g. a legacy-format cohort)
can build the `slides` list themselves and reuse the identical page layout.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hest_valis import config as hvconfig
from hest_valis import provenance as hvprov
from hest_valis import xenium as hvxen

PLOT_CAP = 60000        # max points drawn per layer (metrics use ALL points; plotting is subsampled)


def disposition(chosen, r, u, th):
    """(label, color) from the real acceptance gate; rescue/coarse are flagged separately."""
    if chosen is None:
        return ("QC-FAILED", "#b0b0b0")
    if chosen not in hvprov.ASIS:                       # coarse / rescued (variant shown in protocol col)
        return ("rescued", "#3b6fb0")
    if hvprov.gate(r, u, th):
        return ("accepted", "#2e8b57")
    return ("manual-review", "#d98a30")


def _subsample(xy, cap=PLOT_CAP, seed=0):
    if xy is None or len(xy) <= cap:
        return xy
    idx = np.random.default_rng(seed).choice(len(xy), cap, replace=False)
    return xy[idx]


def render_cohort_pdf(slides, out_pdf, thresholds):
    """slides: list of dicts {sample_id, xen(Nx2), he(Mx2), protocol, density_r, median_um, note?}.
    Writes <out_pdf>: a summary page + one overlay page per slide. Returns the page count."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    th = thresholds
    rows = []
    for s in slides:
        lbl, col = disposition(s.get("protocol"), s.get("density_r"), s.get("median_um"), th)
        rows.append({**s, "_disp": lbl, "_col": col})

    n_acc = sum(r["_disp"] == "accepted" for r in rows)
    n_rev = sum(r["_disp"] == "manual-review" for r in rows)
    n_res = sum(r["_disp"].startswith("rescued") for r in rows)

    with PdfPages(out_pdf) as pdf:
        # ---------- page 1: summary table ----------
        fig, ax = plt.subplots(figsize=(8.5, 11)); ax.axis("off")
        ax.set_title("Cohort registration QC -- selected protocol per slide", fontsize=14, weight="bold", pad=14)
        ax.text(0.5, 0.955,
                f"{len(rows)} slides   |   accepted {n_acc}   manual-review {n_rev}   rescued {n_res}"
                f"      gate: density_r >= {th['density_r_accept']} AND median_um <= {th['median_um_accept']}",
                ha="center", va="top", transform=ax.transAxes, fontsize=9.5)
        header = ["slide", "protocol", "density_r", "median_um", "disposition"]
        cells = [[r["sample_id"], str(r.get("protocol")),
                  ("%.3f" % r["density_r"]) if r.get("density_r") is not None else "-",
                  ("%.2f" % r["median_um"]) if r.get("median_um") is not None else "-",
                  r["_disp"]] for r in rows]
        tbl = ax.table(cellText=cells, colLabels=header, loc="center", cellLoc="center",
                       bbox=[0.02, 0.02, 0.96, 0.90])
        tbl.auto_set_font_size(False); tbl.set_fontsize(8)
        for j in range(len(header)):
            tbl[0, j].set_facecolor("#333333"); tbl[0, j].set_text_props(color="w", weight="bold")
        for i, r in enumerate(rows, start=1):
            tbl[i, 4].set_facecolor(r["_col"]); tbl[i, 4].set_text_props(color="w")
        pdf.savefig(fig); plt.close(fig)

        # ---------- one overlay page per slide (selected protocol only) ----------
        for r in rows:
            fig, ax = plt.subplots(figsize=(8.5, 8.5))
            xen, he = r.get("xen"), r.get("he")
            if xen is not None and len(xen):
                xs = _subsample(xen); ax.scatter(xs[:, 0], xs[:, 1], s=1.0, c="0.72",
                                                 marker=".", linewidths=0, rasterized=True,
                                                 label=f"Xenium DAPI nuclei (n={len(xen):,})")
            if he is not None and len(he):
                hs = _subsample(he); ax.scatter(hs[:, 0], hs[:, 1], s=1.0, c="#c0392b",
                                                marker=".", linewidths=0, alpha=0.55, rasterized=True,
                                                label=f"H&E nuclei -- {r.get('protocol')} (n={len(he):,})")
            ax.set_aspect("equal"); ax.invert_yaxis()
            ax.set_xlabel("Xenium x (px)"); ax.set_ylabel("Xenium y (px)")
            rr = ("%.3f" % r["density_r"]) if r.get("density_r") is not None else "-"
            uu = ("%.2f" % r["median_um"]) if r.get("median_um") is not None else "-"
            ax.set_title(f"{r['sample_id']}   [{r['_disp']}]\n"
                         f"selected protocol = {r.get('protocol')}   density_r = {rr}   median_um = {uu}",
                         fontsize=11)
            lg = ax.legend(loc="upper right", fontsize=8, markerscale=6, framealpha=0.9)
            if r.get("note"):
                ax.text(0.01, 0.01, r["note"], transform=ax.transAxes, fontsize=7.5,
                        color="0.35", va="bottom")
            pdf.savefig(fig); plt.close(fig)

    return len(rows) + 1


def _load_slide(sample_id, sdir, xenium_cells, pixel_um):
    """Assemble one slide dict from standard pipeline outputs; None if the decision/nuclei are absent."""
    qcp = os.path.join(sdir, "qc.json")
    if not os.path.exists(qcp):
        return None
    dec = (json.load(open(qcp)) or {}).get("decision") or {}
    chosen = dec.get("chosen")
    if chosen is None:
        return {"sample_id": sample_id, "protocol": None, "xen": None, "he": None,
                "density_r": None, "median_um": None, "note": "no decision in qc.json"}
    he_path = os.path.join(sdir, f"he_nuclei_{chosen}.npy")
    if not os.path.exists(he_path):
        return {"sample_id": sample_id, "protocol": chosen, "xen": None, "he": None,
                "density_r": dec.get("sel_density_r"), "median_um": dec.get("sel_median_um"),
                "note": f"he_nuclei_{chosen}.npy not found (rescued slides warp via run_rescue)"}
    he = np.load(he_path)
    xen = hvxen.load_xenium_nuclei(xenium_cells, pixel_um) if os.path.exists(xenium_cells) else None
    return {"sample_id": sample_id, "protocol": chosen, "xen": xen, "he": he,
            "density_r": dec.get("sel_density_r"), "median_um": dec.get("sel_median_um")}


def main(argv=None):
    ap = argparse.ArgumentParser(description="cohort QC overlay PDF (selected protocol per slide)")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args(argv)

    import csv
    cfg = hvconfig.load_config(a.config)
    th = hvconfig.thresholds(cfg)
    pixel_um = th.get("pixel_um", 0.2125)
    out_root = cfg["output_dir"]

    slides = []
    with open(a.samples) as f:
        for row in csv.DictReader(f):
            sid = row["sample_id"]
            s = _load_slide(sid, os.path.join(out_root, sid), row.get("xenium_cells", ""), pixel_um)
            if s is not None:
                slides.append(s)
    if not slides:
        print("no slides with a qc.json decision found under", out_root); return 1

    out_pdf = a.out or os.path.join(out_root, "cohort_qc.pdf")
    n = render_cohort_pdf(slides, out_pdf, th)
    print(f"wrote {out_pdf}  ({n} pages: 1 summary + {len(slides)} slides)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
