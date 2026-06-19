#!/usr/bin/env python
"""Per-cell annotation transfer: tag each Xenium cell with its registered-H&E region.

Uses the chosen-protocol warped nuclei (from run_register / the coarse fallback) to build a
tumor/stroma/background region map and assign every Xenium cell its region.
Writes <out>/<sample>/cell_labels.parquet (cell_id, x_um, y_um, he_region) + a region overlay.
Run in an env with numpy/scipy/pandas/scikit-learn/matplotlib.
  python run_annotate.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hest_valis import config, annotate, xenium

COL = {"tumor": "#c0392b", "stroma": "#2c7fb8", "background": "#dddddd"}


def chosen_nuclei_path(out):
    chosen = json.load(open(os.path.join(out, "qc.json")))["decision"]["chosen"]
    name = {"micro": "he_nuclei_micro.npy", "nomicro": "he_nuclei_nomicro.npy",
            "coarse": "he_nuclei_coarse.npy"}[chosen]
    return os.path.join(out, name), chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    s = config.get_sample(config.load_samples(a.samples), a.sample)
    out = os.path.join(cfg["output_dir"], a.sample)
    um = cfg["pixel_um"]

    nuc_path, chosen = chosen_nuclei_path(out)
    he = np.load(nuc_path)
    cells = xenium.load_xenium_cells(s["xenium_cells"], um, in_um=cfg["centroids_in_um"])
    labels = annotate.assign_cells(he, cells, um)
    labels.to_parquet(os.path.join(out, "cell_labels.parquet"))

    vc = labels.he_region.value_counts()
    nbx = int(max(he[:, 0].max(), cells.x_px.max()) / (50.0 / um)) + 2
    nby = int(max(he[:, 1].max(), cells.y_px.max()) / (50.0 / um)) + 2
    region, _ = annotate.region_map(he, um, nbx, nby)
    rgb = np.zeros((nby, nbx, 3))
    for k, c in COL.items():
        rgb[region == k] = tuple(int(c[i:i + 2], 16) / 255 for i in (1, 3, 5))
    fig, ax = plt.subplots(figsize=(4, 5))
    ax.imshow(rgb, origin="upper"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{a.sample} ({chosen})\n{vc.get('tumor',0):,} tumor / {vc.get('stroma',0):,} stroma cells", fontsize=10)
    fig.savefig(os.path.join(out, "region_overlay.png"), dpi=110, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"[{a.sample}] {len(labels):,} cells labeled ({chosen}) -> "
          f"tumor {vc.get('tumor',0):,} stroma {vc.get('stroma',0):,} background {vc.get('background',0):,}", flush=True)


if __name__ == "__main__":
    main()
