#!/usr/bin/env python
"""Step 3/3 (QC env): compute concordance QC for each warped variant and pick the protocol.

Writes <out>/<sample>/qc.json with both variants' metrics, the chosen protocol, and the rule.
Run in an env with numpy/scipy/pandas/tifffile/zarr.
  python run_qc.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import json
import os
import numpy as np
from hest_valis import config, concordance, select, xenium


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--no-occupancy", action="store_true", help="skip the tissue-occupancy check")
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    s = config.get_sample(config.load_samples(a.samples), a.sample)
    out = os.path.join(cfg["output_dir"], a.sample)
    um = cfg["pixel_um"]

    xen = xenium.load_xenium_nuclei(s["xenium_cells"], um, in_um=cfg["centroids_in_um"])
    mask, mstep = (None, None)
    if not a.no_occupancy:
        try:
            mask, mstep = xenium.tissue_mask_from_dapi(s["dapi_path"])
        except Exception as e:
            print(f"[{a.sample}] occupancy mask skipped ({e})", flush=True)

    metrics = {}
    for variant in ("micro", "nomicro"):
        p = os.path.join(out, f"he_nuclei_{variant}.npy")
        if not os.path.exists(p):
            metrics[variant] = None
            continue
        he = np.load(p)
        metrics[variant] = concordance.compute_qc(
            he, xen, um, tissue_mask=mask,
            mask_pixel_um=(um * mstep) if mask is not None else None)

    decision = select.choose(metrics["micro"], metrics["nomicro"])
    result = {"sample_id": a.sample, "metrics": metrics, "decision": decision}
    json.dump(result, open(os.path.join(out, "qc.json"), "w"), indent=2)
    d = decision
    print(f"[{a.sample}] chosen={d['chosen']} ({d['rule']}) "
          f"median={d['sel_median_um']}um density_r={d['sel_density_r']}", flush=True)


if __name__ == "__main__":
    main()
