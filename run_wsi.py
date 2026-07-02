#!/usr/bin/env python
"""Warp the chosen-protocol H&E image into the Xenium frame -> a registered OME-TIFF.

Reads the chosen protocol from <out>/<sample>/qc.json (or pass --micro 0/1 to override).
Writes, into <out>/<sample>/registered/:
  - aligned_fullres_HE.ome.tiff   the Xenium-frame H&E image (the deliverable), and
  - he_nuclei_registered.npy      the H&E nuclei warped by the SAME registration -- a consistent
                                  (image, nuclei) pair for the downstream cell / mask step.
Run in the PATCHED valis_hest_env (env/setup.md) or valis deadlocks at the COLLECTING step.
  python run_wsi.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import glob
import json
import os
import numpy as np
from hest_valis import config, registration


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--micro", type=int, default=-1, help="force 1/0; default = read qc.json")
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    s = config.get_sample(config.load_samples(a.samples), a.sample)
    out = os.path.join(cfg["output_dir"], a.sample)
    reg_out = os.path.join(out, "registered")
    nuc_out = os.path.join(reg_out, "he_nuclei_registered.npy")   # step-1 input: nuclei in Xenium frame
    img_ok = any(registration.ometiff_pages(f) > 0
                 for f in glob.glob(os.path.join(reg_out, "*.ome.tif*")))
    # a present-but-truncated OME-TIFF (0 pages) does NOT count as done -> re-warp
    if img_ok and os.path.exists(nuc_out):
        print(f"[{a.sample}] WSI + registered nuclei exist, skip", flush=True)
        return

    micro = a.micro
    if micro < 0:
        with open(os.path.join(out, "qc.json")) as _f:
            chosen = json.load(_f)["decision"]["chosen"]
        if chosen in ("coarse", "rescued"):
            print(f"[{a.sample}] chosen={chosen}: this slide was rescued from a failed "
                  f"registration, so warping the ORIGINAL H&E here would be wrong. Generate its "
                  f"WSI with `run_rescue.py ... --warp-image` (uses the pre-rotated H&E). "
                  f"Skipping run_wsi. See README 'Self-healing'.", flush=True)
            return
        micro = 1 if chosen == "micro" else 0

    reg = registration.register_slide(s["he_path"], s["dapi_path"], out, micro=bool(micro))
    if not img_ok:
        registration.warp_image(reg, reg_out, level=0, non_rigid=True)
    # step-1 input: warp the H&E nuclei with THIS registration so the (image, nuclei) pair the
    # downstream cell/mask step consumes is consistent. Cheap point-warp, saved beside the OME-TIFF.
    os.makedirs(reg_out, exist_ok=True)
    nuc_src = os.path.join(out, "he_nuclei.npy")
    if os.path.exists(nuc_src):
        wn = registration.warp_points(reg, np.load(nuc_src).astype(float), non_rigid=True)
        np.save(nuc_out, wn)
        print(f"[{a.sample}] registered nuclei -> {nuc_out} ({len(wn)} pts)", flush=True)
    else:
        print(f"[{a.sample}] WARN: {nuc_src} missing (run_segment first); no registered nuclei", flush=True)
    print(f"[{a.sample}] WSI written (micro={micro}) -> {reg_out}", flush=True)
    registration.shutdown()


if __name__ == "__main__":
    main()
