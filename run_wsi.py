#!/usr/bin/env python
"""Warp the chosen-protocol H&E image into the Xenium frame -> a registered OME-TIFF.

Reads the chosen protocol from <out>/<sample>/qc.json (or pass --micro 0/1 to override).
Writes <out>/<sample>/registered/aligned_fullres_HE.ome.tiff (Xenium-frame H&E).
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
    if any(registration.ometiff_pages(f) > 0
           for f in glob.glob(os.path.join(reg_out, "*.ome.tif*"))):
        print(f"[{a.sample}] WSI exists, skip", flush=True)
        return
    # a present-but-truncated OME-TIFF (0 pages) does NOT count as done -> fall through and re-warp

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
    registration.warp_image(reg, reg_out, level=0, non_rigid=True)
    print(f"[{a.sample}] WSI written (micro={micro}) -> {reg_out}", flush=True)
    registration.shutdown()


if __name__ == "__main__":
    main()
