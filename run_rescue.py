#!/usr/bin/env python
"""Rescue a coarse-flagged slide into a FULL registration (valis env).

When run_qc.py selected 'coarse' (a gross mis-orientation VALIS could not recover), this finishes
the job automatically: pick the cardinal rotation, losslessly pre-rotate the H&E, re-register with
VALIS + micro, re-QC, and -- if it beats the coarse density-r -- adopt it so the slide gets a
proper micro / no-micro registration AND a warpable image.

  python run_rescue.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import json
import os
import numpy as np
from hest_valis import config, registration, coarse_align, concordance, xenium


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--warp-image", action="store_true", help="also warp the WSI if rescue wins")
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    s = config.get_sample(config.load_samples(a.samples), a.sample)
    out = os.path.join(cfg["output_dir"], a.sample)
    um = cfg["pixel_um"]
    with open(os.path.join(out, "qc.json")) as _f:
        qc = json.load(_f)
    if qc["decision"]["chosen"] != "coarse":
        print(f"[{a.sample}] not coarse-flagged; nothing to rescue"); return
    coarse_r = qc["decision"]["sel_density_r"]

    xen = xenium.load_xenium_nuclei(s["xenium_cells"], um, in_um=cfg["centroids_in_um"])
    he_src = np.load(os.path.join(out, "he_nuclei.npy")).astype(float)
    scale = xenium.he_pixel_um(s["he_path"], cfg.get("he_pixel_um") or xenium.HE_FALLBACK_MPP) / um

    import pyvips
    W, H = (lambda im: (im.width, im.height))(pyvips.Image.new_from_file(s["he_path"], access="sequential"))
    rot, _ = coarse_align.cardinal_rotation(he_src, xen, scale, (W, H), dapi_um=um)
    pre, _ = registration.prerotate_he(s["he_path"], rot, os.path.join(out, "he_prerot.tiff"))
    print(f"[{a.sample}] pre-rotated {rot} deg; re-registering with micro", flush=True)

    reg = registration.register_slide(pre, s["dapi_path"], out, micro=True)
    he_pre = coarse_align.rotate_points(he_src, rot, W, H).astype(float)
    aligned = registration.warp_points(reg, he_pre, non_rigid=True)
    m = concordance.compute_qc(aligned, xen, um)
    print(f"[{a.sample}] rescue density_r = {m['density_r']} (coarse was {coarse_r})", flush=True)

    if m["density_r"] > coarse_r:
        np.save(os.path.join(out, "he_nuclei_rescued.npy"), aligned)
        qc["metrics"]["rescued"] = m
        qc["decision"] = {"chosen": "rescued", "rule": "prerotate_reregister",
                          "sel_median_um": round(m["nucleus_coincidence"]["median_um"], 3),
                          "sel_density_r": round(m["density_r"], 3),
                          "prerotate_deg": rot}
        with open(os.path.join(out, "qc.json"), "w") as _f:
            json.dump(qc, _f, indent=2)
        if a.warp_image:
            registration.warp_image(reg, os.path.join(out, "registered"), level=0, non_rigid=True)
            print(f"[{a.sample}] rescued WSI written", flush=True)
        print(f"[{a.sample}] RESCUED -> adopted (density_r {m['density_r']})", flush=True)
    else:
        print(f"[{a.sample}] rescue did not beat coarse; keeping coarse", flush=True)
    registration.shutdown()


if __name__ == "__main__":
    main()
