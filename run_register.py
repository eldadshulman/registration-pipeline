#!/usr/bin/env python
"""Step 2/3 (valis env): register H&E onto DAPI and warp the H&E nuclei into the Xenium frame
for BOTH protocols from a single registration.

  register()        -> warp nuclei (non-rigid)  = no-micro variant
  register_micro()  -> warp nuclei (non-rigid)  = micro variant   (skipped on failure)

Writes <out>/<sample>/he_nuclei_{nomicro,micro}.npy (Xenium-frame pixels).
Run in the patched valis_hest_env (see env/setup.md).
  python run_register.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import os
import numpy as np
from hest_valis import config, registration


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    s = config.get_sample(config.load_samples(a.samples), a.sample)
    out = os.path.join(cfg["output_dir"], a.sample)
    nuc = np.load(os.path.join(out, "he_nuclei.npy")).astype(float)

    reg = registration.register_slide(s["he_path"], s["dapi_path"], out, micro=False)
    np.save(os.path.join(out, "he_nuclei_nomicro.npy"),
            registration.warp_points(reg, nuc, non_rigid=True))
    print(f"[{a.sample}] no-micro nuclei warped", flush=True)

    try:
        reg.register_micro(
            max_non_rigid_registration_dim_px=registration.MICRO_MAX_DIM_PX,
            align_to_reference=True,
            brightfield_processing_cls=registration.preprocessing.HEDeconvolution,
            reference_img_f=s["dapi_path"],
        )
        np.save(os.path.join(out, "he_nuclei_micro.npy"),
                registration.warp_points(reg, nuc, non_rigid=True))
        print(f"[{a.sample}] micro nuclei warped", flush=True)
    except Exception as e:
        print(f"[{a.sample}] MICRO FAILED ({type(e).__name__}: {str(e)[:120]}) -> no-micro only", flush=True)
    registration.shutdown()


if __name__ == "__main__":
    main()
