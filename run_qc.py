#!/usr/bin/env python
"""Step 3/3 (QC env): compute concordance QC for each warped variant and pick the protocol.

Self-healing: if the best variant still has a negative / near-zero density-r (the signature of a
gross mis-orientation that VALIS could not recover), run the coarse rotation/flip search as a
fallback and, if it does better, select it. This would have auto-rescued a 270-deg-rotated slide.

Writes <out>/<sample>/qc.json with all variants' metrics, the chosen protocol, and the rule.
  python run_qc.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import json
import os
import numpy as np
from hest_valis import config, concordance, select, xenium, coarse_align

COARSE_TRIGGER = 0.10  # if the selected density-r is below this, try the coarse fallback


def _sanitize(obj):
    """Recursively coerce NaN/inf to None so the dict is JSON-serialisable."""
    import math
    if isinstance(obj, float):
        return None if not math.isfinite(obj) else obj
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    return obj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--no-occupancy", action="store_true")
    ap.add_argument("--no-coarse-fallback", action="store_true")
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
    qc = lambda he: concordance.compute_qc(he, xen, um, tissue_mask=mask,
                                           mask_pixel_um=(um * mstep) if mask is not None else None)

    metrics = {}
    for variant in ("micro", "nomicro"):
        p = os.path.join(out, f"he_nuclei_{variant}.npy")
        metrics[variant] = qc(np.load(p)) if os.path.exists(p) else None

    decision = select.choose(metrics["micro"], metrics["nomicro"])

    # self-healing fallback on a failed registration (negative / near-zero density-r)
    sel_r = decision["sel_density_r"]
    if (not a.no_coarse_fallback and sel_r is not None and sel_r < COARSE_TRIGGER
            and os.path.exists(os.path.join(out, "he_nuclei.npy"))):
        scale = xenium.he_pixel_um(s["he_path"], cfg.get("he_pixel_um") or xenium.HE_FALLBACK_MPP) / um
        he_src = np.load(os.path.join(out, "he_nuclei.npy")).astype(float)
        aligned, params, _ = coarse_align.coarse_align(he_src, xen, scale, dapi_um=um)
        mc = qc(aligned)
        metrics["coarse"] = {**mc, "coarse_params": params}
        print(f"[{a.sample}] coarse fallback: density_r {mc['density_r']} (was {sel_r}) "
              f"angle={params['angle']} flip={params['flip']}", flush=True)
        if mc["density_r"] > sel_r:
            np.save(os.path.join(out, "he_nuclei_coarse.npy"), aligned)
            decision = {"chosen": "coarse", "rule": "coarse_rescue_negative_density_r",
                        "sel_median_um": round(mc["nucleus_coincidence"]["median_um"], 3),
                        "sel_density_r": round(mc["density_r"], 3)}

    payload = _sanitize({"sample_id": a.sample, "metrics": metrics, "decision": decision})
    with open(os.path.join(out, "qc.json"), "w") as _f:
        json.dump(payload, _f, indent=2)
    d = decision
    print(f"[{a.sample}] chosen={d['chosen']} ({d['rule']}) "
          f"median={d['sel_median_um']}um density_r={d['sel_density_r']}", flush=True)


if __name__ == "__main__":
    main()
