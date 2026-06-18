#!/usr/bin/env python
"""Aggregate per-sample qc.json into one decision table + a WSI manifest.

  python run_select.py --samples samples.csv --config config.json
Writes <output_dir>/per_slide_decision.csv and <output_dir>/wsi_manifest.csv
"""
import argparse
import csv
import json
import os
from hest_valis import config


def g(m, *keys, default=None):
    for k in keys:
        if m is None:
            return default
        m = m.get(k) if isinstance(m, dict) else None
    return m if m is not None else default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    samples = config.load_samples(a.samples)
    od = cfg["output_dir"]

    rows, manifest = [], []
    for s in samples:
        sid = s["sample_id"]
        qp = os.path.join(od, sid, "qc.json")
        if not os.path.exists(qp):
            rows.append({"sample_id": sid, "status": "no_qc"})
            continue
        q = json.load(open(qp)); d = q["decision"]; M = q["metrics"]
        rows.append({
            "sample_id": sid, "status": "ok",
            "micro_med": g(M, "micro", "nucleus_coincidence", "median_um"),
            "micro_r": g(M, "micro", "density_r"),
            "nomicro_med": g(M, "nomicro", "nucleus_coincidence", "median_um"),
            "nomicro_r": g(M, "nomicro", "density_r"),
            "chosen": d["chosen"], "rule": d["rule"],
            "sel_median_um": d["sel_median_um"], "sel_density_r": d["sel_density_r"]})
        manifest.append({"sample_id": sid, "micro": 1 if d["chosen"] == "micro" else 0})

    cols = ["sample_id", "status", "micro_med", "micro_r", "nomicro_med", "nomicro_r",
            "chosen", "rule", "sel_median_um", "sel_density_r"]
    with open(os.path.join(od, "per_slide_decision.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with open(os.path.join(od, "wsi_manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "micro"]); w.writeheader(); w.writerows(manifest)

    ok = [r for r in rows if r["status"] == "ok"]
    nmic = sum(1 for r in ok if r["chosen"] == "micro")
    print(f"decided {len(ok)}/{len(rows)} | {nmic} micro / {len(ok)-nmic} no-micro "
          f"-> {od}/per_slide_decision.csv, wsi_manifest.csv")


if __name__ == "__main__":
    main()
