#!/usr/bin/env python
"""Step 1/3 (StarDist env): detect H&E nuclei -> <out>/<sample>/he_nuclei.npy

Run in an env that has stardist + tensorflow. GPU optional but much faster.
  python run_segment.py --samples samples.csv --config config.json --sample <id>
"""
import argparse
import os
import numpy as np
from hest_valis import config, segment


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    s = config.get_sample(config.load_samples(a.samples), a.sample)
    out = os.path.join(cfg["output_dir"], a.sample)
    os.makedirs(out, exist_ok=True)
    npy = os.path.join(out, "he_nuclei.npy")
    if os.path.exists(npy):
        print(f"[{a.sample}] nuclei exist, skip", flush=True)
        return
    nuc = segment.segment_he(s["he_path"])
    np.save(npy, nuc)
    print(f"[{a.sample}] {len(nuc):,} H&E nuclei -> {npy}", flush=True)


if __name__ == "__main__":
    main()
