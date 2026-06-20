#!/usr/bin/env python
"""Cohort provenance / acceptance table -> <output_dir>/provenance.csv (+ triage summary).

Reads every <output_dir>/<sample>/qc.json -- the SAME source the cohort report reads -- and
gates each slide on the SAME thresholds the report draws (config.thresholds), so the move table
and the cohort report can never disagree. CPU-only, no re-registration.

  python run_provenance.py --samples samples.csv --config config.json
  python run_provenance.py --config config.json            # batch column -> UNKNOWN

samples.csv may carry an optional `batch` column (scanner-run key) for the per-batch
orientation-consistency check; it is not required.
"""
import argparse
import csv
import os

from hest_valis import config, provenance


def _batch_map(samples_path):
    """sample_id -> batch from samples.csv (tolerant: no required-column validation)."""
    out = {}
    if samples_path and os.path.exists(samples_path):
        with open(samples_path) as f:
            for r in csv.DictReader(f):
                sid = r.get("sample_id")
                if sid:
                    out[sid] = (r.get("batch") or "UNKNOWN").strip() or "UNKNOWN"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="")
    ap.add_argument("--config", default="")
    a = ap.parse_args()
    cfg = config.load_config(a.config)
    th = config.thresholds(cfg)
    csv_path, rows = provenance.write_provenance(
        cfg["output_dir"], th, batch_by_sample=_batch_map(a.samples))
    print(f"wrote {csv_path}")
    print(provenance.summarize(rows, th))


if __name__ == "__main__":
    main()
