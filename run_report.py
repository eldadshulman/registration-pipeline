#!/usr/bin/env python
"""Render the alignment QC report (no GPU, no re-registration; plotting only).

Per slide  -> <out>/<sample>/report.pdf (+ report.png for the review notebook)
Cohort     -> <out>/cohort_report.pdf (the triage page)

It reads what run_qc / run_wsi / run_annotate already produced (qc.json, the chosen
he_nuclei_<variant>.npy, the Xenium centroids, and -- if present -- the warped raster and
cell_labels.parquet). It never re-registers and only reloads a WSI when the warped raster
already exists on disk.

  python run_report.py --samples samples.csv --config config.json --sample <id>
  python run_report.py --samples samples.csv --config config.json --cohort
"""
import argparse
import os

from hest_valis import config, report, xenium


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", help="render one slide's report.pdf")
    ap.add_argument("--cohort", action="store_true", help="render the cohort triage page")
    ap.add_argument("--no-png", action="store_true", help="skip the per-slide report.png")
    a = ap.parse_args()
    if not a.sample and not a.cohort:
        ap.error("pass --sample <id> and/or --cohort")
    cfg = config.load_config(a.config)

    if a.sample:
        s = config.get_sample(config.load_samples(a.samples), a.sample)
        out = os.path.join(cfg["output_dir"], a.sample)
        um = cfg["pixel_um"]
        xen = xenium.load_xenium_nuclei(s["xenium_cells"], um, in_um=cfg["centroids_in_um"])
        pdf = report.render_sample(out, xen, um, sample_id=a.sample,
                                   dapi_path=s.get("dapi_path"), save_png=not a.no_png)
        print(f"[{a.sample}] report -> {pdf}", flush=True)

    if a.cohort:
        pdf = report.render_cohort(cfg["output_dir"])
        print(f"[cohort] report -> {pdf}", flush=True)


if __name__ == "__main__":
    main()
