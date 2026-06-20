"""Per-slide provenance / acceptance table -- the audit + move record for a cohort.

ONE triage source: this reads the SAME `qc.json` the cohort report reads (it does NOT read a
separate rescue-metrics file), and gates on the SAME thresholds the report draws its lines from
(`config.thresholds`). So the cohort report and the move-provenance can never tell you different
things about the same slide.

For each `<output_dir>/<sample>/qc.json`:
  - final outcome  = decision.chosen + decision.sel_density_r / sel_median_um
                     (already incorporates the rescue: run_rescue writes chosen="rescued",
                      metrics["rescued"], decision.sel_* = post-rescue accepted values)
  - acceptance gate: density_r >= density_r_accept AND median_um <= median_um_accept
  - audit columns reconstruct as-is (micro/nomicro) vs rescued/coarse from the metrics dict so you
    can see how far the rescue moved each slide, plus per-batch orientation consistency.

A slide that fails the gate is flagged "manual" -- NOT auto-accepted.

build_rows()/write_provenance() are the API; run_provenance.py is the cohort driver.
"""
import csv
import glob
import json
import os
import collections

NAN = float("nan")

ASIS = ("micro", "nomicro")          # the as-is VALIS registration variants
RESCUE = ("rescued", "coarse")       # orientation-rescue variants (rescued preferred over coarse)


def gate(density_r, median_um, th):
    """Acceptance gate: both metrics must pass, using the single threshold set `th`."""
    return (density_r is not None and density_r == density_r and density_r >= th["density_r_accept"]
            and median_um is not None and median_um == median_um and median_um <= th["median_um_accept"])


def _variant_rum(metrics, variant):
    """(density_r, median_um) for a metrics variant, or (None, None) if absent / no_nuclei."""
    m = (metrics or {}).get(variant)
    if not m or m.get("status") == "no_nuclei":
        return None, None
    return m.get("density_r"), (m.get("nucleus_coincidence") or {}).get("median_um")


def _asis_best(metrics):
    """Best as-is (micro/nomicro) variant by median_um -> (r, um, variant) or (None, None, None).
    Lower median_um wins, matching select.py's primary metric."""
    cands = [(um, r, v) for v in ASIS for (r, um) in [_variant_rum(metrics, v)] if um is not None]
    if not cands:
        return None, None, None
    um, r, v = min(cands)
    return r, um, v


def _rescue_variant(metrics):
    """The orientation-rescue variant present (rescued > coarse) -> (r, um, variant) or Nones."""
    for v in RESCUE:
        r, um = _variant_rum(metrics, v)
        if um is not None:
            return r, um, v
    return None, None, None


def _reason(accepted, rescued_chosen, triggered, small_delta, chosen,
            pre_r, pre_um, final_r, final_um, has_rescue, orient):
    if accepted and not rescued_chosen and not triggered:
        return "accepted: as-is good"
    if accepted and rescued_chosen:
        return f"accepted: {chosen} ({orient})" + (" [SMALL DELTA -> eyes]" if small_delta else "")
    if accepted and not rescued_chosen and triggered:
        return "accepted: as-is (rescue not better)"
    if not triggered and not accepted:
        return f"manual: mediocre as-is r={_n(pre_r)} um={_n(pre_um)}"
    if triggered and not has_rescue and not rescued_chosen:
        return "manual: triggered but no rescue variant (not run / failed)"
    return f"manual: failed gate r={_n(final_r)} um={_n(final_um)}"


def _n(v, nd=2):
    return "?" if v is None or v != v else f"{v:.{nd}f}"


def build_rows(output_dir, th, batch_by_sample=None):
    """One row per <output_dir>/<sample>/qc.json. Pure (no file writes). Returns list[dict]."""
    batch_by_sample = batch_by_sample or {}
    rows = []
    for qp in sorted(glob.glob(os.path.join(output_dir, "*", "qc.json"))):
        sample = os.path.basename(os.path.dirname(qp))
        try:
            with open(qp) as f:
                qc = json.load(f)
        except Exception:
            continue
        metrics = qc.get("metrics", {}) or {}
        dec = qc.get("decision", {}) or {}
        chosen = dec.get("chosen")
        final_r, final_um = dec.get("sel_density_r"), dec.get("sel_median_um")

        pre_r, pre_um, _pre_v = _asis_best(metrics)
        post_r, post_um, _post_v = _rescue_variant(metrics)
        has_rescue = post_um is not None
        rescued_chosen = chosen in ("coarse", "rescued")
        delta_r = (post_r - pre_r) if (post_r is not None and pre_r is not None) else NAN

        triggered = (rescued_chosen
                     or pre_r is None
                     or (pre_r == pre_r and pre_r < th["rescue_trigger_r"]))
        accepted = gate(final_r, final_um, th)
        small_delta = bool(rescued_chosen and delta_r == delta_r and delta_r < th["rescue_delta_min"])

        deg = dec.get("prerotate_deg")
        k90 = (int(deg) // 90) if isinstance(deg, (int, float)) else ""
        flip = dec.get("prerotate_flip", "")          # set once the flip-capable rescue lands
        orient = f"k90={k90}" + (f" flip={flip}" if flip != "" else "")

        reason = _reason(accepted, rescued_chosen, triggered, small_delta, chosen,
                         pre_r, pre_um, final_r, final_um, has_rescue, orient)

        nuc = f"he_nuclei_{chosen}.npy" if chosen else ""
        wsi = sorted(glob.glob(os.path.join(os.path.dirname(qp), "registered", "*.ome.tif*")))
        rows.append({
            "slide": sample, "batch": batch_by_sample.get(sample, "UNKNOWN"), "chosen": chosen,
            "pre_r": _r(pre_r), "post_r": _r(post_r),
            "pre_median_um": _r(pre_um, 2), "post_median_um": _r(post_um, 2),
            "recovered_flip": flip, "recovered_k90": k90,
            "delta_r": _r(delta_r), "triggered": triggered, "accepted": accepted,
            "reason": reason, "orientation_outlier": "",
            "final_r": _r(final_r), "final_median_um": _r(final_um, 2),
            "small_delta_flag": small_delta,
            "chosen_nuclei": os.path.join(sample, nuc) if nuc else "",
            "registered_wsi": os.path.relpath(wsi[0], output_dir) if wsi else "",
        })

    _flag_orientation_outliers(rows)
    return rows


def _r(v, nd=3):
    return round(v, nd) if isinstance(v, (int, float)) and v == v else ""


def _flag_orientation_outliers(rows):
    """Within a batch, flag rescued slides whose recovered orientation disagrees with the
    batch majority (only when the batch has >= 3 rescued slides)."""
    by_batch = collections.defaultdict(list)
    for r in rows:
        if r["chosen"] in ("coarse", "rescued") and r["recovered_k90"] != "":
            by_batch[r["batch"]].append(r)
    for grp in by_batch.values():
        if len(grp) < 3:
            continue
        votes = collections.Counter((g["recovered_flip"], g["recovered_k90"]) for g in grp)
        maj, _ = votes.most_common(1)[0]
        for g in grp:
            if (g["recovered_flip"], g["recovered_k90"]) != maj:
                g["orientation_outlier"] = (f"batch_majority={maj} "
                                            f"this=({g['recovered_flip']},{g['recovered_k90']})")


COLUMNS = ["slide", "batch", "chosen", "pre_r", "post_r", "pre_median_um", "post_median_um",
           "recovered_flip", "recovered_k90", "delta_r", "triggered", "accepted", "reason",
           "orientation_outlier", "final_r", "final_median_um", "small_delta_flag",
           "chosen_nuclei", "registered_wsi"]


def write_provenance(output_dir, th, batch_by_sample=None, out_csv=None):
    """Write <output_dir>/provenance.csv. Returns (csv_path, rows)."""
    rows = build_rows(output_dir, th, batch_by_sample)
    out_csv = out_csv or os.path.join(output_dir, "provenance.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    return out_csv, rows


def summarize(rows, th):
    """Human triage summary string (accepted vs manual + flags). Mirrors the report's verdicts."""
    acc = sum(r["accepted"] for r in rows)
    out = [f"{len(rows)} slides | accepted={acc} | manual-review={len(rows) - acc}",
           f"thresholds: gate r>={th['density_r_accept']} & um<={th['median_um_accept']}, "
           f"trigger r<{th['rescue_trigger_r']}, delta_min={th['rescue_delta_min']}"]
    manual = [r for r in rows if not r["accepted"]]
    if manual:
        out.append("MANUAL REVIEW:")
        out += [f"  {r['slide']:8s} [{r['batch']}] {r['reason']}" for r in manual]
    outl = [r for r in rows if r["orientation_outlier"]]
    if outl:
        out.append(f"ORIENTATION OUTLIERS ({len(outl)}):")
        out += [f"  {r['slide']:8s} {r['orientation_outlier']}" for r in outl]
    sd = [r["slide"] for r in rows if r["small_delta_flag"]]
    if sd:
        out.append(f"SMALL-DELTA (accepted but small density-r jump -> eyes): {sd}")
    return "\n".join(out)
