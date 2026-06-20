"""Tests for hest_valis.provenance: the acceptance gate + audit table read the SAME qc.json the
report reads and gate on the SAME config thresholds, so the cohort report and the move-provenance
cannot disagree. Pins that invariant plus the per-slide verdict logic.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hest_valis import config, provenance, report

TH = config.thresholds(config.DEFAULTS)


def _metric(density_r, median_um, collapse=0.5):
    """A compute_qc-shaped metric dict with the fields provenance/report read."""
    return {"density_r": density_r, "n_he": 100, "n_xenium": 100,
            "nucleus_coincidence": {"median_um": median_um, "n_matched": 80, "frac_matched": 0.8},
            "negative_control_100um": {"density_collapse": collapse}}


def _write(out, sample, metrics, chosen, sel_r, sel_um, extra_decision=None):
    sd = os.path.join(out, sample)
    os.makedirs(sd, exist_ok=True)
    dec = {"chosen": chosen, "rule": "test", "sel_density_r": sel_r, "sel_median_um": sel_um}
    if extra_decision:
        dec.update(extra_decision)
    with open(os.path.join(sd, "qc.json"), "w") as f:
        json.dump({"sample_id": sample, "metrics": metrics, "decision": dec}, f)
    return sd


# --------------------------------------------------------------------------------------
# the core invariant: ONE set of thresholds for report line AND provenance gate
# --------------------------------------------------------------------------------------
def test_report_line_equals_provenance_gate():
    assert report.QC_DENSITY_R == TH["density_r_accept"]
    assert report.QC_MEDIAN_UM == TH["median_um_accept"]


def test_gate_boundaries():
    ra, ua = TH["density_r_accept"], TH["median_um_accept"]
    assert provenance.gate(ra, ua, TH)               # exactly on the line passes (>=, <=)
    assert not provenance.gate(ra - 0.01, ua, TH)    # just below density_r fails
    assert not provenance.gate(ra, ua + 0.01, TH)    # just above median_um fails
    assert not provenance.gate(None, ua, TH)
    assert not provenance.gate(float("nan"), ua, TH)


# --------------------------------------------------------------------------------------
# per-slide verdicts
# --------------------------------------------------------------------------------------
def test_as_is_good_accepted(tmp_path):
    out = str(tmp_path)
    _write(out, "S1", {"nomicro": _metric(0.85, 3.0)}, "nomicro", 0.85, 3.0)
    rows = provenance.build_rows(out, TH)
    r = rows[0]
    assert r["accepted"] and not r["triggered"]
    assert r["reason"] == "accepted: as-is good"
    assert r["chosen"] == "nomicro"


def test_rescued_accepted_with_big_delta(tmp_path):
    out = str(tmp_path)
    # as-is failed (low r), rescued recovered well -> accepted, large delta (no small-delta flag)
    _write(out, "S2", {"nomicro": _metric(0.15, 12.0), "rescued": _metric(0.88, 2.5)},
           "rescued", 0.88, 2.5, extra_decision={"prerotate_deg": 270})
    rows = provenance.build_rows(out, TH)
    r = rows[0]
    assert r["accepted"] and r["triggered"]
    assert r["chosen"] == "rescued" and r["recovered_k90"] == 3
    assert not r["small_delta_flag"]
    assert "accepted: rescued" in r["reason"]


def test_rescued_small_delta_flagged(tmp_path):
    out = str(tmp_path)
    # rescued passes the gate but barely moved vs as-is -> small-delta flag (eyes)
    _write(out, "S3", {"nomicro": _metric(0.72, 4.5), "rescued": _metric(0.80, 4.0)},
           "rescued", 0.80, 4.0, extra_decision={"prerotate_deg": 0})
    rows = provenance.build_rows(out, TH)
    r = rows[0]
    assert r["accepted"] and r["small_delta_flag"]      # delta_r = 0.08 < rescue_delta_min (0.20)
    assert "SMALL DELTA" in r["reason"]


def test_failed_gate_manual(tmp_path):
    out = str(tmp_path)
    _write(out, "S4", {"nomicro": _metric(0.30, 8.0)}, "nomicro", 0.30, 8.0)
    rows = provenance.build_rows(out, TH)
    r = rows[0]
    assert not r["accepted"]
    assert r["reason"].startswith("manual:")


def test_write_and_summarize(tmp_path):
    out = str(tmp_path)
    _write(out, "S1", {"nomicro": _metric(0.85, 3.0)}, "nomicro", 0.85, 3.0)
    _write(out, "S4", {"nomicro": _metric(0.30, 8.0)}, "nomicro", 0.30, 8.0)
    csv_path, rows = provenance.write_provenance(out, TH)
    assert os.path.exists(csv_path)
    text = provenance.summarize(rows, TH)
    assert "accepted=1" in text and "manual-review=1" in text
    assert "MANUAL REVIEW" in text


def test_batch_orientation_outlier(tmp_path):
    out = str(tmp_path)
    # 3 rescued slides in one batch: two agree on k90=3, one disagrees -> flagged
    for s, deg in [("A", 270), ("B", 270), ("C", 90)]:
        _write(out, s, {"nomicro": _metric(0.1, 15.0), "rescued": _metric(0.85, 3.0)},
               "rescued", 0.85, 3.0, extra_decision={"prerotate_deg": deg})
    batch = {"A": "run1", "B": "run1", "C": "run1"}
    rows = provenance.build_rows(out, TH, batch_by_sample=batch)
    by = {r["slide"]: r for r in rows}
    assert by["C"]["orientation_outlier"] and not by["A"]["orientation_outlier"]
