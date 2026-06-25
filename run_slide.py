#!/usr/bin/env python
"""Thin single-slide driver:  register -> QC -> selection -> (optional) WSI warp.

A LOCAL, one-slide convenience for demos/debugging. The documented PRODUCTION path stays the
SLURM arrays (`slurm/*.sbatch` + `run_select.py`). This wrapper does NOT re-implement any QC
thresholds, protocol-selection, or rescue logic: it invokes the existing entry points and
consumes their authoritative machine-readable decision (`qc.json`), and imports the real
acceptance gate (`provenance.gate`, `config.thresholds`) and eligibility sets
(`provenance.ASIS`/`RESCUE`) for terminology.

Stages and their artifacts (under `<output_dir>/<sample>/`):
  segment   run_segment.py  -> he_nuclei.npy                 (auto-run only if missing)
  register  run_register.py -> he_nuclei_{nomicro,micro}.npy
  qc        run_qc.py       -> qc.json   (metrics + decision{chosen,rule,sel_*}; coarse self-heal)
  [select]  consumed from qc.json's `decision` -- no separate call for one slide
            (run_select.py is the COHORT aggregator; here the per-slide decision is already in qc.json)
  wsi       run_wsi.py      -> registered/aligned_fullres_HE.ome.tiff   (only with --warp AND eligible)

Per-stage Python envs differ (StarDist / patched valis / QC). Set $STARDIST_PY, $VALIS_PY,
$QC_PY; each falls back to the current interpreter if unset.

Exit codes:
  0  completed through selection (and the warp, if --warp + eligible)
  1  a stage command failed (implementation error -> stop)
  2  REVIEW_REQUIRED  (selected but below the accept gate, OR an orientation rescue -> run_rescue)
  3  QC_FAILED        (no usable registration decision)
"""
import argparse
import glob
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hest_valis import config as hvconfig
from hest_valis import provenance as hvprov          # gate() + ASIS/RESCUE sets (authoritative)

HERE = os.path.dirname(os.path.abspath(__file__))
STAGES = ("segment", "register", "qc", "wsi")


def _py(var):
    """Per-stage interpreter from env (StarDist/valis/QC), else the current one."""
    return os.environ.get(var) or sys.executable


def run_stage(cmd):
    """Invoke one pipeline entry point as a subprocess; return its exit code. Patched in tests."""
    print("    $ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd).returncode


def _valid_ometiff(reg_dir):
    """True if a non-truncated registered OME-TIFF already exists (pages > 0). tifffile is imported
    only when there is a file to check, so the wrapper/tests stay importable without it."""
    files = glob.glob(os.path.join(reg_dir, "*.ome.tif*"))
    if not files:
        return False
    import tifffile
    for p in files:
        try:
            with tifffile.TiffFile(p) as tf:
                if len(tf.pages) > 0:
                    return True
        except Exception:
            pass
    return False


def classify(qc, th):
    """Map the authoritative qc.json decision -> (status, protocol, eligible_for_warp, reason, code).

    No thresholds/rules are redefined here: `chosen` and `sel_*` come from qc.json (run_qc ->
    select.choose), the gate is `provenance.gate`, and eligibility uses `provenance.ASIS`."""
    dec = (qc or {}).get("decision") or {}
    chosen = dec.get("chosen")
    r, u = dec.get("sel_density_r"), dec.get("sel_median_um")
    if chosen is None:
        return ("QC_FAILED", None, False,
                "no protocol chosen (no nuclei / registration produced no decision)", 3)
    passed = hvprov.gate(r, u, th)
    if chosen not in hvprov.ASIS:                       # coarse/rescued: a valid outcome, but NOT a run_wsi target
        return ("REVIEW_REQUIRED", chosen, False,
                f"orientation-rescue ('{chosen}'): warp via `run_rescue.py --warp-image`, not run_wsi", 2)
    if not passed:                                      # selected as-is but under the accept gate
        return ("REVIEW_REQUIRED", chosen, False,
                f"selected '{chosen}' but below accept gate "
                f"(density_r>={th['density_r_accept']} AND median_um<={th['median_um_accept']})", 2)
    return ("ELIGIBLE_FOR_WARP", chosen, True, f"'{chosen}' passes the accept gate", 0)


def _plan(out, samples, cfgpath, sample):
    """[(stage, output_exists, command)] for the auto/register/qc stages."""
    common = ["--samples", samples, "--config", cfgpath, "--sample", sample]
    return [
        ("segment",  os.path.exists(os.path.join(out, "he_nuclei.npy")),
         [_py("STARDIST_PY"), os.path.join(HERE, "run_segment.py"), *common]),
        ("register", os.path.exists(os.path.join(out, "he_nuclei_nomicro.npy")),
         [_py("VALIS_PY"), os.path.join(HERE, "run_register.py"), *common]),
        ("qc",       os.path.exists(os.path.join(out, "qc.json")),
         [_py("QC_PY"), os.path.join(HERE, "run_qc.py"), *common]),
    ]


def main(argv=None):
    ap = argparse.ArgumentParser(description="single-slide register->QC->select->(optional)warp")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--warp", action="store_true",
                    help="actually run the expensive full-res WSI warp (only if eligible)")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and intended actions; run nothing")
    ap.add_argument("--resume", action="store_true", help="skip stages whose (valid) outputs already exist")
    ap.add_argument("--force-stage", action="append", default=[], choices=STAGES,
                    help="re-run this stage even if its output exists (repeatable)")
    a = ap.parse_args(argv)

    cfg = hvconfig.load_config(a.config)
    th = hvconfig.thresholds(cfg)
    out = os.path.join(cfg["output_dir"], a.sample)
    qcp = os.path.join(out, "qc.json")
    regdir = os.path.join(out, "registered")
    common = ["--samples", a.samples, "--config", a.config, "--sample", a.sample]
    plan = _plan(out, a.samples, a.config, a.sample)

    def will_run(stage, done):
        if stage in a.force_stage:
            return True
        return not (a.resume and done)

    print(f"=== single-slide: {a.sample}  ->  {out} ===")

    # ---- dry-run: print the plan + intended warp decision, execute nothing ----
    if a.dry_run:
        print("DRY-RUN (no commands executed):")
        for st, done, cmd in plan:
            print(f"  {st:9} {'RUN ' if will_run(st, done) else 'skip'}  ({os.path.basename(cmd[1])})")
        if os.path.exists(qcp):
            status, prot, elig, reason, _ = classify(json.load(open(qcp)), th)
            print(f"  selection  status={status}  protocol={prot}  reason={reason}")
            warp_intent = ("WOULD WARP" if (elig and a.warp)
                           else "would NOT warp (no --warp)" if elig
                           else f"would NOT warp (not eligible: {status})")
            print(f"  wsi        {warp_intent}")
        else:
            print("  selection/wsi  (qc.json not present yet -> would run register+qc first)")
        return 0

    os.makedirs(out, exist_ok=True)
    stage_status = {}

    # ---- segment -> register -> qc (stop on the first failure) ----
    for st, done, cmd in plan:
        if not will_run(st, done):
            stage_status[st] = "skipped (resume)"
            print(f"[{st}] skip -- output exists (--resume)")
            continue
        print(f"[{st}] run")
        rc = run_stage(cmd)
        if rc != 0:
            stage_status[st] = f"FAILED rc={rc}"
            print(f"[{st}] FAILED (rc={rc}) -> stopping (later stages not run)")
            _summary(a.sample, stage_status, None, "STAGE_FAILED", None, False,
                     f"{st} exited rc={rc}", "not run", qcp, "", th)
            return 1
        stage_status[st] = "ran"

    # ---- selection: consume the authoritative qc.json decision ----
    if not os.path.exists(qcp):
        print("qc.json missing after the QC stage -> cannot decide"); return 1
    qc = json.load(open(qcp))
    status, protocol, eligible, reason, code = classify(qc, th)
    dec = qc.get("decision", {})

    # ---- WSI warp: only with --warp AND eligible; never silently overwrite ----
    force_wsi = "wsi" in a.force_stage
    wsi_out = ""
    if not eligible:
        wsi_state = f"skipped (not eligible: {status})"
    elif not a.warp:
        wsi_state = "skipped (no --warp; eligible -> would warp)"
    elif a.resume and not force_wsi and _valid_ometiff(regdir):
        wsi_state = "skipped (--resume; valid registered image already present)"
    else:
        print("[wsi] run")
        micro = "1" if protocol == "micro" else "0"
        rc = run_stage([_py("VALIS_PY"), os.path.join(HERE, "run_wsi.py"), *common, "--micro", micro])
        if rc != 0:
            stage_status["wsi"] = f"FAILED rc={rc}"
            print(f"[wsi] FAILED (rc={rc}) -> stopping")
            _summary(a.sample, stage_status, protocol, status, dec, eligible, reason,
                     f"FAILED rc={rc}", qcp, "", th)
            return 1
        wsi_state = "ran"
    if eligible and a.warp:
        w = glob.glob(os.path.join(regdir, "*.ome.tif*"))
        wsi_out = w[0] if w else ""
    stage_status["wsi"] = wsi_state

    _summary(a.sample, stage_status, protocol, status, dec, eligible, reason, wsi_state, qcp, wsi_out, th)

    # exit: 0 if we reached an eligible selection (warp is opt-in); else the documented decision code
    return 0 if status == "ELIGIBLE_FOR_WARP" else code


def _summary(sample, stage_status, protocol, status, dec, eligible, reason, wsi_state, qcp, wsi_out, th):
    dec = dec or {}
    print("\n================ STAGE SUMMARY ================")
    print(f"  sample             {sample}")
    print(f"  stages             " + ", ".join(f"{k}={v}" for k, v in stage_status.items()))
    print(f"  status             {status}")
    print(f"  selected protocol  {protocol}   (rule: {dec.get('rule')})")
    print(f"  QC                 density_r={dec.get('sel_density_r')}  median_um={dec.get('sel_median_um')}"
          f"   [gate: r>={th['density_r_accept']} AND um<={th['median_um_accept']}]")
    print(f"  reason             {reason}")
    print(f"  WSI warp           {wsi_state}")
    print(f"  outputs            qc.json = {qcp}")
    if wsi_out:
        print(f"                     registered = {wsi_out}")
    print("==============================================")


if __name__ == "__main__":
    sys.exit(main())
