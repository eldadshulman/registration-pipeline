"""Tests for run_slide.py -- the thin single-slide wrapper.

Stage commands are MOCKED (run_stage is patched), so nothing real runs. We verify the
orchestration contract:
  - WSI is not invoked after a QC failure,
  - WSI is not invoked without --warp,
  - WSI IS invoked after an eligible selection with --warp (and the protocol is passed through),
  - a stage failure stops subsequent stages,
  - --resume skips completed stages.
"""
import json
import os
import sys
from unittest import mock

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_slide


# --------------------------------------------------------------------------- helpers
def _setup(tmp_path, qc=None, artifacts=()):
    """Build config.json + samples.csv under tmp; optionally seed artifacts + a qc.json."""
    out_root = tmp_path / "out"
    sample = "S1"
    sdir = out_root / sample
    sdir.mkdir(parents=True, exist_ok=True)
    (tmp_path / "config.json").write_text(json.dumps({"output_dir": str(out_root)}))
    (tmp_path / "samples.csv").write_text(
        "sample_id,he_path,dapi_path,xenium_cells\nS1,a.svs,b.tif,c.parquet\n")
    for f in artifacts:
        (sdir / f).write_text("x")
    if qc is not None:
        (sdir / "qc.json").write_text(json.dumps(qc))
    return str(tmp_path / "config.json"), str(tmp_path / "samples.csv"), sample


def _qc(chosen, r, u, rule="lower_um"):
    return {"decision": {"chosen": chosen, "rule": rule, "sel_density_r": r, "sel_median_um": u},
            "metrics": {}}


def _recorder():
    calls = []
    def fake(cmd):
        calls.append(list(cmd))
        return 0
    return calls, fake


def _has(calls, token):
    return [c for c in calls if any(token in part for part in c)]


# --------------------------------------------------------------------------- tests
def test_wsi_not_called_after_qc_failure(tmp_path):
    cfg, samples, s = _setup(tmp_path, qc=_qc(None, None, None),
                             artifacts=("he_nuclei.npy", "he_nuclei_nomicro.npy"))
    calls, fake = _recorder()
    with mock.patch.object(run_slide, "run_stage", side_effect=fake):
        rc = run_slide.main(["--samples", samples, "--config", cfg, "--sample", s, "--warp", "--resume"])
    assert not _has(calls, "run_wsi.py")          # QC_FAILED -> never warp, even with --warp
    assert rc == 3                                  # documented QC_FAILED exit code


def test_wsi_not_called_without_warp(tmp_path):
    cfg, samples, s = _setup(tmp_path, qc=_qc("micro", 0.92, 3.0),
                             artifacts=("he_nuclei.npy", "he_nuclei_nomicro.npy"))
    calls, fake = _recorder()
    with mock.patch.object(run_slide, "run_stage", side_effect=fake):
        rc = run_slide.main(["--samples", samples, "--config", cfg, "--sample", s, "--resume"])
    assert not _has(calls, "run_wsi.py")          # eligible, but no --warp -> stop after selection
    assert rc == 0


@pytest.mark.parametrize("chosen,micro_flag", [("micro", "1"), ("nomicro", "0")])
def test_wsi_called_with_warp_and_protocol_passthrough(tmp_path, chosen, micro_flag):
    cfg, samples, s = _setup(tmp_path, qc=_qc(chosen, 0.92, 3.0),
                             artifacts=("he_nuclei.npy", "he_nuclei_nomicro.npy"))
    calls, fake = _recorder()
    with mock.patch.object(run_slide, "run_stage", side_effect=fake):
        rc = run_slide.main(["--samples", samples, "--config", cfg, "--sample", s, "--resume", "--warp"])
    wsi = _has(calls, "run_wsi.py")
    assert wsi, "eligible slide with --warp must invoke run_wsi"
    cmd = wsi[0]
    assert "--micro" in cmd and cmd[cmd.index("--micro") + 1] == micro_flag   # protocol passed through
    assert rc == 0


def test_stage_failure_stops_execution(tmp_path):
    cfg, samples, s = _setup(tmp_path)            # no artifacts -> segment + register will run
    calls = []
    def fail_on_register(cmd):
        calls.append(list(cmd))
        return 1 if any("run_register.py" in p for p in cmd) else 0
    with mock.patch.object(run_slide, "run_stage", side_effect=fail_on_register):
        rc = run_slide.main(["--samples", samples, "--config", cfg, "--sample", s, "--warp"])
    assert _has(calls, "run_register.py")         # register attempted
    assert not _has(calls, "run_qc.py")           # qc not reached
    assert not _has(calls, "run_wsi.py")          # wsi not reached
    assert rc == 1                                  # stage-failure exit code


def test_resume_skips_completed_stages(tmp_path):
    cfg, samples, s = _setup(tmp_path, qc=_qc("micro", 0.92, 3.0),
                             artifacts=("he_nuclei.npy", "he_nuclei_nomicro.npy"))
    calls, fake = _recorder()
    with mock.patch.object(run_slide, "run_stage", side_effect=fake):
        rc = run_slide.main(["--samples", samples, "--config", cfg, "--sample", s, "--resume"])
    assert calls == []                             # segment/register/qc all skipped (outputs exist)
    assert rc == 0


def test_below_gate_is_review_required_not_warped(tmp_path):
    # selected as-is but under the accept gate -> REVIEW_REQUIRED, not warped even with --warp
    cfg, samples, s = _setup(tmp_path, qc=_qc("micro", 0.50, 3.0),
                             artifacts=("he_nuclei.npy", "he_nuclei_nomicro.npy"))
    calls, fake = _recorder()
    with mock.patch.object(run_slide, "run_stage", side_effect=fake):
        rc = run_slide.main(["--samples", samples, "--config", cfg, "--sample", s, "--resume", "--warp"])
    assert not _has(calls, "run_wsi.py")
    assert rc == 2                                  # documented REVIEW_REQUIRED exit code
