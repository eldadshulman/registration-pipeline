"""Tests for run_cohort.sh -- the blocking whole-cohort SLURM driver.

`sbatch` (on PATH) and the inline python interpreter ($SELECT_PY) are replaced with stubs that log
their argv, so nothing real is submitted. We assert the orchestration contract:
  - stages run in order: qc_array -> run_select -> wsi_array -> report_array -> run_provenance,
  - array sizes are correct (qc/report = #samples, wsi = #wsi_manifest rows),
  - --dry-run submits nothing,
  - --strict aborts the cohort on an array-task failure; the default warns and continues,
  - --skip-wsi drops the warp stage.
"""
import os
import stat
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
WRAPPER = os.path.join(HERE, "..", "run_cohort.sh")

SBATCH_STUB = """#!/usr/bin/env bash
echo "sbatch $*" >> "$RUN_LOG"
if [ -n "${SBATCH_FAIL_ON:-}" ] && [[ "$*" == *"$SBATCH_FAIL_ON"* ]]; then exit 1; fi
exit 0
"""

# stub interpreter: `-c ...` -> print the fake output_dir; otherwise log "PY <script> <args>"
SELECT_STUB = """#!/usr/bin/env bash
if [ "${1:-}" = "-c" ]; then echo "$FAKE_OUT"; exit 0; fi
echo "PY $(basename "${1:-}") ${@:2}" >> "$RUN_LOG"
exit 0
"""


def _exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def _run(tmp_path, extra_args=(), sbatch_fail_on="", nsamp=3, nwsi=2):
    bindir = tmp_path / "bin"; bindir.mkdir()
    _exe(bindir / "sbatch", SBATCH_STUB)
    select_py = _exe(tmp_path / "select_stub.sh", SELECT_STUB)

    out = tmp_path / "out"; out.mkdir()
    (out / "wsi_manifest.csv").write_text("sample_id\n" + "".join(f"S{i}\n" for i in range(nwsi)))

    samples = tmp_path / "samples.csv"
    samples.write_text("sample_id,he_path,dapi_path,xenium_cells\n"
                       + "".join(f"S{i},a,b,c\n" for i in range(nsamp)))
    config = tmp_path / "config.json"; config.write_text("{}")
    log = tmp_path / "calls.log"

    env = dict(os.environ, PATH=f"{bindir}:{os.environ['PATH']}", RUN_LOG=str(log),
               FAKE_OUT=str(out), SELECT_PY=select_py, SBATCH_FAIL_ON=sbatch_fail_on)
    r = subprocess.run(["bash", WRAPPER, "--samples", str(samples), "--config", str(config), *extra_args],
                       env=env, capture_output=True, text=True)
    lines = log.read_text().splitlines() if log.exists() else []
    return r.returncode, lines


def test_full_cohort_order_and_array_sizes(tmp_path):
    rc, lines = _run(tmp_path)
    assert rc == 0
    assert len(lines) == 5
    assert "--array=0-2" in lines[0] and "qc_array.sbatch" in lines[0]      # 3 samples
    assert lines[1].startswith("PY run_select.py")
    assert "--array=0-1" in lines[2] and "wsi_array.sbatch" in lines[2]     # 2 manifest rows
    assert "--array=0-2" in lines[3] and "report_array.sbatch" in lines[3]  # 3 samples
    assert lines[4].startswith("PY run_provenance.py")


def test_dry_run_submits_nothing(tmp_path):
    rc, lines = _run(tmp_path, extra_args=("--dry-run",))
    assert rc == 0
    assert lines == []


def test_strict_aborts_on_array_task_failure(tmp_path):
    rc, lines = _run(tmp_path, extra_args=("--strict",), sbatch_fail_on="wsi_array.sbatch")
    joined = "\n".join(lines)
    assert "wsi_array.sbatch" in joined          # wsi was attempted
    assert "report_array.sbatch" not in joined   # and the cohort aborted before report
    assert "run_provenance.py" not in joined
    assert rc == 1


def test_default_continues_past_array_task_failure(tmp_path):
    rc, lines = _run(tmp_path, sbatch_fail_on="wsi_array.sbatch")   # no --strict
    joined = "\n".join(lines)
    assert "report_array.sbatch" in joined       # aggregators/report still run
    assert "run_provenance.py" in joined
    assert rc == 0


def test_skip_wsi_drops_the_warp_stage(tmp_path):
    rc, lines = _run(tmp_path, extra_args=("--skip-wsi",))
    joined = "\n".join(lines)
    assert "wsi_array.sbatch" not in joined
    assert "qc_array.sbatch" in joined and "report_array.sbatch" in joined
    assert "PY run_provenance.py" in joined
    assert rc == 0
