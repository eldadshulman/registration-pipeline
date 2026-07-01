"""Tests for run_slide.sh -- the transparent shell chain (segment -> register -> qc).

We replace the three per-stage interpreters with a stub that just logs its argv (and can be told
to fail on a chosen stage), then assert the wrapper's contract -- i.e. that it behaves EXACTLY
like running the three commands by hand with `&&`:
  - all three stages run, in order, each receiving the forwarded CLI args,
  - a stage failure stops the chain (later stages are not run) and its exit code propagates,
  - no extra logic (no qc.json read-back / remapped exit codes): the stub stands in for run_qc.py.
"""
import os
import stat
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
WRAPPER = os.path.join(HERE, "..", "run_slide.sh")

STUB = """#!/usr/bin/env bash
# stub interpreter: log "<script-basename> <forwarded-args...>"; fail if $1 matches $FAIL_ON
echo "$(basename "$1") ${@:2}" >> "$RUN_LOG"
if [ -n "$FAIL_ON" ] && [[ "$1" == *"$FAIL_ON"* ]]; then exit 7; fi
exit 0
"""


def _stub(tmp_path):
    p = tmp_path / "stub.sh"
    p.write_text(STUB)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return str(p)


def _run(tmp_path, fail_on=""):
    log = tmp_path / "calls.log"
    stub = _stub(tmp_path)
    env = dict(os.environ, RUN_LOG=str(log), FAIL_ON=fail_on,
               STARDIST_PY=stub, VALIS_PY=stub, QC_PY=stub)
    r = subprocess.run(
        ["bash", WRAPPER, "--samples", "samples.csv", "--config", "config.json", "--sample", "S1"],
        env=env, capture_output=True, text=True)
    lines = log.read_text().splitlines() if log.exists() else []
    return r.returncode, lines


def test_runs_all_three_in_order_with_forwarded_args(tmp_path):
    rc, lines = _run(tmp_path)
    assert rc == 0
    assert [l.split()[0] for l in lines] == ["run_segment.py", "run_register.py", "run_qc.py"]
    for l in lines:                                   # every stage got the same forwarded args
        assert "--samples samples.csv --config config.json --sample S1" in l


def test_stage_failure_stops_chain_and_propagates_code(tmp_path):
    rc, lines = _run(tmp_path, fail_on="run_register.py")
    assert [l.split()[0] for l in lines] == ["run_segment.py", "run_register.py"]   # qc not reached
    assert rc == 7                                                                   # code propagates
