#!/usr/bin/env bash
# run_slide.sh -- transparent shell wrapper: run segment -> register -> qc in order,
# chained with && (stop on the first failure), each stage in its own conda env.
#
# This is EXACTLY equivalent to running the three commands by hand:
#     stardist_env/bin/python  run_segment.py  --samples ... --config ... --sample ...  && \
#     valis_env/bin/python      run_register.py --samples ... --config ... --sample ...  && \
#     qc_env/bin/python         run_qc.py       --samples ... --config ... --sample ...
#
# It adds NO extra behaviour: no qc.json read-back, no remapped exit codes, no selection
# summary, no WSI warp. The script's exit status IS the failing stage's status (or 0 if all
# three succeed). Protocol selection still happens inside run_qc.py and is written to qc.json,
# exactly as in the manual chain.
#
# For the decision-aware version (stage summary, exit code 0/2/3 by QC outcome, --resume,
# --dry-run, and the optional --warp), use run_slide.py instead.
#
# Usage:
#     ./run_slide.sh --samples samples.csv --config config.json --sample 0101165
#   (every argument is forwarded verbatim to all three stages, just like the manual chain.)
#
# Per-stage interpreters (override via env; defaults are the sibling envs next to this script):
#     STARDIST_PY   default: <script dir>/stardist_env/bin/python
#     VALIS_PY      default: <script dir>/valis_env/bin/python
#     QC_PY         default: <script dir>/qc_env/bin/python

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STARDIST_PY="${STARDIST_PY:-$HERE/stardist_env/bin/python}"
VALIS_PY="${VALIS_PY:-$HERE/valis_env/bin/python}"
QC_PY="${QC_PY:-$HERE/qc_env/bin/python}"

echo "+ $STARDIST_PY $HERE/run_segment.py  $*" >&2
"$STARDIST_PY" "$HERE/run_segment.py"  "$@" && \
echo "+ $VALIS_PY $HERE/run_register.py $*" >&2 && \
"$VALIS_PY"    "$HERE/run_register.py" "$@" && \
echo "+ $QC_PY $HERE/run_qc.py       $*" >&2 && \
"$QC_PY"       "$HERE/run_qc.py"       "$@"
