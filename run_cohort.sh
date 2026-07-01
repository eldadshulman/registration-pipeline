#!/usr/bin/env bash
# run_cohort.sh -- drive the whole-cohort SLURM pipeline end to end, in order, BLOCKING on each
# stage before the next starts. This is the cluster analog of run_slide.sh, but it cannot be a
# plain `&&` chain: the stages are SLURM array jobs (sbatch submits and returns immediately), and
# the wsi array's size is only known after run_select writes wsi_manifest.csv mid-pipeline. So each
# SLURM step is launched with `sbatch --wait` (blocks until the whole array finishes), and the wsi
# array is sized right before it is submitted.
#
# Stages (documented order):
#   1. qc_array       (GPU array, 1 task/sample)   segment -> register -> QC + per-slide selection
#   2. run_select     (CPU)                        aggregate -> per_slide_decision.csv + wsi_manifest.csv
#   3. wsi_array      (CPU array, sized to manifest) warp chosen-protocol H&E -> registered OME-TIFF
#   4. report_array   (CPU array, 1 task/sample)   per-slide + cohort QC report
#   5. run_provenance (CPU)                        acceptance gate + audit table
#
# run_wsi is slow (hours per slide), so this driver blocks for a long time. Run it under tmux/screen
# or `nohup`, or submit it as a tiny long-lived defq job -- it uses almost nothing itself (it waits).
#
# FAILURE POLICY: individual slide (array-task) failures do NOT abort the cohort by default --
# run_select and run_provenance aggregate whatever completed and report the misses. Pass --strict to
# abort if any array task fails. A failure of an inline aggregation step (run_select / run_provenance)
# always stops the pipeline.
#
# Usage:
#   ./run_cohort.sh [--samples samples.csv] [--config config.json]
#                   [--dry-run] [--strict] [--skip-wsi] [--skip-report]
#
# The array TASKS read the samples.csv/config.json baked into slurm/*.sbatch (their EDIT-THESE
# block); keep --samples/--config here pointing at the SAME files (the defaults do). SELECT_PY is the
# interpreter for the two inline python steps + the output_dir lookup (default: qc_env next to this).
# Make sure your SLURM client is on PATH first (on Cedars: module load slurm/slurm-compbio/23.11.10).

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"; export PYTHONPATH="$HERE:${PYTHONPATH:-}"

SAMPLES="$HERE/samples.csv"; CONFIG="$HERE/config.json"
DRYRUN=0; STRICT=0; SKIP_WSI=0; SKIP_REPORT=0
SELECT_PY="${SELECT_PY:-$HERE/qc_env/bin/python}"

while [ $# -gt 0 ]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2;;
    --config)  CONFIG="$2";  shift 2;;
    --dry-run) DRYRUN=1; shift;;
    --strict)  STRICT=1; shift;;
    --skip-wsi) SKIP_WSI=1; shift;;
    --skip-report) SKIP_REPORT=1; shift;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not on PATH (load your SLURM client first)"; exit 2; }
[ -f "$SAMPLES" ] || { echo "ERROR: samples file not found: $SAMPLES"; exit 2; }

nrows() { if [ -f "$1" ]; then local n; n=$(tail -n +2 "$1" | grep -c .); echo "${n:-0}"; else echo 0; fi; }

NSAMP=$(nrows "$SAMPLES")
[ "$NSAMP" -ge 1 ] || { echo "ERROR: no data rows in $SAMPLES"; exit 2; }

# output_dir (where run_select writes wsi_manifest.csv) straight from the pipeline config loader
OUT=$("$SELECT_PY" -c "from hest_valis import config; print(config.load_config('$CONFIG')['output_dir'])" 2>/dev/null) \
  || { echo "ERROR: could not read output_dir from $CONFIG via $SELECT_PY"; exit 2; }

echo "=== cohort run ==="
echo "  samples      $SAMPLES  ($NSAMP samples)"
echo "  config       $CONFIG"
echo "  output_dir   $OUT"
echo "  select_py    $SELECT_PY"
echo "  task failure $([ "$STRICT" = 1 ] && echo 'abort (--strict)' || echo 'warn + continue (aggregators report misses)')"
[ "$DRYRUN" = 1 ] && echo "  MODE         DRY-RUN (nothing submitted)"

# submit one array with --wait; honor --strict on task failure
run_array() {  # label  sbatch-relpath  N
  local label="$1" f="$2" n="$3" rc
  if [ "$n" -lt 1 ]; then echo "[$label] 0 tasks -> skip"; return 0; fi
  echo "[$label] sbatch --wait --array=0-$((n-1)) $f"
  [ "$DRYRUN" = 1 ] && return 0
  if sbatch --wait --array=0-$((n-1)) "$HERE/$f"; then
    echo "[$label] all $n tasks OK"
  else
    rc=$?
    [ "$STRICT" = 1 ] && { echo "[$label] task failure (rc=$rc) -> abort (--strict)"; exit 1; }
    echo "[$label] WARNING: >=1 task failed (rc=$rc) -> continuing (aggregators report misses)"
  fi
}

run_py() {  # label  python-args...
  local label="$1"; shift
  echo "[$label] $SELECT_PY $*"
  [ "$DRYRUN" = 1 ] && return 0
  "$SELECT_PY" "$@" || { echo "[$label] FAILED -> abort"; exit 1; }
}

# 1) QC array: segment -> register -> qc + per-slide selection
run_array qc_array slurm/qc_array.sbatch "$NSAMP"

# 2) aggregate the per-slide decisions -> wsi_manifest.csv (feeds the wsi array)
run_py run_select "$HERE/run_select.py" --samples "$SAMPLES" --config "$CONFIG"

# 3) WSI array, sized to the manifest run_select just wrote (coarse/rescued are already excluded)
if [ "$SKIP_WSI" = 1 ]; then
  echo "[wsi_array] skipped (--skip-wsi)"
elif [ "$DRYRUN" = 1 ]; then
  echo "[wsi_array] would size --array to rows of $OUT/wsi_manifest.csv (unknown until run_select runs)"
else
  run_array wsi_array slurm/wsi_array.sbatch "$(nrows "$OUT/wsi_manifest.csv")"
fi

# 4) per-slide + cohort QC report (after wsi so raster/annotation panels are populated)
if [ "$SKIP_REPORT" = 1 ]; then
  echo "[report_array] skipped (--skip-report)"
else
  run_array report_array slurm/report_array.sbatch "$NSAMP"
fi

# 5) acceptance gate + audit table (accepted vs manual-review)
run_py run_provenance "$HERE/run_provenance.py" --samples "$SAMPLES" --config "$CONFIG"

echo "=== cohort run complete ==="
echo "  decisions    $OUT/per_slide_decision.csv"
echo "  wsi manifest $OUT/wsi_manifest.csv"
echo "  (acceptance gate / audit printed by run_provenance above)"
