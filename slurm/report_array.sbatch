#!/bin/bash
# Report array: one task per sample -> render <out>/<sample>/report.pdf (+ report.png).
# CPU-only, no GPU -- this is just plotting. Run AFTER qc_array (needs qc.json + warped nuclei);
# the raster/annotation panels also appear if run_wsi / run_annotate have run.
#   sbatch --array=0-$(( $(tail -n +2 samples.csv | wc -l) - 1 )) slurm/report_array.sbatch
# Then the cohort triage page (one quick CPU job, or run locally):
#   REPORT_PY run_report.py --samples samples.csv --config config.json --cohort
#SBATCH --job-name=hv_report
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --output=logs/report_%A_%a.log
set -uo pipefail

# ---- EDIT THESE ----
PIPE=/common/ruppinelab/Eldad/st_ref_v2/hest_valis_pipeline   # this repo
SAMPLES=$PIPE/samples.csv
CONFIG=$PIPE/config.json
REPORT_PY=/path/to/qc_env/bin/python   # env with numpy/scipy/matplotlib/pandas/tifffile/zarr
#                                        (SimpleITK optional -> enables the checkerboard panel)
# --------------------

mkdir -p "$PIPE/logs"
cd "$PIPE"; export PYTHONPATH="$PIPE:${PYTHONPATH:-}"
export MPLBACKEND=Agg
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$SAMPLES" | cut -d, -f1)
echo "[$(date +%T)] $(hostname) report sample=$SAMPLE"

"$REPORT_PY" run_report.py --samples "$SAMPLES" --config "$CONFIG" --sample "$SAMPLE" \
  || { echo "REPORT FAIL $SAMPLE"; exit 1; }
echo "[$(date +%T)] DONE $SAMPLE"

# Last task also renders the cohort triage page (every per-slide qc.json exists by now).
NSAMP=$(( $(tail -n +2 "$SAMPLES" | wc -l) ))
if [ "$SLURM_ARRAY_TASK_ID" -eq "$(( NSAMP - 1 ))" ]; then
  echo "[$(date +%T)] rendering cohort_report.pdf"
  "$REPORT_PY" run_report.py --samples "$SAMPLES" --config "$CONFIG" --cohort \
    || echo "COHORT REPORT FAIL (non-fatal)"
fi
