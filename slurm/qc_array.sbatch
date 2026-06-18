#!/bin/bash
# QC array: one task per sample -> segment (StarDist) -> register+warp nuclei (valis) -> QC.
# Submit with --array sized to your samples.csv (rows after the header), e.g. 19 samples:
#   sbatch --array=0-18 slurm/qc_array.sbatch
#SBATCH --job-name=hv_qc
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=logs/qc_%A_%a.log
set -uo pipefail

# ---- EDIT THESE ----
PIPE=/common/ruppinelab/Eldad/st_ref_v2/hest_valis_pipeline   # this repo
SAMPLES=$PIPE/samples.csv
CONFIG=$PIPE/config.json
STARDIST_PY=/path/to/stardist_env/bin/python                  # env with stardist+tensorflow
VALIS_PY=/path/to/valis_hest_env/bin/python                   # patched valis env (env/setup.md)
QC_PY=/path/to/qc_env/bin/python                              # env with numpy/scipy/pandas/tifffile/zarr
VALIS_ENV=/path/to/valis_hest_env                             # for LD_LIBRARY_PATH / JAVA_HOME
# --------------------

mkdir -p "$PIPE/logs"
cd "$PIPE"; export PYTHONPATH="$PIPE:${PYTHONPATH:-}"
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$SAMPLES" | cut -d, -f1)
echo "[$(date +%T)] $(hostname) sample=$SAMPLE"

echo "[$(date +%T)] 1/3 segment (StarDist)"
"$STARDIST_PY" run_segment.py --samples "$SAMPLES" --config "$CONFIG" --sample "$SAMPLE" \
  || { echo "SEGMENT FAIL $SAMPLE"; exit 1; }

echo "[$(date +%T)] 2/3 register + warp nuclei (valis, CPU)"
export LD_LIBRARY_PATH="$VALIS_ENV/lib"; export PATH="$VALIS_ENV/bin:$PATH"; export JAVA_HOME="$VALIS_ENV"
export CUDA_VISIBLE_DEVICES=""
"$VALIS_PY" run_register.py --samples "$SAMPLES" --config "$CONFIG" --sample "$SAMPLE" \
  || { echo "REGISTER FAIL $SAMPLE"; exit 1; }

echo "[$(date +%T)] 3/3 concordance QC + protocol selection"
"$QC_PY" run_qc.py --samples "$SAMPLES" --config "$CONFIG" --sample "$SAMPLE" \
  || { echo "QC FAIL $SAMPLE"; exit 1; }
echo "[$(date +%T)] DONE $SAMPLE"
