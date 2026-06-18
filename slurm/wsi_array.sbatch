#!/bin/bash
# WSI array: warp the chosen-protocol H&E image for each sample -> registered OME-TIFF.
# Run AFTER run_select.py has written wsi_manifest.csv. Size --array to the manifest rows.
#   sbatch --array=0-18 slurm/wsi_array.sbatch
#
# NOTE: the image warp is CPU-only and slow (serial tile read; see env/setup.md). Give it a
# GENEROUS time limit. It does NOT need a GPU.
#SBATCH --job-name=hv_wsi
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=14:00:00
#SBATCH --requeue
#SBATCH --output=logs/wsi_%A_%a.log
set -uo pipefail

# ---- EDIT THESE ----
PIPE=/common/ruppinelab/Eldad/st_ref_v2/hest_valis_pipeline
SAMPLES=$PIPE/samples.csv
CONFIG=$PIPE/config.json
VALIS_PY=/path/to/valis_hest_env/bin/python   # PATCHED valis env (env/setup.md)
VALIS_ENV=/path/to/valis_hest_env
# --------------------

mkdir -p "$PIPE/logs"
cd "$PIPE"; export PYTHONPATH="$PIPE:${PYTHONPATH:-}"
OUT=$("$VALIS_PY" -c "from hest_valis import config; print(config.load_config('$CONFIG')['output_dir'])")
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$OUT/wsi_manifest.csv" | cut -d, -f1)
echo "[$(date +%T)] $(hostname) sample=$SAMPLE"

export LD_LIBRARY_PATH="$VALIS_ENV/lib"; export PATH="$VALIS_ENV/bin:$PATH"; export JAVA_HOME="$VALIS_ENV"
export CUDA_VISIBLE_DEVICES=""
"$VALIS_PY" run_wsi.py --samples "$SAMPLES" --config "$CONFIG" --sample "$SAMPLE" \
  || { echo "WSI FAIL $SAMPLE"; exit 1; }
echo "[$(date +%T)] DONE $SAMPLE"
