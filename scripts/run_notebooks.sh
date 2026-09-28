#!/usr/bin/env bash
# Execute the analysis notebooks in place on the analysis machine (kernel: conda env `sc`). Usage: sh scripts/run_notebooks.sh [notebook ...]
set -u
cd "$(dirname "$0")/.." || exit 1
PY=${OLIGOMETAB_PYTHON:-$HOME/miniconda3/envs/sc/bin}
set -a; [ -f .env ] && . ./.env; set +a
mkdir -p logs
NBS=${@:-notebooks/analysis/analysis_spatial_metabolism.ipynb notebooks/analysis/analysis_public_datasets_metabolism.ipynb}
for nb in $NBS; do
  name=$(basename "$nb" .ipynb)
  echo "[$(date '+%F %T')] start $name" | tee -a logs/run.log
  ( cd "$(dirname "$nb")" && "$PY/jupyter" nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=sc "$(basename "$nb")" ) > "logs/$name.log" 2>&1
  echo "[$(date '+%F %T')] end $name (exit $?)" | tee -a logs/run.log
done
