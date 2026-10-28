#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/compbio/projects/FFPE_filtering"
CONDA_BIN="/gpfs/commons/groups/nygcfaculty/kancero/jupyterhub/alt/miniconda3/bin/conda"
OUTPUT_DIR="$PROJECT_ROOT/prediction/pickles/retrained_shared_bins"

eval "$("$CONDA_BIN" shell.bash hook)"
set +u
conda activate pandasDash
set -u

which python
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/fontcache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

python "$PROJECT_ROOT/repos/fifa/src/retrain_shared_bins.py" \
    --data-root "$PROJECT_ROOT/prediction" \
    --output-dir "$OUTPUT_DIR" \
    "$@"
