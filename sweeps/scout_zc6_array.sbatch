#!/bin/bash
#SBATCH --job-name=zc6_scout
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --cpus-per-task=32
#SBATCH --time=4:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

n_simuls=(78 180)
n_simul=${n_simuls[$SLURM_ARRAY_TASK_ID]}

uv run python sweeps/scout_zc6_fer_targets.py "${n_simul}"
