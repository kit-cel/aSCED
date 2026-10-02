#!/bin/bash
#SBATCH --job-name=asced_sched_study
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err

# Manifest-driven array runner for the scheduling/iteration-count study
# (sweeps/scheduling_study.py), modeled on sweeps/run_manifest_array.sbatch.
# Submit with the manifest path as $1 and an explicit --array range matching
# that manifest's row count, e.g.:
#   sbatch --array=1-126%10 sweeps/scheduling_study_array.sbatch sweeps/scheduling_manifest.csv
#
# One task per manifest row: (variant, n_simul, snr, scheduling_config,
# max_iterations) for sweeps/scheduling_study.py.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

MANIFEST="$1"
mkdir -p sweeps/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

# +1 to skip the header line; SLURM array tasks are 1-indexed to match.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")

IFS=',' read -r variant n_simul snr scheduling_config max_iterations <<< "${LINE}"

echo "Task ${SLURM_ARRAY_TASK_ID} (${MANIFEST}): variant=${variant} n_simul=${n_simul} snr=${snr} scheduling_config=${scheduling_config} max_iterations=${max_iterations}"

uv run python sweeps/scheduling_study.py \
    "${variant}" "${n_simul}" "${snr}" "${snr}" "${scheduling_config}" "${max_iterations}"
