#!/bin/bash
#SBATCH --job-name=asced_alpha_sweep
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err

# Manifest-driven array runner for the norm_factor (NMSA normalization
# constant) sweep. Submit with an --array range matching the manifest's row
# count, e.g.:
#   sbatch --array=1-189%10 sweeps/alpha_sweep_array.sbatch sweeps/alpha_sweep_manifest.csv
#
# One task per manifest row: (variant, n_simul, snr, norm_factor) for
# sweeps/alpha_sweep.py.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

MANIFEST="${1:-sweeps/alpha_sweep_manifest.csv}"
mkdir -p sweeps/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

# +1 to skip the header line; SLURM array tasks are 1-indexed to match.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")

IFS=',' read -r variant n_simul snr norm_factor <<< "${LINE}"

echo "Task ${SLURM_ARRAY_TASK_ID} (${MANIFEST}): variant=${variant} n_simul=${n_simul} snr=${snr} norm_factor=${norm_factor}"

uv run python sweeps/alpha_sweep.py \
    "${variant}" "${n_simul}" "${snr}" "${snr}" "${norm_factor}"
