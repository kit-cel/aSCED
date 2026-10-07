#!/bin/bash
#SBATCH --job-name=asced_zc6_greedy_prod
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --cpus-per-task=32
#SBATCH --time=08:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err

# Manifest-driven array runner for the Zc=6 greedy-block-selection production
# validation sweep (reproduce_fig5_scatter_plot_zc6_5G_LDPC_greedy.py).
# Row schema: (decoder_variant, n_simul, snr_start, target_fer).
#
# Submit with the manifest path as $1 and an explicit --array range matching
# that manifest's row count, e.g.:
#   sbatch --array=1-70%10 sweeps/run_scatter_zc6_greedy_array.sbatch sweeps/scatter_zc6_greedy_manifest.csv

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

MANIFEST="$1"
mkdir -p sweeps/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

# +1 to skip the header line; SLURM array tasks are 1-indexed to match.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")

IFS=',' read -r decoder_variant n_simul snr_start target_fer <<< "${LINE}"

echo "Task ${SLURM_ARRAY_TASK_ID} (${MANIFEST}): decoder_variant=${decoder_variant} n_simul=${n_simul} snr_start=${snr_start} target_fer=${target_fer}"

uv run python reproduce_fig5_scatter_plot_zc6_5G_LDPC_greedy.py \
    "${decoder_variant}" "${n_simul}" "${snr_start}" "${target_fer}"
