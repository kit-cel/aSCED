#!/bin/bash
#SBATCH --job-name=asced_greedy
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err

# Manifest-driven array runner for the greedy RL/QC-block selection search
# (sweeps/greedy_block_search.py). Modeled on sweeps/run_manifest_array.sbatch,
# but the row schema differs (ensemble_name, fixed_blocks, candidate_block, snr)
# since the greedy search's CLI differs from the main reproduce script.
#
# Submit with the manifest path as $1 and an explicit --array range matching
# that manifest's row count, e.g.:
#   sbatch --array=1-68%34 sweeps/run_greedy_search_array.sbatch sweeps/greedy_round_manifest.csv

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

MANIFEST="$1"
mkdir -p sweeps/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

# +1 to skip the header line; SLURM array tasks are 1-indexed to match.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")

IFS=',' read -r ensemble_name fixed_blocks candidate_block snr <<< "${LINE}"

echo "Task ${SLURM_ARRAY_TASK_ID} (${MANIFEST}): ensemble=${ensemble_name} fixed_blocks=${fixed_blocks} candidate_block=${candidate_block} snr=${snr}"

uv run python sweeps/greedy_block_search.py \
    "${ensemble_name}" "${fixed_blocks}" "${candidate_block}" "${snr}"
