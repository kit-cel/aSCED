#!/bin/bash
#SBATCH --job-name=asced_seq_rl
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err

# Generic manifest-driven array runner. Submit with the manifest path as $1
# and an explicit --array range matching that manifest's row count, e.g.:
#   sbatch --array=1-56%10 sweeps/run_manifest_array.sbatch sweeps/parallelism_ordering_manifest.csv
#
# One task per manifest row: (variant, n_simul, snr, selector,
# members_per_group, target_num_converged) for
# reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

MANIFEST="$1"
mkdir -p sweeps/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

# +1 to skip the header line; SLURM array tasks are 1-indexed to match.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")

IFS=',' read -r variant n_simul snr selector members_per_group target_num_converged <<< "${LINE}"

echo "Task ${SLURM_ARRAY_TASK_ID} (${MANIFEST}): variant=${variant} n_simul=${n_simul} snr=${snr} selector=${selector} members_per_group=${members_per_group} target_num_converged=${target_num_converged}"

uv run python reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py \
    "${variant}" "${n_simul}" "${snr}" "${snr}" "${selector}" "${members_per_group}" "${target_num_converged}"
