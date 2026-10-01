#!/bin/bash
#SBATCH --job-name=asced_seq_rl
#SBATCH --partition=all
#SBATCH --account=cel-assis
#SBATCH --array=1-70%10
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=sweeps/logs/%x_%A_%a.out
#SBATCH --error=sweeps/logs/%x_%A_%a.err

# One task per row of sweeps/sequential_rl_manifest.csv (see
# generate_sequential_rl_manifest.py). Each task runs one (variant, selector,
# snr) config of reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

mkdir -p sweeps/logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate asced

MANIFEST="sweeps/sequential_rl_manifest.csv"
# +1 to skip the header line; SLURM array tasks are 1-indexed to match.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")

IFS=',' read -r variant n_simul snr selector members_per_group target_num_converged <<< "${LINE}"

echo "Task ${SLURM_ARRAY_TASK_ID}: variant=${variant} n_simul=${n_simul} snr=${snr} selector=${selector} members_per_group=${members_per_group} target_num_converged=${target_num_converged}"

uv run python reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py \
    "${variant}" "${n_simul}" "${snr}" "${snr}" "${selector}" "${members_per_group}" "${target_num_converged}"
