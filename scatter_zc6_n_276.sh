#!/bin/bash

#SBATCH --job-name=asced
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --time=2-00:00:00

# Tell OpenMP how many threads to use
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


source /home/pj9034/miniconda3/etc/profile.d/conda.sh    # adjust path if needed
conda activate asced


n_simul=276
snr_start=2.6
target_fer=1e-3

variants=(
    nmsa
    12_split0
    24_split0
    48_split0
    12_split1
    24_split1
    48_split1
    16
    48
    64
    64_nosplit
    128
)

for e in "${variants[@]}"; do
    echo "Running $e (n_simul=${n_simul}, SNR=${snr_start}:${target_fer})"
    uv run reproduce_zc6_asced_5G_LDPC.py \
        "$e" \
        "$n_simul" \
        "$snr_start" \
        "$target_fer"
done
