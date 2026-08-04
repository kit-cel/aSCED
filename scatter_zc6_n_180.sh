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


n_simul=180
snr_start=2.5
target_fer=1e-3

variants=(
    # nmsa
    # 12_split1 # splitting pattern [1,2,3,4,5]
    # 24_split1
    # 48_split1
    # 60_split1
    # 96_split1
    # 12_split2 # splitting pattern [2,4]
    # 24_split2
    # 48_split2
    # 96_split2
    # 16_split3 # splitting pattern [3]
    # 48_split3
    # 64_split3
    # 96_split3
    # 128_split3
    # 64_nosplit # no splitting
    # 128_nosplit
    192_nosplit
)

for e in "${variants[@]}"; do
    echo "Running $e (n_simul=${n_simul}, SNR=${snr_start}:${target_fer})"
    uv run reproduce_scatter_plot_zc6_5G_LDPC.py \
        "$e" \
        "$n_simul" \
        "$snr_start" \
        "$target_fer"
done
