#!/bin/bash

#SBATCH --job-name=asced
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --time=1-00:00:00

# Tell OpenMP how many threads to use
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source /home/pj9034/miniconda3/etc/profile.d/conda.sh    # adjust path if needed
conda activate asced


n_simul=132
snr_start=1.0
snr_end=5


variants=(
    #nmsa
    #asced22
    #asced44
    #asced88
    #asced24
    asced48
    #asced96_4batch
    #asced192_8batch
    #asced96
    #asced192
    #asced288
    asced384
    #asced2048
)

for e in "${variants[@]}"; do
    echo "Running $e (n_simul=${n_simul}, SNR=${snr_start}:${snr_end})"
    uv run reproduce_zc11_asced_5G_LDPC.py \
        "$e" \
        "$n_simul" \
        "$snr_start" \
        "$snr_end"
done
