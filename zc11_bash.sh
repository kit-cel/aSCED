#!/bin/bash

#!/bin/bash
#SBATCH --job-name=asced
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=256



source /home/pj9034/miniconda3/etc/profile.d/conda.sh    # adjust path if needed
conda activate asced

n_simul=143
snr_start=2.0
snr_end=4.5

variants=(
    nmsa
    asced22
    asced44
    asced88
    asced24
    asced48
    asced96_4batch
    asced96
    asced192
    asced2048
)

for e in "${variants[@]}"; do
    echo "Running $e (n_simul=${n_simul}, SNR=${snr_start}:${snr_end})"
    uv run reproduce_zc11_asced_5G_LDPC.py \
        "$e" \
        "$n_simul" \
        "$snr_start" \
        "$snr_end"
done
