#!/bin/bash

#!/bin/bash
#SBATCH --job-name=asced
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=256



source /home/pj9034/miniconda3/etc/profile.d/conda.sh    # adjust path if needed
conda activate asced


n_simul=78
snr_start=2.0
snr_end=4.5

variants=(
    #nmsa
    #12_split0
    #24_split0
    #48_split0
    #12_split1
    #24_split1
    #48_split1
    #16
    #48
    64
    64_nosplit
    128
)

for e in "${variants[@]}"; do
    echo "Running $e (n_simul=${n_simul}, SNR=${snr_start}:${snr_end})"
    uv run reproduce_zc6_asced_5G_LDPC.py \
        "$e" \
        "$n_simul" \
        "$snr_start" \
        "$snr_end"
done
