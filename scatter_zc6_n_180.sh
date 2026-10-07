#!/bin/bash

#SBATCH --array=0-37
#SBATCH --cpus-per-task=50
#SBATCH --ntasks=1
#SBATCH --time=0-04:00:00

# Tell OpenMP how many threads to use
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK



source /home/pj9034/miniconda3/etc/profile.d/conda.sh    # adjust path if needed
conda activate asced


n_simul=180
snr_start=2.5
target_fer=1e-3

variants=(
    nmsa
    aed
    2_split1_subsplit1 # splitting pattern [1,2,3,4,5]
    4_split1_subsplit2
    6_split1_subsplit3
    8_split1_subsplit4
    10_split1_subsplit5       
    12_split1 
    24_split1
    36_split1
    48_split1
    60_split1
    72_split1
    84_split1
    96_split1
    4_split2_subsplit1 # splitting pattern [2,4]
    8_split2_subsplit2
    12_split2 
    24_split2
    36_split2
    48_split2
    60_split2
    72_split2
    84_split2
    96_split2
    8_split3_subsplit1 # splitting pattern [3]
    16_split3 
    32_split3 
    48_split3
    64_split3
    80_split3
    96_split3
    112_split3
    128_split3
    64_nosplit # no splitting
    128_nosplit
    192_nosplit
)

variant=${variants[$SLURM_ARRAY_TASK_ID]}

echo "Running ${variant} (task ${SLURM_ARRAY_TASK_ID})"

uv run reproduce_scatter_plot_zc6_5G_LDPC.py \
    "$variant" \
    "$n_simul" \
    "$snr_start" \
    "$target_fer"