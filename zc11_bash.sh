#!/bin/bash

#!/bin/bash
#SBATCH --job-name=asced
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=

conda activate asced


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

for e in "${variants[@]}"
do
    echo "Running $e"
    uv run reproduce_zc11_asced_5G_LDPC.py "$e"
done