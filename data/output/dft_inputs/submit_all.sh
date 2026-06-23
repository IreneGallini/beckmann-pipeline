#!/bin/bash
# Run from the dft_inputs/ directory on the HPC cluster:
#   bash submit_all.sh
for dir in */; do
    name="${dir%/}"
    cd "$dir" && sbatch "${name}_submit.sh" && cd ..
done
