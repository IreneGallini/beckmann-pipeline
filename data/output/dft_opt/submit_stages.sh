#!/bin/bash
# Two-stage submission for DFT test set.
# Run from the dft_opt_test/ directory on the HPC cluster.
#
# Stage 1: submit all optimisations
for dir in */; do
    name="${dir%/}"
    cd "$dir" && sbatch "${name}_opt_submit.sh" && cd ..
done
#
# Stage 2: after ALL Stage 1 jobs finish, submit NBO single-points
# (the _nbo.gjf reads the .chk from Stage 1, so Stage 1 must complete first)
for dir in */; do
    name="${dir%/}"
    cd "$dir" && sbatch "${name}_nbo_submit.sh" && cd ..
done
