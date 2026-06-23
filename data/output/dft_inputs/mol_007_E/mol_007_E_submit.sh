#!/bin/bash
#SBATCH --job-name=mol_007_E
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --output=mol_007_E_%j.out
#SBATCH --error=mol_007_E_%j.err
# Adjust the module name and add --partition / --account as needed:
module load gaussian/16

g16 < mol_007_E.gjf > mol_007_E.log
