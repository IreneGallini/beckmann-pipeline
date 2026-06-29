#!/bin/bash
#SBATCH --job-name=mol_023_Z
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --output=mol_023_Z_%j.out
#SBATCH --error=mol_023_Z_%j.err
# Adjust the module name and add --partition / --account as needed:
module load gaussian/16

g16 < mol_023_Z.gjf > mol_023_Z.log
