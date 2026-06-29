#!/bin/bash
#SBATCH --job-name=mol_020_E_nbo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --output=mol_020_E_nbo_%j.out
#SBATCH --error=mol_020_E_nbo_%j.err
# Adjust module name and --partition / --account for your cluster:
module load gaussian/16

g16 < mol_020_E_nbo.gjf > mol_020_E_nbo.log
