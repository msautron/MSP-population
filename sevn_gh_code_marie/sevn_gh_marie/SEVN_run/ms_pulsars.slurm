#!/bin/bash
#SBATCH --job-name=WD_pulsars_study
#SBATCH --partition=Long
#SBATCH --time infinite
#SBATCH --nodes=1
#SBATCH --exclude=obas-hpc-nr6
#SBATCH --mem-per-cpu=1gb
#SBATCH --output=logs/log-%j.out
#SBATCH --error=logs/log-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=matteo.sautron@astro.unistra.fr

NB_SIM=$1

./birth_and_evol 10000 "$NB_SIM"
