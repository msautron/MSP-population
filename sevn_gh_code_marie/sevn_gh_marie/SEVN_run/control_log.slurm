#!/bin/bash
#SBATCH --job-name=controls
#SBATCH --partition=Long
#SBATCH --time infinite
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --output=logs/log-%j.out
#SBATCH --error=logs/log-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=matteo.sautron@astro.unistra.fr

while true; do
	for file in logs/*.out;do
		[ -e "$file" ] || continue
		> "$file"
	done
	sleep 1800
done
