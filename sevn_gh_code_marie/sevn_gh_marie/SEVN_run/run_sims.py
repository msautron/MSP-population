import os

for i in range(1,101):
    os.system(f'sbatch ms_pulsars.slurm {i}')
