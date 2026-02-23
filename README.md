# User Guide

This github contains:
- A code allowing to produce a MSP population.
- A code allowing to produce pulsars obtained through AIC of WD.

SEVN (N-body code) is in this repository, however M.Sautron modified some files, especially concerning the evolution of NS. 

## Step 1 - Generate and evolve the population 

You will need to use an HPC to have this code run in less than 24 hours. Go to the folder `sevn_gh/SEVN_run` and compile SEVN:
```bash
bash compile.sh
```
Once SEVN is compiled, you need to compile the C code allowing to generate each binary within the listbin.dat files: 
```bash
bash compile_be.sh
```
Check the file `ms_pulsars.slurm` to check how many binaries you want to simulate, modify the number after `./birth_and_evol` if you want to modify the number of simulated binaries. Then, to choose the number of different simulation that you are doing, go to `run_sims.py`. 

In the files `birth_evolution_pop.c` and `birth_pulsars.c` you can modify the initial distribution related to the pulsars and the binary parameters. `run.sh` and `run2.sh` allowed to run an instance of sevn, for companion with a mass > 2M and < 2M respectively. The code related to the evolution of the NS is in `sevn_gh/src/star/remnant.cpp`
and the code related to the binary evolution (that I modified) is in `sevn_gh/src/binary/Processes.cpp`. Especially, the irradiation from the pulsar onto the companion is in the function related to the RLO evolution (evaluated at each step !).  

Once your simulations are finished you will obtain 5 different types of files: `data_pop_msXX.txt`, `data_beta_initXX.txt`, `data_betaXX.txt`, `data_submsXX.txt` and `data_taudXX.txt`. Use the command: 
```bash
cat data_pop_ms*.txt > data_pop_ms_all.txt 
```
in order to put all of the data in one single file for each type of data (do it for data_beta_init, data_beta etc...). You can change every **inf** in the generated `data_pop_ms_all.txt` file with this command: 
```bash
sed -i 's/inf/1.000000e+55/g' data_pop_ms_all.txt
```
Otherwise, you can check any issues with this command (display lines with more than 20 columns):
```bash
awk '{if (NF > 20) print NR ": " $0}' data_pop_ms_all.txt 
```
You are now all set to simulate the detection ! 

## Step 2 - Detection and plots

You will need a GPU for this part. Go to the folder `ms_pulsars_SEVN_evol`. You will need to go to the drive to get the temperature map of the sky obtained by Remazeilles 2014 at 408 MHz: https://drive.google.com/drive/folders/1g5Z820uUpipQSCmaCKylrhwtkVzcgv_w?usp=sharing 

All of your `data*all.txt` put them in this folder. Then you can execute: 
```bash
bash run_ms_pop.sh 
```
It allows to run the detection pipeline (radio + gamma), evolving spatially the pulsars in the Galaxy. Then it makes the plot thanks to the file `plots.py`. 
`get_birth_info.cu` -> Allows to get the info from the data_pop_ms_all.txt and other files. 

`macro.h` -> Contains the declaration of most variables and lists.

`initialize.cu` -> Declare also some variables + allows to finalize the declaration of the tables.

`main.cu` -> Spatial evolution, made with GPU is here, allows to choose the order of use of each function. Free the memory at the end.

`cn.h`, `dmdtau.cu`, `dora.cu`, `fermibubble.cu`, `frb_d.cu`, `galcen.cu`, `gum.cu`, `lmc.cu`, `localbubble.cu`, `ne_crd.cu`, `nps.cu`, `smc.cu`, `spiral.cu`, `spiral.txt`, `thick.cu`, `thin.cu`, `ymw16`, `ymw16a`, `ymw16.cu`, `ymw16_ne`, `ymw16_ne.cu`, `ymw16par.cu` and `ymw16par.txt` -> files related to ISM dispersion of the width of the radio pulse profile. Adapted to GPU from Yao et al. (2017)

`galac_pot.cu` -> Contains function for initial speed.

`3PC_Catalog_20230803.xls` and `3PC_SensitivityMap_20230629.fits` -> Contains the 3PC catalog and Fermi/LAT sensitivity map

`AIC_info.txt`, `data_ATNF_fastfermipmps.txt`, `gl_gb_ATNF.txt`, `med_comp_mass_ATNF.txt`, `mpulsar_obs.txt`, `w10_atnf.txt` and `ATNF_canonical_pulsars.txt` -> First file = result of pulsars formed through AIC that became recycled. Others = ATNF data files, should be updated if used long after 02/2026 and if you add surveys. 

`detection.cu` -> Handles the detection in radio + gamma. Most important functions are there.


`facteur_omega.dat` -> File obtained from J.Pétri allowing to obtain the anisotropy factor from simulation (cf J.Pétri (2024)). 

`get_temp.py` -> Python file allowing to compute the temperature of the sky at a position in the sky. 

`plots.py` -> Makes the plots. Look here if you want to know which file is needed to obtain plots.

`sensitivity_3PC.py` -> Compute sensitivity in gamma thanks to the fits map. 

## Simulate a populatin of WD-MS Star

You will have to go to the folder `sevn_gh_code_marie`. Here do the same as step 1. Only `data_pop_ms_all.txt` matters. 

Then for step 2, modify `get_birth_info.cu` in consequence to use only `data_pop_ms_all.txt`. You will have to generate a random age for each pulsar also in `get_birth_info.cu`

