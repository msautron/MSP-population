#!/bin/sh
start_time=$(date +%s)
#Npulsars_simulated=$1
#bash /home/matteo.sautron/Documents/SEVN/sevn/SEVN_run/compile_be.sh
#/home/matteo.sautron/Documents/SEVN/sevn/SEVN_run/birth_and_evol $Npulsars_simulated
nvcc -g -O0 -rdc=true -fmad=false -Xcompiler -fno-fast-math --ptxas-options=-v main.cu get_birth_info.cu initialize.cu detection.cu galac_pot.cu ymw16.a -o msPop -lm -lgsl -lgslcblas -lc -lc -lcudadevrt
./msPop
sed -i 's/-nan/0.000000e+00/g' wr.txt
sed -i 's/nan/0.000000e+00/g' wint.txt
sed -i 's/-nan/0.000000e+00/g' xi_rho_data.txt
python3 plots.py
end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))
hours=$(( elapsed_time / 3600 ))
minutes=$(( (elapsed_time % 3600) / 60 ))
seconds=$(( elapsed_time % 60 ))
echo "Elapsed time : $(printf "%02d" $hours):$(printf "%02d" $minutes):$(printf "%02d" $seconds)"
