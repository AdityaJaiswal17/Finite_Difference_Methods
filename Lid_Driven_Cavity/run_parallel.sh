#!/bin/bash

rm -rf parallel_output
export OMP_NUM_THREADS=4;
g++ -O3 lid_driven_cavity_parallel.cpp -fopenmp -lm;
./a.out;
gnuplot -e "set pm3d map; set size ratio 1.0; sp 'output_parallel.dat' u 1:2:7 with image; pause -1"
mkdir parallel_output
mv output_parallel.dat U-Re_x=0.5_parallel.dat V-Re_x=0.5_parallel.dat ~/Desktop/CFD/FVM_FDM/Lid_driven_cavity/parallel_output





