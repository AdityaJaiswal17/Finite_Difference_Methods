#!/bin/bash

rm -rf serial_output
g++ -O3 lid_driven_cavity_serial.cpp -lm;
./a.out;
gnuplot -e "set pm3d map; set size ratio 1.0; sp 'output_serial.dat' u 1:2:7 w image; pause -1"
mkdir serial_output
mv output_serial.dat U-Re_x=0.5_serial.dat V-Re_x=0.5_serial.dat ~/Desktop/CFD/FVM_FDM/Lid_driven_cavity/serial_output





