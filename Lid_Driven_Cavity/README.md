Steady state Lid driven Cavity finite difference code.
To run, run the script files "run_serial.sh" or "run_parallel.sh" note that you will have to change "OMP_NUM_THREADS" values in the parallel bash script according your machine specifications.
A gnuplot surface plot is also created automatically after you run the script.

Parameters (U_inf is the lid velocity, length, KinVisc is the kinetic viscosity) can be changed in the code accordingly to adjust the reynolds number value.

A validation graph is also presented in which simulation results are compared with the experimental values obtained at RE=100 which was performed by Ghia-Ghia et al.
