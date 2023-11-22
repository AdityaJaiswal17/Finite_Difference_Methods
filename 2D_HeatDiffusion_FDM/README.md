Finite difference 2D Steady Heat diffusion code written in C++ using openMP for parallelization. 
The code has a serial segment as well as a parallel segment which is implemented using openMP, an "output.dat" file is created after you run the program. 
First column in the output.dat file gives you the domain points in X-axis(eg i=1,2,3....LX), Second column gives you the domain points as well but in Y-axis, Third column gives youtemperature that is evaluated computationally.
Two output files are created that are "output_parallel.dat" and "output_serial.dat".

Boundary conditions can be set below the comment in the program given as "//initialization" (line 55 & line 163). 
The results can be viewed using gnuplot or any other softwares of users choice.

Syntax to view a surface map on gnuplot;

set pm3d map;
set size ratio 1.0; //(ratio as per needed)
sp 'output_serial.dat' u 1:2:3 with image;


To compile use the following command;

export OMP_NUM_THREADS= ("value of the threads you can set according to your machine specifications")
; g++ 2D_heat_diff.cpp -fopenmp

NOTE : fopenmp flag is just to compile the openMP part in the code.
