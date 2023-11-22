Finite difference 1D Steady Heat diffusion code written in C++ using openMP for parallelization.
The code has a serial segment as well as a parallel segment which is implemented using openMP, an "output.dat" file is created after you run the program. 
First column in the output.dat file gives you the domain points (eg i=1,2,3....LX), Second column gives you the temperature that is evaluated parallely and finally, the third column gives you the temperature evaluated serially. The last column gives you the difference between two evaluated temperatures (serially and parallely).

Boundary conditions can be set below the comment in the program given as "//initialization".
The results can be viewed using gnuplot or any other softwares of users choice.


