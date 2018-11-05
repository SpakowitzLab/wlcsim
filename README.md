# WormLike Chain SIMulator

This directory contains the code to run the simulations referred to in “Bottom-up modeling of chromatin segregation due to epigenetic modifications.”  The included methylation input sequence used to generate Figures 1,2,3,5. 

Direct questions to:
Quinn MacPherson, qmac@stanford.edu or quinnmacp@gmail.com

Software Requirements:
Open MPI with mpifort compiler
Make

To compile:
$ make

To run a multi-thread job (15 different HP1 concentrations):
$ ./aMPIrunwlcsim.sh

To change settings edit src/defines.inc, to change input methylation sequence
edit input/meth, to edit code edit src/

Output is found in data/ and is of the form r110v9 where 110 refers to save point
and 9 refers to which HP1 chemical potential.  Fromat is of the form
x y z

