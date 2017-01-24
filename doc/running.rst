.. _running:

Running Simulations
###################

It is recommended that you always use the scan_wlcsim.py script to run your
simulations. This will make it easy to leverage the :py:mod:`wlcsim`
functionality in your analysis, and keep track of lots of easy-to-forget things
for you, like saving the version of the code used, the parameters used, random
seeds, etc.

To use scan_wlcsim.py, simply fill out the input file, input/input, with the
default parameters you want to use, then follow the instructions in
scan_wlcsim.py to fill out the parameters that you want to scan, how many
repeats of each parameter, and what name you want to give to the simulation run.


