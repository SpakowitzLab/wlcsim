## Chromatin Simulator

This code is designed to simulation chromatin as described in our manuscript submitted to Nucleic Acids Research in 2019.  It built upon the code described in “Bottom-up modeling of chromatin segregation due to epigenetic modifications,” PNAS 2018.
The code is a joint effort of the members of the Spakowitz lab with particular contributions made by Quinn MacPherson for this version.

## Downloading

Clone from github repository SpakowitzLab/wlcsim or download using github’s download feature.  Make sure you clone and checkout the branch ``MS2019_NAR``

## Requirements

Running this code with require phyton3 for the Gillespie algorithm and initial LEF setup.  A few standard python packages will be required such as ``numpy``, ``random``, ``bisect``, ``contextlib``, and ``sys``.  These can be installed with ``pip`` if not already available on your system.
To run the fortran code you will need ``gfortran`` version 2003 or newer and ``make``.
You will need a few gigabites of disk space for output data and a few weeks of cpu time (depending on the speed of your processor).  RAM requirements are minimal.  
Built in visualization (optional) requires the freely available software package ``pymol``.

## To Run

Run the following script (you may need to change its permissions to executable if it is not already on your system).
``$ ./runwlcsim.sh``

## Results

The results will be saved in ``data/`` and are numbered in time order.  ``data/r110`` contains the 110th positions values of nucleosome beads with a single line corresponding to a single bead.  Units are in 28.7nm, the discretization length of the simulation.
Results are saved as time goes so a structure which may not be fully equilibrated will be generated in less time.  For reasonable results wait for save points greater than about 40.

## Changing the Settings

The manuscript submitted to NAR also refers to a simulation in which 150 beads are randomly attached to the boundary.  To run this simulation instead copy the relevant input file shown below before running the program.
``$ cp input/example_defines/defines_only150.inc src/defines.inc``

## Visualization 

We have included a visualization tool which requires ``pymol``.  To visualize the final save point:
``$ cd visualization``
``$ python3 chromosomesLoc.py``

