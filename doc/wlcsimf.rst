.. _wlcsimf:

wlcsim Fortran Codebase
#######################

.. f:program:: wlcsim

    Use a universal discretization scheme to simulate from WLCs through
    Gaussian chains.

How to use
==========

The parameters of the simulation are specified at compile time in the large
`src/defines.inc` file. The comments in this file should describe the parameters
well enough for one to piece together how to use them. Any parameters that
aren't relevant to a given simulation can typically be ignored and left to their
default values.

Some parameters that MUST be set are given default values that prevent the code
from compiling. This is on purpose so that the code is not accidentally run with
something arbitrary for these values (like the length of the chain, the
persistence length, etc).

To run the codebase, simply specify the parameters in `src/defines.inc` then use
`make run` to build and run the code with the given parameters. The code will
write its output to the `data` folder in the current working directory. To run
the code many times (to get statistics for a given set of parameters, or to scan
over many parameters) use the script `scan_wlcsim.py` to make one new directory
per run of the code, and automatically save a bunch of useful information about
each run.

Structure of Codebase
=====================

The main entry point of the program is pretty useless as a starting point to
understand the codebase, and simply calls one of several "versions" of the
program built from our Brownian Dynamics/Monte Carlo API.

.. f:autosrcfile:: wlcsim.f03

The actual hard work is done in the `wlcsim_*` files

.. f:autosrcfile:: wlcsim_bruno.f03

.. f:autosrcfile:: wlcsim_quinn.f03

The main entry points for our API are the `MCsim.f03` and `BDsim.f03` routines,
which allow you to do Monte Carlo moves and Brownian Dynamics timesteps
(respectively) given a set of energies (or forces, respectively) specified by
the input in `defines.inc`.

`MCsim.f03` calls different subroutines that define possible Monte Carlo moves
(like `MC_reptationMove.f03`) and then checks for whether the move should
succeed by summing the various energies that are turned on (like `MC_wlc.f90`).

.. f:autosrcfile:: MC_reptationMove.f03

.. f:autosrcfile:: MC_wlc.f90

.. f:autosrcfile:: adapt.f90

