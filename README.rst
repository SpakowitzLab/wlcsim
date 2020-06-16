WormLike Chain SIMulator
========================

.. image:: https://travis-ci.com/SpakowitzLab/wlcsim.svg?branch=master
    :target: https://travis-ci.com/SpakowitzLab/wlcsim

.. image:: https://readthedocs.org/projects/wlcsim/badge/?version=latest
    :target: https://wlcsim.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

What wlcsim is
--------------

This is a project started by the Spakowitz lab to make varioius polymer physics
simulations / calulations.  The bulk of the codebalse is built around a FORTRAN
code that is designed to efficiently simulate various systems of wormlike
chain polymer(s) using various coarse-grainings where applicable.  This code
performs Monte Carlo or Brownian dynamics simulations.

The wlcsim projects also contains various related code in Python  These include
a more modern Browning dynamics simulation as well as a number of analytical
calculations.  For more details on these see :ref:`wlcsimpy`.

The remainder of this section focusus on the FORTRAN code.
The features of this codebase are described in :ref:`features`.

More details on how the code is structured see :ref:`wlcsimf`

Setting up a simulation
-----------------------

To define the system you would like to simulate, set the approparte values
``src\defines.inc``.  Discriptions of each parameter are found along with their
definitions in ``src\defines.inc``.  In pracatice, the best approach is often to
start from examples provided in ``input/example_defines/``.

Some parameters that MUST be set are given default values that prevent the code
from compiling. This is on purpose so that the code is not accidentally run with
something arbitrary for these values (like the length of the chain, the
persistence length, etc).

For tips on setting up and running simulations see :ref:`tips`.

To Run
------


Simply typing ``make`` in the top level directory will build the simulator
from source. The executable created (``wlcsim.exe``) will data from the ``input/``
directoyr and write its output to the ``data`` directory.
To force a rerun without having to manually delete all the old output files, you
can also simply type ``make run`` at any time.

By default, specifying multiple polymers just simulates them in parallel in the same
reaction volume, no interactions are assumed.

To scan parameters, the Python script ``scan_wlcsim.py`` should be used. It takes
care of saving the current git commit\_hash, all inputs, etc. into a unique
directory, and preventing race conditions even on shared filesystems, among
other things.

To perform parallel tempering using MPI for multiprocessing using 10 threads
first compile using ``make`` then type ``mpirun -np 10 wlcsim.exe``.  For more
details on parallel tempering see :ref:`parallel_temp`.

Output
------

There are several ways to easily visualize simulation output. There are PyMol
scripts in the ``vizualization`` directory, ``python -m wlcsim.plot_wlcsim``
from the repo's top level directory will launch a GUI designed to visualize BD
simulations, and one can of course simply use the output in the ``data``
directory, which contains rank two arrays of shape
``num_beads*num_polymers-by-3``, with one file per time point.

For more details see :ref:`output`.


Disclaimer
----------

This codebase is internal to the Spakowitz lab and is not guaranteed to be
bug-free at any point. For battle-tested versions of our software, please see
the links in the relevant papers.
