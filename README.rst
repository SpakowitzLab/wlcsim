WormLike Chain SIMulator
========================

.. image:: https://travis-ci.com/SpakowitzLab/wlcsim.svg?branch=master
    :target: https://travis-ci.com/SpakowitzLab/wlcsim

.. image:: https://readthedocs.org/projects/wlcsim/badge/?version=latest
    :target: https://wlcsim.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

What is *wlcsim*?
-------------------

A project started by the Spakowitz lab for carrying out various polymer physics
calculations, especially multi-scale, coarse-grained simulation and theory
relating to semiflexible polymers. The library has been applied largely to
simulate DNA, and to compare results from polymer field theory to measurements
of biological polymer systems, but our `universal coarse graining procedure
<https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.013304>`_ and our
`field <https://pubs.acs.org/doi/abs/10.1021/acs.macromol.5b02639>`_ `theoretic
<https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.067802>`_
`results <https://pubs.acs.org/doi/abs/10.1021/acsmacrolett.7b00638>`_ should be
broadly applicable to any semiflexible polymer system.

For example, combining our coarse graining procedure with `field theoretic
simulations
<https://pubs.rsc.org/ko/content/articlelanding/2017/sm/c7sm00164a/unauth#!divAbstract>`_,
we were `able to simulate <https://www.pnas.org/content/115/50/12739>`_ the
phase segregation of an entire chromosome due to H3K9 methylation.


There are two largely independent codebases that are each called *wlcsim*.
One is a Fortran program that implements our universal coarse graining procedure
to allow Monte Carlo and Brownian dynamics simulations of semiflexible polymers
with high discretization lengths (i.e. with a very small number of beads). The
Monte Carlo routines in this codebase are highly optimized, and integrated with
our field theoretic simulations. For details see :ref:`wlcsimf`.

There is an associated Python package, :mod:`wlcsim`, which can be used to help
process the output of our Fortran code, but also contains much easier to use
Monte Carlo and Brownian dynamics routines. For more details on this package,
see the :ref:`wlcsimpy` docs.

The remainder of this section focusus on the FORTRAN code. The features of this
codebase are described in :ref:`features`.

Setting up a Fortran simulation
-------------------------------

To define the system you would like to simulate, set the approparte values
``src/defines.inc``.  Discriptions of each parameter are found along with their
definitions in ``src/defines.inc``.  In pracatice, the best approach is often to
start from examples provided in ``input/example_defines/``.

Some parameters that MUST be set are given default values that prevent the code
from compiling. This is on purpose so that the code is not accidentally run with
something arbitrary for these values (like the length of the chain, the
persistence length, etc).

For tips on setting up and running simulations see :ref:`tips`.

Running the Fortran code
------------------------

Simply typing ``make`` in the top level directory will build the simulator from
source. The executable created (``wlcsim.exe``) will data from the ``input/``
directoyr and write its output to the ``data`` directory.  To force a rerun
without having to manually delete all the old output files, you can also simply
type ``make run`` at any time.

By default, specifying multiple polymers just simulates them in parallel in the
same reaction volume, no interactions are assumed.

To scan parameters, the Python script ``scan_wlcsim.py`` should be used. It
takes care of saving the current git commit\_hash, all inputs, etc. into a
unique directory, and preventing race conditions even on shared filesystems,
among other things.

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
