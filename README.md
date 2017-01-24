# WormLike Chain SIMulator

[![Build Status](https://travis-ci.org/brunobeltran/wlcsim.svg?branch=master)](https://travis-ci.org/brunobeltran/wlcsim)

This code is designed to efficiently simulate the wormlike chain polymer model
using various coarse-grainings where applicable.

For very stiff polymers, the usual wormlike chain is simulated.

For relatively more flexible polymers, the "stretchable, shearable" chain is
used.

For *VERY* stretchable polymers, a purely Gaussian chain is used.

## To Run

Simply typing ``make`` in the top level directory will build the simulator
from source. The executable created (``wlcsim.exe``) will use the parameters in
the file ``input/input`` and write its output to the ``data`` directory.
Descriptions of the available parameters can be found at their definitions and
where they are read in ``src/wlcsim/params.f03``. Example input files are
usually more useful, and can be found in the ``input`` directory.
To force a rerun without having to manually delete all the old output files, you
can also simply type ``make run`` at any time.

There are several ways to easily visualize simulation output. There are PyMol
scripts in the ``vizualization`` directory, ``python -m wlcsim.plot_wlcsim``
from the repo's top level directory will launch a GUI designed to visualize BD
simulations, and one can of course simply use the output in the ``data``
directory, which contains rank two arrays of shape
``num_beads*num_polymers-by-3``, with one file per time point. By default,
specifying multiple polymers just simulates them in parallel in the same
reaction volume, no interactions are assumed.

To scan parameters, the Python script ``scan_wlcsim.py`` should be used. It takes
care of saving the current git commit\_hash, all inputs, etc. into a unique
directory, and preventing race conditions even on shared filesystems, among
other things.

## Disclaimer

This codebase is internal to the Spakowitz lab and is not guaranteed to be
bug-free at any point. For battle-tested versions of our software, please see
the links in the relevant papers.
