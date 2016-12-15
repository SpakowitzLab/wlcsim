WormLike Chain SIMulator
========================

This code is designed to efficiently simulate the wormlike chain polymer model
using various coarse-grainings where applicable.

For very stiff polymers, the usual wormlike chain is simulated.

For relatively more flexible polymers, the "stretchable, shearable" chain is
used.

For *VERY* stretchable polymers, a purely Gaussian chain is used.

To Run
------

Simply typing ``make`` in the top level directory will build the simulator
from source. The executable created (``wlcsim.exe``) will use the parameters in
the file ``input/input`` and write its output to the ``data`` directory.
To force a rerun without having to manually delete all the old output files, you
can also simply type ``make run`` at any time.

The output can be vizualized using the PyMol scripts in the ``vizualization``
directory or by hand using the output in the ``data`` directory, which contains
rank two arrays of shape ``num_beads*num_polymers-by-3``, with one file per time
point. Specifying multiple polymers just simulates them in parallel in the same
reaction volume, no interactions are assumed.

To scan parameters, the Python script ``scan_wlcsim.py`` should be used. It takes
care of saving the current git commit\_hash, all inputs, etc. into a unique
directory, preventing race conditions even on shared filesystems.

Disclaimer
----------

This codebase is internal to the Spakowitz lab and is not guaranteed to be
bug-free at any point. For battle-tested versions of our software, please see
the links in the relevant papers.
