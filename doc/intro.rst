.. _intro:

Introducing the WormLike Chain SIMulator
########################################

This code is designed to efficiently simulate the wormlike chain polymer model
using various coarse-grainings where applicable to decrease run time.

For very stiff polymers, the usual wormlike chain is simulated.

For relatively more flexible polymers, the "stretchable, shearable" chain is
used.

For *VERY* flexible polymers, a purely Gaussian chain is used.

.. _quickstart:

Quickstart
----------

Simply typing ``make`` in the top level directory will build the simulator
from source. The executable created (``wlcsim.exe``) will use the parameters in
the file ``input/input`` and write its output to the ``data`` directory.
To force a rerun without having to manually delete all the old output files, you
can also simply type ``make run`` at any time.

Several tools for vizualizing the output are included. For simply plotting
polymer configurations, see the Python modules described in :doc:`plotting`.
Plotting can also be done by hand using the output in the ``data`` directory,
which contains rank two arrays of shape ``num_beads*num_polymers-by-3``, with
one file per time point. Specifying multiple polymers just simulates them in
parallel in the same reaction volume, no interactions are assumed. For more
complex analysis, see :doc:`analyzing`.

To scan parameters, the Python script ``scan_wlcsim.py`` should be used. It takes
care of saving the current git commit\_hash, all inputs, etc. into a unique
directory, preventing race conditions even on shared filesystems. For more
details, see :doc:`running`.

.. _disclaimer:

Disclaimer
----------

This codebase is internal to the Spakowitz lab and is not guaranteed to be
bug-free at any point. For battle-tested versions of our software, please see
the links in the relevant papers.

