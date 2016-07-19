# WormLike Chain SIMulator

This code is designed to efficiently simulate the wormlike chain polymer model
using various coarse-grainings where applicable.

For very stiff polymers, the usual wormlike chain is simulated.

For relatively more flexible polymers, the "stretchable, shearable" chain is
used.

For *VERY* stretchable polymers, a purely Gaussian chain is used.

## To Run

Simply typing `make run` in the top level directory will build the simulator
from source and run it with the parameters in the file `input/input`.

The output can be vizualized using the PyMol scripts in the `vizualization`
directory or by hand using the output in the `data` directory, which contains
rank two arrays of shape `num_beads*num_polymers-by-3`, with one file per time
point. Specifying multiple polymers just simulates them in parallel in the same
reaction volume, no interactions are assumed.

To scan parameters, the Perl scripts `run_parameter.pl` or
`run_parameter_series.pl` can be used. Care must be taken to match the
parameters in these scripts with the parameters of the simulation. To run these
in parallel over many cores, simply call them via the `run-parallel.sh` script.

## Disclaimer

This codebase is internal to the Spakowitz lab and is not guaranteed to be
bug-free at any point. For battle-tested versions of our software, please see
the links in the relevant papers.
