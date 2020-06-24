"""
Compute results related to (twistable, stretchable, shearable) wormlike chains.

This package can be largely separated into

1. Modules for processing simulation output from our FORTRAN codebase
(:mod:`wlcsim.data`, :mod:`wlcsim.input`)
2. Modules for generating polymer chains at equilibrium (both directly in
:mod:`wlcsim.chains` and using Monte Carlo :mod:`wlcsim.mc`).
3. Modules for simulating diffusing polymers (:mod:`wlcsim.bd`).
4. Modules for plotting these various chains and
their properties (:mod:`wlcsim.plot` and :mod:`wlcsim.plot_wlcsim`).

"""
from . import bd
from . import mc
from . import chains
from . import utils
from . import analytical
from . import special

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
