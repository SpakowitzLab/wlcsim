"""Compute results related to the (stretchable-shearable) worm-like chain

This module can be largely separated into modules for processing simulation
output from our FORTRAN codebase (:mod:`wlcsim.data`, :mod:`wlcsim.input`),
modules for generating polymer chains at equilibrium directly
(:mod:`wlcsim.chains`), modules for sampling polymer chains from equilibrium
using Monte Carlo (:mod:`wlcsim.mc`), modules for simulating diffusing polymer
chains (:mod:`wlcsim.bd`), and modules for plotting these various chains and
their properties (:mod:`wlcsim.plot` and :mod:`wlcsim.plot_wlcsim`).
"""

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from . import bd
from . import mc
from . import chains
from . import utils
from . import analytical
