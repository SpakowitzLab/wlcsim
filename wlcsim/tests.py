"""
wlcsim.tests contains functions for validating the output of wlcsim.

For now, most of these are just wrappers to the functions used in Quinn and
Shifan's simulations (SpakowitzLab/MCparallelMPI/utility/*)
"""

from pymatbridge import Matlab
import os
import .input as winput
import .data as wdata

def test_init_end_to_end_distance(run_dir):
    """
    Plot the end-to-end distance distribution for a set of runs' r0 versus the
    theoretical distribution.
    """
    if not os.path.isdir(os.path.join(run_dir, 'data', 'end_to_end_r0')):
        wdata.aggregate_r0_group_NP_L0(run_dir)
    m = Matlab()
    m.start()
    m.run_code("help pcalc")
    m.stop()

