"""Various routines that are too slow to run in real time are tabulated here.

Especially of note are Lena's ssWLC parameters code, which nobody has access to
anymore to the best of my knowledge (except Andy, who doesn't know how to use
it) and Tom's Rouse velocity cross-correlation code, which required an entire
rotation student to even run back when it did exist. I think it's saved under
some folder called ParConV4, but I haven't put in the time to try to figure out
if it's usable."""
import numpy as np

from pathlib import Path

# precomupted velocity cross-correlation for rouse polymer
# from Lampo et al, was pre-computed on a grid...
deltas_ = np.linspace(-3, 3, 25) # placeholders for corresponding values in logspace
alphas_ = np.linspace(0.25, 1, 31) # steps of 0.025
tOverDeltas_ = np.linspace(0, 5, 501) # steps of 0.01
vvcf_table_ = np.reshape(np.loadtxt(Path(__file__).parent / Path('vvcf_table.csv'),
        delimiter=','), (31, 25, 501))
def calc_vel_corr_fixed_(tOverDelta, deltaOverTDeltaN, alpha):
    # this performs interpolation in logspace for "delta"/deltaOverTDeltaN
    deltaOverTDeltaN = np.log10(deltaOverTDeltaN)
    return scipy.interpolate.interpn((alphas_, deltas_, tOverDeltas_),
                                     vvcf_table_,
                                     (alpha, deltaOverTDeltaN, tOverDelta))
calc_vel_corr_fixed_.vvcf = None
frac_vel_corr = np.vectorize(calc_vel_corr_fixed_)

