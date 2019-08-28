from bruno_util.runge_kutta import *

from numba import jit
import numpy as np

from pathlib import Path

def rouse(N, L, b, D, t, x0=None):
    r"""Simulate a Rouse polymer made of N beads free in solution.

    The polymer beads can be described by the discrete Rouse equations

    .. math::
        \xi \frac{dx(i, t)}{dt} = - k (x(i, t) - x(i+1, t)) - k (x(i, t) - x(i-1, t)) + R(t)


    where :math:`\xi = k_BT/D`, :math:`k = 3k_BT/b^2`, :math:`b` is the Kuhn
    length of the polymer, :math:`D` is the self-diffusion coefficient of a
    bead, and :math:`R(t)/\xi` is a delta-correlated stationary Gaussian
    process with mean zero and :math:`\langle R(t) R(t')/\xi^2 \rangle =
    2DI\delta(t-t')`.

    Notice that in practice, :math:`k/\xi = 3D/b^2`, so we do not need to
    include mass units (i.e. there's no dependence on :math:`k_BT`).
    """
    k_over_xi = 3*D/b**2
    def rouse_f(x, t):
        dx = np.diff(x, axis=0)
        zero = np.array([[0,0,0]])
        return -k_over_xi*(np.concatenate([zero, dx]) + np.concatenate([dx, zero]))
    if x0 is None:
        L0 = L/(N-1) # length per bead
        Lb = L0/b # kuhn lengths per bead
        x0 = b*np.sqrt(Lb)*np.cumsum(np.random.normal(size=(N,3)), axis=0)
    X = rk4_thermal_lena(rouse_f, D, t, x0)
    return X

def fast_rouse(N, L, b, D, t):
    return jit_rouse(N, L, b, D, t)
