"""
Routines for setting the initial locations of the beads for BD.

Functions here are largely used in order to make sure that the initial
configuration that is used to perform BD is valid.

.. note::

    Writing a correct ``init_*`` function can be arbitrarily difficult due to
    compounding factors such as confinement, bead connectivity, inter-polymer
    interactions, and so forth. If your system is complicated enough, consider
    using an existing function here and simply adding your energies to the
    Monte-Carlo code in order to generate approximately equilibrium initial
    configurations using MC.

This file also contains routines to directly initialize some simple types of
chains extremely near (or at) equilibrium to reduce the number of
initialization BD or MC steps required to equilibrate the polymer.
Alternatively, some of these chains are as far as possible from equilibrium
(such as `init_straight_x`) in order to ensure that e.g. the end-to-end
distance reaching its equilibrium distribution is representative of the time it
takes for the chain itself to reach equilibrium.


.. note::

    If you know the analog of the "Rouse" time for your system (i.e. the
    time-scale of the largest mode), this is typically an easier and  more
    principled way to determine how long to run initialization BD than doing a
    huge number of simulations and trying to measure how long it takes for
    something like the end-to-end distribution to equilibrate.
"""
import numpy as np
from numba import jit


def init_straight_x(N, L):
    """Return straight line from 0 along positive x-axis."""
    x0 = np.zeros((N, 3))
    x0[:, 0] = np.linspace(0, L, N)
    return x0


def init_circle(N, L):
    """Return a circle with circumference L in the x-y plane."""
    # want circumference L
    R = L/(2*np.pi)
    r = np.zeros((N, 3))
    r[:, 0] = R*np.cos(np.linspace(0, 2*np.pi, N))
    r[:, 1] = R*np.sin(np.linspace(0, 2*np.pi, N))
    u = np.zeros((N, 3))
    u[:, 0] = -np.sin(np.linspace(0, 2*np.pi, N))
    u[:, 1] = np.cos(np.linspace(0, 2*np.pi, N))
    return r, u


@jit(nopython=True)
def init_linear_rouse_conf(N, bhat, rx, ry, rz):
    """Ad-hoc method for Rouse polymers in elliptical confinement.

    Will not necessarily be exactly at equilibrium, but pretty close.

    Parameters
    ----------
    bhat : float
        b*np.sqrt(Lb), the std dev of the distance between beads
    rx, ry, rz : float
        Principle axes of confinement ellipse

    Returns
    -------
    (N, 3) array of float
        The initial locations of the beads.
    """
    x0 = np.zeros((N, 3))
    # x0 = np.cumsum(x0, axis=0)
    for i in range(1, N):
        # 1/sqrt(3) since generating per-coordinate
        x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
        while x0[i, 0]**2/rx**2 + x0[i, 1]**2/ry**2 + x0[i, 2]**2/rz**2 > 1:
            x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
    return x0


@jit(nopython=True)
def init_homolog_rouse_conf(N, N_tot, loop_list, bhat, rx, ry, rz):
    """VERY ad-hoc, probably not near equilibrium, or even in confinement."""
    # first initialize the first chain
    x0 = np.zeros((N_tot, 3))
    x0[:N, :] = init_linear_rouse_conf(N, bhat, rx, ry, rz)
    # now initialize the loops
    num_loops, _ = loop_list.shape
    for j in range(1, num_loops-1):
        k1l, k1r, k2l, k2r = loop_list[j]
        # need to initialize beads [k2l,..,k2r] so that they form a Brownian
        # bridge between the (already determined) positions of k1l and k1r
        # first, build a regular brownian walk starting at k1l
        num_beads = k1r - k1l - 1
        for i in range(num_beads):
            i1 = k2l + i - 1 if i > 0 else k1l
            i2 = k2l + i
            x0[i2] = x0[i1] + bhat/np.sqrt(3)*np.random.randn(3)
            while 1 < x0[i2, 0]**2/rx**2 + x0[i2, 1]**2/ry**2 \
                    + x0[i2, 2]**2/rz**2:
                x0[i2] = x0[i1] + bhat/np.sqrt(3)*np.random.randn(3)
        # then subtract off the correct amount from each step to make the
        # brownian bridge. This guy is no longer guaranteed to be in the
        # confinement....but prevents being catastrophically far outside
        bridge = x0[k1r] - x0[k1l]  # target end-to-end vector
        # center bridge temporarily about zero
        for jj in range(k2l, k2r+1):
            x0[jj] -= x0[k1l]
        # jit >:( x0[k2l:k2r+1] = x0[k2l:k2r+1] - x0[k1l][None,:]
        # now fix end-to-end vector
        for i in range(num_beads):
            # fraction of the way along the bridge
            f = (i + 1)/(num_beads + 1)
            x0[k2l+i] = x0[k2l+i]*(1 - f) + f*bridge
        # re-center bridge around x0[k1l]
        for jj in range(k2l, k2r+1):
            x0[jj] += x0[k1l]
        # x0[k2l:k2r+1] = x0[k2l:k2r+1] + x0[k1l][None,:]
    # finally, initialize the free ends
    if num_loops == 1:
        x0[N:, :] = init_linear_rouse_conf(N, bhat, rx, ry, rz)
        return x0
    # "left" free end must be built backwards from first connection point
    k1l, k1r, k2l, k2r = loop_list[0]
    num_beads = k1r - k1l - 1  # or k2r - k2l + 1
    for i in range(num_beads):
        i1 = k2r - i + 1 if i > 0 else k1r
        i2 = k2r - i
        x0[i2] = x0[i1] + bhat/np.sqrt(3)*np.random.randn(3)
        while x0[i2, 0]**2/rx**2 + x0[i2, 1]**2/ry**2 + x0[i2, 2]**2/rz**2 > 1:
            x0[i2] = x0[i1] + bhat/np.sqrt(3)*np.random.randn(3)
    # right free end can be built as normal with no bridge stuff
    k1l, k1r, k2l, k2r = loop_list[-1]
    num_beads = k1r - k1l - 1  # or k2r - k2l + 1
    for i in range(num_beads):
        i1 = k2l + i - 1 if i > 0 else k1l
        i2 = k2l + i
        x0[i2] = x0[i1] + bhat/np.sqrt(3)*np.random.randn(3)
        while x0[i2, 0]**2/rx**2 + x0[i2, 1]**2/ry**2 + x0[i2, 2]**2/rz**2 > 1:
            x0[i2] = x0[i1] + bhat/np.sqrt(3)*np.random.randn(3)
    return x0
