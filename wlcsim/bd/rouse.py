from ..plot import PolymerViewer
from bruno_util.runge_kutta import *

from numba import jit
import numpy as np

from pathlib import Path

# for testing
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt

def recommended_dt(b, D):
    r"""Recommended "dt" for use with rouse*jit family of functions.

    Currently set to :math:`\frac{1}{10}\frac{b^2}{6D}`.

    See the rouse_jit docstring for source of this time scale."""
    return (1/10)*b*b/6/D

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
        return -k_over_xi*(np.concatenate([zero, dx]) + np.concatenate([dx, zero]))
    if x0 is None:
        L0 = L/(N-1) # length per bead
        Lb = L0/b # kuhn lengths per bead
        x0 = b*np.sqrt(Lb)*np.cumsum(np.random.normal(size=(N,3)), axis=0)
    X = rk4_thermal_lena(rouse_f, D, t, x0)
    return X

def test_rouse_msd_line_approx():
    """Figure out what Dapp is exactly using our analytical result.

    if we use msd_approx(t) = 3*bhat*np.sqrt(Dhat*t)/np.sqrt(3)*1.283791(6)
    then |msd(t) - msd_approx(t)|/msd(t) = np.sqrt(2)/N*np.power(t, -1/2)
    msd_approx(t) > msd(t) in this range, so that means that
    msd_approx(t)/(np.sqrt(2)/N*np.power(t, -1/2) + 1) = msd(t)

    if we redefine msd_approx with this correction, the new relative error is
    about np.sqrt(2)/N/100*t**(-1/2)? gotta see if this carries over to other chain
    parameters though...it does not...
    this time msd(t) is bigger, so
    msd_approx = msd_approx(t)/(1 - np.sqrt(2)/N/100*t**(-1/2))

    okay actually looks like extra factor is additive?

    for N = 1e8+1; L = 25; b = 2; D = 1;, we have
    msd_approx(t) - msd(t) = 1.01321183e-07

    for N = 1e8+1; L = 174; b = 150; D = 166;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-05

    for N = 1e8+1; L = 174; b = 15; D = 166;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-06

    for N = 1e8+1; L = 17.4; b = 15; D = 166;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-07

    for N = 1e8+1; L = 17.4; b = 15; D = 16.6;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-07 (no change)

    for N = 1e7+1; L = 17.4; b = 15; D = 16.6;, we have
    msd_approx(t) - msd(t) = 5.2889657(3)e-06 (no change)

    so the answer is like
    3*bhat*np.sqrt(Dhat*t)/np.sqrt(3)*1.283791(6) - 0.202642385398*b*L/N
    """
    N = 1e8+1; L = 174; b = 150; D = 166; dt = 1e-2; Nt = 1e6; Nt_save = 1e4;
    t = np.arange(0, Nt*dt, dt); t_save = t[::int(Nt/Nt_save)];
    Nhat = L/b; L0 = L/(N-1); Dhat = D*(N)/Nhat; bhat = np.sqrt(L0*b)
    plt.plot(tmsd, np.abs(msd_rouse - 3*bhat*np.sqrt(Dhat*tmsd)/np.sqrt(3)*1.12837916)/msd_rouse, 'k')
    plt.yscale('log')
    plt.xscale('log')

@jit(nopython=True)
def jit_rouse(N, L, b, D, t, t_save=None):
    r"""faster version of wlcsim.bd.rouse.rouse using jit
    #TODO: change to input units of "D" = 3*bhat*np.sqrt(Dhat)/np.sqrt(3)*1.28379161(5)
    which is the "apparent" diffusion coefficient of the Rouse-part of the MSD.

    N=101,L=100,b=1,D=1 takes about 3.5min to run when
    t=np.linspace(0, 1e5, 1e7+1)
    adding t_save does not slow function down

    our srk1 scheme is accurate as long as :math:`\Delta t` is less than the
    transition time where the MSD goes from high k behavior (t^1) to rouse
    behavior (t^(1/2)).  This is exactly the time required to diffuse a Kuhn
    length, so we just need :math:`\Delta t < \frac{b^2}{6D}` in 3D. the
    "crossover" from fast-k to rouse-like behavior takes about one order of
    magnitude in time, but adding more than that doesn't seem to make the MSD
    any more accurate, so we suggest setting :math:`\Delta t =
    \frac{1}{10}\frac{b^2}{6D}`.

    the number of orders of magnitude of "rouse" scaling the simulation will
    capture is exactly dictated by the ratio between this time scale and the
    rouse relaxation time (so like ?N**2?)

    recall (doi & edwards, eq 4.25) that the first mode's relaxation time is
    :math:`\tau_1 = \frac{\xi N^2}{k \pi^2 }`.
    and the :math:`p`th mode is :math:`\tau_p = \tau_1/p^2` (this is the
    exponential falloff rate of the :math:`p`th mode's correlation function).
    """
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position, sqrt(3) since generating per-coordinate
    x0 = bhat/np.sqrt(3)*np.random.randn(N, 3)
    # x0 = np.cumsum(x0, axis=0)
    for i in range(1,N):
        x0[i] = x0[i-1] + x0[i]
    if t_save is None:
        t_save = t
    x = np.zeros(t_save.shape + x0.shape)
    dts = np.diff(t)
    # -1 or 1, p=1/2
    S = 2*(np.random.rand(len(t)) < 0.5) - 1
    save_i = 0
    if 0 == t_save[save_i]:
        x[0] = x0
        save_i += 1
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = dts[i-1]
        t0 = t[i-1]
        dW = np.random.randn(*x0.shape)
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S[i])
        # estimate for slope at interval start
        f = np.zeros(x0.shape)
        for j in range(1,N):
            for n in range(3):
                f[j,n] += -k_over_xi*(x0[j,n] - x0[j-1,n])
                f[j-1,n] += -k_over_xi*(x0[j-1,n] - x0[j,n])
        K1 = f + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S[i])
        # estimate for slope at interval end
        x1 = x0 + h*K1
        f = np.zeros(x0.shape)
        for j in range(1,N):
            for n in range(3):
                f[j,n] += -k_over_xi*(x1[j,n] - x1[j-1,n])
                f[j-1,n] += -k_over_xi*(x1[j-1,n] - x1[j,n])
        K2 = f + Fbrown
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x

@jit(nopython=True)
def jit_rouse_confined(N, L, b, D, Aex, rx, ry, rz, t, t_save):
    """ adds an elliptical confinement. energy is like cubed of distance
    outside of ellipsoid, pointing normally back in. times some factor Aex
    controlling its strength """
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    x0 = np.zeros((N, 3))
    # x0 = np.cumsum(x0, axis=0)
    for i in range(1,N):
        x0[i] = x0[i-1] + b*np.sqrt(Lb)*np.random.randn(3)
        while x0[i,0]**2/rx**2 + x0[i,1]**2/ry**2 + x0[i,2]**2/rz**2 > 1:
            x0[i] = x0[i-1] + b*np.sqrt(Lb)*np.random.randn(3)
    if t_save is None:
        t_save = t
    x = np.zeros(t_save.shape + x0.shape)
    dts = np.diff(t)
    # -1 or 1, p=1/2
    S = 2*(np.random.rand(len(t)) < 0.5) - 1
    save_i = 0
    if 0 == t_save[save_i]:
        x[0] = x0
        save_i += 1
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = dts[i-1]
        t0 = t[i-1]
        dW = np.random.randn(*x0.shape)
# D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S[i])
        # estimate for slope at interval start
        f = np.zeros(x0.shape)
        j = 0
        conf = x0[j,0]**2/rx**2 + x0[j,1]**2/ry**2 + x0[j,2]**2/rz**2
        if conf > 1:
            conf_u = np.array([-x0[j,0]/rx**2, -x0[j,1]/ry**2, -x0[j,2]/rz**2])
            conf_u = conf_u/np.linalg.norm(conf_u)
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3) # Steph: np.power(np.sqrt(conf) - 1, 3)
        for j in range(1,N):
            conf = x0[j,0]**2/rx**2 + x0[j,1]**2/ry**2 + x0[j,2]**2/rz**2
            if conf > 1:
                conf_u = np.array([-x0[j,0]/rx**2, -x0[j,1]/ry**2, -x0[j,2]/rz**2])
                conf_u = conf_u/np.linalg.norm(conf_u)
                f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3) # Steph: np.power(np.sqrt(conf) - 1, 3)
            for n in range(3):
                f[j,n] += -k_over_xi*(x0[j,n] - x0[j-1,n])
                f[j-1,n] += -k_over_xi*(x0[j-1,n] - x0[j,n])
        K1 = f + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S[i])
        # estimate for slope at interval end
        x1 = x0 + h*K1
        f = np.zeros(x1.shape)
        j = 0
        conf = x1[j,0]**2/rx**2 + x1[j,1]**2/ry**2 + x1[j,2]**2/rz**2
        if conf > 1:
            conf_u = np.array([-x1[j,0]/rx**2, -x1[j,1]/ry**2, -x1[j,2]/rz**2])
            conf_u = conf_u/np.linalg.norm(conf_u)
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3) # Steph: np.power(np.sqrt(conf) - 1, 3)
        for j in range(1,N):
            conf = x1[j,0]**2/rx**2 + x1[j,1]**2/ry**2 + x1[j,2]**2/rz**2
            if conf > 1:
                conf_u = np.array([-x1[j,0]/rx**2, -x1[j,1]/ry**2, -x1[j,2]/rz**2])
                conf_u = conf_u/np.linalg.norm(conf_u)
                f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3) # Steph: np.power(np.sqrt(conf) - 1, 3)
            for n in range(3):
                f[j,n] += -k_over_xi*(x1[j,n] - x1[j-1,n])
                f[j-1,n] += -k_over_xi*(x1[j-1,n] - x1[j,n])
        K2 = f + Fbrown
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x

@jit(nopython=True)
def _init_in_confinement(N, bLb, rx, ry, rz):
    """ad-hoc method to get close to equilibrium before starting

    Parameters
    ----------
    bLb : float
        b*np.sqrt(Lb), the std dev of the distance between beads
    rx, ry, rz : float
        principle axes of confinement ellipse

    """
    x0 = np.zeros((N, 3))
    # x0 = np.cumsum(x0, axis=0)
    for i in range(1,N):
        x0[i] = x0[i-1] + bLb*np.random.randn(3)
        while x0[i,0]**2/rx**2 + x0[i,1]**2/ry**2 + x0[i,2]**2/rz**2 > 1:
            x0[i] = x0[i-1] + bLb*np.random.randn(3)
    return x0

@jit(nopython=True)
def f_conf(x0, Aex, rx, ry, rz):
    """compute soft (cubic) force due to elliptical confinement"""
    N, _ = x0.shape
    f = np.zeros(x0.shape)
    for j in range(N):
        conf = x0[j,0]**2/rx**2 + x0[j,1]**2/ry**2 + x0[j,2]**2/rz**2
        if conf > 1:
            conf_u = np.array([-x0[j,0]/rx**2, -x0[j,1]/ry**2, -x0[j,2]/rz**2])
            conf_u = conf_u/np.linalg.norm(conf_u)
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3) # Steph: np.power(np.sqrt(conf) - 1, 3)
    return f

@jit(nopython=True)
def f_elas(x0, k_over_xi):
    """compute spring forces on single, linear rouse polymer"""
    N, _ = x0.shape
    f = np.zeros(x0.shape)
    for j in range(1,N):
        for n in range(3):
            f[j,n] += -k_over_xi*(x0[j,n] - x0[j-1,n])
            f[j-1,n] += -k_over_xi*(x0[j-1,n] - x0[j,n])
    return f

@jit(nopython=True)
def jit_rouse_confinement_clean(N, L, b, D, Aex, rx, ry, rz, t, t_save):
    """Unfortunately, by "cleaning up" this function, we also make it 2x slower"""
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    x0 = _init_in_confinement(N, b*np.sqrt(Lb), rx, ry, rz)
    # pre-alloc output
    if t_save is None:
        t_save = t
    x = np.zeros(t_save.shape + x0.shape)
    # setup for saving only requested time points
    save_i = 0
    if 0 == t_save[save_i]:
        x[0] = x0
        save_i += 1
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = t[i] - t[i-1]
        t0 = t[i-1]
        dW = np.random.randn(*x0.shape)
        # -1 or 1, p=1/2
        S = 2*(np.random.rand() < 0.5) - 1
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S)
        # estimate for slope at interval start
        K1 = f_conf(x0, Aex, rx, ry, rz) + f_elas(x0, k_over_xi) + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S)
        x1 = x0 + h*K1
        # estimate for slope at interval end
        K2 = f_conf(x1, Aex, rx, ry, rz) + f_elas(x1, k_over_xi) + Fbrown
        # average the slope estimates
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x

def _homolog_points_to_loops_list(N, homolog_points):
    """Create loops list for jit_rouse_linked.

    Parameters
    ----------
    N : int
        number of beads in each polymer if they were unbound
    homolog_points : (L,) array_like
        list of loci indices (beads) that are homologous paired.

    Returns
    -------
    N_tot : int
        final size of the x0 array for the simulation.
    loop_list : (L,4) np.array
        data structure used by jit_rouse_linked to correctly compute the spring
        forces in the network.

    Notes
    -----
    You can use the output like

    >>> k1l, k1r, k2l, k2r = loop_list[i]

    k1l, k1r denote the connected beads, so the array contains a "loop"
    representing the beads of the second chain from k1l+1 to k1r-1
    at array locations [k2l, k2r]. k1l==-1 and k1r==N denote the free ends.
    """
    homolog_points = np.sort(np.array(homolog_points))
    num_points = len(homolog_points)
    num_loops = num_points + 1 # includes two end non-"loops" (aka free ends)
    loop_list = np.zeros((num_loops, 4))
    if num_loops == 1:
        return np.array([[-1,N,N,2*N-1]])
    # keep track of length of "x0"
    curr_n = N
    # there's at least one connection point, thanks to num_loops==1 check above
    for i in range(0,num_points+1):
        k1l = homolog_points[i-1] if i > 0 else -1
        k1r = homolog_points[i] if i < num_points else N
        k2l = curr_n
        beads_in_loop = k1r - k1l - 1
        curr_n += beads_in_loop
        k2r = curr_n - 1
        loop_list[i,:] = np.array([k1l, k1r, k2l, k2r])
    return curr_n, loop_list.astype(int)

@jit(nopython=True)
def f_elas_homolog(x0, N, loop_list, k_over_xi):
    """compute spring forces on two linear rouse polymers hooked together
    (homologously) at beads specified by loop_list. While each polymer is of length
    N beads, the beads that are hooked together move together identically, so
    you have x0.shape[0] == 2*N - len(loop_list)"""
    f = np.zeros(x0.shape)
    num_loops, _ = loop_list.shape
    # first get forces on "first" polymer
    for j in range(1,N):
        for n in range(3):
            f[j,n] += -k_over_xi*(x0[j,n] - x0[j-1,n])
            f[j-1,n] += -k_over_xi*(x0[j-1,n] - x0[j,n])
    for i in range(num_loops):
        # k1l, k1r denote the connected beads, so the array contains a "loop"
        # representing the beads of the second chain from k1l+1 to k1r-1
        # at array locations [k2l, k2r].
        k1l, k1r, k2l, k2r = loop_list[i]
        # k1l == -1 is flag for "end" loop
        if k1l >= 0:
            f[k2l] += -k_over_xi*(x0[k2l] - x0[k1l])
        # k1r == N is flat for "end" loop
        if k1r < N:
            f[k2r] += -k_over_xi*(x0[k2r] - x0[k1r])
        N_loop = k1r - k1l + 1
        # analagous to loop above, since end-bead forces already done
        for k in range(1,N_loop):
            j = k2l+k
            f[j] += -k_over_xi*(x0[j] - x0[j-1])
            f[j-1] += -k_over_xi*(x0[j-1] - x0[j])
    return f

@jit(nopython=True)
def _init_homologs(N, N_tot, loop_list, bLb, rx, ry, rz):
    """VERY ad-hoc, probably not that near equilibrium, or even strictly in the
    confinement"""
    # first initialize the first chain
    x0 = np.zeros((N_tot,3))
    x0[:N,:] = _init_in_confinement(N, bLb, rx, ry, rz)
    # now initialize the loops
    num_loops, _ = loop_list.shape
    for j in range(1,num_loops-1):
        k1l, k1r, k2l, k2r = loop_list[j]
        # need to initialize beads [k2l,..,k2r] so that they form a Brownian
        # bridge between the (already determined) positions of k1l and k1r
        # first, build a regular brownian walk starting at k1l
        num_beads = k1r - k1l - 1
        for i in range(num_beads):
            i1 = k2l + i - 1 if i > 0 else k1l
            i2 = k2l + i
            x0[i2] = x0[i1] + bLb*np.random.randn(3)
            while x0[i2,0]**2/rx**2 + x0[i2,1]**2/ry**2 + x0[i2,2]**2/rz**2 > 1:
                x0[i2] = x0[i1] + bLb*np.random.randn(3)
        # then subtract off the correct amount from each step to make the
        # brownian bridge. This guy is no longer guaranteed to be in the
        # confinement....but prevents being catastrophically far outside
        bridge = x0[k1r] - x0[k1l] # target end-to-end vector
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
        x0[N:,:] = _init_in_confinement(N, bLb, rx, ry, rz)
    # "left" free end must be built backwards from first connection point
    k1l, k1r, k2l, k2r = loop_list[0]
    num_beads = k1r - k1l - 1 # or k2r - k2l + 1
    for i in range(num_beads):
        i1 = k2r - i + 1 if i > 0 else k1r
        i2 = k2r - i
        x0[i2] = x0[i1] + bLb*np.random.randn(3)
        while x0[i2,0]**2/rx**2 + x0[i2,1]**2/ry**2 + x0[i2,2]**2/rz**2 > 1:
            x0[i2] = x0[i1] + bLb*np.random.randn(3)
    # right free end can be built as normal with no bridge stuff
    k1l, k1r, k2l, k2r = loop_list[-1]
    num_beads = k1r - k1l - 1 # or k2r - k2l + 1
    for i in range(num_beads):
        i1 = k2l + i - 1 if i > 0 else k1l
        i2 = k2l + i
        x0[i2] = x0[i1] + bLb*np.random.randn(3)
        while x0[i2,0]**2/rx**2 + x0[i2,1]**2/ry**2 + x0[i2,2]**2/rz**2 > 1:
            x0[i2] = x0[i1] + bLb*np.random.randn(3)
    return x0

def split_homologs_X(X, N, loop_list):
    N_t, N_tot, _ = X.shape
    X1 = X[:,:N,:].copy()
    X2 = np.zeros_like(X1)
    num_loops, _ = loop_list.shape
    if num_loops == 1:
        X2 = X[:,N:,:]
        return X1, X2
    curr_i2 = 0
    curr_i0 = N
    for i in range(num_loops):
        k1l, k1r, k2l, k2r = loop_list[i]
        num_beads = k1r - k1l - 1
        X2[:,curr_i2:curr_i2+num_beads] = X[:,curr_i0:curr_i0+num_beads]
        curr_i2 += num_beads
        if i < num_loops - 1:
            X2[:,curr_i2] = X[:,k1r]
            curr_i2 += 1
        curr_i0 += num_beads
    return X1, X2

def split_homolog_x(x0, N, loop_list):
    N_tot, _ = x0.shape
    x1 = x0[:N,:].copy()
    x2 = np.zeros((N,3))
    num_loops, _ = loop_list.shape
    if num_loops == 1:
        x2 = x0[N:,:]
        return x1, x2
    curr_i2 = 0
    curr_i0 = N
    for i in range(num_loops):
        k1l, k1r, k2l, k2r = loop_list[i]
        num_beads = k1r - k1l - 1
        x2[curr_i2:curr_i2+num_beads] = x0[curr_i0:curr_i0+num_beads]
        curr_i2 += num_beads
        if i < num_loops - 1:
            x2[curr_i2] = x0[k1r]
            curr_i2 += 1
        curr_i0 += num_beads
    return x1, x2

def rouse_homologs(N, FP, L, b, D, Aex, rx, ry, rz, t, t_save=None):
    """2 polymers + allow some fraction to be homolog paired"""
    homolog_points = np.where(np.random.rand(N) < FP)[0]
    N_tot, loop_list = _homolog_points_to_loops_list(N, homolog_points)
    if t_save is None:
        t_save = t
    return loop_list, jit_rouse_homologs(N, N_tot, loop_list, L, b, D, Aex, rx, ry, rz, t, t_save)

@jit(nopython=True)
def jit_rouse_homologs(N, N_tot, loop_list, L, b, D, Aex, rx, ry, rz, t, t_save):
    """do rouse_homologs given specific homolog pairs"""
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    x0 = _init_homologs(N, N_tot, loop_list, b*np.sqrt(Lb), rx, ry, rz)
    # pre-alloc output
    if t_save is None:
        t_save = t
    x = np.zeros(t_save.shape + x0.shape)
    # setup for saving only requested time points
    save_i = 0
    if 0 == t_save[save_i]:
        x[0] = x0
        save_i += 1
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = t[i] - t[i-1]
        t0 = t[i-1]
        dW = np.random.randn(*x0.shape)
        # -1 or 1, p=1/2
        S = 2*(np.random.rand() < 0.5) - 1
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S)
        # estimate for slope at interval start
        K1 = f_conf(x0, Aex, rx, ry, rz) + f_elas_homolog(x0, N, loop_list, k_over_xi) + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S)
        x1 = x0 + h*K1
        # estimate for slope at interval end
        K2 = f_conf(x1, Aex, rx, ry, rz) + f_elas_homolog(x1, N, loop_list, k_over_xi) + Fbrown
        # average the slope estimates
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x

def plot_homologs(X1, X2, loop_list, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    points = ax.scatter(*X1[loop_list[:-1,1],:].T, c='r', s=50)
    ax.plot(*X1.T)
    ax.plot(*X2.T)
    return points

class HomologViewer(PolymerViewer):

    def __init__(self, r, N, loop_list):
        # will keep track of which *index* in the time axis should be used to
        # draw the current scene
        self.t = 0
        self.N = N
        self.loop_list = loop_list
        self.points = None
        # define the actual axes layout
        self.fig = plt.figure()
        # area to plot polymers
        self.ax3d = plt.axes([0.05, 0.15, 0.9, 0.8], facecolor='w', projection='3d')
        self.r = r
        # area to plot "slider" for selecting what time to plot
        self.ax = plt.axes([0.1, 0.025, 0.8, 0.05], facecolor='lightgoldenrodyellow')
        self.update_ax_limits()
        # slider to control what "time" is plotted, (0,1) is rescaled total
        # simulation time units
        self.slider_t = Slider(self.ax, 't', 0, 1, valinit=0)
        # handler to check if the "time" slider has been moved
        self.slider_t.on_changed(self.update_t)
        plt.show()

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, new_r):
        self._r = new_r
        self._x1, self._x2 = split_homologs_X(self._r, self.N, self.loop_list)
        self._num_time_points, self._num_beads, self._d = self._r.shape
        if hasattr(self, 'ax3d'):
            self.update_ax_limits()
            self.update_drawing()

    # should never be called before a Sim has been loaded successfully
    def update_drawing(self):
        num_polymers = 2
        # if num_polymers != len(self.ax3d.lines):
        self.ax3d.lines = []
        self.ax3d.plot(*self._x1[self.t].T, c='b')
        self.ax3d.plot(*self._x2[self.t].T, c='g')
        if self.points is not None:
            self.points.remove()
            del self.points
        self.points = self.ax3d.scatter(*self._x1[self.t,self.loop_list[:-1,1],:].T, c='r', s=50)
        # for i, r in enumerate([self._x1, self._x2]):
        #     x, y, z = (self.r[self.t,:,0], self.r[self.t,:,1],
        #                self.r[self.t,:,2])
        #     self.ax3d.lines[i].set_data(x, y)
        #     self.ax3d.lines[i].set_3d_properties(z)
        #     self.points.remove()
        #     del self.points
        #     self.points = self.ax3d.scatter(*self._x1[self.t,self.loop_list[:-1,1],:].T, c='r', s=50)
        self.fig.canvas.draw_idle()

    def update_ax_limits(self):
        self.ax3d.set_xlim((np.nanmin(self.r[:,:,0]), np.nanmax(self.r[:,:,0])))
        self.ax3d.set_ylim((np.nanmin(self.r[:,:,1]), np.nanmax(self.r[:,:,1])))
        self.ax3d.set_zlim((np.nanmin(self.r[:,:,2]), np.nanmax(self.r[:,:,2])))

    def update_t(self, val):
        if self.r is None:
            return
        # calculate value on slider
        new_t = int(round(val*self._num_time_points))
        # bound it to 0, num_time_points-1
        new_t = max(0, min(self._num_time_points-1, new_t))
        self.t = new_t
        self.update_drawing()

