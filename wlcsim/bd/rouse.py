r"""Simulate Rouse polymers


Notes
-----

There are various ways to parameterize Rouse polymers. In this module, we use
the convention that `N` is the number of beads of a Rouse polymer (contrary to
e.g. Doi & Edwards, where `N` is the number of Kuhn lengths in the polymer). We
instead use `Nhat` for the number of Kuhn lengths in the polymer.

.. warning:

    This is different from the convention in :mod:`wlcsim.analytical.rouse`,
    where :math:`N` is the number of Kuhn lengths in the polymer, to match the
    usual conventions in e.g. Doi & Edwards.

`D` is the "diffusivity of a Kuhn length", i.e. `kbT/xi`, where `xi` is the
dynamic viscosity of the medium in question, as typically found in the Langevin
equation for a Rouse polymer

    .. math::

        \xi \frac{d}{dt} \vec{r}(n, t) = k \frac{d^2}{dn^2} \vec{r}(n, t) + f^{(B)}

This means that to simulate a Rouse polymer diffusing in a medium with
viscosity `xi`, the diffusion coefficient of each bead should be set to `Dhat =
D (N/Nhat)`. Since :math:`\xi` is in units of "viscosity per Kuhn length" and
`(N/Nhat)` is in units of "number of beads over number of Kuhn lengths", this
can be thought of as changing units from "viscosity per Kuhn length" to
"viscosity per bead".

Some books use `b` for the mean (squared) bond length between beads in the
discrete gaussian chain. We instead use `b` to be the real Kuhn length of the
polymer, so that the mean squared bond length `bhat**2` is instead given by
`L0*b`, where `L0` is L0 = L/(N-1) is the amount of polymer "length"
represented by the space between two beads.

The units chosen here make the most sense if you don't actually consider the
polymer to be a "real" Rouse polymer, but instead an e.g. semiflexible chain
with a real "length" whose long distance statistics are being captured by the
simulation.

To compare these results to the :mod:`rouse.analytical.rouse` module, use

>>> plt.plot(t_save, sim_msd)
>>> plt.plot(t_save, wlcsim.analytical.rouse.rouse_mid_msd(t_save, b, Nhat, D, num_modes=int(N/2)))
>>> plt.plot(t_save, 6*Dhat*t_save)
>>> plt.plot(t_save, 6*(Dhat/N)*t_save) # or 6*(D/Nhat)*t_save
>>> # constant prefactor determined empirically...
>>> # otherwise, ~6*Dhat*np.sqrt(t_R * t_save)/N, where t_R is the terminal
>>> # relaxation time of the polymer, t_R = b**2 * N**2 / Dhat
>>> plt.plot(t_save, 1.9544100*bhat*np.sqrt(t_save/Dhat))

where cutting off the number of modes corresponds to ensuring that the rouse
behavior only continues down to the length scale of a single "bead", thus
matching the simulation at arbitrarily short time scales. (Notice that we abide
by the warning above. Our `Nhat` is :mod:`wlcsim.analytical.rouse`'s `N`).

Example
-------

For example, if you want to simulate a megabase of DNA, which has a real linear
length of about 1e6 megabase/base * 0.33 nm/base ~ 3e5 nm, and a Kuhn length of
1e2 nm, then

"""
from ..plot import PolymerViewer
from .runge_kutta import *

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

def measured_D_to_rouse(Dapp, d, N, bhat=None, regime='rouse'):
    r"""Get the full-polymer diffusion coefficient from the "apparent" D

    In general, a discrete Rouse polymer's MSD will have three regimes. On time
    scales long enough that the whole polymer diffuses as a large effective
    particle, the MSD is of the form :math:`6 D/\hat{N} t`, where, as is true
    throughout this module, :math:`D` is the diffusivity of a Kuhn length, and
    :math:`\hat{N}` is the number of Kuhn lengths in the polymer.

    This is true down to the full chain's relaxation time (the relaxation time
    of the 0th Rouse mode, also known as the "Rouse time") :math:`t_R = N^2 b^2
    / D`. For times shorter than this, the MSD will scale as :math:`t^{1/2}`.
    Imposing continuity of the MSD, this means that the MSD will behave as
    :math:`\kappa_0 D/\hat{N} (t_R t)^{1/2}`, where :math:`\kappa_0` is a
    constant that I'm too lazy to compute using :math:`\lim_{t\to 0}` of the
    analytical Rouse MSD result, so I just determine it empirically. We rewrite
    this MSD as :math:`\kappa b D^{-1/2} t^{1/2}`, and find by comparing to the
    analytical Rouse theory that :math:`\kappa = 1.9544100(4)`.

    Eventually (at extremely short times), most "real" polymers will eventually
    revert to a MSD that scales as :math:`t^1`. This cross-over time/length
    scale defines a "number of segments" :math:`\tilde{N}`, where the
    diffusivity of a segment of length :math:`\tilde{L} = L/(\tilde{N} - 1)`
    matches the time it takes stress to propagate a distance :math:`\tilde{L}`
    along the polymer. Said in other words, for polymer lengths smaller than
    :math:`\tilde{L}`, the diffusivity of the segment outruns the stress
    communication time between two neighboring segments.

    The diffusivity at these extremely short times will look like :math:`6 D
    (\tilde{N} / \hat{N}) t`. In order to simulate this behavior exactly, one
    can simply use a polymer with :math:`\tilde{N}` beads, then all three
    length scales of the MSD behavior will match.

    This three-regime MSD behavior can be very easily visualized in log-log
    space, where it is simply the continuous function defined by the following
    three lines. From shortest to longest time scales, :math:`\log{t}
    + \log{6D(\tilde{N}/\hat{N})}`, :math:`(1/2)\log{t} + \log{\kappa b
    D^{1/2}}`, and :math:`\log{t} + \log{6D\hat{N}}`. The lines
    (in log-log space), have slopes 1, 1/2, and 1, respectively, and their
    "offsets" (y-intercept terms with the log removed) are typically referred
    to as :math:`D_\text{app}`.

    The simulations in this module use the diffusivity of a single Kuhn length
    as input (i.e. plain old :math:`D` is expected as input) but typically,
    measurements of diffusing polymers are done below the Rouse time, and for
    long enough polymers (like a mammalian chromosome), it may be difficult to
    impossible to measure unconfined diffusion at time scales beyond the Rouse
    time. This means you often can't just look at the diffusivity of the whole
    polymer and multiply by :math:`\hat{N}` to get the diffusivity of a Kuhn
    length. And since there's not really a principled way in general to guess
    :math:`\tilde{N}`, you usually can't use the short-time behavior to get
    :math:`D` either, even supposing you manage to measure short enough times
    to see the cross-over back to :math:`t^1` scaling. This means that usually
    the only way one can extract :math:`D` is by measuring
    :math:`D_\text{app}`, since we can use the fact that
    :math:`D_\text{app} = \kappa b D^{1/2}` (the :math:`\kappa` value
    quoted above only works for 3-d motion. I haven't bothered to check if it
    scales linearly with the number of dimensions or like :math:`\sqrt{d}` for
    :math:`d`-dimensional motion, but this shouldn't be too hard if you are
    interested).

    In short, this function gives you :math:`D` given :math:`D_\text{app}`.

    Parameters
    ----------
    D_app : float
        The measured :math:`D_\text{app} based on the y-intercept of the MSD in
        log-log space. This can be computed as `np.exp(scipy.stats.linregress(
        np.log(t), np.log(msd))[1])`, assuming that `t` is purely in the regime
        where the MSD scales like :math:`t^{1/2}`.
    Nhat : int
        number of beads in the polymer
    d : int
        number of dimensions of the diffusion (e.g. 3 in 3-d space)

    Returns
    -------
    D_rouse : float
        the diffusivity of a Kuhn length of the polymer

    Notes
    -----
    What follows is an alternate way of "thinking" about the three regimes of
    the MSD, from a simulation-centric perspective, that I typed up back when I
    first wrote this function, but has been obsoleted by the description at the
    start of this docstring.

    In general, a discrete Rouse polymer's MSD will have three regimes. On time
    scales short enough that chain connectivity is not relevant, the MSD is of
    the form :math:`6 \hat{D} t`.

    As soon as chain connectivity begins to affect the dynamics, the form of
    the MSD will change to :math:`\kappa \hat{b} \hat{D}^{1/2} t^{1/2}`, where
    we have determined the constant :math:`\kappa` empirically (using our exact
    analytical theory) to be approximately 1.9544100(4). Note that (up to a
    constant factor), this is just :math:`6 \hat{D} (t_R t)^{1/2}`, where
    :math:`t_R` is the relaxation time (aka Rouse time) of the polymer, given
    by :math:`\hat{b}^2 \tilde{N}^2 / \hat{D}`.

    This is because the Rouse behavior only continues up to the relaxation time
    :math:`t = t_R`, where the MSD transitions into the form :math:`6
    \frac{\hat{D}}{\tilde{N}} t`. If you ask what "line" with slope 1/2 in log-log
    space intersects :math:`6 \frac{\hat{D}}{\tilde{N}} t` at :math:`t = t_R`, you get
    the above form for the sub-Rouse time MSD.

    Here, :math:`\tilde{N}` is distinguished from :math:`N` only as a
    formality. For a simulation of a discrete Rouse polymer with a fixed number
    of beads, :math:`N = \tilde{N}` is exactly the number of beads. For a
    "real" polymer, it is simply a numerical parameter that describes how long
    the Rouse behavior lasts (in the limiit of short times). For a "true" Rouse
    polymer (the fractal object) N is infinity and there is no time scale on
    which the MSD transitions back into :math:`t^1` behavior. But, for any real
    polymer, this will "almost always" eventually happen.
    """
    if d != 3:
        raise ValueError("We've only calculated kappa for d=3")
    kappa = 1.9544100
    return (Dapp/(kappa*b))**2


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
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    def rouse_f(x, t):
        dx = np.diff(x, axis=0)
        return -k_over_xi*(np.concatenate([zero, dx]) + np.concatenate([dx, zero]))
    if x0 is None:
        L0 = L/(N-1) # length per bead
        Lb = L0/b # kuhn lengths per bead
        x0 = b*np.sqrt(Lb)*np.cumsum(np.random.normal(size=(N,3)), axis=0)
    X = rk4_thermal_lena(rouse_f, Dhat, t, x0)
    return X

@jit(nopython=True)
def jit_rouse(N, L, b, D, t, t_save=None):
    r"""faster version of wlcsim.bd.rouse.rouse using jit

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
        # 1/sqrt(3) since generating per-coordinate
        x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
        while x0[i,0]**2/rx**2 + x0[i,1]**2/ry**2 + x0[i,2]**2/rz**2 > 1:
            x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
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
def _init_in_confinement(N, bhat, rx, ry, rz):
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
        # 1/sqrt(3) since generating per-coordinate
        x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
        while x0[i,0]**2/rx**2 + x0[i,1]**2/ry**2 + x0[i,2]**2/rz**2 > 1:
            x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
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
def f_tether(x0, Aex, rx, ry, rz):
    """compute soft (cubic) force towards an elliptical confinement from the
    inside. to tether to a surface, use f_conf + f_tether"""
    N, _ = x0.shape
    f = np.zeros(x0.shape)
    for j in range(N):
        conf = x0[j,0]**2/rx**2 + x0[j,1]**2/ry**2 + x0[j,2]**2/rz**2
        if conf < 1:
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
    """Unfortunately, by "cleaning up" this function, we also make it 2x slower

    #worthit
    """
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    x0 = _init_in_confinement(N, bhat, rx, ry, rz)
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

def homolog_points_to_loops_list(N, homolog_points):
    r"""Create loops list for jit_rouse_linked.

    Parameters
    ----------
    N : int
        number of beads in each polymer if they were unbound
    homolog_points : (L,) array_like
        list of loci indices :math:`{l_i}_1^L, l_i\in[1,N]` (beads) that are
        homologously paired.

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

    where "k[i][l,r]" refers to the bead demarcating the [l]eft or [r]ight end
    of a given "homologous loop" on each of the polymers polymer
    :math:`i\in[1,2]`.

    Namely, the connectivity structure of our polymer "network" is that for
    each `k1l, k1r, k2l, k2r in loop_list`, the beads `k1l` thru `k1r` are
    linked sequentially (as are `k2l` thru `k2r`), but the beads `k1[l,r]` are
    linked to the beads `k2[l,r]`, respectively. This network structure can
    represent arbitrary "homologously linked" polymers, i.e. those of the form

    .. code-block:: text

        --  ---  -------  --  -
          \/   \/       \/  \/
          /\   /\       /\  /\
        --  ---  -------  --  -

    Notice in the example below that `loop_list[1:,1]` and `loop_list[:-1,0]`
    are both the same as `homolog_list`.

    Example
    -------
    Here's an example for a typical case. With 101 beads per polymer (so 100
    inter-bead segments), and one tenth of the beads bound to each other, we
    might get the following Poisson-space homolog points:

    >>> N = 101; FP = 0.1

    >>> homolog_points = np.where(np.random.rand(int(N)) < FP)[0]

    >>> homolog_points
    array([ 5, 18, 27, 59, 68, 82, 90, 94])

    These eight homologous junctions would lead to a simulation array length of
    `2*N - len(homolog_points) = 202 - 8 = 194`. So `N_tot` will be 194, and
    the `loop_list` will look like

    >>> N_tot, loop_list = rouse.homolog_points_to_loops_list(N, homolog_points)

    >>> loop_list
    array([[ -1,   5, 101, 105],
            [  5,  18, 106, 117],
            [ 18,  27, 118, 125],
            [ 27,  59, 126, 156],
            [ 59,  68, 157, 164],
            [ 68,  82, 165, 177],
            [ 82,  90, 178, 184],
            [ 90,  94, 185, 187],
            [ 94, 101, 188, 193]])

    This can be read as follows. The first four beads of each polymer are
    "free", so first row is saying that beads 0 thru 4 correspond to unlinked
    beads in polymer 1. bead 5 is the linked bead (has half the diffusivity).
    since the polymer is of length 101, if there were no connections, the first
    polymer would extend from bead 0 to bead 100, and the second from bead 101
    to bead 201. but since bead "5" is linked to (what would be) bead 106, that
    bead is omitted from the second polymer, and beads 101-105 are the "free
    end" of the second polymer.

    The second row corresponds to the first actual "loop". The two tethering
    points are beads 5 and 18 of the first polymer. What would have been beads
    106 and 119 (that's beads 5 and 18 of the second polymer) are omitted from
    the list of simulated beads. So bead 106 is linked to bead 5 and bead 117
    is linked to bead 18, and beads 106-117 are linked sequentially.

    Similarly, beads 18 and 118 (and 27 and 125) are linked, with beads 118-125
    being linked sequentially.
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
    you have `x0.shape[0] == 2*N - len(loop_list)`"""
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
    """Make the (Nt,Ntot,3) output of rouse_homologs to (2,Nt,N,3)."""
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
    """Make an (N_tot,3) array from rouse_homologs into (2,N,3)."""
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

def homo_sim_loc(i, N, loop_list):
    """Return actual indices for bead "i" of the second polymer in the
    trucated simulation array.

    Parameters
    ----------
    i : (M,) int, array_like
        beads index is requested for
    N : int
        actual number of beads in each polymer
    loop_list : (~N*FP, 4) int
        connectivity of the homologous polymers being simulated

    Returns
    -------
    iloc : (M,) Optional<int>, array_like
        index in simulation array (or None if the bead is a paired bead)
        for each bead in `i`.
    """
    iloc = np.array(i).copy()
    homolog_points = loop_list[1:,0]
    is_homolog = np.isin(iloc, homolog_points)
    if np.any(is_homolog):
        iloc = iloc.astype(object)
        iloc[is_homolog] = None
    for ii in np.where(~is_homolog)[0]:
        iloc[ii] = N + i[ii] - np.sum(homolog_points < i[ii])
    return iloc

def rouse_homologs(N, FP, L, b, D, Aex, rx, ry, rz, t, t_save=None,
                   tether_list=None):
    """BD simulation of two homologous yeast chromosomes in meiosis

    An arbitrary elliptical-shaped confinement can be chosen, and in order to
    simulate synaptonemal formation in prophase, some fraction `FP` of the
    "homologous" beads of the polymer can be "rigidly tethered" to each other.
    The simulation treats these pairs as being one bead (with half the diffusion
    coefficient). The elastic force can then be communicated between the two
    polymers via the connecting loci.

    Depending on what part of prophase is being simulated, either (or both) of
    the telomeres or centromeres can be tethered to the confinement by
    including them in `tether_list`. Currently the confinement force and the
    tethering force are both cubic w.r.t. the distance from the confining
    ellipse with the same strength (`Aex`), and have the form described in
    :func:`f_conf` (and :func:`f_tether`).


    Parameters
    ----------
    N : int
        Number of beads to use in the discretized Rouse chain
    FP : float
        `FP`:math:`\in[0,1] is fraction of beads in the chain (out of `N`) that
        will be tethered to each other
    L : float
        Length of the chain (in desired length units)
    b : float
        Kuhn length of the chain (in same length units as L)
    D : float
        The diffusivity of a Kuhn length of the chain. This can often difficult
        to measure directly, but see the documentation for the function
        :func:`measured_D_to_rouse` for instructions on how to compute this
        value given clean measurements of short-time regular diffusion or
        Rouse-like MSD behavior
    Aex : float
        multiplicative prefactor controlling strength of confinment and
        tethering forces
    rx, ry, rz : float
        radius of confinement in x, y, z direction(s) respectively
    t : (M,) float, array_like
        The BD integrator will step time from t[0] to t[1] to t[2] to ...
    t_save : (m,) float, array_like
        `t_save`:math:`\subset``t` will be the time steps where the simulation
        output is saved. Default is `t_save = t`.
    tether_list : int, array_like
        list of bead indices that will be tethered to the confinement surface,
        for now, both polymers must be tethered at the same beads (but this
        would not be hard to change).

    Returns
    -------
    tether_list : List<int>
        Indices (into compressed `Ntot` output) that were tethered to the
        confinement surface
    loop_list : (~N*FP, 4) int
        Output of homolog_points_to_loops_list used to specify which points are
        tethered
    X : (len(t_save), Ntot, 3)
        output of BD simulation at each time in t_save. The number of output
        beads Ntot <= 2*N excludes the beads on the second polymer whose
        positions are determined by the fact that they are tethered to the
        first polymer. Use :func:`split_homologs_X` to reshape this output into
        shape `(2, len(t_save), N, 3)`.
    """
    if t_save is None:
        t_save = t

    homolog_points = np.where(np.random.rand(int(N)) < FP)[0]
    N_tot, loop_list = homolog_points_to_loops_list(N, homolog_points)

    # homolog points aren't in the simulation array for the second polymer
    tethered_p2 = set(tether_list) - set(loop_list[1:,0])
    # calculate their actual positions in the simulation array
    tether_list = list(tether_list) + list(sim_loc(tethered_p2, N, loop_list))
    tether_list = np.sort(np.array(tether_list))

    x = _jit_rouse_homologs(N, N_tot, tether_list, loop_list, L, b, D, Aex, rx, ry, rz, t, t_save)
    return tether_list, loop_list, x

@jit(nopython=True)
def _jit_rouse_homologs(N, N_tot, tether_list, loop_list, L, b, D, Aex, rx, ry, rz, t, t_save):
    """"Inner loop" of rouse_homologs."""
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
    # extract from loop_list for repeated use in loop
    homolog_points = loop_list[1:,0]
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
        # tethered beads have half the diffusivity
        Fbrown[homolog_points] = Fbrown[homolog_points]/2
        # estimate for slope at interval start
        K1 = f_conf(x0, Aex, rx, ry, rz) + f_elas_homolog(x0, N, loop_list, k_over_xi) + Fbrown
        K1[tether_list] += f_tether(x0[tether_list], Aex, rx, ry, rz)
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S)
        # tethered beads have half the diffusivity
        Fbrown[homolog_points] = Fbrown[homolog_points]/2
        x1 = x0 + h*K1
        # estimate for slope at interval end
        K2 = f_conf(x1, Aex, rx, ry, rz) + f_elas_homolog(x1, N, loop_list, k_over_xi) + Fbrown
        K2[tether_list] += f_tether(x1[tether_list], Aex, rx, ry, rz)
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

