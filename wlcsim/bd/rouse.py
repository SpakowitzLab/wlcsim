r"""
Simulate Rouse polymers.

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

        \xi \frac{d}{dt} \vec{r}(n, t)
        = k \frac{d^2}{dn^2} \vec{r}(n, t) + f^{(B)}

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
>>> plt.plot(t_save, wlcsim.analytical.rouse.rouse_mid_msd(t_save, b, Nhat, D,
...                                                        num_modes=int(N/2)))
>>> plt.plot(t_save, 6*Dhat*t_save)
>>> plt.plot(t_save, 6*(Dhat/N)*t_save)  # or 6*(D/Nhat)*t_save
>>> # constant prefactor determined empirically...
>>> # otherwise, ~6*Dhat*np.sqrt(t_R * t_save)/N, where t_R is the terminal
>>> # relaxation time of the polymer, t_R = b**2 * N**2 / Dhat
>>> plt.plot(t_save, 1.9544100*b*np.sqrt(D)*np.sqrt(t_save))

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
from .runge_kutta import rk4_thermal_lena
from .init_beads import init_linear_rouse_conf
from .forces import f_conf_ellipse, f_elas_linear_rouse

from numba import njit
import numpy as np


def recommended_dt(N, L, b, D):
    r"""
    Recommended "dt" for use with rouse*jit family of functions.

    Currently set to :math:`\frac{1}{10}\frac{b^2}{6D}`.

    See the :func:`jit_rouse` docstring for source of this time scale.
    """
    Nhat = L/b
    L0 = L/(N-1)
    Dhat = D*N/Nhat
    bhat = np.sqrt(L0*b)
    return (1/10)*bhat**2/(6*Dhat)


def measured_D_to_rouse(Dapp, d, N, bhat=None, regime='rouse'):
    r"""
    Get the full-polymer diffusion coefficient from the "apparent" D.

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
    :math:`t = t_R`, where the MSD transitions into the form
    :math:`6 \frac{\hat{D}}{\tilde{N}} t`. If you ask what "line" with slope
    1/2 in log-log space intersects :math:`6 \frac{\hat{D}}{\tilde{N}} t` at
    :math:`t = t_R`, you get the above form for the sub-Rouse time MSD.

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
    return (Dapp/(kappa*bhat))**2


def with_integrator(N, L, b, D, t, x0=None, integrator=rk4_thermal_lena):
    r"""
    Simulate a Rouse polymer made of N beads free in solution.

    Parameters
    ----------
    N : float
        Number of beads in the chain.
    L : float
        Length of chain.
    b : float
        Kuhn length of the chain (same units as *L*).
    D : float
        Diffusion coefficient. (Units of ``length**2/time``). In order to
        compute *D* from a single-locus MSD, use `measured_D_to_rouse`.
    t : (Nt,) array_like of float
        Time points to use for stepping the integrator. None of our current
        integrators is currently capable of taking time steps larger than the
        time scale of the highest Rouse mode of the finite chain (use
        `~.bd.recommended_dt` to compute the optimal *dt*).
    x0 : (N, 3) array_like of float, optional
        The initial locations of each bead. If not specified, defaults to
        initializing from the free-draining equilibrium, with the first bead at
        the origin.
    integrator : Callable[[ForceFunc, float, Times, Positions], Positions]
        Either `~.runge_kutta.rk4_thermal_lena` or `~.runge_kutta.srk1_roberts`

    Returns
    -------
    (Nt, N, 3) array_like of float
        The positions of the *N* beads at each of the *Nt* time points.

    Notes
    -----
    The polymer beads can be described by the discrete Rouse equations

    .. math::
        \xi \frac{dx(i, t)}{dt} = - k (x(i, t) - x(i+1, t))
                                  - k (x(i, t) - x(i-1, t))
                                  + R(t)

    where :math:`\xi = k_BT/D`, :math:`k = 3k_BT/b^2`, :math:`b` is the Kuhn
    length of the polymer, :math:`D` is the self-diffusion coefficient of a
    bead, and :math:`R(t)/\xi` is a delta-correlated stationary Gaussian
    process with mean zero and :math:`\langle R(t) R(t')/\xi^2 \rangle =
    2DI\delta(t-t')`.

    Notice that in practice, :math:`k/\xi = 3D/b^2`, so we do not need to
    include mass units (i.e. there's no dependence on :math:`k_BT`).
    """
    # derived parameters
    L0 = L/(N-1)  # length per bead
    bhat = np.sqrt(L0*b)  # mean squared bond length of discrete gaussian chain
    Nhat = L/b  # number of Kuhn lengths in chain
    Dhat = D*N/Nhat  # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    if x0 is None:
        x0 = bhat/np.sqrt(3)*np.cumsum(np.random.normal(size=(N, 3)), axis=0)
    # zero = np.array([[0, 0, 0]])

    @njit
    def rouse_f(x, t):
        N, _ = x.shape
        f = np.zeros(x.shape)
        for j in range(1, N):
            for n in range(3):
                f[j, n] += -k_over_xi*(x[j, n] - x[j-1, n])
                f[j-1, n] += -k_over_xi*(x[j-1, n] - x[j, n])
        return f
        # dx = np.diff(x, axis=0)
        # return -k_over_xi*(
        #     np.concatenate([zero, dx]) + np.concatenate([dx, zero])
        # )

    X = integrator(rouse_f, Dhat, t, x0)
    return X


@njit
def jit_srk1(N, L, b, D, t, t_save=None):
    r"""
    Faster version of `wlcsim.bd.rouse.rouse` using jit.

    ``N=101,L=100,b=1,D=1`` takes about 3.5min to run when
    ``t=np.linspace(0, 1e5, 1e7+1)``
    adding ``t_save`` does not slow function down

    Our srk1 scheme is accurate as long as :math:`\Delta t` is less than the
    transition time where the MSD goes from high k behavior (:math:`t^1`) to
    Rouse behavior (:math:`t^{1/2}`).  This is exactly the time required to
    diffuse a Kuhn length, so we just need
    :math:`\Delta t < \frac{\hat{b}^2}{6\hat{D}}` in 3D. the "crossover" from
    fast-k to rouse-like behavior takes about one order of magnitude in time,
    but adding more than that doesn't seem to make the MSD any more accurate,
    so we suggest setting :math:`\Delta t` to one tenth of the bound above.

    the number of orders of magnitude of "rouse" scaling the simulation will
    capture is exactly dictated by the ratio between this time scale and the
    rouse relaxation time (so like ``N**2``?)

    recall (doi & edwards, eq 4.25) that the first mode's relaxation time is
    :math:`\tau_1 = \frac{\xi N^2}{k \pi^2 }`.
    and the :math:`p`\th mode is :math:`\tau_p = \tau_1/p^2` (this is the
    exponential falloff rate of the :math:`p`\th mode's correlation function).
    """
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1)  # length per bead
    bhat = np.sqrt(L0*b)  # mean squared bond length of discrete gaussian chain
    Nhat = L/b  # number of Kuhn lengths in chain
    Dhat = D*N/Nhat  # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position, sqrt(3) since generating per-coordinate
    x0 = bhat/np.sqrt(3)*np.random.randn(N, 3)
    # for jit, we unroll ``x0 = np.cumsum(x0, axis=0)``
    for i in range(1, N):
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
        dW = np.random.randn(*x0.shape)
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S[i])
        # estimate for slope at interval start
        f = np.zeros(x0.shape)
        for j in range(1, N):
            for n in range(3):
                f[j, n] += -k_over_xi*(x0[j, n] - x0[j-1, n])
                f[j-1, n] += -k_over_xi*(x0[j-1, n] - x0[j, n])
        K1 = f + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S[i])
        # estimate for slope at interval end
        x1 = x0 + h*K1
        f = np.zeros(x0.shape)
        for j in range(1, N):
            for n in range(3):
                f[j, n] += -k_over_xi*(x1[j, n] - x1[j-1, n])
                f[j-1, n] += -k_over_xi*(x1[j-1, n] - x1[j, n])
        K2 = f + Fbrown
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x


@njit
def jit_confined_srk1(N, L, b, D, Aex, rx, ry, rz, t, t_save):
    """
    Add an elliptical confinement.

    Energy is like cubed of distance
    outside of ellipsoid, pointing normally back in. times some factor Aex
    controlling its strength.
    """
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1)  # length per bead
    bhat = np.sqrt(L0*b)  # mean squared bond length of discrete gaussian chain
    Nhat = L/b  # number of Kuhn lengths in chain
    Dhat = D*N/Nhat  # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    x0 = np.zeros((N, 3))
    # x0 = np.cumsum(x0, axis=0)
    for i in range(1, N):
        # 1/sqrt(3) since generating per-coordinate
        x0[i] = x0[i-1] + bhat/np.sqrt(3)*np.random.randn(3)
        while x0[i, 0]**2/rx**2 + x0[i, 1]**2/ry**2 + x0[i, 2]**2/rz**2 > 1:
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
        dW = np.random.randn(*x0.shape)
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S[i])
        # estimate for slope at interval start
        f = np.zeros(x0.shape)
        j = 0
        conf = x0[j, 0]**2/rx**2 + x0[j, 1]**2/ry**2 + x0[j, 2]**2/rz**2
        if conf > 1:
            conf_u = np.array([
                -x0[j, 0]/rx**2, -x0[j, 1]/ry**2, -x0[j, 2]/rz**2
            ])
            conf_u = conf_u/np.linalg.norm(conf_u)
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3)
        for j in range(1, N):
            conf = x0[j, 0]**2/rx**2 + x0[j, 1]**2/ry**2 + x0[j, 2]**2/rz**2
            if conf > 1:
                conf_u = np.array([
                    -x0[j, 0]/rx**2, -x0[j, 1]/ry**2, -x0[j, 2]/rz**2
                ])
                conf_u = conf_u/np.linalg.norm(conf_u)
                f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3)
            for n in range(3):
                f[j, n] += -k_over_xi*(x0[j, n] - x0[j-1, n])
                f[j-1, n] += -k_over_xi*(x0[j-1, n] - x0[j, n])
        K1 = f + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S[i])
        # estimate for slope at interval end
        x1 = x0 + h*K1
        f = np.zeros(x1.shape)
        j = 0
        conf = x1[j, 0]**2/rx**2 + x1[j, 1]**2/ry**2 + x1[j, 2]**2/rz**2
        if conf > 1:
            conf_u = np.array([
                -x1[j, 0]/rx**2, -x1[j, 1]/ry**2, -x1[j, 2]/rz**2
            ])
            conf_u = conf_u/np.linalg.norm(conf_u)
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3)
        for j in range(1, N):
            conf = x1[j, 0]**2/rx**2 + x1[j, 1]**2/ry**2 + x1[j, 2]**2/rz**2
            if conf > 1:
                conf_u = np.array([
                    -x1[j, 0]/rx**2, -x1[j, 1]/ry**2, -x1[j, 2]/rz**2
                ])
                conf_u = conf_u/np.linalg.norm(conf_u)
                f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3)
            for n in range(3):
                f[j, n] += -k_over_xi*(x1[j, n] - x1[j-1, n])
                f[j-1, n] += -k_over_xi*(x1[j-1, n] - x1[j, n])
        K2 = f + Fbrown
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x


@njit
def jit_confinement_clean(N, L, b, D, Aex, rx, ry, rz, t, t_save):
    """
    Unfortunately, by "cleaning up" this function, we also make it 2x slower.

    #worthit
    """
    f_c = f_conf_ellipse
    f_e = f_elas_linear_rouse
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1)  # length per bead
    bhat = np.sqrt(L0*b)  # mean squared bond length of discrete gaussian chain
    Nhat = L/b  # number of Kuhn lengths in chain
    Dhat = D*N/Nhat  # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    x0 = init_linear_rouse_conf(N, bhat, rx, ry, rz)
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
    # in order to make it easy to pull into a new function later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = t[i] - t[i-1]
        dW = np.random.randn(*x0.shape)
        # -1 or 1, p=1/2
        S = 2*(np.random.rand() < 0.5) - 1
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S)
        # estimate for slope at interval start
        K1 = f_c(x0, Aex, rx, ry, rz) + f_e(x0, k_over_xi) + Fbrown
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S)
        x1 = x0 + h*K1
        # estimate for slope at interval end
        K2 = f_c(x1, Aex, rx, ry, rz) + f_e(x1, k_over_xi) + Fbrown
        # average the slope estimates
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x
