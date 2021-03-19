"""
Some spare stochastic integrators.

Not currently used because if you pass a function as an argument you can't
@numba.jit that.

For rouse chains, Bruno (that's me) has tested the srk1 integrator and found it
to work extremely well (see notes in wlcsim.bd.rouse on suggested dt).

For WLC/ssWLC, Lena found it useful to use a scheme that was higher
deterministic order in time, to resolve the high elastic energies involved.
Bruno has not yet tested this explicitly.

You can find Lena's algorithm (hold brownian force constant over an rk4 step)
below. It is unlikely to be strongly convergent except if you subsample below
the actual desired time resolution.
"""
from numba import jit
import numpy as np


def rk4_thermal_lena(f, D, t, x0):
    r"""
    Integrate using the method of Steph Weber (Phys. Rev. E., 82, 011913).

    Solves the equation

    .. math::

        x'(t) = f(x(t), t) + \xi(t),

    where :math:`\xi` represents thermal forces, and each dimension of ``x0``
    has diffusivity :math:`D`.

    ``x0`` is :math:`x(t=0)`, and :math:`f: \mathbb{R}^n \times \mathbb{R} ->
    \mathbb{R}^n`
    """
    t = np.array(t)
    x0 = np.array(x0)
    x = np.zeros(t.shape + x0.shape)
    dts = np.diff(t)
    x[0] = x0
    dxdt = np.zeros((4,) + x0.shape)  # one for each RK step
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = dts[i-1]
        t0 = t[i-1]
        x0 = x[i-1]
        x_est = x0
        Fbrown = np.sqrt(2*D/(t[i]-t[i-1]))*np.random.normal(size=x0.shape)
        dxdt[0] = f(x0, t0) + Fbrown  # slope at beginning of time step
        x_est = x0 + dxdt[0]*(h/2)  # first estimate at midpoint
        dxdt[1] = f(x_est, t0 + (h/2)) + Fbrown  # estimated slope at midpoint
        x_est = x0 + dxdt[1]*(h/2)  # second estimate at midpoint
        dxdt[2] = f(x_est, t0 + (h/2)) + Fbrown  # second slope at midpoint
        x_est = x0 + dxdt[2]*h  # first estimate at next time point
        dxdt[3] = f(x_est, t0 + h) + Fbrown  # slope at end of time step
        # final estimate is weighted average of slope estimates
        x[i] = x0 + h*(dxdt[0] + 2*dxdt[1] + 2*dxdt[2] + dxdt[3])/6
    return x


def euler_maruyama_with_torque(ft, D_r, D_u, t, x0, u0, *args, **kwargs):
    r"""
    The most well-known stochastic integrator, this time with Torques!

    Parameters
    ----------
    ft : Func[(Array[Pos], Array[Orientation]), (Array[Force], Array[Torq])]
        The function which take the locations and orientations of each bead and
        returns the force and torques on each bead.
    D_r : float
        Translational diffusivity (same as *D* in `euler_maruyama`).
    D_u : float
        Rotational diffusivity (if unsure, use 1/xiu as returned by
        `tabulation.dsswlc_params`).
    t : (N, ) array-like of float
        Times over which to intergrate.
    x0 : (M, 3) array-like of float
        Initial position of beads.
    u0 : (Q, 3) array-like of float
        Initial orientations.

    Returns
    -------
    x : (N, M, 3)
        Position of all beads at each time point requested.
    u : (N, Q, 3)
        Orientations at each time point requested.
    """
    t = np.array(t)
    x0 = np.array(x0)
    x = np.zeros(t.shape + x0.shape)
    u = np.zeros(t.shape + u0.shape)
    dts = np.diff(t)
    x[0] = x0
    u[0] = u0
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = dts[i-1]
        x0 = x[i-1]
        u0 = u[i-1]
        Fbrown = np.sqrt(2*D_r/(t[i]-t[i-1]))*np.random.normal(size=x0.shape)
        Tbrown = np.sqrt(2*D_u/(t[i]-t[i-1]))*np.random.normal(size=u0.shape)
        F, T = ft(x0, u0, *args, **kwargs)
        x[i] = x0 + h*(Fbrown + F)
        u[i] = u0 + h*np.cross(Tbrown + F, u0)
        u[i] /= np.linalg.norm(u[i])
    return x, u


def srk1_roberts(f, D, t, x0):
    r"""
    From wiki, from A. J. Roberts. Modify the improved Euler scheme to
    integrate stochastic differential equations. [1], Oct 2012.

    If we have an Ito SDE given by

    .. math::

        d\vec{X} = \vec{a}(t, \vec{X}) + \vec{b}(t, \vec{X}) dW

    then

    .. math::

        \vec{K}_1 = h \vec{a}(t_k, \vec{X}_k) + (\Delta W_k - S_k\sqrt{h})
            \vec{b}(t_k, \vec{X}_k)

    .. math::

        \vec{K}_2 = h \vec{a}(t_{k+1}, \vec{X}_k + \vec{K}_1)
            + (\Delta W_k - S_k\sqrt{h}) \vec{b}(t_{k+1},\vec{X}_k + \vec{K}_1)

    .. math::

        \vec{X}_{k+1} = \vec{X}_k + \frac{1}{2}(\vec{K}_1 + \vec{K}_2)

    where :math:`\Delta W_k = \sqrt{h} Z_k` for a normal random :math:`Z_k \sim
    N(0,1)`, and :math:`S_k=\pm1`, with the sign chosen uniformly at random
    each time.
    """
    t = np.array(t)
    x0 = np.array(x0)
    x = np.zeros(t.shape + x0.shape)
    dts = np.diff(t)
    # -1 or 1, p=1/2
    S = 2*(np.random.random_sample(size=t.shape) < 0.5) - 1
    x[0] = x0
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = dts[i-1]
        t0 = t[i-1]
        x0 = x[i-1]
        dW = np.random.normal(size=x0.shape)
# D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*D/h)*(dW - S[i])
        # estimate for slope at interval start
        K1 = f(x0, t0) + Fbrown
        Fbrown = np.sqrt(2*D/h)*(dW + S[i])
        # estimate for slope at interval end
        K2 = f(x0+h*K1, t0+h) + Fbrown
        x[i] = x0 + h * (K1 + K2)/2
    return x


# simple test case
def ou(x0, t, k_over_xi, D, method=rk4_thermal_lena):
    r"""
    Test simulation: the Ornstein-Uhlenbeck process.

    Simulates a standard Ornstein-Uhlenbeck process with ``theta=k_over_xi``
    and ``sigma^2/2=D``. This is the appropriate parameterization to match the
    limiting case of ``N=2`` in `wlcsim.bd.rouse.with_integrator`.
    """
    def f(x, t):
        return -k_over_xi*x
    return method(f, D=D, t=t, x0=x0)


@jit(nopython=True)
def _get_scalar_corr(X):
    "fast correlation calculation for testing"
    num_samples, num_t = X.shape
    corr = np.zeros((num_t,))
    count = np.zeros((num_t,))
    for i in range(num_samples):
        for j in range(num_t):
            for k in range(j, num_t):
                corr[k-j] = corr[k-j] + X[i, k]*X[i, j]
                count[k-j] = count[k-j] + 1
    return corr, count


@jit(nopython=True)
def _get_vector_corr(X):
    "fast correlation calculation for testing"
    num_samples, num_t, d = X.shape
    corr = np.zeros((num_t,))
    count = np.zeros((num_t,))
    for i in range(num_samples):
        for j in range(num_t):
            for k in range(j, num_t):
                corr[k-j] = corr[k-j] + X[i, k]@X[i, j]
                count[k-j] = count[k-j] + 1
    return corr, count


@jit(nopython=True)
def _get_bead_msd(X, k=None):
    """center bead by default

    for 1e4-long time arrays, this takes ~10-30s on my laptop"""
    num_t, num_beads, d = X.shape
    if k is None:
        k = max(0, num_beads/2 - 1)
    k = int(k)
    ta_msd = np.zeros((num_t,))
    count = np.zeros((num_t,))
    for i in range(num_t):
        for j in range(i, num_t):
            ta_msd[j-i] += (X[j, k] - X[i, k])@(X[j, k] - X[i, k])
            count[j-i] += 1
    return ta_msd, count


@jit(nopython=True)
def _msd(x):
    result = np.zeros_like(x)
    for delta in range(1, len(x)):
        for i in range(delta, len(x)):
            result[delta] += (x[i] - x[i-delta])**2
        result[delta] = result[delta] / (len(x) - delta)
    return result


def test_ou_autocorr(method=srk1_roberts):
    """Test an individual integrator visually using an OU process."""
    import scipy.stats
    import matplotlib.pyplot as plt
    k = 2
    xi = 4
    kbT = 1
    D = kbT/xi
    k_over_xi = k/xi
    x0 = scipy.stats.norm.rvs(scale=np.sqrt(D/k_over_xi), size=(10_000,))
    t = np.linspace(0, 1e2, 1e3+1)
    X = ou(x0, t, k_over_xi, D, method=method)
    assert(np.abs(np.var(X) - D/k_over_xi)/(D/k_over_xi) < 0.1)
    corr, count = _get_scalar_corr(X.T)
    err = corr/count - (kbT/k)*np.exp(-k_over_xi*t)
    plt.figure()
    plt.plot(t, err)
    plt.figure()
    plt.hist(X[-100:].flatten(), bins=100, density=1)
    x = np.linspace(-3, 3, 100)
    plt.plot(x, scipy.stats.norm(scale=np.sqrt(D/k_over_xi)).pdf(x))
    return X
