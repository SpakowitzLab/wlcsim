"""Rouse polymer, analytical results.

Notes
-----
There are two parameterizations of the "Rouse" polymer that are commonly used,
and they use the same variable name for two different things.

In one, N is the number of Kuhn lengths, and in the other, N is the number of
beads, each of which can represent an arbitrary number of Kuhn lengths.
"""
from bruno_util.mittag_leffler import ml as mittag_leffler

import numpy as np
import scipy
from scipy.special import gamma
from scipy.signal import savgol_filter, savgol_coeffs
from numba import jit
import mpmath

from functools import lru_cache
from pathlib import Path
import os


@jit
def rouse_mode(p, n, N=1):
    """Eigenbasis for Rouse model.

    Indexed by p, depends only on position n/N along the polymer of length N.
    N=1 by default.

    Weber, Phys Rev E, 2010 (Eq 14)"""
    p = np.atleast_1d(p)
    phi = np.sqrt(2)*np.cos(p*np.pi*n/N)
    phi[p == 0] = 1
    return phi

@jit(nopython=True)
def rouse_mode_coef(p, b, N, kbT=1):
    """k_p: Weber Phys Rev E 2010, after Eq. 18."""
    # alternate: k*pi**2/N * p**2, i.e. k = 3kbT/b**2
    return 3*np.pi**2*kbT/(N*b**2)*p**2

@jit(nopython=True)
def kp_over_kbt(p : float, b : float, N : float):
    """k_p/(k_B T) : "non-dimensionalized" k_p is all that's needed for most
    formulas, e.g. MSD."""
    return (3*np.pi*np.pi)/(N*b*b) * p*p

@jit(nopython=True)
def rouse_mid_msd(t, b, N, D, num_modes=1000):
    """
    modified from Weber Phys Rev E 2010, Eq. 24.
    """
    rouse_corr = np.zeros_like(t)
    for p in range(1, num_modes+1):
        # k2p = rouse_mode_coef(2*p, b, N, kbT)
        # rouse_corr += 12*kbT/k2p*(1 - np.exp(-k2p*t/(N*xi)))
        k2p_norm = kp_over_kbt(2*p, b, N)
        rouse_corr += (1/k2p_norm)*(1 - np.exp(-k2p_norm*(D/N)*t))
    # return rouse_corr + 6*kbT/xi/N*t
    return 12*rouse_corr + 6*D*t/N

def rouse_large_cvv_g(t, delta, deltaN, b, D):
    """Cvv^delta(t) for infinite polymer.

    Lampo, BPJ, 2016 Eq. 16."""
    # k = 3*kbT/b**2
    ndmap = lambda G, arr: np.array(list(map(G, arr)))
    G = lambda x: float(mpmath.meijerg([[],[3/2]], [[0,1/2],[]], x))
    # we can use the fact that :math:`k/\xi = 3D/b^2` to replace
    # gtmd = ndmap(G, np.power(deltaN, 2)*xi/(4*k*np.abs(t-delta)))
    # gtpd = ndmap(G, np.power(deltaN, 2)*xi/(4*k*np.abs(t+delta)))
    # gt = ndmap(G, np.power(deltaN, 2)*xi/(4*k*np.abs(t)))
    # with the same formulas in terms of "D"
    gtmd = ndmap(G, np.power(deltaN, 2)*b*b/(4*3*D*np.abs(t-delta)))
    gtpd = ndmap(G, np.power(deltaN, 2)*b*b/(4*3*D*np.abs(t+delta)))
    gt = ndmap(G, np.power(deltaN, 2)*b*b/(4*3*D*np.abs(t)))
    # tricky conversion from Tom's formula
    # :math:`k_BT / \sqrt{\xi*k} = \frac{k_BT}{\xi} / \sqrt{k / \xi}`
    # and we know :math:`k/\xi = \frac{3D}{b^2}`.
    # so when the dust settles,
    # :math:`\frac{3k_BT}{\delta^2\sqrt{\xi k}} = \frac{b\sqrt{3D}}{\delta^2}
    # so instead of the expression from Tom's paper:
    # return 3*kbT/(np.power(delta, 2)*np.sqrt(xi*k)) * (
    #     np.power(np.abs(t - delta), 1/2)*gtmd
    #   + np.power(np.abs(t + delta), 1/2)*gtpd
    #   - 2*np.power(np.abs(t), 1/2)*gt
    # )
    # we can simply return
    return b*np.sqrt(3*D)*np.power(delta, -2) * (
        np.power(np.abs(t - delta), 1/2)*gtmd
      + np.power(np.abs(t + delta), 1/2)*gtpd
      - 2*np.power(np.abs(t), 1/2)*gt
    )

mod_file = os.path.abspath(__file__)
mod_path = os.path.dirname(mod_file)
def end2end_distance(r, lp, N, L):
    """For now, always returns values for r = np.linspace(0, 1, 50001).
    TODO: fix this

    The end to end distance for a polymer with given parameters.
    end2end_distance(r: array_like(n), lp: float, N: int, L: float):
        -> array_like(n)
    Params:
        r - the values at which to evaluate the end-to-end probability
            distribution
        lp - persistence length of polymer
        N - number of beads in polymer
        L - polymer length
    Output:
        (x,g)
        x = np.linspace(0, 1, 5001)
        g = P(|R| = x | lp, N, L)

    Uses the gaussian chain whenever applicable, ssWLC tabulated values
    otherwise. If you request parameters that require a WLC end-to-end
    distance, the function will ValueError."""
    Delta = L/(lp*(N-1)) # as in wlcsim code
    # WLC case
    actual_r = np.linspace(0, 1, 5001)
    if Delta < 0.01: # as in wlcsim code
        ValueError('end2end_distance: doesn\'t know how to compute WLC case!')
    # ssWLC case
    elif Delta < 10: # as in wlcsim code
        Eps = N/(2*lp*L) # as in wlcsim code
        file_num = round(100*Eps) # index into file-wise tabulated values
        file_name = os.path.join('pdata', 'out' + str(file_num) + '.txt')
        file_name = os.path.join(mod_path, file_name)
        G = np.loadtxt(file_name)
        return (actual_r, G)
    # GC case
    else:
        return (actual_r, end2end_distance_gauss(actual_r, b=2*lp, N=N, L=L))

def end2end_distance_gauss(r, b, N, L):
    """ in each dimension... ? seems to be off by a factor of 3 from the
    simulation...."""
    r2 = np.power(r, 2)
    return 3.0*r2*np.sqrt(6/np.pi)*np.power(N/(b*L), 1.5) \
            *np.exp(-(3/2)*(N/(b*L))*r2)

def test_rouse_msd_line_approx():
    """Figure out what Dapp is exactly using our analytical result.

    if we use msd_approx(t) = 3*bhat*np.sqrt(Dhat*t)/np.sqrt(3)*1.1283791(6)
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
    3*bhat*np.sqrt(Dhat*t)/np.sqrt(3)*1.1283791615 - 0.202642385398*b*L/N
    """
    N = 1e8+1; L = 174; b = 150; D = 166; dt = 1e-2; Nt = 1e6; Nt_save = 1e4;
    t = np.arange(0, Nt*dt, dt); t_save = t[::int(Nt/Nt_save)];
    Nhat = L/b; L0 = L/(N-1); Dhat = D*(N)/Nhat; bhat = np.sqrt(L0*b)
    plt.plot(tmsd, np.abs(msd_rouse - 3*bhat*np.sqrt(Dhat*tmsd)/np.sqrt(3)*1.12837916)/msd_rouse, 'k')
    plt.yscale('log')
    plt.xscale('log')

