"""Fractional (Rouse) polymer. A Rouse polymer in a viscoelastic medium.

For the case alpha=1, these formulas should reduce to the analogous ones in
:mod:`wlcsim.analytical.rouse`."""
from .rouse import kp_over_kbt, rouse_mode, linear_mid_msd
from ..tabulation import frac_vel_corr

from bruno_util.mittag_leffler import ml as mittag_leffler
import numpy as np

from pathlib import Path

def frac_msd(t, alpha, D):
    """MSD of fractionally diffusing free particle.

    Weber, Phys Rev E, 2010 (Eq 10)"""
    return 3*D*np.sin(alpha*np.pi)*t**alpha \
        /(np.pi*(1-alpha/2)*(1-alpha)*alpha)

def frac_cv(t, alpha, D):
    """Velocity autocorrelation of a fractionally-diffusing particle.
    Weber, Phys Rev E, 2010 (Eq 32)"""
    return -(3*D)*np.sin(alpha*np.pi)/(np.pi*(2-alpha))*np.abs(np.power(t, alpha-2))

def frac_discrete_cv(t, delta, alpha, D):
    """Discrete velocity autocorrelation of a fractionally-diffusing particle.
    Weber, Phys Rev E, 2010 (Eq 33)"""
    t = np.atleast_1d(t)
    delta = np.atleast_1d(delta)
    # too many divide by t's to rewrite to avoid these warnings in less than
    # 5min, not worth it.
    with np.errstate(divide='ignore', invalid='ignore'):
        eta = delta/t
        t = t + np.zeros_like(eta) # fix to full size if not already
        delta = delta + np.zeros_like(eta) # fix to full size if not already
        cv_delta_t = frac_cv(t, alpha, D)/(eta*eta*alpha*(1 - alpha))
        cv_delta_t[eta<=1] = cv_delta_t[eta<=1] \
            *(2 - np.power(1 - eta[eta<=1], alpha) - np.power(1 + eta[eta<=1], alpha))
        cv_delta_t[eta>1] = cv_delta_t[eta>1] \
            *(2 + np.power(eta[eta>1] - 1, alpha) - np.power(1 + eta[eta>1], alpha)) \
            + frac_msd(delta[eta>1] - t[eta>1], alpha, D)/delta[eta>1]/delta[eta>1]
        cv_delta_t[t == 0] = frac_msd(delta[t == 0], alpha, D)/np.power(delta[t == 0], 2)
    return cv_delta_t

def frac_discrete_cv_normalized(t, delta, alpha):
    """Normalized discrete velocity autocorrelation of a fractionally-diffusing
    particle. Should be equivalent to

        frac_discrete_cv(t, delta, 1, 1)/frac_discrete_cv(0, delta, 1, 1)

    Lampo, BPJ, 2016 (Eq 5)"""
    return (np.power(np.abs(t - delta), alpha)
        - 2*np.power(np.abs(t), alpha)
        + np.power(np.abs(t + delta), alpha)
        )/(2*np.power(delta, alpha))

def terminal_relaxation_time(alpha, b, N, D):
    return np.power(N*N*b*b/D, 1/alpha)

def frac_rouse_mode_corr(p, t, alpha, b, N, D):
    """Weber Phys Rev E 2010, Eq. 21."""
    # equivalent to the following
    # kp = rouse_mode_coef(p, b, N, kbT)
    # return (3*kbT/kp)*mittag_leffler((-kp*t**alpha)/(N*xi*gamma(3-alpha)), alpha, 1)
    kp_norm = kp_over_kbt(p, b, N)
    # notice kp_over_kbt*D == kp/xi
    return 3*kp_norm*mittag_leffler(-(kp_norm*D)/(N*gamma(3-alpha)) * t**alpha, alpha, 1)

def frac_rouse_mid_msd(t, alpha, b, N, D, num_modes=1000):
    """
    Weber Phys Rev E 2010, Eq. 24.
    """
    if alpha == 1:
        return linear_mid_msd(t, b, N, D, num_modes)
    rouse_corr = 0
    for p in range(1, num_modes+1):
        # k2p = rouse_mode_coef(2*p, b, N, kbT)
        # rouse_corr += 12*kbT/k2p*(1 - mittag_leffler(-k2p*t**alpha/(N*xi*gamma(3-alpha)),
        #         alpha, 1))
        k2p_norm = kp_over_kbt(2*p, b, N)
        rouse_corr += (1/k2p_norm)*(1 - mittag_leffler(
                -(k2p_norm*D)/(N*gamma(3-alpha)) * t**alpha, alpha, 1
        ))
    # return rouse_corr + frac_msd(t, alpha, kbT, xi)/N
    return 12*rouse_corr + frac_msd(t, alpha, D)/N

def rouse_cv_mid(t, alpha, b, N, D, min_modes=1000):
    """Velocity autocorrelation of midpoint of a rouse polymer.

    Weber Phys Rev E 2010, Eq. 33."""
    t = np.atleast_1d(t)
    psum = np.zeros_like(t)
    for p in range(1, min_modes+1):
        # k2p = rouse_mode_coef(2*p, b, N, kbT)
        k2p_norm = kp_over_kbt(2*p, b, N)
        # coef = -k2p/(N*xi*gamma(3-alpha))
        coef = -(k2p_norm*D)/(N*gamma(3-alpha))
        psum += mittag_leffler(coef*np.power(t, alpha), alpha, alpha-1)
    gam = 1 # why Steph, why...
    return frac_cv(t, alpha, D)*(1 + 2*gam*(alpha - 1)*psum)

def tR(alpha, b, N, D):
    """Lampo et al, BPJ, 2016, eq 8"""
    return np.power(N*N*b*b/D, 1/alpha)

def tDeltaN(n1, n2, alpha, b, D):
    """Lampo et al, BPJ, 2016, eq 11"""
    delN = np.abs(n2 - n1)
    return np.power(delN*delN*b*b/D, 1/alpha)

def rouse_cvv_ep(t, delta, p, alpha, b, N, D):
    """Term in parenthesis in Lampo, BPJ, 2016 Eq. 10."""
    t = np.atleast_1d(t)
    delta = np.atleast_1d(delta)
    p = np.atleast_1d(p)
    tpdelta = np.power(np.abs(t + delta), alpha)
    tmdelta = np.power(np.abs(t - delta), alpha)
    ta = np.power(np.abs(t), alpha)
    # kp = rouse_mode_coef(p, b, N, kbT)
    kp_norm = kp_over_kbt(p, b, N)
    # z = -kp/(N*xi*gamma(3 - alpha))
    z = -(kp_norm*D)/(N*gamma(3-alpha))
    return -mittag_leffler(z*tpdelta, alpha, 1) \
        + 2*mittag_leffler(z*ta, alpha, 1) \
        - mittag_leffler(z*tmdelta, alpha, 1)


def frac_rouse_cvv(t, delta, n1, n2, alpha, b, N, D, min_modes=500,
              rtol=1e-5, atol=1e-8, force_convergence=True):
    """Velocity cross-correlation of two points on fractional Rouse polymer

    rtol/atol specify when to stop adding rouse modes. a particle (t,delta)
    pair is considered to have converged when np.isclose returns true give
    rtol/atol and p is even (p odd contributes almost nothing)

    Lampo, BPJ, 2016 Eq. 10."""
    t = np.atleast_1d(t)
    delta = np.atleast_1d(delta)
    tpdelta = np.power(np.abs(t + delta), alpha)
    tmdelta = np.power(np.abs(t - delta), alpha)
    ta = np.power(np.abs(t), alpha)
    # center of mass velocity correlation
    c0 = 3*D/(delta*delta*N*gamma(3-alpha)*gamma(1+alpha)) \
            * (tpdelta + tmdelta - 2*ta)
    # correction factor due to being on a polymer is an expansion in terms of
    # the normal modes of the polymer
    rouse_corr = np.zeros_like(tpdelta)
    # because the rouse modes are cos(p*pi*(n/N)), then the number of p's to
    # include per pass before checking for convergence again should be
    theta_eps = min(n1/N/2, n2/N/2)
    chunk_size = int(np.ceil(2*np.pi/theta_eps))
    # first include the first min_modes for all values of t
    for p in range(1, min_modes+1):
        # kp = rouse_mode_coef(p, b, N, kbT)
        kp_norm = kp_over_kbt(p, b, N)
        # z = -kp/(N*xi*gamma(3-alpha))
        z = -(kp_norm*D)/(N*gamma(3-alpha))
        pdiff = 3/kp_norm*rouse_mode(p, n1, N)*rouse_mode(p, n2, N) \
                *(2*mittag_leffler(z*ta, alpha, 1)
                - mittag_leffler(z*tpdelta, alpha, 1)
                - mittag_leffler(z*tmdelta, alpha, 1))
        rouse_corr = rouse_corr + pdiff
    counts = min_modes*np.ones_like(rouse_corr)
    if not force_convergence:
        return c0 + rouse_corr
    # now include more terms in chuck_size number of p's at a time until
    # convergence.
    p += 1
    old_rouse_corr = rouse_corr.copy() # force alloc
    tol_i = np.ones_like(rouse_corr).astype(bool) # always at least one correction term
    while np.any(tol_i):
        tpd = tpdelta[tol_i]
        tmd = tmdelta[tol_i]
        tai = ta[tol_i]
        old_rouse_corr[:] = rouse_corr # force memmove
        for i in range(chunk_size):
            # kp = rouse_mode_coef(p, b, N, kbT)
            kp_norm = kp_over_kbt(p, b, N)
            # z = -kp/(N*xi*gamma(3-alpha))
            z = -(kp_norm*D)/(N*gamma(3-alpha))
            pdiff = 3/kp_norm*rouse_mode(p, n1, N)*rouse_mode(p, n2, N) \
                    *(2*mittag_leffler(z*tai, alpha, 1)
                        - mittag_leffler(z*tpd, alpha, 1)
                        - mittag_leffler(z*tmd, alpha, 1))
            rouse_corr[tol_i] = rouse_corr[tol_i] + pdiff
            counts[tol_i] += 1
            p += 1
        tol_i = ~np.isclose(rouse_corr, old_rouse_corr, rtol, atol)
    return c0 + rouse_corr, counts

def rouse_nondim_(t, delta, n1, n2, alpha, b, N, D):
    """uses parameters defined by Lampo et al, BPJ, 2016, eq 12"""
    tOverDelta = t/delta
    deltaOverTDeltaN = delta/tR(alpha, b, N, D) \
            *np.power(np.abs(n2-n1)/N, -2/alpha)
    return tOverDelta, deltaOverTDeltaN, alpha
rouse_nondim = np.vectorize(rouse_nondim_)

def un_rouse_nondim(tOverDelta, deltaOverTDeltaN, alpha, delN=0.001):
    """Takes a requested choice of parameters to compare to Tom's calculated
    values and generates a full parameter list that satisfies those
    requirements.

    uses approximation defined by Lampo et al, BPJ, 2016, eq 12

    In Tom's paper, delN is non-dimensinoalized away into deltaOverTDeltaN, but
    that only works in the case where the chain is infinitely long. Otherwise,
    where exactly we choose n1,n2 will affect the velocity correlation because
    of the effects of the chain ends. While other choices are truly arbitrary
    (e.g. D, b, N, etc), delN can be used to specify the size of the segment
    (relative to the whole chain) used to compute the cross correlation."""
    # get arbitrary choices out of the way first
    N = 1
    b = 1
    D = 1
    n1 = N/2
    n2 = N/2 + delN
    # now the computed quantities
    tr = tR(alpha, b, N, D)
    delta = deltaOverTDeltaN*tr*np.power(delN/N, 2/alpha)
    t = tOverDelta*delta
    return t, delta, n1, n2, alpha, b, N, D

def vc(t, delta, beta):
    """velocity correlation of locus on rouse polymer. beta = alpha/2."""
    return ( np.power(np.abs(t - delta), beta)
           + np.power(np.abs(t + delta), beta)
           - 2*np.power(np.abs(t), beta)
           )/( 2*np.power(delta, beta) )

def vcp(t, delta, beta):
    return ( np.power(np.abs(t + delta), beta)
           - np.power(np.abs(t), beta) )*np.power(delta, beta)

def vvc_unscaled_theory(t, delta, beta, A, tDeltaN):
    """velocity cross correlation of two points on rouse polymer."""
    return 2*(vc(t*delta, delta, beta) - frac_vel_corr(t, delta/tDeltaN, 2*beta))

def vvc_rescaled_theory(t, delta, beta, A, tDeltaN):
    """velocity cross correlation of two points on rouse polymer."""
    return 2*A*np.power(delta, beta)*(vc(t*delta, delta, beta) - frac_vel_corr(t, delta/tDeltaN, 2*beta))

def vvc_normalized_theory(t, delta, beta, A, tDeltaN):
    return vvc_unscaled_theory(t, delta, beta, A, tDeltaN) \
         / vvc_unscaled_theory(0, delta, beta, A, tDeltaN)

