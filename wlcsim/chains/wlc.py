"""Generate equilibrium conformations of Wormlike Chains (WLCs) exactly."""
from ..utils.rotation import *

import numpy as np
from numpy.random import random_sample as urand
import scipy
import scipy.optimize
from scipy.special import erfi
from numba import jit

from functools import partial
import multiprocessing
from multiprocessing import Pool

proc_max = multiprocessing.cpu_count()

def phi_cdf(phi, eps):
    """The cumulative distribution of a given azimuthal angle :math:`phi` given
    a distance between neighboring beads of :math:`\Epsilon = ds/l_p`."""
    f0, fpi = phi_Z_0_1_(eps)
    return (phi_indef_(phi, eps) - f0)/(fpi - f0)

def phi_indef_(phi, eps):
    """proportional to the indefinite integral of sin(phi)exp(-eps*phi**2/2)"""
    return np.real(
            erfi((1 + 1j*eps*phi)/np.sqrt(2*eps))
          + erfi((1 - 1j*eps*phi)/np.sqrt(2*eps))
    )

def phi_Z_0_1_(eps):
    """return F(0), F(pi), so that (phi_indef_(theta)-F(0))/(F(pi)-F(0)) is a
    CDF."""
    return phi_indef_(0, eps), phi_indef_(np.pi, eps)

def phi_indef_2_(phi, eps):
    """alternate representation of the indefinite integral of
    sin(phi)exp(-phi**2/(2*sigma)), different notation?"""
    return np.real(np.sqrt(np.pi*eps/2)/2*np.exp(-eps/2)*(
        - erfi((eps - 1j*phi)/np.sqrt(2*eps))
        - erfi((eps + 1j*phi)/np.sqrt(2*eps))
    ))

def phi_z_2_(eps):
    """normalization factor. phi_indef_(pi)-phi_indef_(0)"""
    return np.real(np.sqrt(np.pi*eps/2)/2*np.exp(-eps/2)*(
        - erfi((eps - 1j*np.pi)/np.sqrt(2*eps))
        + 2*erfi(np.sqrt(eps/2))
        - erfi((eps + 1j*np.pi)/np.sqrt(2*eps))
    ))

def phi_cdf_2(phi, eps):
    """alternate representation of The cumulative distribution of a given
    azimuthal angle :math:`phi` given a distance between neighboring beads of
    :math:`\Epsilon = ds/l_p`."""
    return (phi_indef_(phi, eps) - phi_indef_(0, eps))/phi_z_(eps)

def phi_prob(phi, eps):
    """The probability distribution of a given azimuthal angle :math:`phi`
    given a distance between neighboring beads of :math:`\Epsilon = ds/l_p`."""
    return np.sin(phi)*np.exp(-np.power(phi, 2)/2/eps)/phi_z_(eps)

def dphi_prob(phi, eps):
    """Derivative of phi_prob w.r.t. phi."""
    dphi = np.cos(phi)*np.exp(-np.power(phi, 2)/2/eps) \
         - np.sin(phi)*np.exp(-np.power(phi, 2)/2/eps)*2*phi/2/eps
    return dphi/phi_z_(eps)


def Fzero_(phi, eps, xi, f0, fpi):
    return (phi_indef_(phi, eps) - f0)/(fpi - f0) - xi
def inv_sample_(xi, eps, f0, fpi):
    return scipy.optimize.toms748(Fzero_, 0, np.pi, k=2, args=(eps, xi, f0, fpi))

def dwlc_k(N, L, lp, n_proc=proc_max-1):
    r"""Draw discrete wormlike chain steps, piecewise constant curvature.

    The energy of the wormlike chain is proportional to the curvature of the
    chain squared. We make the assumption that the curvature is constant, and
    that we are approximating the curve via its osculating circle at each
    point.

    For the discrete WLC, we want the next bead to be on the sphere centered at
    the previous bead with radius equal to the length of chain between the
    beads (we call this eps).

    To decide where on that sphere the next bead should be placed, we write the
    curve as r(s) and define the azimuthal angle between the two beads to be
    phi(ds) = arccos(u(s) . u(s+ds)), where ds is the path length between the
    beads. Assuming constant curvature between the two beads, we have simply
    that phi(ds) = kappa*ds. Thus, the probability of a given azimuthal angle
    is just

    .. math::

        P(\phi | \theta) \propto e^{-\beta \frac{l_p}{2} \phi / ds}

    or just a Gaussian, with variance given by :math:`ds/l_p`, the number of
    persistence lengths between the beads.

    Recall the shape of the sphere means that correctly integrating over the
    uniform distribution in the spherical angle leaves us with

    .. math::

        P(\phi) \propto \sin(\phi) e^{-\beta \frac{l_p}{2} \phi / ds}

    the normalization factor to make this a probability distribution for a
    given :math:`ds/l_p` can be calculating analytically (restricting
    explicitly to :math:`\phi\in[0,\pi]` makes the formula uglier but
    WolframAlpha can still do it), yielding an exact probability distribution

    .. math::
        P(\phi) = \sqrt{\frac{2}{\pi (ds/l_p)}} \frac{-1}{e^{ds/2l_p}
        \erf{\sqrt{ds/2l_p}}} \sin(\phi) e^{phi^2/2(ds/l_p)}

    which we can use to inverse-transform sample phi, then draw theta
    uniformly.

    Given theta, phi, we can then get the next tangent vector (and thus next
    location, since for discrete WLC we define u = diff(r)) from the current
    one.

    Notes
    -----
    We using the formalism above, we can get an estimate for what values of
    :math:`ds/l_p` are reasonable to use. For example, for eps=0.1,
    we will typically see values of phi less than pi/3. for eps=0.01, we see
    phi less than phi/10 or so.
    """
    l0 = L/(N-1)
    eps = lp/l0 # ds

    r = np.zeros((N,3))
    Omega = np.identity(3) # n, b, u vectors, arranged in columns
    r[1] = l0*Omega[:,2] # initial in z-direction

    # pre-generate the random numbers for speed
    thetas = 2*np.pi*urand(size=(N-2,))
    xis = urand(size=(N-2,))
    f0, fpi = phi_Z_0_1_(eps)
    with Pool(n_proc) as p:
        phis = p.map(partial(inv_sample_, eps=eps, f0=f0, fpi=fpi), xis)

    for j in range(2, N):
        # first rotate about tangent vector by theta, then "bend" by rotating
        # about the new "normal" vector
        Omega = Omega @ Rz(thetas[j-2]) @ Ry(phis[j-2])

        # third column of orientation matrix is new tangent vector
        r[j] = r[j-1] + l0*Omega[:,2]

    return r

def dwlc(N, L, lp):
    pass

def wlc_init(N, L, lp, lt=0):
    fortran_code = """
subroutine wlc_init(R, U, NB, EPS, l0, rand_stat)
        ! takes R(3,NB) with R(:,1) preset and makes a WLC given EPS
    use mersenne_twister, only : random_number, random_stat
    use params, only : dp, pi

    implicit none

    integer, intent(in) :: NB
    real(dp), intent(in) :: EPS, l0
    real(dp), intent(inout) :: R(3,NB), U(3,NB)
    type(random_stat), intent(inout) :: rand_stat

    integer J
    real(dp) N1(3), N2(3), z, theta
    real(dp) urand(3)

    do J = 2,NB

         call random_number(urand,rand_stat)
         theta = urand(1)*2.0_dp*pi
         z = (1.0_dp/EPS)*log(2.0_dp*sinh(EPS)*urand(2)+exp(-EPS))

         N1 = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
         N1 = N1 - dot_product(N1, U(:,J-1))*U(:,J-1)
         N1 = N1/norm2(N1)

         N2 = cross(N1, U(:,J-1))
         N2 = N2/norm2(N2)

         U(:,J) = sqrt(1-z*z)*(cos(theta)*N1 + sin(theta)*N2 + z*U(:,J-1))
         U(:,J) = U(:,J)/norm2(U(:,J))

         if (WLC_P__LOCAL_TWIST) then
             print*, "wlc chain initialization is not implimented for local twist"
             stop
         endif

         R(:,J) = R(:,J-1) + l0*U(:,J)
     enddo
end subroutine wlc_init
    """
    l0 = L/(N-1)
    eps = l0/lp

    r = np.zeros((N,3))
    u = np.zeros((N,3))
    u[:,2] = 1
    n = np.zeros((N,3))
    n[:,1] = 1
    for j in range(1, N):
        # first get an arbitrary basis n1, n2 to the cotangent plane to u[j-1]
        n1 = np.array([0, 0, 1])
        n1 = n1 - (n1@u[j-1])*u[j-1]
        mag = np.linalg.norm(n1)
        # if we accidentally chose vector parallel to u[j-1]
        if np.isclose(mag, 0):
            n1 = np.array([0, 1, 0])
            n1 = n1 - (n1@u[j-1])*u[j-1]
            mag = np.linalg.norm(n1)
        n1 = n1/mag
        n2 = np.cross(n1, u[j-1])
        n2 = n2/np.linalg.norm(n2)

        # now at a random angle in cotangent plane
        theta = 2*np.pi*urand(1)
        # WLC formula for the magnitude of u[j]@u[j-1]
        z = (1/eps)*np.log(2*urand()*np.sinh(eps)+np.exp(-eps))
        u[j] = np.sqrt(1-z*z)*(np.cos(theta)*n1 + np.sin(theta)*n2 + z*u[j-1])
        u[j] = u[j]/np.linalg.norm(u[j])

        r[j] = r[j-1] + l0*u[j]
    return r, u, n


def effective_wormlike_chain_init(N, L, lp, lt):
    """Generate a wormlike chain with potentially less beads than is wise for
    direct sampling of WLC with 'quadratic' approximation used by wlc
    function, by "regridding" more finely then only saving the beads that were
    originally requested."""

    fortran_code = """subroutine effective_wormlike_chain_init(R, U, NT, wlc_p, rand_stat)
    use mersenne_twister
    use params, only : wlcsim_params, dp, max_wlc_l0, maxWlcDelta
    use vector_utils, only: randomUnitVec
    use polydispersity, only: length_of_chain
    implicit none
    integer, intent(in) :: nt
    type(wlcsim_params), intent(in) :: wlc_p
    real(dp), intent(out) :: R(3,nt), U(3,nt)
    type(random_stat), intent(inout) :: rand_stat

    integer IB, NgB, i, j
    real(dp) l0, eps
    real(dp) urand(3)
    real(dp), dimension(:,:), allocatable :: tmpR, tmpU

    IB = 1
    if (wlc_p%DEL > max_wlc_l0) then
        NgB = ceiling(wlc_p%DEL/max_wlc_l0) + 1
    else
        NgB = 1 + 1
    endif
    allocate(tmpR(3,NgB))
    allocate(tmpU(3,NgB))
    l0 = wlc_p%DEL/(NgB - 1)
    EPS = WLC_P__LP/l0 ! bending rigidity for wormlike chain
    do I = 1,WLC_P__NP
        ! uniformly first bead inside box
        call random_number(urand,rand_stat)
        R(1,IB) = urand(1)*WLC_P__LBOX_X
        R(2,IB) = urand(2)*WLC_P__LBOX_Y
        R(3,IB) = urand(3)*WLC_P__LBOX_Z
        ! uniformly from unit sphere first tan vec
        call randomUnitVec(U(:,IB),rand_stat)
        IB = IB + 1
        do J = 2,length_of_chain(I)
            tmpR(:,1) = R(:,IB-1)
            tmpU(:,1) = U(:,IB-1)
            call wlc_init(tmpR, tmpU, NgB, EPS, l0, rand_stat)
            R(:,IB) = tmpR(:,NgB)
            U(:,IB) = tmpU(:,NgB)
            IB = IB + 1
        enddo
    enddo
    deallocate(tmpR)
    deallocate(tmpU)
end subroutine effective_wormlike_chain_init
"""

