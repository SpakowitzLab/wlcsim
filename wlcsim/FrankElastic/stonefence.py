from numpy import sqrt
import numpy as np
#from util import sphinx_compat_jit as jit
from numba import jit

ORDER_L=50

@jit
def alpha(l,m):
    return sqrt((3*(l-m)*(l+m))/(4*np.pi*(2*l-1)*(2*l+1)))
    
@jit
def alpha_plus(l,m):
    return sqrt((3*(l+m)*(l+m+1))/(8*np.pi*(2*l-1)*(2*l+1)))

@jit
def Alm(l,m):
    return alpha(l,m)*alpha(l-1,m)/alpha(2,0)

@jit
def Blm(l,m):
    return (alpha(l+1,m)*alpha(l+1,m) 
            - alpha(1,0)*sqrt(0.25/np.pi) +
            alpha(l,m)*alpha(l,m))/alpha(2,0)

@jit
def PgammaB_vec(m, p, gamma):
    """ P - \gamma \\beta
        
    Returns:
        vector with index ell
    """
    PgammaB = np.zeros(ORDER_L, np.complex128)
    for ell in range(abs(m),ORDER_L):
        PgammaB[ell] = p+ell*(ell+1) - gamma*Blm(ell,m)
    return PgammaB

@jit
def Alm_vec(m):
    Am_vec = np.zeros(ORDER_L, np.complex128)
    for ell in range(abs(m)+2,ORDER_L):
        Am_vec[ell] = Alm(ell,m)
    return Am_vec

@jit 
def Wplus_vec(m, gamma, p, Am, PgammaB):
    Wplus = np.zeros(ORDER_L, np.complex128)
    for ell in (ORDER_L-1,ORDER_L-2):
        Wplus[ell] = 1.0/PgammaB[ell]
    for ell in range(ORDER_L-3,abs(m)-1,-1):
        Wplus[ell] = 1.0/(PgammaB[ell] - (gamma*Am[ell+2])**2*Wplus[ell+2])
    return Wplus

@jit 
def Wminus_vec(m, gamma, p, Am, PgammaB):
    Wminus = np.zeros(ORDER_L, np.complex128)
    for ell in (abs(m),abs(m)+1):
        Wminus[ell] = 1.0/PgammaB[ell]
    for ell in range(abs(m)+2,ORDER_L):
        Wminus[ell] = 1.0/(PgammaB[ell] - (gamma*Am[ell])**2*Wminus[ell-2])
    return Wminus


@jit
def Gmll_matrix(Wplus, Wminus, Am, PgammaB, gamma, m):
    """"Matrox of propagators between starting and ending l value.

    Args:
        Wplus (numpy array): Result of Wplus_vec for same m, p
        Wminus (numpy array): Reult of Wminus_vec for same m, p
        Am (numpy array): Result of Am_vec for same m
        PgammaB (numpy array): Result of PgammaB_vec for same m, p, gamma
        gamma (float): alignment strength, in kT's per Kuhn length
        m (int): z component of agular momentum quantum number

    Returns:
        An ORDER_L x ORDER_L numpy matrix with propagators that use Maier-Saupe
        steps to get from l0 to lf.
    """
    Wpm = np.zeros(ORDER_L,np.complex128)
    absm = abs(m)
    Wpm[absm:] = (Wplus[absm:]*Wminus[absm:])\
                 /(Wminus[absm:] 
                   - PgammaB[absm:]*Wplus[absm:]*Wminus[absm:] + Wplus[absm:])
    Gmll = np.zeros((ORDER_L,ORDER_L), np.complex128)
    for ell in range(abs(m),ORDER_L):
        Gmll[ell,ell] = Wpm[ell]
        for lf in range(ell+2,ORDER_L,2):
            Gmll[ell, lf] = Gmll[ell, lf-2]*Wplus[lf]*Am[lf]*gamma
            Gmll[lf, ell] = Gmll[ell, lf] # Must be symmetric
        # the following loop is an alturnative the using the symmetric propory
        #for lf in range(ell-2,-1,-2):
        #    Gmll[ell, lf] = Gmll[ell, lf+2]*Wminus[lf]*Am[lf+2]*gamma
    return Gmll

def precalculate_data(p, gamma, m_values=[0]):
    """Precalculate W_plus, W_minus, W_pm, and G_m_ll

    Args:
        p (complex): laplace conjugate of path length
        gamma (real): aligning l=2 (Maier-Saupe) field strength
        m_values (list): list of integer m values to precalculate for
    """
    
    Wps = {}
    Wms = {}
    Gmlls = {}
    for m in m_values:
        Am = Alm_vec(m)
        PgammaB = PgammaB_vec(m, p, gamma)
        Wplus = Wplus_vec(m, gamma, p, Am, PgammaB)
        Wminus = Wminus_vec(m, gamma, p, Am, PgammaB)
        Gmll = Gmll_matrix(Wplus, Wminus, Am, PgammaB, gamma, m)

        Wps[m]=Wplus
        Wms[m]=Wminus
        Gmlls[m] = Gmll

    return {"Wplus":Wps, "Wminus":Wms, "Gmll":Gmlls, "ms":m_values, "p":p,
            "gamma":gamma}

