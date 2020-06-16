from numpy import sqrt
from numpy import pi
import numpy as np
from .stonefence import ORDER_L


def a(ell, m):
    if ((3*(ell-m)*(ell+m))/(4*pi*(2*ell-1)*(2*ell+1))<0):
        print(ell, m)
        raise ValueError("negative value")

    return np.sqrt((3*(ell-m)*(ell+m))/(4*pi*(2*ell-1)*(2*ell+1)))

def apm(ell, m, pm):
    if 3*(ell+m)*(ell+m+pm)/(8*pi*(2*ell-1)*(2*ell+1)) < 0:
        print(ell, m, pm)
        raise ValueError("negative value")
    return np.sqrt(3*(ell+m)*(ell+m+pm)/(8*pi*(2*ell-1)*(2*ell+1)))
def ap(ell, m):
    return apm(ell, m, 1)
def am(ell, m):
    return apm(ell, m, -1)

Y00 = 1.0/sqrt(4*pi)





def YYYreal(m1,L,M,m2, ORDER_L):
    '''
    YYYreal_{l_1,l_2} = \\int d\\vec{u} Y_{l_1, m1} Y_{L,M} Y_{l2,m2}

    Not able to handle all combiniations of inputs.   Can only handle
    Y_{l_1,0} Y_{2,M} Y_{l_2,M}
    or
    Y_{l_1,M} Y_{2,M} Y_{l_2,0}
    or
    Y_{l_1,m_1} Y_{1,M} Y{l_2,m_2}
    or
    Y_{l_1,m_1} Y_{0,0} Y{l_2,m_2}

    '''
    
    assert type(ORDER_L) == type(20), "ORDER_L must be an integer"
    assert abs(M) <= L, "l must be >= |m|"

    if m1 == 0 and L==2 and M==m2:
        return YYY_0mm(M, ORDER_L)

    if L==1:
        return YYY_l1(m1, M, m2, ORDER_L)

    if m2 ==0 and L==2 and M==m1:
        return np.transpose( YYY_0mm(M, ORDER_L) )

    if L==0 and M==0:
        if m1 == m2:
            return np.identity(ORDER_L)/np.sqrt(np.pi*4)
        else:
            return np.zeros((ORDER_L,ORDER_L))
    
    print(m1,L,M,m2)
    raise ValueError("Case not handeled")


# -----  Store precalculate Vaues ------
saveYYY = {}
def YYY(m1, L, M, m2, ORDER_L=ORDER_L):
    """
    Same as YYYreal
    """
    key = (m1, L,M, m2)
    if key in saveYYY:
        return saveYYY.get(key)
    else:
        value = YYYreal(m1, L, M, m2, ORDER_L)
        saveYYY[key] = value
        return value

def YYY_0mm(m, ORDER_L):
    """
    Y_{ell1,0} Y_{2,m} Y_{ell2,m}
    """
    out = np.zeros((ORDER_L, ORDER_L))
    if m==0:
        for l1 in range(0, ORDER_L-2):
            out[l1, l1+2] = a(l1+1,0)*a(l1+2,0)/a(2,0)
        for l1 in range(0, ORDER_L):
            out[l1,l1] = (a(l1+1,0)**2 - a(1,0)*Y00 + a(l1,0)**2)/a(2,0)
        for l1 in range(2, ORDER_L):
            out[l1, l1-2] = a(l1,0)*a(l1-1,0)/a(2,0)

    elif abs(m) == 1:
        for l1 in range(0, ORDER_L-2):
            out[l1, l1+2] = ap(l1+1,0)*a(l1+2,1)/a(2,1) 
        for l1 in range(abs(m), ORDER_L):
            out[l1, l1] = (ap(l1+1,0)*a(l1+1,1)-am(l1,0)*a(l1,1))/a(2,1)
        for l1 in range(abs(m)+2, ORDER_L):
            out[l1, l1-2] = -am(l1,0)*a(l1-1,1)/a(2,1)
    elif abs(m) == 2:
        for l1 in range(0, ORDER_L-2):
            out[l1, l1+2] = ap(l1+1,0)*ap(l1+2,1)/ap(2,1)
        for l1 in range(abs(m), ORDER_L):
            out[l1, l1] = -(ap(l1+1,0)*am(l1+1,-1)+am(l1,0)*ap(l1,1))/ap(2,1)
        for l1 in range(abs(m)+2, ORDER_L):
            out[l1, l1-2] = am(l1,0)*am(l1-1,-1)/ap(2,1)
    else:
        raise ValueError("m out or range")
    return out
    

def f(m):
    if m<0:
        return -1.0/sqrt(2)
    if m==0:
        return 0.0
    if m==1:
        return 1.0
    if m>1:
        return 1.0/sqrt(2)

def fzz(m):
    if m<0:
        return -1.0/sqrt(2)
    if m==0:
        return 0.0
    if m==1:
        return 0.0
    if m>1:
        return 1.0/sqrt(2)

def ftt(m):
    if m<-1:
        return -1/sqrt(2)
    if m==-1:
        return -1.0
    if m==0:
        return 1.0
    if m>0:
        return 1.0/sqrt(2)

def YYY_l1(m1,m,m2,ORDER_L):
    """
    Y_{l1,m1} Y_{1,m} Y{l2,m2}
    """
    out = np.zeros((ORDER_L, ORDER_L))
    if m == 0:
        if m1 != m2:
            return out
        for l1 in range(0,ORDER_L-1):
            if l1<abs(m1) or l1+1<abs(m2):
                continue
            out[l1,l1+1] = a(l1+1,m1)
        for l1 in range(1,ORDER_L):
            if l1<abs(m1) or l1-1<abs(m2):
                continue
            out[l1,l1-1] = a(l1,m1)
    if m == 1:
        if m2 == m1+1:
            for l1 in range(0,ORDER_L-1):
                if l1<abs(m1) or l1+1<abs(m2):
                    continue
                out[l1,l1+1] = ap(l1+1,m1)*f(m1+1)
            for l1 in range(1,ORDER_L):
                if l1<abs(m1) or l1-1<abs(m2):
                    continue
                out[l1,l1-1] = -ap(l1,-m1-1)*f(m1+1)
        if m2 == m1-1:
            for l1 in range(0,ORDER_L-1):
                if l1<abs(m1) or l1+1<abs(m2):
                    continue
                out[l1,l1+1] = -ap(l1+1,-m1)*f(m1)
            for l1 in range(1,ORDER_L):
                if l1<abs(m1) or l1-1<abs(m2):
                    continue
                out[l1,l1-1] = ap(l1,m1-1)*f(m1)
    if m == -1:
        if m2 == -m1+1:
            for l1 in range(0,ORDER_L-1):
                if l1<abs(m1) or l1+1<abs(m2):
                    continue
                out[l1,l1+1] = ap(l1+1,-m1)*fzz(m1)
            for l1 in range(1,ORDER_L):
                if l1<abs(m1) or l1-1<abs(m2):
                    continue
                out[l1,l1-1] = -ap(l1,m1-1)*fzz(m1)
        if m2 == -m1-1:
            for l1 in range(0,ORDER_L-1):
                if l1<abs(m1) or l1+1<abs(m2):
                    continue
                out[l1,l1+1] = ap(l1+1,m1)*ftt(m1)
            for l1 in range(1,ORDER_L):
                if l1<abs(m1) or l1-1<abs(m2):
                    continue
                out[l1,l1-1] = -ap(l1,-m1-1)*ftt(m1)
    return out


#raise Exception('exit')
# ------------------  Testing ----------------------------------

from  scipy.special import sph_harm as Y
from scipy.integrate import dblquad as intt

def ind_sphere(fun):
    """Integrate function on sphere.

    Args:
        fun (callable) : fun(theta, phi)
    """
    def refun(phi, theta):
        return fun(theta, phi)*np.sin(theta)
    def gfun(x):
        return 0.0
    def hfun(x):
        return 2*np.pi
    return intt(refun, 0, np.pi, gfun, hfun)[0]

def YR(m, ell, phi, theta):
    '''
    Real spherical harmonic function.
    '''
    if m<0:
        out = (1.0j/sqrt(2))*(Y(m, ell, phi, theta)
                               - (-1)**m * Y(-m, ell, phi, theta))
    if m==0:
        out = Y(m, ell, phi, theta)
    if m>0:
        out = (1.0/sqrt(2))*(Y(-m, ell, phi, theta) 
                             + (-1)**m * Y(m, ell, phi, theta))

    if np.isnan(np.real(out)):
        import pdb
        pdb.set_trace()

    return np.real(out)

def intYYYnum(l1, m1, L, M, l2, m2):
    '''
    Numerical integral over 3 spherical harmics
    '''
    def num_YYY(theta, phi):
        return YR(m1, l1, phi, theta)*\
               YR(M, L, phi, theta)*\
               YR(m2, l2, phi, theta)
    return ind_sphere(num_YYY)



def testYYYreal():
    ORDER_L=20
    values = {}
    L = 1
    for M in [-1,0,1]:
        elmax=5
        mtrx = YYYreal(0, L, M, M, ORDER_L)
        for l1 in range(0,elmax):
            for l2 in range(abs(M),elmax):
                values[(l1,0,L,M,l2,M)] = mtrx[l1,l2]
                # compair mtrx[l1,l2] with Y_{l1,0} Y_{L,M} Y_{l2,M}
    
    L = 2
    for M in [-2,-1,0,1,2]:
        m1=0
        m2=M
        mtrx = YYYreal(m1,L,M,m2,ORDER_L)
        for l1 in range(abs(m1),elmax):
            for l2 in range(abs(m2),elmax):
                values[(l1,m1,L,M,l2,m2)] = mtrx[l1,l2]
                #compair mtrx[l1,l2] with Y_{l1,m1} Y_{L,M} Y_{l2,m2}

    for key, value in values.items():
        print("-------------------")
        print(key)
        num = intYYYnum(*key)
        diff = abs(num - value)
        print(value)
        print(num)
        if diff > 0.00001 or np.isnan(diff):
            print("diff")
            print(diff)
            print("num")
            print(num)
            print("value")
            print(value)
            raise ValueError("Doesn't Match")
#testYYYreal()
