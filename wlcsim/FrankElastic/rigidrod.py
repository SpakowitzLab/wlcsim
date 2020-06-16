import numpy as np
import pdb
from scipy.linalg import expm
import matplotlib
import matplotlib.pyplot as plt

from .yyyreal import YYY
def alpha(l,m):
    return np.sqrt((3*(l-m)*(l+m))/(4*np.pi*(2*l-1)*(2*l+1)))
Y00 = 1.0/np.sqrt(4*np.pi)

lmax = 50


def getM2(m):
    M = np.zeros((lmax+1, lmax+1))

    for l in range(abs(m)+2, lmax+1):
        M[l, l-2] = alpha(l, m)*alpha(l-1, m)/alpha(2, 0)
    for l in range(abs(m), lmax+1):
        M[l, l] = (alpha(l+1, m)**2 - alpha(1, 0)*Y00 +alpha(l, m)**2)/alpha(2, 0)
    for l in range(abs(m), lmax+1-2):
        M[l, l+2] = alpha(l+1, m)*alpha(l+2, m)/alpha(2, 0)

    return M

M2= {}
for m in [-2, -1, 0, 1, 2]:
    M2[m] = getM2(m)


def rodFun(m_1, m_2, m_k_1, m_k_2, eM):
    const = (2*np.pi/3.0)*(eM[0][0,2]**2)/(eM[m_1][2,2]*eM[m_2][2,2]*eM[0][0,0])
    summ = 0.0

    # Calculate M's allowable by selection rules
    if m_k_1==0:
        Mset = {m_1}
    elif m_k_1==-1:
        Mset = {-m_1+1, -m_1-1}
    elif m_k_1==1:
        Mset = {m_1+1, m_1-1}
    if m_k_2==0:
        Mset = Mset.intersection({m_2})
    elif m_k_2==-1:
        Mset = Mset.intersection({-m_2-1,-m_2+1})
    elif m_k_2==1:
        Mset = Mset.intersection({m_2-1, m_2+1})

    for M in Mset:
        summ = summ + (YYY(m_1, 1, m_k_1, M, ORDER_L=lmax+1)@eM[M]
                       @YYY(M, 1, m_k_2, m_2, ORDER_L=lmax+1))[2,2]
    return const*summ





basis = [(-1,-1), (-1,1), (0,-1), (0,1), (1,-1), (1,1)] # (m_k_1, m_1)

def showMatrix(gamma=1.2345):
    out = np.zeros((6,6))
    np.set_printoptions(precision=5)
    np.set_printoptions(linewidth=200)

    eM = {}
    for m in [-2,-1,0,1,2]:
        eM[m] = expm( gamma*M2[m])
    print("eM[-2]")
    print(eM[-2][0:6,0:6])
    print("eM[-1]")
    print(eM[-1][0:6,0:6])
    print("eM[1]")
    print(eM[1][0:6,0:6])
    print("eM[2]")
    print(eM[2][0:6,0:6])

    for ii, (m_k_1, m_1) in enumerate(basis):
        for jj, (m_k_2, m_2) in enumerate(basis):
            out[ii,jj] = rodFun(m_1, m_2, m_k_1, m_k_2, eM)

    print("Frank Matrix")
    print(out)
#showMatrix()

def RR_datapoint(gamma):
    """Takes what is considered gamma*N for WLC
    """
    eM = {}
    for m in [-2,-1,0,1,2]:
        eM[m] = expm( gamma*M2[m])

    K_bend = rodFun(1, 1, 0, 0, eM)
    K_twist = rodFun(-1, -1, 1, 1, eM)
    K_splay = rodFun(1, 1, 1, 1, eM)
    a_prime = 15*gamma*eM[0][0,0]/(4*np.sqrt(np.pi)*eM[0][0,2])
    return {'aLAf':a_prime, "K'bend":K_bend, "K'twist":K_twist, "K'splay":K_splay}

def rodFrank(gammas=None, plot=False): 
    if gammas is None:
        npts=200
        gammas = np.linspace(0.0,50.0,npts)
    else:
        npts=len(gammas)
    K_bend = np.zeros(npts)
    K_twist = np.zeros(npts)
    K_splay = np.zeros(npts)
    a_prime = np.zeros(npts)
    
    for ii in range(npts):
        eM = {}
        for m in [-2,-1,0,1,2]:
            eM[m] = expm( gammas[ii]*M2[m])

        K_bend[ii] = rodFun(1, 1, 0, 0, eM)
        K_twist[ii] = rodFun(-1, -1, 1, 1, eM)
        K_splay[ii] = rodFun(1, 1, 1, 1, eM)
        
        a_prime[ii] = 15*gammas[ii]*eM[0][0,0]/(4*np.sqrt(np.pi)*eM[0][0,2])
   
    out = {"K_bend":K_bend, "K_twist":K_twist, "K_splay":K_splay}
    temp=[]
    for idx, val in enumerate(a_prime):
        if np.isnan(val):
            continue
        temp.append((val,idx))
    val, idx = min(temp)

    out["a_min"] = val
    out["idx_a_min"] = idx
    out["gammas*N"] = gammas
    out["aLAf"] = a_prime


    if plot:
        matplotlib.rcParams['text.usetex'] = True
        font = {'size':16}
        matplotlib.rc('font', **font)
        fig, ax = plt.subplots(tight_layout=True)
        ax.plot(gammas, K_bend, label="$K_{bend}$")
        ax.plot(gammas, K_twist, label="$K_{twist}$")
        ax.plot(gammas, K_splay, label="$K_{splay}$")
        plt.xlabel(r'$ \gamma $')
        plt.ylabel(r'$ K\left(\frac{A}{L}\frac{V}{nAL}\right)$ ')
        #plt.ylabel("K'")
        plt.legend()
        plt.show()

        plt.plot(gammas, a_prime)
        plt.title(r"$a\left(\frac{nA^2 L^2}{V}\right)$")
        plt.xlabel(r'$\gamma$')
        plt.tight_layout()
        plt.show()

        fig, ax = plt.subplots(tight_layout=True)
        ax.plot(a_prime, K_bend, label="$K_{bend}$")
        ax.plot(a_prime, K_twist, label="$K_{twist}$")
        ax.plot(a_prime, K_splay, label="$K_{splay}$")
        plt.xlabel(r"$a\left(\frac{nA^2 L^2}{V}\right)$")
        plt.ylabel(r'$ K\left(\frac{A}{L}\frac{V}{nAL}\right)$ ')
        #plt.ylabel("K'")
        plt.legend()
        plt.tight_layout()
        plt.show()
    else:
        return out

#rodFrank(plot=True)
    

