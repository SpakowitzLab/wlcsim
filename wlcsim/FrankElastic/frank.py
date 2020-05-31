import numpy as np
import pdb

from .stonefence import precalculate_data as calc_G_set
from .stonefence import ORDER_L
from .yyyreal import YYY
from .invlaplace import invLaplace, make_path, plot_invLaplace_path, invLaplace_path, plot_path
from .rigidrod import rodFrank
from .rigidrod import RR_datapoint

def get_first_pole(gamma, tol=0.000001):
    '''Search the real axes for polls of G[0][0,0]. Returns highest real value.

    Args:
        gamma (float): Field strength.
        tol (float): absolute tolerance of result
    '''
    nninitial=50
    if gamma<=0:
        return 0.0
    if gamma<200:
        p_initial = np.linspace(-0.01, gamma/1.99, nninitial)
    else:
        p_initial = np.linspace(-0.01, gamma, nninitial)
    G_initial = np.zeros(nninitial)
    nums = np.zeros(nninitial)
    for ii, p in enumerate(p_initial):
        G = calc_G_set(p, gamma, m_values=[0])['Gmll']
        nums[ii] = np.real(G[0][0,0])
    i_above = nninitial-1
    if nums[i_above]<0.0:
        print("gamma="+str(gamma))
        raise ValueError("Out of range! Search larger area")
    for ii in range(nninitial-1):
        if nums[i_above-1] >= 0.0:
            i_above=i_above-1
        else:
            break
    lower = p_initial[i_above-1]
    upper = p_initial[i_above]
    while True:
        if upper-lower<tol:
            return (upper+lower)*0.5
        p = (upper+lower)*0.5
        G = calc_G_set(p, gamma, m_values=[0])['Gmll']
        test = np.real(G[0][0,0])
        if test>=0:
            upper=p
        else:
            lower=p
            
            
e0 = np.zeros(ORDER_L); e0[0] = 1.0 # Unit vector
def GJG(G):
    J1 = YYY(0,2,0,0, ORDER_L=ORDER_L)
    return e0@G[0]@J1@G[0]@e0

def GJGJG(G, l1, l2, m1, m2):
    if m1 != m2:
        return 0.0
    J1 = YYY(0, l1, m1, m1, ORDER_L=ORDER_L)
    J2 = YYY(m1, l2, m2, 0, ORDER_L=ORDER_L)
    return e0@G[0]@J1@G[m1]@J2@G[0]@e0

def GJGJGJG(G, l1, l2, m1, m2, mkmk):
    (mk1, mk2) = mkmk
    J1 = YYY(0, l1, m1, m1, ORDER_L=ORDER_L)
    J4 = YYY(m2, l2, m2, 0, ORDER_L=ORDER_L)

    out=0.0+0.0j
    for M in {m1+mk1, m1-mk1, mk1-m1}:
        if M not in {m2+mk2, m2-mk2, mk2-m2}:
            continue
        J2 = YYY(m1, 1, mk1, M, ORDER_L=ORDER_L)
        J3 = YYY(M, 1, mk2, m2, ORDER_L=ORDER_L)
        out = out + e0@G[0]@J1@G[m1]@J2@G[M]@J3@G[m2]@J4@G[0]@e0
    return out

def get_a(gamma, N, gammaNmax=None):
    """
    Maier-Saupe parameter need to generate field strength gamma.

    Args:
        gamma (real): Field strength in kT's per Kuhn length of chain
        N (real): How many Kuhn lengths long the polymers are

    Note: when N=0.0, gamma is assume to be gamma*N, that is kT's per chain.
    """
    if N==0.0:
        return {'aLAf': RR_datapoint(gamma)['aLAf']}
    args={'factor':400,
          'scale':3.0/(N),
          'pole_offset':get_first_pole(gamma)}
    [p_values, cutts] = make_path(**args)
    if gammaNmax is not None and gamma*N > gammaNmax and gamma>gammaNmax:
        return {"K'splay":np.NaN, "K'twist":np.NaN, "K'bend":np.NaN}

    conditions = {}
    conditions["G000"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJG"] = np.zeros(len(p_values), dtype=complex)

    # Do stone fence calculations
    for ii, p in enumerate(p_values):
        G = calc_G_set(p, gamma, m_values=[0])["Gmll"]
        conditions["G000"][ii] = G[0][0,0]
        conditions["GJG"][ii] = GJG(G)

    # Do inverse laplace transforms
    data={}
    for key, values in conditions.items():
        reduction=(2.0*args["scale"]+args["pole_offset"])*N
        data[key]=np.real(invLaplace_path(p_values, cutts, values, N,
                                          reduction=reduction))
    
    data["aLAf"] = gamma*(N**2)*15*data["G000"]\
                   /(8*np.pi*data["GJG"])

    data["a2lpAf"] = data["aLAf"]/N
    return data

def get_Frank_values(gamma, N, gammaNmax=None, laplace_args={}):
    """Returns dictionary of Frank Elastic Values for specified system.

    Nondimensionalized Frank elastc constant K'=K*A/(phi_00*L).
    Maier saupe parameter ether nondimensionalized by L*A*phi_00 or
    2*l_p*A*phi_00 for aLAf and a2lpAf respectively.
    Correlation functions S are also returned.

    Args:
        gamma (real): Field strength in kT's per Kuhn length of chain
        N (real): How many Kuhn lengths long the polymers are
        gammaNmax (real): Cutoff to prevent wasting time.
        laplace_args (dict): Arguements to pass to inverse laplace calculation.


    Note: when N=0.0, gamma is assume to be gamma*N, that is kT's per chain.
    """
    if N==0.0:
        return RR_datapoint(gamma)
    args={'factor':400,
          'scale':3.0/(N),
          'pole_offset':get_first_pole(gamma)}
    for key, value in laplace_args.items():
        args[key]=value
        
    [p_values, cutts] = make_path(**args)
    if gammaNmax is not None and gamma*N > gammaNmax and gamma>gammaNmax:
        return {"K'splay":np.NaN, "K'twist":np.NaN, "K'bend":np.NaN}

    conditions = {}
    conditions["G000"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJG"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_00"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_11"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_couple"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_density"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_-1-1"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_22"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJG_-2-2"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJGJG_11_mk00"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJGJG_11_mk11"] = np.zeros(len(p_values), dtype=complex)
    conditions["GJGJGJG_-1-1_mk11"] = np.zeros(len(p_values), dtype=complex)

    # Do stone fence calculations
    for ii, p in enumerate(p_values):
        G = calc_G_set(p, gamma, m_values=[-3,-2,-1,0,1,2,3])["Gmll"]
        conditions["G000"][ii] = G[0][0,0]
        conditions["GJG"][ii] = GJG(G)
        conditions["GJGJG_00"][ii] = GJGJG(G, 2, 2, 0, 0)
        conditions["GJGJG_11"][ii] = GJGJG(G, 2, 2, 1, 1)
        conditions["GJGJG_couple"][ii] = GJGJG(G, 0, 2, 0, 0)
        conditions["GJGJG_density"][ii] = GJGJG(G, 0, 0, 0, 0)
        conditions["GJGJG_-1-1"][ii] = GJGJG(G, 2, 2, -1, -1)
        conditions["GJGJG_22"][ii] = GJGJG(G, 2, 2, 2, 2)
        conditions["GJGJG_-2-2"][ii] = GJGJG(G, 2, 2, -2, -2)
        conditions["GJGJGJG_11_mk00"][ii] = GJGJGJG(G, 2, 2, 1, 1, (0,0))
        conditions["GJGJGJG_11_mk11"][ii] = GJGJGJG(G, 2, 2, 1, 1, (1,1))
        conditions["GJGJGJG_-1-1_mk11"][ii] = GJGJGJG(G, 2, 2, -1, -1, (1,1))

    # Do inverse laplace transforms
    data={}
    for key, values in conditions.items():
        reduction=(2.0*args["scale"]+args["pole_offset"])*N
        data[key]=np.real(invLaplace_path(p_values, cutts, values, N,
                                          reduction=reduction))
    
    for (m1m2, mkmk) in [("11","11"), ("11","00"), ("-1-1", "11")]:
        data["ddk"+mkmk+" S'22"+m1m2]=4*(4*np.pi/5)*(2/(N**2))\
                                     *data["GJGJGJG_"+m1m2+"_mk"+mkmk]\
                                     /data["G000"]
    for m1m2 in ["11", "-1-1","00","22","-2-2"]:
        data["S'22"+m1m2] = (4*np.pi/5)*(2/(N**2))\
                             *data["GJGJG_"+m1m2]\
                             /data["G000"]
    data["S'0200"] = (4*np.pi/np.sqrt(5))*(2/(N**2))*data["GJGJG_couple"]\
                     /data["G000"]
    data["S'0000"] = (4*np.pi)*(2/(N**2))*data["GJGJG_density"]\
                     /data["G000"]

    data["phi_20/phi_00"] = (N**-1)*np.sqrt(4*np.pi/5.0)*\
                            data["GJG"]/data["G000"]

    #data["K'splay 2"] = 4*np.pi*(data["GJG"]**2)*data["GJGJGJG_11_mk11"]\
    #                  /((N*data["GJGJG_11"])**2 * data["G000"])
    #data["K'twist 2"] = 4*np.pi*(data["GJG"]**2)*data["GJGJGJG_-1-1_mk11"]\
    #                  /((N*data["GJGJG_-1-1"])**2 * data["G000"])
    #data["K'bend 2"] = 4*np.pi*(data["GJG"]**2)*data["GJGJGJG_11_mk00"]\
    #                  /((N*data["GJGJG_11"])**2 * data["G000"])

    data["K'splay"] = (2*np.pi/(N**2))*(data["phi_20/phi_00"]**2)*\
                        (data["S'2211"]**-2)*data["ddk11 S'2211"]
    data["K'twist"] = (2*np.pi/(N**2))*(data["phi_20/phi_00"]**2)*\
                        (data["S'22-1-1"]**-2)*data["ddk11 S'22-1-1"]
    data["K'bend"] = (2*np.pi/(N**2))*(data["phi_20/phi_00"]**2)*\
                        (data["S'2211"]**-2)*data["ddk00 S'2211"]

    data["aLAf"] = gamma*(N**2)*15*data["G000"]\
                   /(8*np.pi*data["GJG"])

    data["a2lpAf"] = data["aLAf"]/N

    #print("splay %f = %f"%(data["K'splay"],data["K'splay 2"]))
    #print("twist %f = %f"%(data["K'twist"],data["K'twist 2"]))
    #print("bend %f = %f"%(data["K'bend"],data["K'splay 2"]))
    return data


def test_Frank_accuracy(N=10, gamma=1000):
    """
    Compair result for different laplace_args.

    Looks really good at N=10, gamma=25
    To within ~0.3% for N=100, gamma=50. Increasing factor is impoartant.
    At N=100 gamma=100 the factor=400 level gets off by a few percent.
    At N=10, gamma=100 it is accurate to may descimal points
    At N=10, gamma=1000 there is a several persent error
    """
    
    args={'factor':400, 'maxp':6000, 'nwing':100}
    data = get_Frank_values(gamma, N, laplace_args=args)
    print(args)
    print([data["K'splay"], data["K'twist"], data["K'bend"]])

    args={'factor':400, 'maxp':24000, 'nwing':400}
    data = get_Frank_values(gamma, N, laplace_args=args)
    print(args)
    print([data["K'splay"], data["K'twist"], data["K'bend"]])

    args={'factor':1200, 'maxp':6000, 'nwing':100}
    data = get_Frank_values(gamma, N, laplace_args=args)
    print(args)
    print([data["K'splay"], data["K'twist"], data["K'bend"]])
#test_Frank_accuracy()


