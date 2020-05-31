import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from multiprocessing import Pool
import pdb
import seaborn as sns
import pickle

from .stonefence import ORDER_L
from .rigidrod import rodFrank
from .frank import *

def multi_process_fun(inputs):
    "pass through for muiltiproccessing"
    return get_Frank_values(*inputs, gammaNmax=10**3)
    
def plot_Frank(log_gamma=True, MultiProcess=True, npts=24, N=1.0,
               gammaN=False, maxgamma=35, nthreads=25):
    """Calculate Frank Elastic data for range of gamma values.

    Args:
        log_gamma (bool): Space gamma logrithmeicly
        MultiProcess (bool): Use multiple threads
        npts (int): number gamma values
        N (real): Kuhn lenths in polymer
        gammaN (bool): Devide gamma by N.  I.e. nondimensionalized by polyer length
        nthreads (int): Number of threads to use
    """
    if log_gamma:
        gammas = np.logspace(0.0,3.5,npts)
    else:
        gammas = np.linspace(0,maxgamma,npts)
    if gammaN:
        gammas=gammas/N
   
    inputs = [(gamma, N) for gamma in gammas]
    if MultiProcess:
        #if __name__ == '__main__':
        with Pool(nthreads) as my_pool:
            results = my_pool.map(multi_process_fun, inputs)
    else:
        results = []
        for inputt in inputs:
            print("Done with %d of %d"%(ii, len(gammas)))
            results.append(get_Frank_values(*inputt))
    
    # Now we need to reorganize the output for plotting
    toplot={}
    for key in results[0].keys():
        toplot[key] = np.zeros(npts)
        for ii, gamma in enumerate(gammas):
            if key not in results[ii].keys():
                toplot[key][ii]=np.NaN
                continue
            if key not in toplot.keys():
                toplot[key]=np.zeros(npts)*np.NaN
            try:
                toplot[key][ii] = results[ii][key]
            except:
                pdb.set_trace()

    toplot["gammas"]=gammas
    return toplot

def get_data_for_Splot(name, nthreads=25, N=0.05, npts=45):
    """Get data for plot_S

    Args:
        name (string): Where to save result.
        nthreads (int): Number of computational threads
        N (real): Kuhn lenths in polymer
        npts (int): number gamma values
    """
    gammaN=False
    data = plot_Frank(N=N, log_gamma=True, MultiProcess=True,
                      npts=npts, gammaN=gammaN, nthreads=nthreads)
    data["npts"]=45
    data["gammaN"]=gammaN
    data["N"]=N
    pickle.dump(data,open(name,"wb"))


def plot_S(name):
    """Plot S over diffetn gamma values.
    Args:
        name (string): File name to load data from.
    """
    data = pickle.load(open(name,"rb"))
    for (l1, l2, m1, m2) in [(0,0,0,0), (2,2,0,0), (0,2,0,0), (0,2,0,0), (2,2,-1,-1),
                             (2,2,1,1), (2,2,-2,-2), (2,2,2,2)]:
        if m1<0:
            linetype = "o-"
        else:
            linetype = ".-"
        label="l1=%d, l2=%d, m1=%d, m2=%d"%(l1, l2, m1, m2)
        llmm = "S'"+str(l1)+str(l2)+str(m1)+str(m2)
        plt.plot(data["gammas"], data[llmm], linetype, label=label)
    plt.xscale("log")
    plt.xlim([1,1000])
    plt.ylim([0,1.05])
    plt.xlabel("gamma")
    plt.legend()
    plt.savefig("Splot.pdf")
    plt.show()

def get_Multi_frank_data(log_gamma=False, npts=48,
                     gammaN=False, maxgamma=35, saveas="Frank.p",
                     nthreads=25, Nset=None):
    """Calculate Frank elastic constants for a range of gamma and N.

    Args:
        log_gamma (bool): Space gamma logrithmeicly
        npts (int): number gamma values
        N (real): Kuhn lenths in polymer
        saveas= Where to save results.
        gammaN (bool): Devide gamma by N.  I.e. nondimensionalized by polyer length
        nthreads (int): Number of threads to use
        Nset (iterable): Where to save results
    """
    if Nset is None:
        Nset = [0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0]
    FofN={}

    for N in Nset:
        FofN[N] = plot_Frank(N=N, log_gamma=log_gamma,
                             MultiProcess=True, npts=npts, gammaN=gammaN,
                            maxgamma=maxgamma, nthreads=nthreads)
    FofN["ORDER_L"]=ORDER_L
    FofN["npts"]=npts
    FofN["log_gamma"]=log_gamma
    FofN["Nset"]=Nset
    FofN["maxgamma"]=maxgamma
    FofN["gammaN"]=gammaN
    pickle.dump(FofN, open(saveas,"wb"))

def Multi_frank_plot(load=None, logy=False):
    """Plot data from get_Multi_frank_data

    Args:
        load (string): File to load
        logy (bool): Log y axis
    """
    FofN=pickle.load(open(load,"rb"))
    Nset = FofN["Nset"]
    log_gamma=FofN['log_gamma']
    gammaN=FofN['gammaN']
    maxgamma=FofN['maxgamma']
    if load in {"./Frank_as_gammaN.p","./Frank_as_gamma.p"}:
        for N in Nset:
            if type(N) != type(0.1):
                continue
            FofN[N]["K'splay"]=FofN[N]["K'splay"]/2.0
            FofN[N]["K'twist"]=FofN[N]["K'twist"]/2.0
            FofN[N]["K'bend"]=FofN[N]["K'bend"]/2.0
    
    for N in Nset:
        if type(N) != type(0.1):
            continue
        a_values = FofN[N]["aLAf"]
        temp=[]
        for idx, val in enumerate(a_values):
            if np.isnan(val):
                continue
            temp.append((val,idx))
        val, idx = min(temp)

        FofN[N]["a_min"] = val
        FofN[N]["idx_a_min"] = idx
    colors = sns.color_palette("bright",len(Nset))

    if (gammaN):
        # ----------------------------------------------------------#
        #                       Low N
        # ----------------------------------------------------------#
        rod_data = rodFrank()

        idx = rod_data["idx_a_min"]
        plt.plot(rod_data["gammas*N"][idx:], rod_data["K_splay"][idx:], "--", color="blue")
        plt.plot(rod_data["gammas*N"][idx:], rod_data["K_bend"][idx:], "--", color="red")
        plt.plot(rod_data["gammas*N"][idx:], rod_data["K_twist"][idx:], "--", color="green")
        for N in Nset:
            if N==0.0316:
                alpha=1.0
                labs=["splay","bend","twist"]
            else:
                alpha=0.3
                labs=[None, None, None]
            if type(N) != type(0.1):
                continue
            gammas=FofN[N]["gammas"]
            idx = FofN[N]["idx_a_min"]
            plt.plot(gammas[idx:]*N, FofN[N]["K'splay"][idx:], color="blue", alpha=alpha,
                     label=labs[0])
            plt.plot(gammas[idx:]*N, FofN[N]["K'bend"][idx:], color="red", alpha=alpha,
                     label=labs[1])
            plt.plot(gammas[idx:]*N, FofN[N]["K'twist"][idx:], color="green", alpha=alpha,
                     label=labs[2])
        if log_gamma:
            plt.xscale("log")
        if logy:
            plt.yscale("log")
        plt.xlabel("gamma*N")
        plt.ylabel("K'")
        plt.title("Low N")
        plt.legend()
        plt.show()

        idx = rod_data["idx_a_min"]
        plt.plot(rod_data["aLAf"][idx:], rod_data["K_splay"][idx:], "--", color="blue")
        plt.plot(rod_data["aLAf"][idx:], rod_data["K_bend"][idx:], "--", color="red")
        plt.plot(rod_data["aLAf"][idx:], rod_data["K_twist"][idx:], "--", color="green")

        for N in Nset:
            #if N==0.0316:
            if N==1.0:
                alpha=0.5
                labs=["splay","bend","twist"]
            else:
                alpha=0.5
                labs=[None, None, None]
            if type(N) != type(0.1):
                continue
            a_values = FofN[N]["aLAf"]
            idx = FofN[N]["idx_a_min"]
            plt.plot(a_values[idx:], FofN[N]["K'splay"][idx:], color="blue", alpha=alpha,
                     label=labs[0])
            plt.plot(a_values[idx:], FofN[N]["K'bend"][idx:], color="red", alpha=alpha,
                     label=labs[1])
            plt.plot(a_values[idx:], FofN[N]["K'twist"][idx:], color="green", alpha=alpha,
                     label=labs[2])
        if log_gamma:
            plt.xscale("log")
        if logy:
            plt.yscale("log")
        plt.xlabel("a(LAf)")
        plt.ylabel("K'")
        plt.title("Low N")
        plt.xlim([5,30])
        plt.ylim([0,2.5])
        plt.legend()
        plt.savefig("Low_N_Kplot.pdf")
        plt.show()
        
        plt.plot(rod_data["aLAf"], rod_data["gammas*N"], "--", color="black",
                 label=["N=0"])
        for iN, N in enumerate(Nset):
            if type(N) != type(0.1):
                continue
            if N > 5:
                continue
            gammas=FofN[N]["gammas"]
            a_values = FofN[N]["aLAf"]
            plt.plot(a_values, gammas*N, label=str(N), color=colors[iN])
        if log_gamma:
            plt.xscale("log")
            plt.yscale("log")
        plt.ylabel("gamma*N")
        plt.xlabel("a(LAf)")
        plt.xlim([0,70])
        plt.ylim([0,35])
        plt.legend()
        plt.savefig("gamma_lowN.pdf")
        plt.show()
    else:
        # ----------------------------------------------------------#
        #                       High N
        # ----------------------------------------------------------#
        for N in Nset:
            #if N==100.:
            if N==1.0:
                alpha=1.0
                labs=["splay","bend","twist"]
            else:
                alpha=0.3
                labs=[None, None, None]
            if type(N) != type(0.1):
                continue
            gammas=FofN[N]["gammas"]
            plt.plot(gammas, FofN[N]["K'splay"]*N, color="blue", alpha=alpha,
                     label=labs[0])
            plt.plot(gammas, FofN[N]["K'bend"]*N, color="red", alpha=alpha,
                     label=labs[1])
            plt.plot(gammas, FofN[N]["K'twist"]*N, color="green", alpha=alpha,
                     label=labs[2])
        if log_gamma:
            plt.xscale("log")
        if logy:
            plt.yscale("log")
        plt.xlabel("gamma")
        plt.ylabel("K'*N")
        plt.title("High N")
        plt.xlim([0,100])
        plt.ylim([0,5])
        plt.legend()
        plt.show()

        for N in Nset:
            #if N==100.:
            if N==1.0:
                alpha=1.0
                labs=["splay","bend","twist"]
            else:
                alpha=0.3
                labs=[None, None, None]
            if type(N) != type(0.1):
                continue
            a_values = FofN[N]["a2lpAf"]
            idx = FofN[N]["idx_a_min"]
            plt.plot(a_values[idx:], FofN[N]["K'splay"][idx:]*N, color="blue", alpha=alpha,
                     label=labs[0])
            plt.plot(a_values[idx:], FofN[N]["K'bend"][idx:]*N, color="red", alpha=alpha,
                     label=labs[1])
            plt.plot(a_values[idx:], FofN[N]["K'twist"][idx:]*N, color="green", alpha=alpha,
                     label=labs[2])
        if log_gamma:
            plt.xscale("log")
        if logy:
            plt.yscale("log")
        plt.xlabel("a(2lpAf)")
        plt.ylabel("K'*N")
        plt.title("High N")
        plt.xlim([0,100])
        plt.ylim([0,5])
        plt.legend()
        plt.show()

        for iN, N in enumerate(Nset):
            if type(N) != type(0.1):
                continue
            if N<0.3:
                continue
            gammas=FofN[N]["gammas"]
            a_values = FofN[N]["a2lpAf"]
            plt.plot(a_values, gammas, label=str(N), color=colors[iN])
        if log_gamma:
            plt.xscale("log")
            plt.yscale("log")
        plt.legend()
        plt.ylabel("gamma")
        plt.xlabel("a(2lpAf)")
        plt.xlim([18,40])
        plt.ylim([0.0,30])
        plt.savefig("gamma_highN.pdf")
        plt.show()
        
def seperate_bend_twist_splay(load, gammaN=False, logy=False):
    """Plot data from get_Multi_frank_data seperatly

    Args:
        load (string): File to load
        logy (bool): Log y axis
    """
    font={'size':12}
    matplotlib.rc('font',**font)
    FofN=pickle.load(open(load,"rb"))
    Nset = FofN["Nset"]
    log_gamma = FofN["log_gamma"]
    if load in {"./Frank_as_gammaN.p","./Frank_as_gamma.p"}:
        for N in Nset:
            if type(N) != type(0.1):
                continue
            FofN[N]["K'splay"]=FofN[N]["K'splay"]/2.0
            FofN[N]["K'twist"]=FofN[N]["K'twist"]/2.0
            FofN[N]["K'bend"]=FofN[N]["K'bend"]/2.0
    for N in Nset:
        if type(N) != type(0.1):
            continue
        a_values = FofN[N]["aLAf"]
        temp=[]
        for idx, val in enumerate(a_values):
            if np.isnan(val):
                continue
            temp.append((val,idx))
        val, idx = min(temp)

        FofN[N]["a_min"] = val
        FofN[N]["idx_a_min"] = idx
    colors = sns.color_palette("bright",len(Nset))
    fig, axs=plt.subplots(3)
    fig.set_figheight(8)
    fig.set_figwidth(5)
    for ii, mode in enumerate(["K'splay","K'bend","K'twist"]):
        for iN, N in enumerate(list(Nset)[::-1]):
            iN=len(Nset)-1-iN
            alpha=1.0
            if type(N) != type(0.1):
                continue
            if N<0.3:
                continue
            a_values = FofN[N]["a2lpAf"]
            idx = FofN[N]["idx_a_min"]
            if mode=="K'splay":
                if N>0.9:
                    label=str(N)
                else:
                    label=None
            elif mode=="K'twist":
                if N<4.0:
                    label=str(N)
                else:
                    label=None
            else:
                label=None
            axs[ii].plot(a_values[idx:], FofN[N][mode][idx:]*N, 
                     color=colors[iN], alpha=alpha,
                     label=label)
        if log_gamma:
            axs[ii].xscale("log")
        if logy:
            axs[ii].yscale("log")
        #axs[ii].set_xlabel("a(2lpAf)")
        axs[ii].set_ylabel("K'*N")
        axs[ii].set_xlim([17,107])
        if mode=="K'splay":
            axs[ii].set_ylim([0,45])
            axs[ii].legend(loc=1, bbox_to_anchor=(0.97,0.97))
        elif mode=="K'bend":
            axs[ii].set_ylim([0,0.8])
        elif mode=="K'twist":
            axs[ii].set_ylim([0,0.04])
            axs[ii].legend(loc=4)
        #axs[ii].set_title(mode)
    plt.tight_layout()
    fig.savefig("Subplots.pdf")
    plt.show()

def find_gamma(N, aset = None, gammaN = True, gamma_TOL=10**-6, a_TOL=10**-6):
    """Frank Elastic data for given N and a

    Args:
        N (float): Number of Kuhn length
        aset (itrable): aLAf or a2lpAf values
        gammaN (bool): If ture use aLAf, else use a2lpAf, amoung other things
    """
    if gammaN:
        if aset is None:
            aset = [20, 30, 40, 50, 60]
        if N==0:
            #in this case get_a will take gammaN
            gamma_min=2.0
            gamma_max=80.0
        else:
            gamma_min = 2.0/N
            gamma_max = 80.0/N
        astring = "aLAf"
    else:
        if aset is None:
            aset = [20, 25, 30, 35, 40]
        gamma_min = 5.0
        gamma_max = 25.0
        astring="a2lpAf"
    gamma_max2 = gamma_max
    
    pairs = {}
    pairs[gamma_max] = get_a(gamma_max, N)[astring]
    pairs[gamma_min] = get_a(gamma_min, N)[astring]


    while True:
        if gamma_max-gamma_min<gamma_TOL:
            break
        if abs(pairs[gamma_max]-pairs[gamma_min])<a_TOL:
            break
    
        gamma_middle = (gamma_max+gamma_min)*0.5
        pairs[gamma_middle] = get_a(gamma_middle, N)[astring]
        #print("gamma %f, a %f"%(gamma_middle,pairs[gamma_middle]))
        if pairs[gamma_max]>pairs[gamma_min]:
            gamma_max = gamma_middle
        else:
            gamma_min = gamma_middle
    
    out={"a_min":pairs[gamma_max], "gamma_at_amin":gamma_max}
    out["aset"]=aset
    out["frank_data"] = []
    out["gammas"] = []
    out["N"]=N

    for a_desire in aset:
        #print("a desire %f"%(a_desire))
        gamma_max = gamma_max2
        gamma_min = out["gamma_at_amin"]
    
        for (gamma, aa) in pairs.items():
            if gamma<=gamma_min:
                continue
            if gamma>=gamma_max:
                continue
            if aa>a_desire:
                gamma_max = gamma
            elif aa<a_desire:
                gamma_min = gamma
        
        while True:
            if a_desire < out["a_min"]:
                gamma_best=None
                break
            a_at_min = pairs[gamma_min]
            a_at_max = pairs[gamma_max]
            if gamma_max-gamma_min<0 or a_at_max-a_at_min<0:
                print("N"+str(N))
                pdb.set_trace()
                raise ValueError("max<min")
            if gamma_max-gamma_min<gamma_TOL or a_at_max-a_at_min<a_TOL:
                gamma_best=(gamma_max+gamma_min)*0.5
                break
            
            gamma_middle = gamma_min + (gamma_max-gamma_min)\
                                       *(a_desire-a_at_min)\
                                       /(a_at_max-a_at_min)

            a_middle = get_a(gamma_middle, N)[astring]
            #print("gamma %f, a %f"%(gamma_middle, a_middle))
            pairs[gamma_middle] = a_middle
            if abs(a_middle-a_desire)<a_TOL:
                gamma_best = gamma_middle
                break
            if a_middle<a_desire:
                gamma_min = gamma_middle
            else:
                gamma_max = gamma_middle 

        out["gammas"].append(gamma_best)
        if gamma_best is None:
            #print("gamma_best = None")
            out["frank_data"].append(None)
        else:
            #print("gamma_best %f"%(gamma_best))
            out["frank_data"].append(get_Frank_values(gamma_best, N))

    return out

def get_Nsweep_data(saveas="Nsweep.p", MultiProcess=True, nthreads=25, Ns=None):
    """Calculate Frank Elastic Constants for sweep over N


    Args:
        saveas= Where to save results.
        nthreads (int): Number of threads to use
        MultiProcess (bool): Use many threads
        Ns (iterable): N values
    """
    if Ns is None:
        npts = 80
        #Ns = list( 3.1*np.logspace(-2,0,npts) )
        Ns = list( np.linspace(0.02,3.0,npts) )
    data = {"Ns":Ns}
    
    if MultiProcess:
        #if __name__ == '__main__':
        with Pool(nthreads) as my_pool:
            results = my_pool.map(find_gamma, Ns)
         
        data["results"] = results
        for ii, N in enumerate(Ns):
            data[N] = results[ii]
    else:
        for N in Ns:
            print("----- Working on N=%f ----- "%(N))
            data[N] = find_gamma(N)

    pickle.dump(data, open(saveas,"wb"))


def plotNsweep(load):
    """Plot results from get_Nsweep_data

    Args:
        load (string): File to load from.
    """
    #data = pickle.load(open("Nsweep2.p","rb"))
    data = pickle.load(open(load,"rb"))
    dataRR = find_gamma(0.0)["frank_data"]

    Ns = np.array(data["Ns"])
    npts = len(Ns)
    bend = np.zeros(npts)
    twist = np.zeros(npts)
    splay = np.zeros(npts)
    aset = data["results"][0]["aset"]
    colors = sns.color_palette("deep",len(aset))
    for jj, aa in enumerate(aset[::-1]):
        ia=len(aset)-1-jj
        color = colors[ia]
        if False:
            if ia == 0:
                labels=["bend/twist", "splay/twist", "bend/splay"]
            else:
                labels=[None, None, None]
        else:
            labels=[str(aa), None, None]

        for ii, N in enumerate(Ns):
            if data["results"][ii]["frank_data"][ia] is None:
                bend[ii]=np.nan
                twist[ii]=np.nan
                splay[ii]=np.nan
                continue
            bend[ii] = data["results"][ii]["frank_data"][ia]["K'bend"]
            twist[ii] = data["results"][ii]["frank_data"][ia]["K'twist"]
            splay[ii] = data["results"][ii]["frank_data"][ia]["K'splay"]

        linestyle = ["-","--",":"]
        #plt.plot(Ns, bend/twist, color=color, linestyle=linestyle[0], label=labels[0])
        #plt.plot(Ns, splay/twist, color=color, linestyle=linestyle[1], label=labels[1])
        #plt.plot(Ns, bend/splay, color=color, linestyle=linestyle[2], label=labels[2])

        #plt.plot(0.0, dataRR[ia]["K'bend"]/dataRR[ia]["K'twist"], "o", color=color)
        #plt.plot(0.0, dataRR[ia]["K'splay"]/dataRR[ia]["K'twist"], "o", color=color)
        xx = np.insert(Ns, 0, 0.0, axis=0) 
        yy = np.insert(bend/twist, 0,
                       dataRR[ia]["K'bend"]/dataRR[ia]["K'twist"], axis=0) 
        plt.plot(xx, yy, color=color, linestyle=linestyle[0], label=labels[0])
        yy = np.insert(splay/twist, 0,
                       dataRR[ia]["K'splay"]/dataRR[ia]["K'twist"], axis=0) 
        plt.plot(xx, yy, color=color, linestyle=linestyle[1], label=labels[1])
    plt.legend()
    #plt.xscale("log")
    plt.yscale("log")
    plt.xlim([0.0,2.0])
    plt.ylim([3.0,100.0])
    plt.plot()
    plt.ylabel("Ratio of Elastic Constants")
    plt.xlabel("N")
    plt.savefig("Nsweep.pdf")
    plt.show()


