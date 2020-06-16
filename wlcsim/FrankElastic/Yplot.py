import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from  scipy.special import sph_harm as Y
import pdb

def plot_surf(fun, title="",Rmax=None, limits=3.5):
    PHI, THETA = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j] #arrays of angular variables
    R = np.abs(fun(PHI, THETA)) #Array with the absolute values of Ylm
    #Now we convert to cartesian coordinates
    # for the 3D representation
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)
    
    
    if Rmax is None:
        N = R/R.max()    # Normalize R for the plot colors to cover the entire range of colormap.
    else:
        N = R/Rmax
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(12,10))
    im = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(N))
    ax.set_title(title, fontsize=20)
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(R)    # Assign the unnormalized data array to the mappable
                      #so that the scale corresponds to the values of R
    fig.colorbar(m, shrink=0.8)

    Xb = limits*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() 
    Yb = limits*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() 
    Zb = limits*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() 
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')
    plt.axis('off')
    plt.savefig(title+".png")
    plt.show()



def YR(m, ell, phi, theta):
    '''
    Real spherical harmonic function.
    '''
    if m<0:
        out = (1.0j/np.sqrt(2))*(Y(m, ell, phi, theta)
                               - (-1)**m * Y(-m, ell, phi, theta))
    if m==0:
        out = Y(m, ell, phi, theta)
    if m>0:
        out = (1.0/np.sqrt(2))*(Y(-m, ell, phi, theta) 
                             + (-1)**m * Y(m, ell, phi, theta))

    if np.isnan(np.real(out)):
        pdb.set_trace()

    return np.real(out)



fields = [[{"ell":2,"m":0,"phi":1.0}]
         , [{"ell":2,"m":0,"phi":1.0},{"ell":2,"m":0,"phi":0.4}]
         , [{"ell":2,"m":0,"phi":1.0},{"ell":2,"m":1,"phi":0.5}]
         , [{"ell":2,"m":0,"phi":1.0},{"ell":2,"m":-1,"phi":0.5}]
         , [{"ell":2,"m":0,"phi":1.0},{"ell":2,"m":2,"phi":0.9}]
         , [{"ell":2,"m":0,"phi":1.0},{"ell":2,"m":-2,"phi":0.9}] ]
titles = ["main", "m0", "m1", "m_minus_1", "m2", "m_minus_2"]
for ii in range(6):
    title=titles[ii]
    field = fields[ii]

    a=1.0
    def fun(phi, theta):
        if hasattr(phi,"__len__"):
            out = phi*0
            for  ii in range(len(phi)):
                out[ii] = fun(phi[ii], theta[ii])
            return out
        temp=0.0
        for X in field:
            temp += a*X["phi"]*YR(X["m"], X["ell"], phi, theta)
        if np.isnan(temp):
           pdb.set_trace() 
        return np.exp(temp)
    
    plot_surf(fun, Rmax=2.3, title=title)
