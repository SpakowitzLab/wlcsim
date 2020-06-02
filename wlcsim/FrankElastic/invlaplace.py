import numpy as np
from scipy.integrate import simps
from matplotlib import pyplot as plt

def invLaplace(p_values, G_values, L, reduction=0.0):
    r"""
    Numerical Inverse Laplace from :math:`G(p) \to G(L)`.

    Parameters
    ----------
    p_values : array_like
        1D complex array of p values (i.e. the contour)
    G_values : array_like
        1D complex array G(p)
    L : float
        Polymer length (In Kuhn lengths)
    reduction : float
        multiply answer by exp(-reduction) for better numerics
    """
    x_values = np.imag(p_values)
    # print(x_values[0:10])
    y_values = np.array(G_values)*np.exp(p_values*L-reduction)
    return simps(y_values, x_values)/(2*np.pi)

def concentrate_high(low, high, halfpoints):
    if low>high:
        return np.flip(concentrate_high(high, low, halfpoints),0)
    mid = 0.9*high+0.1*low
    upps = np.linspace(mid, high, halfpoints)
    lows = np.linspace(low, mid, halfpoints)
    return np.concatenate((lows, upps[1:]))


def make_path(factor=100, scale=1.0, lambda_offset=0.1, width=2.0, depth=2.0,
              nwing=100, nside=30, nback=20, maxp=6000, pole_offset=0.0):
    """
    Path for complex integral that draws a half rectangle about the origin.

    Parameters
    ----------
    factor : int
        Number of points to caltulate (sort of)
    scale : float
        Scale of path.  Recomended: around 10.0/N
    lambda_offset : float
        distance from complex axis
    width : float
        width of rectangle
    depth : float
        depth of rectangle
    nwing : int
        number of points in wing before ``factor``
    nside : int
        number of points to a side before ``factor``
    back : int
        number of points on the back before ``factor``
    maxp : float
        approximation for infinity
    pole_offset : float
        add to depth
    """
    assert depth >= lambda_offset
    nside=nside*factor
    nback=nback*factor
    nwing=nwing*factor
    depth=depth*scale+pole_offset
    lambda_offset=lambda_offset*scale
    maxp=maxp*scale
    width=width*scale


    lower = np.linspace(-maxp, -width, nwing)*1j+lambda_offset

    low_side = -width*1j+np.linspace(lambda_offset,depth,nside)
    #low_side = -width*1j+concentrate_high(lambda_offset,depth,int(nside/2))

    back = np.linspace(-width, width, nback)*1j+depth

    top_side = width*1j+np.linspace(depth, lambda_offset, nside)
    #top_side = width*1j+concentrate_high(depth, lambda_offset, int(nside/2))

    upper = np.linspace(width, maxp, nwing)*1j+lambda_offset

    p_values = np.concatenate((lower, low_side, back, top_side, upper))

    cutts = {'lower':[0, nwing],
             'low_side':[nwing, nwing+nside],
             'back':[nwing+nside, nwing+nside+nback],
             'top_side':[nwing+nside+nback, nwing+2*nside+nback],
             'upper':[nwing+2*nside+nback, 2*nwing+2*nside+nback]}

    return [p_values, cutts]

def plot_path(path):
    plt.plot(np.real(path), np.imag(path),'.-')
    plt.axis('equal')
    plt.show()

def invLaplace_path(path, cutts, G_values, L, reduction=0.0):
    '''Inverse laplace transform based on path.

    Args:
        path (ndarray): complex path values
        cutts (dict): specified by path
        G_values
    '''
    total=0.0+0.0j
    for key, cutt in cutts.items():
        p_values=path[cutt[0]:cutt[1]]
        y_values=G_values[cutt[0]:cutt[1]]*np.exp(p_values*L-reduction)

        if key == 'low_side':
            x_values=np.real(p_values)
            sgn=1.0/(1j*2*np.pi)
        elif key == 'top_side':
            # Need to reverse to give increasing x_values
            x_values=-np.real(p_values)
            sgn=-1.0/(1j*2*np.pi)
        else:
            x_values=np.imag(p_values)
            sgn=1.0/(2*np.pi) # the 1j/1j cancels

        total = total + sgn*simps(y_values, x_values)

    return total

def plot_int_path(G_values, cutts):
    for key, cutt in cutts.items():
        plt.plot(np.arange(cutt[0],cut[1]), np.real(G_values[cutt[0]:cutt[1]]),
                 label=key)
    plt.show()

def plot_invLaplace_path(path, cutts, G_values, L, title="Inverse Laplace"):
    x_total=0.0
    for key, cutt in cutts.items():
        p_values=path[cutt[0]:cutt[1]]
        y_values=G_values[cutt[0]:cutt[1]]*np.exp(p_values*L)

        if key == 'low_side':
            x_values=np.real(p_values)
            sgn=1.0/(1j*2*np.pi)
        elif key == 'top_side':
            # Need to reverse to give increasing x_values
            x_values=-np.real(p_values)
            sgn=-1.0/(1j*2*np.pi)
        else:
            x_values=np.imag(p_values)
            sgn=1.0/(2*np.pi) # the 1j/1j cancels

        x_values = x_values-x_values[0]+x_total
        x_total=x_values[-1]
        plt.plot(x_values, np.real(y_values), '.-')
        plt.plot(x_values, np.real(G_values[cutt[0]:cutt[1]]), '-')
    plt.title(title)
    plt.show()


def test_inv_laplace():
    p_values=make_path(factor=100)
    nums = p_values*0
    aa = -0.1
    for ii, p in enumerate(p_values):
        nums[ii] = 1.0/(p-aa)

    N_set = {0.05, 0.1, 0.5, 1.0, 10}
    for iN, N in enumerate(N_set):
        print("Correct= %f, Calc=%f"%(np.exp(aa*N),
                                      invLaplace_path(nums, N)))
        plot_invLaplace_path(nums, N, title="N="+str(N))
#test_inv_laplace()

def test_scale():
    N_set = {0.05, 0.1, 0.5, 1.0, 10}
    for iN, N in enumerate(N_set):
        p_values=make_path(factor=100, scale=1.0/N)
        nums = p_values*0
        aa = -0.1
        for ii, p in enumerate(p_values):
            nums[ii] = 1.0/(p-aa)

        print("Correct= %f, Calc=%f"%(np.exp(aa*N),
                                      invLaplace_path(nums, N)))
        plot_invLaplace_path(nums, N, title="N="+str(N))
#test_scale()
