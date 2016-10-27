import numpy as np
import numba

bsave = []
dfsave = []
xsave = []


def init_gaussian_in_sphere_naive(N, k, r):
    """
    Initilize gaussian chain for gaussian_in_spherical_confinement, as
    described by its documentation.
    """
    # distance between two points on a circle of angle $\theta$ apart is just
    # $d(\theta) = 2r\sin(\theta/2)$ by a simple geometric argument.

    # In the rouse model we define the spring constant between beads $k =
    # \frac{3k_{B}T}{b^2}$, so that the characteristic distance between beads
    # will turn out to be $b$.

    # Thus, we can imagine that a reasonable initial set of locations for our
    # beads might be a distance $b$ apart, which we can guarantee by putting
    # them on a circle of radius $r/2$ in the xy-plane, and simply incrementing
    # the angle between adjacent beads as $d^{-1}(b)$
    d = np.sqrt(3/k)
    R = r/2
    dtheta = 2*np.arcsin(d/(2*R))
    theta = dtheta*np.arange(N)
    x0 = np.zeros((N, 3))
    x0[:,0] = (r/2)*np.cos(theta)
    x0[:,1] = (r/2)*np.sin(theta)
    return x0


#@numba.jit
def gaussian_in_spherical_confinement_step(x, N, D, k, t, r):
    """
    Actual simulation code that runs between save interavals for
    gaussian_in_spherical_confinement as described by its documentation.
    """
    dt = np.diff(t)
    for dti in dt:
        dB = np.sqrt(2*dti*D)*np.random.randn(N, 3)
        df = np.zeros_like(x)
        for i in range(N):
            if i == 0:
                df[i,:] = k*(x[i+1,:] - x[i,:])
            elif i == N-1:
                df[i,:] = k*(x[i-1,:] - x[i,:])
            else:
                df[i,:] = k*(x[i+1,:] - 2*x[i,:] + x[i-1,:])
        xsave.append(x)
        dfsave.append(df)
        bsave.append(dB)
        x = x + D*df*dti + dB
        return x

def gaussian_in_spherical_confinement(N, D, k, t, r,
                                      outfile='gaus_spher_sim',
                                      save_interval=1):
    """
    Simulate a gaussian chain with N links with spring constant 'k' through
    discrete times t, within a spherical confinement of radius r, with time
    evolution such that the diffusion coefficient of an individual bead is D.
    Units of $k_{B}T$ are assumed.

    Saves every save_interval values of t to the output file.

    For now, the initialization is a naive laying of the N beads at their
    characteristic distance apart perfectly circularly in a circle in the
    xy-plane of radius r/2. For this reason, 3/k < r/2 is required.
    """
    x = init_gaussian_in_sphere_naive(N, k, r)
    # with open(outfile, 'wb') as f:
        # np.savetxt(f, x)
    with open(outfile + '.out', 'w') as f:
        x.tofile(f)
    Nt = len(t)
    with open(outfile + '.shape', 'wb') as f:
        np.savetxt(f, np.array([N, 3, Nt]), fmt='%d')
    with open(outfile + '.times', 'w') as f:
        t.tofile(f)
    t0 = 0
    tf = min(save_interval, Nt - 1)
    # while we have t's left to iterate through, iterate through them in blocks
    # of size "save_interval"
    while True:
        x = gaussian_in_spherical_confinement_step(x, N, D, k, t[t0:tf+1], r)
        # with open(outfile, 'ab') as f:
            # np.savetxt(f, x)
        with open(outfile + '.out', 'a') as f:
            x.tofile(f)
        t0 = tf
        tf = min(t0 + save_interval, Nt)
        if tf > Nt - 1:
            break







