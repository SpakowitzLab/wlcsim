import numpy as np
import os

mod_file = os.path.abspath(__file__)
mod_path = os.path.dirname(mod_file)

def end2end_distance(r, lp, N, L):
    """For now, always returns values for r = np.linspace(0, 1, 50001).
    TODO: fix this

    The end to end distance for a polymer with given parameters.
    end2end_distance(r: array_like(n), lp: float, N: int, L: float):
        -> array_like(n)
    Params:
        r - the values at which to evaluate the end-to-end probability
            distribution
        lp - persistence length of polymer
        N - number of beads in polymer
        L - polymer length
    Output:
        (x,g)
        x = np.linspace(0, 1, 5001)
        g = P(|R| = x | lp, N, L)

    Uses the gaussian chain whenever applicable, ssWLC tabulated values
    otherwise. If you request parameters that require a WLC end-to-end
    distance, the function will ValueError."""
    Delta = L/(lp*(N-1)) # as in wlcsim code
    # WLC case
    actual_r = np.linspace(0, 1, 5001)
    if Delta < 0.01: # as in wlcsim code
        ValueError('end2end_distance: doesn\'t know how to compute WLC case!')
    # ssWLC case
    elif Delta < 10: # as in wlcsim code
        Eps = N/(2*lp*L) # as in wlcsim code
        file_num = round(100*Eps) # index into file-wise tabulated values
        file_name = os.path.join('pdata', 'out' + str(file_num) + '.txt')
        file_name = os.path.join(mod_path, file_name)
        G = np.loadtxt(file_name)
        return (actual_r, G)
    # GC case
    else:
        return (actual_r, end2end_distance_gauss(actual_r, lp=lp, N=N, L=L))

def end2end_distance_gauss(r, lp, N, L):
    r2 = np.power(r, 2)
    return 3.0*r2*sqrt(6/np.pi)*power(N/(2*lp*L), 1.5) \
            *np.exp(-(3/2)*(N/(2*lp*L))*r2)



