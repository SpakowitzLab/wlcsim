import numpy as np
from . import rouse

def generate_homolog_poisson(lam, L):
    """
    Randomly sample a homologously-linked configuration.

    We model homologous linkages as forming via a Poisson process with rate *lam*.

    Parameters
    ----------
    lam : float
        Poisson parameter, average number of connection.
    L : float
        Length of the chromosome.

    Returns
    -------
    array_like of float
        Location of each linkage.
    """
    num_links = np.random.poisson(lam=lam)
    # old trick for generating poisson samples. uniform random then sort
    linkages = L*np.random.random_sample(size=num_links)
    linkages.sort()
    return linkages

def mscd_homolog_locus(t, loc, linkages, L, D, **kwargs):
    i = np.searchsorted(linkages, loc)
    if i == 0:
        return rouse.linear_mscd(t, D, linkages[0] - loc, 2*linkages[0],
                                 **kwargs)
    elif i == len(linkages):
        return rouse.linear_mscd(t, D, L - linkages[-1],
                                 2*(loc - linkages[-1]), **kwargs)
    else:
        return rouse.ring_mscd(t, D, linkages[i] - loc,
                               2*(linkages[i] - linkages[i - 1]), **kwargs)
