import numpy as np
from . import rouse


def generate_poisson_homologs(mu, chr_size=1):
    r"""
    Generate the number and location of linkage between homologous chromosomes

    Parameters
    ----------
    mu : float
        Average number of (uniformly distributed) linkages between chromosomes.
    chr_size : float
        Size of the chromosome.

    Returns
    -------
    cell : np.array of float
        List of linkage locations between the homologous chromosomes.

    """
    cell = np.sort(np.random.random_sample(size=np.random.poisson(lam=mu)))
    return chr_size*cell


def mscd(t, linkages, label_loc, chr_size, nuc_radius, b, D,
         num_modes=rouse._default_modes):
    r"""
    Calculate the MSCD for the model of linked chromosomes

    Parameters
    ----------
    t : float array
        Time in seconds
    linkages : float array
        List of the link positions between the homologous chromosomes
    label_loc : float
        Location of the fluorescent label along the chromosome (microns)
    chr_size : float
        Length of the chromosome (microns)
    nuc_radius : float
        Radius of the nucleus (microns)
    b : float
        Kuhn length (microns)
    D : float
        Diffusivity (microns ** 2 / second)
    num_modes : int
        Number of normal modes used in the calculation

    Returns
    -------
    mscd_model : float array (size len(t))
        Calculated MSCD (microns ** 2) for the model with defined linkages

    """
    linkages = np.array(linkages) / b
    chr_size = chr_size / b
    label_loc = label_loc / b

    # Evaluate the MSCD if there are no linkages
    if len(linkages) == 0:
        mscd_model = 2 * rouse.linear_mid_msd(t, b, chr_size, D, num_modes)
        return np.minimum(mscd_model, nuc_radius**2)

    # Evaluate the MSCD if there are linkages between the chromosomes
    i = np.searchsorted(linkages, label_loc)
    if i == 0:
        mscd_func = rouse.linear_mscd
        Ndel = linkages[0] - label_loc
        N = 2 * linkages[0]
    elif i == len(linkages):
        mscd_func = rouse.linear_mscd
        Ndel = label_loc - linkages[-1]
        N = 2 * (chr_size - linkages[-1])
    else:
        mscd_func = rouse.ring_mscd
        Ndel = linkages[i] - label_loc
        N = 2*(linkages[i] - linkages[i - 1])

    mscd_model = mscd_func(t, D=D, Ndel=Ndel, N=N, b=b, num_modes=num_modes)

    return np.minimum(nuc_radius**2, mscd_model)


def mscd_plateau(linkages, label_loc, chr_size, nuc_radius, b=1):
    r"""
    Evaluate the plateau values in the MSCD

    Parameters
    ----------
    linkages : float array
        List of the link positions between the homologous chromosomes
    label_loc : float
        Location of the fluorescent label along the chromosome (microns)
    chr_size : float
        Length of the chromosome (microns)
    nuc_radius : float
        Radius of the nucleus (microns)
    b : float
        Kuhn length (microns)

    Returns
    -------
    mscd_plateau : float
        Plateau value of the MSCD in the long-time asymptotic limit
        (microns ** 2).

    """
    chr_size /= b
    label_loc /= b
    linkages = np.array(linkages) / b

    # Evaluate the MSCD if there are no linkages
    if len(linkages) == 0:
        mscd_plateau = nuc_radius ** 2
        return mscd_plateau

    # Evaluate the MSCD if there are linkages between the chromosomes
    i = np.searchsorted(linkages, label_loc)
    if i == 0:
        Ndel = linkages[0] - label_loc
        mscd_plateau = 4 * b ** 2 * Ndel
    elif i == len(linkages):
        Ndel = label_loc - linkages[-1]
        mscd_plateau = 4 * b ** 2 * Ndel
    else:
        Ndel = linkages[i] - label_loc
        N = 2*(linkages[i] - linkages[i - 1])
        mscd_plateau = 2 * b ** 2 / (1 / (2 * Ndel) + 1 / (N - 2 * Ndel))

    return np.minimum(mscd_plateau, nuc_radius**2)
