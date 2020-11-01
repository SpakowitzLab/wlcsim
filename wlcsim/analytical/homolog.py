import numpy as np
import scipy.integrate
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


def get_neighbors(linkages, label_loc):
    i = np.searchsorted(linkages, label_loc)
    neighbors = [None, None]
    if i > 0:
        neighbors[0] = linkages[i - 1]
    if i < len(linkages):
        neighbors[1] = linkages[i]
    return neighbors


def is_linkless(neighbors):
    return neighbors[0] is None and neighbors[1] is None


def is_linear(neighbors):
    # xor
    return (neighbors[0] is None) != (neighbors[1] is None)


def is_ring(neighbors):
    return neighbors[0] is not None and neighbors[1] is not None


def P_is_linkless(mu):
    return np.exp(-mu)


def P_is_ring(mu, l0, L):
    return (1 - np.exp(-(mu / L)*l0)) * (1 - np.exp(-(mu / L)*(L - l0)))


def P_is_linear(mu, l0, L):
    return 1 - P_is_ring(mu, l0, L) - P_is_linkless(mu)


def get_N(neighbors, L):
    if is_linkless(neighbors):
        return None
    neighbors = neighbors.copy()
    if neighbors[0] is None:
        neighbors[0] = 0
    if neighbors[1] is None:
        neighbors[1] = L
    return 2*(neighbors[1] - neighbors[0])


def get_N_ring(neighbors):
    if not is_ring(neighbors):
        return None
    return 2*(neighbors[1] - neighbors[0])


def get_N_linear(neighbors, L):
    if not is_linear(neighbors):
        return None
    return get_N(neighbors, L)


def P_has_left_only(mu, l0, L):
    lam = mu / L
    return (1 - np.exp(-lam*l0)) * np.exp(-lam*(L - l0))


def P_has_right_only(mu, l0, L):
    lam = mu / L
    return (1 - np.exp(-lam*(L - l0))) * np.exp(-lam*l0)


def P_N_linear_left(N, mu, l0, L):
    lam = mu / L
    Z_L = (1 - np.exp(-lam*l0)) / lam
    # (1/2)* from jacobian
    return np.zeros_like(N) + (1/2) * (N <= 2*L) \
        *(N >= 2*(L - l0)) * (1/Z_L)*np.exp(-lam*(N/2 - (L - l0)))


def P_N_linear_right(N, mu, l0, L):
    lam = mu / L
    Z_R = (1 - np.exp(-lam*(L - l0))) / lam
    # (1/2)* from jacobian
    return np.zeros_like(N) + (1/2) * (N <= 2*L) \
        * (N >= 2*l0) * (1/Z_R)*np.exp(-lam*(N/2 - l0))


def P_N_linear(N, mu, l0, L):
    N = np.array(N)
    lam = mu / L
    left = P_has_left_only(mu, l0, L)
    right = P_has_right_only(mu, l0, L)
    return (left*P_N_linear_left(N, mu, l0, L)
        + right*P_N_linear_right(N, mu, l0, L)
    ) / (left + right)


def _inside_interval_weights(N, l0, L):
    line_lengths = (
        (L >= N/2)*(N/2 >= L - l0) * (L - N/2)
      + (L - l0 > N/2)*(N/2 >= l0) * l0
      + (l0 > N/2)*(N/2 >= 0) * N/2
    )
#     total_line_length = np.sqrt(2) * l0 * (L + l0)
    return (0 <= N/2) * (N/2 <= L) * line_lengths


def P_N_ring(N, mu, l0, L):
    N = np.array(N)
    lam = mu / L
    Z_LR = (np.exp(lam*L) - np.exp(lam*l0))*(np.exp(lam*l0) - 1) \
         / (lam**2 * np.exp(lam*(L + l0)))
    return np.zeros_like(N) + (1/2) * (N <= 2*L) *(
        (1/Z_LR)*np.exp(-lam*N/2) * _inside_interval_weights(N, l0, L)
    )


def P_N(N, mu, l0, L):
    N = np.array(N)
    return P_is_linkless(mu)*np.zeros_like(N) \
        + P_is_ring(mu, l0, L)*P_N_ring(N, mu, l0, L) \
        + P_is_linear(mu, l0, L)*P_N_linear(N, mu, l0, L)


def get_Ndel(neighbors, label_loc):
    if neighbors[0] is not None:
        return label_loc - neighbors[0]
    elif neighbors[1] is not None:
        return neighbors[1] - label_loc
    else:
        return None


def get_Ndel_ring(neighbors, label_loc):
    if not is_ring(neighbors):
        return None
    return label_loc - neighbors[0]


def get_Ndel_linear(neighbors, label_loc):
    if not is_linear(neighbors):
        return None
    if neighbors[0] is not None:
        return label_loc - neighbors[0]
    return neighbors[1] - label_loc


def P_Ndel_ring(Ndel, mu, l0, L):
    Ndel = np.array(Ndel)
    lam = mu / L
    # if ring, then it has left-side neighbor
    Z_L = (1 - np.exp(-lam*l0)) / lam
    return np.zeros_like(Ndel) + (Ndel >= 0)*(Ndel <= l0)*(1/Z_L)*np.exp(-lam*Ndel)

def P_Ndel_linear_left(Ndel, mu, l0, L):
    lam = mu / L
    Z_L = (1 - np.exp(-lam*l0)) / lam
    return np.zeros_like(Ndel) + (Ndel >= 0) \
        * (Ndel <= l0) * (1/Z_L)*np.exp(-lam*Ndel)

def P_Ndel_linear_right(Ndel, mu, l0, L):
    lam = mu / L
    Z_R = (1 - np.exp(-lam*(L - l0))) / lam
    return np.zeros_like(Ndel) + (Ndel >= 0) \
        * (Ndel <= L - l0) * (1/Z_R)*np.exp(-lam*Ndel)

def P_Ndel_linear(Ndel, mu, l0, L):
    Ndel = np.array(Ndel)
    lam = mu / L
    left = P_has_left_only(mu, l0, L)
    right = P_has_right_only(mu, l0, L)
    return (
        left*P_Ndel_linear_left(Ndel, mu, l0, L)
        + right*P_Ndel_linear_right(Ndel, mu, l0, L)
    ) / (left + right)

def P_Ndel(Ndel, mu, l0, L):
    Ndel = np.array(Ndel)
    return P_is_linkless(mu)*np.zeros_like(Ndel) \
        + P_is_ring(mu, l0, L)*P_Ndel_ring(Ndel, mu, l0, L) \
        + P_is_linear(mu, l0, L)*P_Ndel_linear(Ndel, mu, l0, L)

def _P_Ndel(Ndel, mu, l0, L):
    """Alternative to P_Ndel for testing."""
    Ndel = np.array(Ndel)
    lam = mu / L
    has_left = (1 - np.exp(-lam*l0))
    Z_L = (1 - np.exp(-lam*l0)) / lam
    Z_R = (1 - np.exp(-lam*(L - l0))) / lam
    Ndel_left = (1/Z_L)*np.exp(-lam*Ndel)
    Ndel_right = (1/Z_R)*np.exp(-lam*Ndel)
    return (Ndel >= 0) * (
        (Ndel <= l0) * has_left * Ndel_left
      + (Ndel <= L - l0) * (1 - has_left) * Ndel_right
    )


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
        return nuc_radius ** 2

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



def P_N_linear_left(N, mu, l0, L):
    lam = mu / L
    Z_L = (1 - np.exp(-lam*l0)) / lam
    # (1/2)* from jacobian
    return np.zeros_like(N) + (1/2) * (N <= 2*L) \
        *(N >= 2*(L - l0)) * (1/Z_L)*np.exp(-lam*(N/2 - (L - l0)))


def P_N_linear_right(N, mu, l0, L):
    lam = mu / L
    Z_R = (1 - np.exp(-lam*(L - l0))) / lam
    # (1/2)* from jacobian
    return np.zeros_like(N) + (1/2) * (N <= 2*L) \
        * (N >= 2*l0) * (1/Z_R)*np.exp(-lam*(N/2 - l0))


def P_N_joint_Ndel_ring(N, Ndel, mu, l0, L):
    """Returns (N, Ndel)."""
    N = np.array(N)
    Ndel = np.array(Ndel)
    lam = mu / L
    return (1/2) * (0 <= N) * (N <= 2*L) \
        * (N/2 - (L - l0) <= Ndel) \
        * (Ndel <= N/2) * (Ndel <= l0) \
        * lam*np.exp(-lam*(N/2 - Ndel)) / (1 - np.exp(-lam*(L - l0))) \
        * lam*np.exp(-lam*Ndel) / (1 - np.exp(-lam*l0))


def mscd_plateau_ensemble(mu, label_loc, L, b, nuc_radius):
    """
    Fully analytical computation of average plateau value over population.

    All parameters should be specified in the same units (of length).

    Parameters
    ----------
    mu : float
        Average number of linkages across entire chromosome.
    label_loc : float
        The location along the chromosome where the label being tracked is.
    L : float
        The full length of each of the two homologous chromosomes.
    b : float
        The Kuhn length of the polymer.
    nuc_radius : float
        The radius of the spherical confinement.

    """
    def lower(Ndel):
        return 2*Ndel
    def upper(Ndel):
        return 2*(Ndel + L - label_loc)
    def ring_plateau_f(N, Ndel, b, nuc_radius, *args):
        return np.minimum(
            nuc_radius**2,
            2 * b**2 / (1 / (2 * Ndel) + 1 / (N - 2*Ndel))
        ) * P_N_joint_Ndel_ring(N, Ndel, *args)
    def right_linear_plateau(Ndel, b, nuc_radius, *args):
        return np.minimum(4 * b**2 * Ndel, nuc_radius**2) \
                * P_Ndel_linear_right(Ndel, *args)
    def left_linear_plateau(Ndel, b, nuc_radius, *args):
        return np.minimum(4 * b**2 * Ndel, nuc_radius**2) \
                * P_Ndel_linear_left(Ndel, *args)
    ring_ave = scipy.integrate.dblquad(
        ring_plateau_f, 0, label_loc, lower, upper,
        args=(b, nuc_radius, mu, label_loc, L)
    )[0]  # throw away error estimate
    left_linear_ave = scipy.integrate.quad(
        left_linear_plateau, 0, label_loc,
        args=(b, nuc_radius, mu, label_loc, L)
    )[0]
    right_linear_ave = scipy.integrate.quad(
        right_linear_plateau, 0, L - label_loc,
        args=(b, nuc_radius, mu, label_loc, L)
    )[0]
    return nuc_radius**2 * P_is_linkless(mu) \
        + ring_ave * P_is_ring(mu, label_loc, L) \
        + left_linear_ave * P_has_left_only(mu, label_loc, L) \
        + right_linear_ave * P_has_right_only(mu, label_loc, L)
