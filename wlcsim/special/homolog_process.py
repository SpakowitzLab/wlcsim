""""Simulation" of two yeast sister chromatids forming a synaptonemal
complex"""
import numpy as np
from numba import jit

from ..analytical.rouse import gaussian_Ploop

UNLOOP = 0
PAIRED = 1
JUNCTION = 2

# @jit(nopython=True)
def Ploop_given_paired(n, starts, ends, b, dN, a):
    r"""Get looping probability for a "homologously paired" pair of Rouse
    chains.

    Parameters
    ----------
    n : int
        number of beads in "state"
    starts/ends : np.array<int>
        output of get_terminatable_loops
    b : float
        Kuhn length of linear chain
    dN : float
        number of kuhn length between each pair of beads
    a : float
        capture radius of looping probability calculation

    Returns
    -------
    ploop : (n,) np.array<float>
        probability that a site will be looped given conformation specified.
        sites not inside the starts/ends dilineators are given 0 probability of
        looping

    Notes
    -----
    Uses the fact that two loci on a ring polymer (Kuhn length `b` and
    separations `N`, `L-N`) can be thought of as being connected by an
    effective linear polymer.

    We simply call $k_1 \propto 1/Nb^2$ and $k_2 \propto 1/(L-N)b^2$ to be the
    spring constants of the two arms of the loop, then they form one combined
    spring additively, with $k_\text{tot} = k_1 + k_2$, meaning that any
    linaer polymer of length $\hat{N}$ and Kuhn length $\hat{b}$ will have the
    same looping probability at its ends as our two loci, as long as

    $$ \hat{N}\hat{b}^2 = \frac{1}{1/Nb^2 + 1/(L-N)b^2}. $$
    """
    ploop = np.zeros((n,))
    if len(starts) == 0:
        return ploop
    # for each loop site, we compute the looping probabilities to its left
    for i, si in enumerate(starts):
        ei = ends[i]
        # the current "loop" is the whole polymer
        if si == -1 and ei == n:
            raise NotImplementedError("there should be at least one paired site")
        # the current "loop" is the free tails at the start of the polymer
        if si == -1:
            ploop[:ei] = gaussian_Ploop(a, 2*dN*(ei - np.arange(ei)), b)
            continue
        # the current "loop" is the free tails at the end of the polymer
        if ei == n:
            tail_len = (n - 1) - si
            ploop[si+1:] = gaussian_Ploop(a,
                    2*dN*(tail_len - np.arange(tail_len)[::-1]), b)
        # the current "loop" is in between si and ei, with total length
        L = 2 * (ei - si)
        # length of the one arm of the loop connecting each homolog pair
        N = 2*np.arange(1, ei - si)
        Nb2_hat = 1 / ( 1/N/b/b + 1/(L-N)/b/b )
        ploop[si+1:ei] = gaussian_Ploop(a, Nb2_hat, 1)
    return ploop

# @jit(nopython=True)
def get_terminatable_loops(state):
    """Take a state array and get the ends of the loops that have yet to
    terminate, so that we can count how many termination reactions are possible
    given the current state."""
    starts = []
    ends = []
    num_beads = len(state)
    # track whether the bead we're at is in the interior of a "loop", i.e. not
    # paired, i.e. whether state[i] == 0
    in_loop = state[0] == UNLOOP
    # track most recent junction, since an unjunctioned pair site can still
    # loop with it
    prev_junc = -1
    # special flag for if the end of the polymer is not looped, treat as though
    # an imaginary "-1" bead is paired. otherwise ensure we start with one more
    # starts than ends
    starts.append(-1 if in_loop else 0)
    # keep track of if we're at loci that have been "zippered" into the
    # synaptonemal complex yet (i.e. those that rest between related "junction"
    # sites (in a string that matches the regex "20*2"))
    in_junction = state[0] == JUNCTION
    for i, si in enumerate(state):
        if i == 0:
            continue
        if not in_junction and si == PAIRED:
            ends.append(i)
            starts.append(i) # now N_ends = N_starts - 1
        if not in_junction and si == JUNCTION:
            ends.append(i) # now N_ends == N_starts
            prev_junc = i
            in_junction = True
        elif not in_junction and si == UNLOOP:
            continue
        elif in_junction and si == PAIRED:
            starts.append(prev_junc)
            ends.append(i)
            starts.append(i) # now N_ends = N_starts - 1
            in_junction = False
        elif in_junction and si == JUNCTION:
            prev_junc = i
            continue
        elif in_junction and si == UNLOOP:
            continue
    if si == PAIRED:
        starts.pop() # unecessarily added a "starts"
    elif in_junction and si == UNLOOP:
        # the end can still loop with the prev junction site
        starts.append(prev_junc)
        # special flag for if the other end of the polymer is not looped, treat
        # as though an imaginary "N+1"th bead is paired
        ends.append(num_beads)
    elif not in_junction and si == UNLOOP:
        ends.append(num_beads)
    elif not in_junction and si == JUNCTION:
        pass
    elif in_junction and si == JUNCTION:
        pass
    return np.array(starts), np.array(ends)

def test_get_terminatable_loops():
    starts, ends = get_terminatable_loops(np.array([0, 0, 1, 1, 0, 1, 2, 2, 0, 0, 1]))
    act_starts = np.array([-1, 2, 3, 5, 7])
    act_ends = np.array([2, 3, 5, 6, 10])
    assert(np.all(act_starts == starts))
    assert(np.all(act_ends == ends))
    starts, ends = get_terminatable_loops(np.array([1, 0, 1, 1, 0, 1, 2, 2, 0, 0, 0]))
    act_starts = np.array([0, 2, 3, 5, 7])
    act_ends = np.array([2, 3, 5, 6, 11])
    assert(np.all(act_starts == starts))
    assert(np.all(act_ends == ends))

# @jit(nopython=True)
def chromatid_gillespie(init_paired, b, N, k_h, k_u, a,
                        tend=np.inf, target_completion=1, k_t=None, R=None,
                        save_times=None):
    """Do Gillespie sim of homolog pairing process.

    Parameters
    ----------
    init_paired : array_like of bool
        Which loci should be initially paired. length of array determines
        number of "beads". for now, at least one locus should be paired.
    b : float
        Kuhn length
    N : float
        total length of polymer in number of Kuhn lengths
    k_h : float
        rate of stable homolog junctions given loop exists
    a : float
        capture radius inside which stable homolog junction can form
    tend : float, optional
        what the cutoff time for the simulation should be. `np.inf` by default,
        where no limit implies `target_completion` will be the only termination
        condition.
    target_completion : float, optional
        $\in[0,1]$ what fraction of beads should be marked as "completed"
        before the simulation ends. since there is not way to undo
        "synaptonemal formation" reactions, it does not make sense to continue
        the simulation after this fraction reaches 1, as nothing can happen. 1
        by default.
    k_t : float, optional
        rate of synaptonemal junction formation (i.e. "termination" reaction),
        equal to k_h by default
    R : float, optional
        the radius of the spherical confinement. only used if there are no
        initial pairings, in order to have a non-zero probability that the
        first reaction occurs. once the first reaction occurs, it is assumed
        that this "random collision" reaction will never happen again

    Notes
    -----
    `state` is an int array of 0's, 1's, and 2's. 0 means unpaired, 1 means
    paired, and 2 means the locus has been sequestered into the synaptonemal
    complex
    """
    if R is not None or not np.any(init_paired):
        raise NotImplementedError("Initial pairing state must be specified.")
    state = np.array(init_paired).astype(int)
    if k_t is None:
        k_t = k_h
    if a <= 0 or tend < 0 or k_h < 0 or k_t < 0 or N <= 0 or b <= 0 \
       or np.any(state < 0) or np.any(state > 2):
        raise ValueError("One of the inputs makes no sense.")
    if save_times is None:
        save_times = [0, np.inf]
    num_beads = len(state)
    dN = N/(num_beads-1)
    t = 0
    save_i = 0
    save_states = []
    # the only loci that can "pair" are the ones that are inside of loops
    # that have not undergone synaptonemal junction formation, so get the
    # boundaries of those loops
    term_starts, term_ends = get_terminatable_loops(state)
    while t < tend and len(term_starts) > 0:
        while t > save_times[save_i]:
            save_states.append(state.copy())
            save_i += 1
        # the synaptonemal junction formation is also a looping process, simply
        # of two loci on exactly opposite ends of the "loop", which is the same
        # as being on a linear polymer with half that separation
        rates_per_term = k_t*gaussian_Ploop(a, dN*(term_ends - term_starts)/2, b)
        # deal with special "end" cases, corresponding to free ends looping
        rates_per_term[term_starts == -1] = k_t*gaussian_Ploop(a, term_ends[term_starts == -1], b)
        rates_per_term[term_ends == num_beads] = k_t*gaussian_Ploop(a, term_starts[term_ends == num_beads], b)
        term_cumulat = np.cumsum(rates_per_term)
        term_tot = term_cumulat[-1]
        # the other two reactions are simple looping and unlooping
        # only non-terminated pairs can unloop
        rates_per_site = (state==1).astype(float)*k_u
        # get looping-weighted reaction rates for all loci that can loop
        rates_per_site += k_h*(state==0)*Ploop_given_paired(num_beads, term_starts, term_ends, b, dN, a)
        site_cumulat = np.cumsum(rates_per_site)
        site_tot = site_cumulat[-1]
        # first select which reaction to do
        if np.random.rand() < term_tot / (term_tot + site_tot):
            ti = np.searchsorted(term_cumulat, term_tot*np.random.rand())
            dt = np.random.exponential(scale=1/term_tot)
            # a termination reaction marks from locus to the next looped site
            # as "done", but make sure to not write out of bounds for two
            # special cases
            si = max(term_starts[ti], 0)
            ei = min(term_ends[ti]+1, num_beads)
            state[si:ei] = 2
        else:
            reacted_site = np.searchsorted(site_cumulat, site_tot*np.random.rand())
            # a looping reaction just switches the state of that locus
            state[reacted_site] = 1 - state[reacted_site]
        # even though we selected the reaction first, we just still step time
        # forward using the total cumulative rate for all reactions, since any
        # of them were in fact "possible"
        dt = np.random.exponential(scale=1/(site_tot+term_tot))
        # update time and our loop list
        t += dt
        term_starts, term_ends = get_terminatable_loops(state)
    if np.inf in save_times:
        save_states.append(state.copy())
    return save_states
