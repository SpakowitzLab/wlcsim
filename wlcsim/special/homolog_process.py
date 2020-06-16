""""Simulation" of two yeast sister chromatids forming a synaptonemal
complex"""
import numpy as np
from numba import jit

from ..analytical.rouse import gaussian_Ploop

UNLOOP = 0
PAIRED = 1
JUNC_LEFT_END = 2
JUNC_RIGHT_END = 4
JUNC_BOTH = 6

@jit(nopython=True)
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
            #TODO: use conf gaus green's function to get actual "free" ploop
            pass # assume floor ploop will be set outside of here
        # the current "loop" is the free tails at the start of the polymer
        elif si == -1:
            ploop[:ei] = gaussian_Ploop(a, 2*dN*(ei - np.arange(ei)), b)
        # the current "loop" is the free tails at the end of the polymer
        elif ei == n:
            tail_len = (n - 1) - si
            ploop[si+1:] = gaussian_Ploop(a,
                    2*dN*(tail_len - np.arange(tail_len)[::-1]), b)
        # the current "loop" is in between si and ei, with total length
        else:
            L = 2*dN*(ei - si)
            # length of the one arm of the loop connecting each homolog pair
            N = 2*dN*np.arange(1, ei - si)
            Nb2_hat = 1 / ( 1/N/b/b + 1/(L-N)/b/b )
            ploop[si+1:ei] = gaussian_Ploop(a, Nb2_hat, 1)
    return ploop

@jit(nopython=True)
def get_terminatable_loops(state):
    """Take a state array and get the ends of the loops that have yet to
    terminate, so that we can count how many termination reactions are possible
    given the current state."""
    starts = []
    ends = []
    num_beads = len(state)
    # keep track of if we're at loci that have been "zippered" into the
    # synaptonemal complex yet (i.e. those that rest between related "junction"
    # sites (in a string that matches the regex "[26]0*[46]"))
    in_junction = state[0] == JUNC_LEFT_END
    if not in_junction:
        # check if the end of the polymer is "free" or in a loop/junction
        in_loop = state[0] == UNLOOP
        # special flag for if the end of the polymer is not looped, treat as though
        # an imaginary "-1" bead is paired. otherwise ensure we start with one more
        # starts than ends
        starts.append(-1 if in_loop else 0) # else PAIRED
    for i, si in enumerate(state):
        if i == 0:
            continue
        elif not in_junction and si == PAIRED:
            ends.append(i)
            starts.append(i) # now N_ends = N_starts - 1
        elif not in_junction and si == JUNC_LEFT_END:
            ends.append(i) # now N_ends == N_starts
            in_junction = True
        elif not in_junction and (si == JUNC_BOTH or si == JUNC_RIGHT_END):
            raise ValueError("How did we end up hitting right end of junction first?")
        elif not in_junction and si == UNLOOP:
            continue
        elif in_junction and si == PAIRED:
            raise ValueError("No loops should form in between junction points.")
        elif in_junction and si == JUNC_LEFT_END:
            raise ValueError("How did we hit left end of junction from inside junction?")
        elif in_junction and si == JUNC_BOTH:
            continue
        elif in_junction and si == JUNC_RIGHT_END:
            starts.append(i) # now N_ends = N_starts - 1
            in_junction = False
        elif in_junction and si == UNLOOP:
            continue
    if si == PAIRED or si == JUNC_RIGHT_END:
        starts.pop() # unecessarily added a "starts"
    elif not in_junction and si == UNLOOP:
        ends.append(num_beads)
    elif in_junction and si == UNLOOP:
        raise ValueError("What is the junction with if the end bead is UNLOOP?")
    # other impossibilities
    # checked above: (not in_junction and si == JUNC_RIGHT_END) or (in_junction
    # and si == PAIRED)
    elif si == JUNC_LEFT_END or si == JUNC_BOTH:
        raise ValueError("Rightmost bead cannot be left end of junction.")
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
                        tend=np.inf, target_completion=1, k_t=None, R=np.inf,
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
    state = np.array(init_paired).astype(int)
    if k_t is None:
        k_t = k_h
    if a <= 0 or tend < 0 or k_h < 0 or k_t < 0 or N <= 0 or b <= 0 \
       or R <= 0 or np.any(state < 0) or np.any(state > 2):
        raise ValueError("One of the inputs makes no sense.")
    if np.isinf(R) and k_u > 0:
        raise ValueError("You probably don't want unbinding reactions if " \
                         "there's no way for a totally unpaired polymer to " \
                         "find a loop (aka R=np.inf")
    if save_times is None:
        save_times = [0, np.inf]
    num_beads = len(state)
    dN = N/(num_beads-1)
    min_ploop = a**3 / R**3 # fraction of volume occupied by interaction sphere
    t = 0
    save_i = 0
    save_states = []
    # the only loci that can "pair" are the ones that are inside of loops
    # that have not undergone synaptonemal junction formation, so get the
    # boundaries of those loops
    term_starts, term_ends = get_terminatable_loops(state)
    while t < tend and len(term_starts) > 0:
        while t >= save_times[save_i]:
            save_states.append(state.copy())
            save_i += 1
        if len(term_ends) > 1:
            # the synaptonemal junction formation is also a looping process, simply
            # of two loci on exactly opposite ends of the "loop", which is the same
            # as being on a linear polymer with half that separation
            ploop = gaussian_Ploop(a, dN*(term_ends - term_starts)/2, b)
            rates_per_term = k_t*np.maximum(min_ploop, ploop)
            # deal with special "end" cases, corresponding to free ends looping
            left_ploop = gaussian_Ploop(a, (num_beads-1) - term_ends[term_starts == -1], b)
            rates_per_term[term_starts == -1] = k_t*np.maximum(min_ploop, left_ploop)
            right_ploop = gaussian_Ploop(a, term_starts[term_ends == num_beads], b)
            rates_per_term[term_ends == num_beads] = k_t*np.maximum(min_ploop, right_ploop)
        else: # special case of full polymer unlooped
            ploop = gaussian_Ploop(a, dN*num_beads, b)
            rates_per_term = k_t*np.maximum(min_ploop, ploop)
        # make array to index into for gillespie reaction "choice"
        term_cumulat = np.cumsum(rates_per_term)
        term_tot = term_cumulat[-1]
        # the other two reactions are simple looping and unlooping
        # only non-terminated pairs can unloop
        rates_per_site = k_u*(state==1).astype(float)
        # get looping-weighted reaction rates for all loci that can loop
        rates_per_site += k_h*(state==0)*Ploop_given_paired(num_beads, term_starts, term_ends, b, dN, a)
        site_cumulat = np.cumsum(rates_per_site)
        site_tot = site_cumulat[-1]
        # first select which reaction to do
        if np.random.rand() < term_tot / (term_tot + site_tot):
            # choose which of the loops will close
            ti = np.searchsorted(term_cumulat, term_tot*np.random.rand())
            # a termination reaction marks from locus to the next looped site
            # as "done", but make sure to not write out of bounds for two
            # special cases
            si = max(term_starts[ti], 0)
            state[si] = (state[si] & ~PAIRED) | JUNC_LEFT_END
            ei = min(term_ends[ti], num_beads-1)
            state[ei] = (state[ei] & ~PAIRED) | JUNC_RIGHT_END
        else:
            # choose which of the sites will react
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
