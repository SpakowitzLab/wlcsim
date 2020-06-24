"""
BD of two homologously-linked linear Rouse polymers.

See `homolog_points_to_loops_list` for a description of how the homologous
network structure is handled.

A specialized
"""
from ..plot import PolymerViewer
from .rouse import recommended_dt
from .init_beads import init_homolog_rouse_conf
from .forces import f_conf_ellipse, f_elas_homolog_rouse, f_tether_ellipse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from numba import jit


# The constraints for dt here are the same as any other Rouse polymer.
recommended_dt = recommended_dt


def points_to_loops_list(N, homolog_points):
    r"""Create loops list for jit_rouse_linked.

    Parameters
    ----------
    N : int
        number of beads in each polymer if they were unbound
    homolog_points : (L,) array_like
        list of loci indices :math:`{l_i}_1^L, l_i\in[1,N]` (beads) that are
        homologously paired.

    Returns
    -------
    N_tot : int
        final size of the x0 array for the simulation.
    loop_list : (L,4) np.array
        data structure used by jit_rouse_linked to correctly compute the spring
        forces in the network.

    Notes
    -----
    You can use the output like

    >>> k1l, k1r, k2l, k2r = loop_list[i]

    where "k[i][l,r]" refers to the bead demarcating the [l]eft or [r]ight end
    of a given "homologous loop" on each of the polymers polymer
    :math:`i\in[1,2]`.

    Namely, the connectivity structure of our polymer "network" is that for
    each `k1l, k1r, k2l, k2r in loop_list`, the beads `k1l` thru `k1r` are
    linked sequentially (as are `k2l` thru `k2r`), but the beads `k1[l,r]` are
    linked to the beads `k2[l,r]`, respectively. This network structure can
    represent arbitrary "homologously linked" polymers, i.e. those of the form

    .. code-block:: text

        --  ---  -------  --  -
          \/   \/       \/  \/
          /\   /\       /\  /\
        --  ---  -------  --  -

    Notice in the example below that `loop_list[1:,1]` and `loop_list[:-1,0]`
    are both the same as `homolog_list`.

    Example
    -------
    Here's an example for a typical case. With 101 beads per polymer (so 100
    inter-bead segments), and one tenth of the beads bound to each other, we
    might get the following Poisson-space homolog points:

    >>> N = 101; FP = 0.1

    >>> homolog_points = np.where(np.random.rand(int(N)) < FP)[0]

    >>> homolog_points
    array([ 5, 18, 27, 59, 68, 82, 90, 94])

    These eight homologous junctions would lead to a simulation array length of
    `2*N - len(homolog_points) = 202 - 8 = 194`. So `N_tot` will be 194, and
    the `loop_list` will look like

    >>> N_tot, loop_list = rouse.homolog_points_to_loops_list(N, homolog_points)

    >>> loop_list
    array([[ -1,   5, 101, 105],
            [  5,  18, 106, 117],
            [ 18,  27, 118, 125],
            [ 27,  59, 126, 156],
            [ 59,  68, 157, 164],
            [ 68,  82, 165, 177],
            [ 82,  90, 178, 184],
            [ 90,  94, 185, 187],
            [ 94, 101, 188, 193]])

    This can be read as follows. The first four beads of each polymer are
    "free", so first row is saying that beads 0 thru 4 correspond to unlinked
    beads in polymer 1. bead 5 is the linked bead (has half the diffusivity).
    since the polymer is of length 101, if there were no connections, the first
    polymer would extend from bead 0 to bead 100, and the second from bead 101
    to bead 201. but since bead "5" is linked to (what would be) bead 106, that
    bead is omitted from the second polymer, and beads 101-105 are the "free
    end" of the second polymer.

    The second row corresponds to the first actual "loop". The two tethering
    points are beads 5 and 18 of the first polymer. What would have been beads
    106 and 119 (that's beads 5 and 18 of the second polymer) are omitted from
    the list of simulated beads. So bead 106 is linked to bead 5 and bead 117
    is linked to bead 18, and beads 106-117 are linked sequentially.

    Similarly, beads 18 and 118 (and 27 and 125) are linked, with beads 118-125
    being linked sequentially.
    """
    homolog_points = np.sort(np.array(homolog_points))
    num_points = len(homolog_points)
    num_loops = num_points + 1 # includes two end non-"loops" (aka free ends)
    loop_list = np.zeros((num_loops, 4))
    if num_loops == 1:
        return 2*int(N), np.array([[-1,N,N,2*N-1]]).astype(int)
    # keep track of length of "x0"
    curr_n = N
    # there's at least one connection point, thanks to num_loops==1 check above
    for i in range(0,num_points+1):
        k1l = homolog_points[i-1] if i > 0 else -1
        k1r = homolog_points[i] if i < num_points else N
        k2l = curr_n
        beads_in_loop = k1r - k1l - 1
        curr_n += beads_in_loop
        k2r = curr_n - 1
        loop_list[i,:] = np.array([k1l, k1r, k2l, k2r])
    return curr_n, loop_list.astype(int)

def split_homologs_X(X, N, loop_list):
    """Make the (Nt,Ntot,3) output of rouse_homologs to (2,Nt,N,3)."""
    N_t, N_tot, _ = X.shape
    X1 = X[:,:N,:].copy()
    X2 = np.zeros_like(X1)
    num_loops, _ = loop_list.shape
    if num_loops == 1:
        X2 = X[:,N:,:]
        return X1, X2
    curr_i2 = 0
    curr_i0 = N
    for i in range(num_loops):
        k1l, k1r, k2l, k2r = loop_list[i]
        num_beads = k1r - k1l - 1
        X2[:,curr_i2:curr_i2+num_beads] = X[:,curr_i0:curr_i0+num_beads]
        curr_i2 += num_beads
        if i < num_loops - 1:
            X2[:,curr_i2] = X[:,k1r]
            curr_i2 += 1
        curr_i0 += num_beads
    return X1, X2

def split_homolog_x(x0, N, loop_list):
    """Make an (N_tot,3) array from rouse_homologs into (2,N,3)."""
    N_tot, _ = x0.shape
    x1 = x0[:N,:].copy()
    x2 = np.zeros((N,3))
    num_loops, _ = loop_list.shape
    if num_loops == 1:
        x2 = x0[N:,:]
        return x1, x2
    curr_i2 = 0
    curr_i0 = N
    for i in range(num_loops):
        k1l, k1r, k2l, k2r = loop_list[i]
        num_beads = k1r - k1l - 1
        x2[curr_i2:curr_i2+num_beads] = x0[curr_i0:curr_i0+num_beads]
        curr_i2 += num_beads
        if i < num_loops - 1:
            x2[curr_i2] = x0[k1r]
            curr_i2 += 1
        curr_i0 += num_beads
    return x1, x2


def sim_loc(i, N, loop_list):
    """
    Return actual indices for bead "i" of the second polymer.

    Useful for indexing directly into the trucated simulation output array.

    Parameters
    ----------
    i : (M,) int, array_like
        beads index is requested for
    N : int
        actual number of beads in each polymer
    loop_list : (~N*FP, 4) int
        connectivity of the homologous polymers being simulated

    Returns
    -------
    iloc : (M,) Optional<int>, array_like
        index in simulation array (or None if the bead is a paired bead)
        for each bead in `i`.
    """
    iloc = np.array(i).copy()
    homolog_points = loop_list[1:,0]
    is_homolog = np.isin(iloc, homolog_points)
    if np.any(is_homolog):
        iloc = iloc.astype(object)
        iloc[is_homolog] = None
    for ii in np.where(~is_homolog)[0]:
        iloc[ii] = N + i[ii] - np.sum(homolog_points < i[ii])
    return iloc


def rouse_homologs(N, FP, L, b, D, Aex, rx, ry, rz, t, t_save=None,
                   tether_list=None):
    r"""
    BD simulation of two homologous yeast chromosomes in meiosis.

    An arbitrary elliptical-shaped confinement can be chosen, and in order to
    simulate synaptonemal formation in prophase, some fraction `FP` of the
    "homologous" beads of the polymer can be "rigidly tethered" to each other.
    The simulation treats these pairs as being one bead (with half the
    diffusion coefficient). The elastic force can then be communicated between
    the two polymers via the connecting loci.

    Depending on what part of prophase is being simulated, either (or both) of
    the telomeres or centromeres can be tethered to the confinement by
    including them in `tether_list`. Currently the confinement force and the
    tethering force are both cubic w.r.t. the distance from the confining
    ellipse with the same strength (`Aex`), and have the form described in
    :func:`f_conf_ellipse` (and :func:`f_tether_ellipse`).


    Parameters
    ----------
    N : int
        Number of beads to use in the discretized Rouse chain
    FP : float
        `FP`:math:`\in[0,1] is fraction of beads in the chain (out of `N`) that
        will be tethered to each other
    L : float
        Length of the chain (in desired length units)
    b : float
        Kuhn length of the chain (in same length units as L)
    D : float
        The diffusivity of a Kuhn length of the chain. This can often difficult
        to measure directly, but see the documentation for the function
        :func:`measured_D_to_rouse` for instructions on how to compute this
        value given clean measurements of short-time regular diffusion or
        Rouse-like MSD behavior
    Aex : float
        multiplicative prefactor controlling strength of confinment and
        tethering forces
    rx, ry, rz : float
        radius of confinement in x, y, z direction(s) respectively
    t : (M,) float, array_like
        The BD integrator will step time from t[0] to t[1] to t[2] to ...
    t_save : (m,) float, array_like
        `t_save`:math:`\subset``t` will be the time steps where the simulation
        output is saved. Default is `t_save = t`.
    tether_list : int, array_like
        list of bead indices that will be tethered to the confinement surface,
        for now, both polymers must be tethered at the same beads (but this
        would not be hard to change).

    Returns
    -------
    tether_list : List<int>
        Indices (into compressed `Ntot` output) that were tethered to the
        confinement surface
    loop_list : (~N*FP, 4) int
        Output of `points_to_loops_list` used to specify which points are
        tethered
    X : (len(t_save), Ntot, 3)
        output of BD simulation at each time in t_save. The number of output
        beads Ntot <= 2*N excludes the beads on the second polymer whose
        positions are determined by the fact that they are tethered to the
        first polymer. Use :func:`split_homologs_X` to reshape this output into
        shape `(2, len(t_save), N, 3)`.
    """
    if t_save is None:
        t_save = t

    homolog_points = np.where(np.random.rand(int(N)) < FP)[0]
    N_tot, loop_list = points_to_loops_list(N, homolog_points)

    # homolog points aren't in the simulation array for the second polymer
    tethered_p2 = set(tether_list) - set(loop_list[1:, 0])
    # calculate their actual positions in the simulation array
    tether_list = list(tether_list) + list(sim_loc(tethered_p2, N, loop_list))
    tether_list = np.sort(np.array(tether_list))

    x = _jit_rouse_homologs(int(N), int(N_tot), tether_list, loop_list, L, b, D, Aex, rx, ry, rz, t, t_save)
    return tether_list, loop_list, x


@jit(nopython=True)
def _jit_rouse_homologs(N, N_tot, tether_list, loop_list, L, b, D, Aex, rx, ry,
                        rz, t, t_save, x0=None):
    """
    "Inner loop" of rouse_homologs.

    Some typing stuff: needs N, N_tot as int, needs tether_list/loop_list as
    dtype == int64, and uses loop_list.shape[0] to get num_loops. This is in
    particular important for empty list cases, where you need to do stuff like

    >>> x = rouse._jit_rouse_homologs(int(N), int(2*N),
    >>>         np.array([]).astype(int),
    >>>         np.array([[]]).T.astype(int),
    >>>         L, b, D, Aex, R, R, R, t, t_save)

    """
    N = int(N)
    N_tot = int(N_tot)
    rtol = 1e-5
    # derived parameters
    L0 = L/(N-1) # length per bead
    bhat = np.sqrt(L0*b) # mean squared bond length of discrete gaussian chain
    Nhat = L/b # number of Kuhn lengths in chain
    Dhat = D*N/Nhat # diffusion coef of a discrete gaussian chain bead
    k_over_xi = 3*Dhat/bhat**2
    # initial position
    if x0 is None:
        x0 = init_homolog_rouse_conf(N, N_tot, loop_list, bhat, rx, ry, rz)
    # pre-alloc output
    if t_save is None:
        t_save = t
    x = np.zeros(t_save.shape + x0.shape)
    # setup for saving only requested time points
    save_i = 0
    if 0 == t_save[save_i]:
        x[0] = x0
        save_i += 1
    # extract from loop_list for repeated use in loop
    homolog_points = loop_list[1:, 0]
    # at each step i, we use data (x,t)[i-1] to create (x,t)[i]
    # in order to make it easy to pull into a new functin later, we'll call
    # t[i-1] "t0", old x (x[i-1]) "x0", and t[i]-t[i-1] "h".
    for i in range(1, len(t)):
        h = t[i] - t[i-1]
        dW = np.random.randn(*x0.shape)
        # -1 or 1, p=1/2
        S = 2*(np.random.rand() < 0.5) - 1
        # D = sigma^2/2 ==> sigma = np.sqrt(2*D)
        Fbrown = np.sqrt(2*Dhat/h)*(dW - S)
        # "homolog paired" beads have half the diffusivity
        Fbrown[homolog_points] = Fbrown[homolog_points]/2
        # estimate for slope at interval start
        K1 = f_conf_ellipse(x0, Aex, rx, ry, rz) + f_elas_homolog_rouse(x0, N, loop_list, k_over_xi) + Fbrown
        K1[tether_list] += f_tether_ellipse(x0[tether_list], Aex, rx, ry, rz)
        Fbrown = np.sqrt(2*Dhat/h)*(dW + S)
        # "homolog paired" beads have half the diffusivity
        Fbrown[homolog_points] = Fbrown[homolog_points]/2
        x1 = x0 + h*K1
        # estimate for slope at interval end
        K2 = f_conf_ellipse(x1, Aex, rx, ry, rz) + f_elas_homolog_rouse(x1, N, loop_list, k_over_xi) + Fbrown
        K2[tether_list] += f_tether_ellipse(x1[tether_list], Aex, rx, ry, rz)
        # average the slope estimates
        x0 = x0 + h * (K1 + K2)/2
        if np.abs(t[i] - t_save[save_i]) < rtol*np.abs(t_save[save_i]):
            x[save_i] = x0
            save_i += 1
    return x


def plot_homologs(X1, X2, loop_list, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    points = ax.scatter(*X1[loop_list[:-1, 1], :].T, c='r', s=50)
    ax.plot(*X1.T)
    ax.plot(*X2.T)
    return points


class HomologViewer(PolymerViewer):

    def __init__(self, r, N, loop_list):
        # will keep track of which *index* in the time axis should be used to
        # draw the current scene
        self.t = 0
        self.N = N
        self.loop_list = loop_list
        self.points = None
        # define the actual axes layout
        self.fig = plt.figure()
        # area to plot polymers
        self.ax3d = plt.axes([0.05, 0.15, 0.9, 0.8], facecolor='w', projection='3d')
        self.r = r
        # area to plot "slider" for selecting what time to plot
        self.ax = plt.axes([0.1, 0.025, 0.8, 0.05], facecolor='lightgoldenrodyellow')
        self.update_ax_limits()
        # slider to control what "time" is plotted, (0,1) is rescaled total
        # simulation time units
        self.slider_t = Slider(self.ax, 't', 0, 1, valinit=0)
        # handler to check if the "time" slider has been moved
        self.slider_t.on_changed(self.update_t)
        plt.show()

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, new_r):
        self._r = new_r
        self._x1, self._x2 = split_homologs_X(self._r, self.N, self.loop_list)
        self._num_time_points, self._num_beads, self._d = self._r.shape
        if hasattr(self, 'ax3d'):
            self.update_ax_limits()
            self.update_drawing()

    # should never be called before a Sim has been loaded successfully
    def update_drawing(self):
        num_polymers = 2
        # if num_polymers != len(self.ax3d.lines):
        self.ax3d.lines = []
        self.ax3d.plot(*self._x1[self.t].T, c='b')
        self.ax3d.plot(*self._x2[self.t].T, c='g')
        if self.points is not None:
            self.points.remove()
            del self.points
        self.points = self.ax3d.scatter(*self._x1[self.t,self.loop_list[:-1,1],:].T, c='r', s=50)
        # for i, r in enumerate([self._x1, self._x2]):
        #     x, y, z = (self.r[self.t,:,0], self.r[self.t,:,1],
        #                self.r[self.t,:,2])
        #     self.ax3d.lines[i].set_data(x, y)
        #     self.ax3d.lines[i].set_3d_properties(z)
        #     self.points.remove()
        #     del self.points
        #     self.points = self.ax3d.scatter(*self._x1[self.t,self.loop_list[:-1,1],:].T, c='r', s=50)
        self.fig.canvas.draw_idle()

    def update_ax_limits(self):
        self.ax3d.set_xlim((np.nanmin(self.r[:,:,0]), np.nanmax(self.r[:,:,0])))
        self.ax3d.set_ylim((np.nanmin(self.r[:,:,1]), np.nanmax(self.r[:,:,1])))
        self.ax3d.set_zlim((np.nanmin(self.r[:,:,2]), np.nanmax(self.r[:,:,2])))

    def update_t(self, val):
        if self.r is None:
            return
        # calculate value on slider
        new_t = int(round(val*self._num_time_points))
        # bound it to 0, num_time_points-1
        new_t = max(0, min(self._num_time_points-1, new_t))
        self.t = new_t
        self.update_drawing()



