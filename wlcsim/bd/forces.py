"""
Force fields used in our BD integrators.

Our BD code can be summarized as typically looking like:

.. code-block:: python

    X = init_function()
    for t_i in requested_times:
        X = integrator_step(f_all, X, prev_t, t_i)
        prev_t = t_i

where ``f_all`` takes ``X`` (size ``N_beads`` by ``N_dimensions``) and ``t``
and returns the total force on each of the beads at time ``t``.

This module holds all the functions that can be used to construct ``f_all``.
Broadly, they are named based on the source of the force:

- ``f_elas_*`` are elastic forces, they define the physics of the polymer being
  simulated. Linear, ring, and different types of network polymers have their
  own independent implementations.
- ``f_conf_*`` are confinement forces.
- ``f_tether_*`` tether individual beads of the polymer to an external
  structure.
- Others may be added in the future.

Notes
-----
Because they are in the "innermost loop", it is important that these functions
are just-in-time compile-able when possible.
"""
import sys
import inspect

import numpy as np
from numba import jit


# Append this to all f_tether_* docstrings to make it clear how they should
# actually be used.
extra_tether_docs = """

    Notes
    -----
    The force computed here is directed at the confinement from the inside
    only. To properly "tether" a bead to a surface, use the sum of the relevant
    ``f_conf + f_tether`` functions.
"""


@jit(nopython=True)
def f_conf_ellipse(x0, Aex, rx, ry, rz):
    """Compute soft (cubic) force due to elliptical confinement."""
    N, _ = x0.shape
    f = np.zeros(x0.shape)
    for j in range(N):
        conf = x0[j, 0]**2/rx**2 + x0[j, 1]**2/ry**2 + x0[j, 2]**2/rz**2
        if conf > 1:
            conf_u = np.array([
                -x0[j, 0]/rx**2, -x0[j, 1]/ry**2, -x0[j, 2]/rz**2
            ])
            conf_u = conf_u/np.linalg.norm(conf_u)
            # Steph's confinement from
            # https://journals.aps.org/pre/abstract/10.1103/PhysRevE.82.011913
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3)
    return f


@jit(nopython=True)
def f_tether_ellipse(x0, Aex, rx, ry, rz):
    """Compute soft (cubic) force towards an elliptical confinement."""
    N, _ = x0.shape
    f = np.zeros(x0.shape)
    for j in range(N):
        conf = x0[j, 0]**2/rx**2 + x0[j, 1]**2/ry**2 + x0[j, 2]**2/rz**2
        if conf < 1:
            conf_u = np.array([
                -x0[j, 0]/rx**2, -x0[j, 1]/ry**2, -x0[j, 2]/rz**2
            ])
            conf_u = conf_u/np.linalg.norm(conf_u)
            # Steph's confinement from
            # https://journals.aps.org/pre/abstract/10.1103/PhysRevE.82.011913
            f[j] += Aex*conf_u*np.power(np.sqrt(conf) - 1, 3)
    return f


@jit(nopython=True)
def f_elas_linear_rouse(x0, k_over_xi):
    """Compute spring forces on single, linear rouse polymer."""
    N, _ = x0.shape
    f = np.zeros(x0.shape)
    for j in range(1, N):
        for n in range(3):
            f[j, n] += -k_over_xi*(x0[j, n] - x0[j-1, n])
            f[j-1, n] += -k_over_xi*(x0[j-1, n] - x0[j, n])
    return f


@jit(nopython=True)
def f_elas_homolog_rouse(x0, N, loop_list, k_over_xi):
    r"""
    Compute spring forces on two "homologously linked linear rouse polymers.

    The two polymers are (homologously) at beads specified by loop_list. While
    each polymer is of length N beads, the beads that are hooked together move
    together identically, so you have `x0.shape[0] == 2*N - len(loop_list)`

    We assume that, per the spec in `~.homolog_points_to_loops_list`, a polymer
    of the shape

    .. code-block:: text

        ++   ^^^   '''''''   ""   =
          \ /   \ /       \ /  \ /
           .     .         .    .
          / \   / \       / \  / \
        --   ---   -------   --   -

    will be laid out in memory like

    .. code-block:: text

        --.---.-------.--.-++^^^'''''''""=

    If we call the lower polymer "polymer 1", and all the beads on the top row
    "polymer 2", then it makes sense to first compute the spring forces between
    adjacent beads in polymer 1 as if there was no "polymer 2". Then, we can
    compute all the forces between beads in "polymer 2" and the "." beads.
    Finally, for each "loop" (represented as symbols of different types above)
    in polymer 2, we can compute the usual linear forces along the polymer.

    """
    f = np.zeros(x0.shape)
    num_loops, _ = loop_list.shape
    # first get forces linearly along "first" polymer
    for j in range(1, N):
        f[j] += -k_over_xi*(x0[j] - x0[j-1])
        f[j-1] += -k_over_xi*(x0[j-1] - x0[j])
    # in what follows, "loop" means the stretches between homologous
    # connections, including from the homologous connection to the end of the
    # polymer
    for i in range(num_loops):
        # k1l, k1r denote the connected beads, so the array contains a "loop"
        # representing the beads of the second chain from k1l+1 to k1r-1
        # at array locations [k2l, k2r].
        k1l, k1r, k2l, k2r = loop_list[i]
        # k1l == -1 is flag for "end" loop (free end)
        # k2l > k2r => no beads in "second chain"
        if k1l >= 0 and k2l <= k2r:
            f[k1l] += -k_over_xi*(x0[k1l] - x0[k2l])
            f[k2l] += -k_over_xi*(x0[k2l] - x0[k1l])
        # k1r == N is flag for "end" loop (free end)
        if k1r < N and k2l <= k2r:
            f[k1r] += -k_over_xi*(x0[k1r] - x0[k2r])
            f[k2r] += -k_over_xi*(x0[k2r] - x0[k1r])
        N_loop = k2r - k2l + 1
        # analagous to loop above, since end-bead forces already done
        for k in range(1, N_loop):
            j = k2l+k
            f[j] += -k_over_xi*(x0[j] - x0[j-1])
            f[j-1] += -k_over_xi*(x0[j-1] - x0[j])
    return f


@jit(nopython=True)
def f_elas_linear_sswlc(x, u, e_b, gam, e_par, e_perp, eta, xi_u, xi_r):
    """Compute spring forces and torques on each bead of dsswlc."""
    N, _ = x.shape
    f = np.zeros(x.shape)
    t = np.zeros(x.shape)
    for i in range(0, N - 1):
        dx = x[i+1] - x[i]
        dx_par = dx @ u[i]
        dx_perp = dx - dx_par*u[i]
        cos_u1_u2 = u[i+1]@u[i]

        Gi = u[i+1] - cos_u1_u2*u[i] - eta*dx_perp
        Fi = -eta*e_b*Gi + e_par*(dx_par - gam)*u[i] + e_perp*dx_perp
        f[i] += Fi
        f[i + 1] -= Fi

        Gi = (u[i+1] - u[i]) - eta*dx_perp
        t[i] += e_b*Gi - eta*e_b*dx_par*Gi + eta*e_b*(1 - cos_u1_u2)*dx \
            - e_par*(dx_par - gam)*dx + e_perp*dx_par*dx_perp
        t[i+1] -= e_b*Gi
    return f, t


for name, obj in inspect.getmembers(sys.modules[__name__]):
    if not inspect.isfunction(obj) or not name.startswith('f_tether'):
        continue
    obj.__doc__ += extra_tether_docs
